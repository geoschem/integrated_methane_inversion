import os
import numpy as np
import xarray as xr
import pandas as pd
import warnings
import subprocess
import re
from scipy.sparse import csr_matrix
from src.inversion_scripts.utils import check_is_OH_element, check_is_BC_element

warnings.filterwarnings("ignore", category=UserWarning, module="xarray")

# common utilities for using different operators
def read_all_geoschem(all_strdate, gc_cache, n_elements, config, build_jacobian=False):
    """
    Call readgeoschem() for multiple dates in a loop.

    Arguments
        all_strdate    [list, str] : Multiple date strings
        gc_cache       [str]       : Path to GEOS-Chem output data
        build_jacobian [log]       : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?

    Returns
        dat            [dict]      : Dictionary of dictionaries. Each sub-dictionary is returned by read_geoschem()
    """

    dat = {}
    for strdate in all_strdate:
        dat[strdate] = read_geoschem(
            strdate, gc_cache, n_elements, config, build_jacobian
        )

    return dat


def read_geoschem(date, gc_cache, n_elements, config, build_jacobian=False):
    """
    Read GEOS-Chem data and save important variables to dictionary.

    Arguments
        date           [str]   : Date of interest, format "YYYYMMDD_HH"
        gc_cache       [str]   : Path to GEOS-Chem output data
        build_jacobian [log]   : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?

    Returns
        dat            [dict]  : Dictionary of important variables from GEOS-Chem:
                                    - CH4
                                    - Latitude
                                    - Longitude
                                    - PEDGE
    """

    UseGCHP = config['UseGCHP']
    # Assemble file paths to GEOS-Chem output collections for input data
    file_species = f"GEOSChem.SpeciesConc.{date}00z.nc4"
    file_pedge = f"GEOSChem.LevelEdgeDiags.{date}00z.nc4"

    # Read lat, lon, CH4 from the SpeciecConc collection
    filename = f"{gc_cache}/{file_species}"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        gc_data = xr.open_dataset(filename)
    # Drop the "anchor" variable if it exists
    if 'anchor' in gc_data:
        gc_data = gc_data.drop_vars('anchor')
    if UseGCHP:
        LON = gc_data["lons"].values
        LAT = gc_data["lats"].values
        CH4 = gc_data["SpeciesConcVV_CH4"].values[0, ...]
        CH4 = CH4 * 1e9  # Convert to ppb
        CH4 = np.einsum("lfyx->fyxl", CH4)
    else:
        LON = gc_data["lon"].values
        LAT = gc_data["lat"].values
        CH4 = gc_data["SpeciesConcVV_CH4"].values[0, :, :, :]
        CH4 = CH4 * 1e9  # Convert to ppb
        CH4 = np.einsum("lij->jil", CH4)
    gc_data.close()

    # Read PEDGE from the LevelEdgeDiags collection
    filename = f"{gc_cache}/{file_pedge}"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        gc_data = xr.open_dataset(filename)
    # Drop the "anchor" variable if it exists
    if 'anchor' in gc_data:
        gc_data = gc_data.drop_vars('anchor')
    if UseGCHP:
        PEDGE = gc_data["Met_PEDGE"].values[0, ...]
        PEDGE = np.einsum("lfyx->fyxl", PEDGE)
    else:
        PEDGE = gc_data["Met_PEDGE"].values[0, :, :, :]
        PEDGE = np.einsum("lij->jil", PEDGE)
    gc_data.close()

    # Store GEOS-Chem data in dictionary
    dat = {}
    dat["lon"] = LON
    dat["lat"] = LAT
    dat["PEDGE"] = PEDGE
    dat["CH4"] = CH4

    # If need to construct Jacobian, read sensitivity data from GEOS-Chem perturbation simulations
    if build_jacobian:

        elements = range(n_elements)
        ntracers = config["NumJacobianTracers"]
        opt_OH = config["OptimizeOH"]
        opt_BC = config["OptimizeBCs"]
        is_Regional = config["isRegional"]

        num_BC = 4
        if is_Regional:
            num_OH = 1
        else:
            num_OH = 2

        n_base_runs = (
            n_elements - int(opt_OH * num_OH) - (int(opt_BC) * num_BC)
        ) / ntracers

        nruns = (
            np.ceil(n_base_runs).astype(int)
            + (int(opt_OH) * num_OH)
            + (int(opt_BC) * num_BC)
        )

        # Dictionary that stores mapping of state vector elements to
        # perturbation simulation numbers
        pert_simulations_dict = {}
        for e in elements:
            # State vector elements are numbered 1..nelements
            sv_elem = e + 1

            is_OH_element = check_is_OH_element(
                sv_elem, n_elements, opt_OH, is_Regional
            )
            is_BC_element = check_is_BC_element(
                sv_elem, n_elements, opt_OH, opt_BC, is_OH_element, is_Regional
            )
            # Determine which run directory to look in
            if is_OH_element:
                if is_Regional:
                    run_number = nruns
                else:
                    num_back = n_elements % sv_elem
                    run_number = nruns - num_back
            elif is_BC_element:
                num_back = n_elements % sv_elem
                run_number = nruns - num_back
            else:
                run_number = np.ceil(sv_elem / ntracers).astype(int)

            run_num = str(run_number).zfill(4)

            # add the element to the dictionary for the relevant simulation number
            if run_num not in pert_simulations_dict:
                pert_simulations_dict[run_num] = [sv_elem]
            else:
                pert_simulations_dict[run_num].append(sv_elem)

        gc_date = pd.to_datetime(date, format="%Y%m%d_%H")
        ds_all = [
            concat_tracers(k, gc_date, config, v, n_elements)
            for k, v in pert_simulations_dict.items()
        ]

        ds_all = [ds.load() for ds in ds_all]

        ds_sensi = xr.concat(ds_all, "element")

        sensitivities = ds_sensi["ch4"].values
        if UseGCHP:
            # Reshape so the data have dimensions (nf, Ydim, Xdim, lev, grid_element)
            sensitivities = np.einsum("klfyx->fyxlk", sensitivities)
        else:
            # Reshape so the data have dimensions (lon, lat, lev, grid_element)
            sensitivities = np.einsum("klji->ijlk", sensitivities)
        dat["jacobian_ch4"] = sensitivities

        # get emis base, which is also BC base
        ds_emis_base = concat_tracers(
            "0001", gc_date, config, [0], n_elements, baserun=True
        )

        if UseGCHP:
            dat["emis_base_ch4"] = np.einsum("klfyx->fyxlk", ds_emis_base["ch4"].values)
        else:
            dat["emis_base_ch4"] = np.einsum("klji->ijlk", ds_emis_base["ch4"].values)

        # get OH base, run RunName_0000
        # it's always here whether OptimizeOH is true or not
        # so we can keep it here for convenience
        ds_oh_base = concat_tracers(
            "0000", gc_date, config, [0], n_elements, baserun=True
        )

        if UseGCHP:
            dat["oh_base_ch4"] = np.einsum("klfyx->fyxlk", ds_oh_base["ch4"].values)
        else:
            dat["oh_base_ch4"] = np.einsum("klji->ijlk", ds_oh_base["ch4"].values)

    return dat


def concat_tracers(run_id, gc_date, config, sv_elems, n_elements, baserun=False):
    """
    Concatenate CH4 tracers from all jacobian GEOS-Chem simulations.
    Tracers are assigned a new dimension: "element"

    Arguments
        run_id     [str]         : ID for Jacobian GEOS-Chem run, e.g. "0001"
        gc_date    [pd.Datetime] : date object, specifies Ymd_h
        config     [dict]        : dictionary of IMI config file
        sv_elems   [list]        : list of state vector element tracers in this simulations
        n_elements [int]         : number of state vector elements in this inversion
        baserun    [bool]        : If True, only the base variable in the simulation will
                                 be opened, and the function will just return this one
                                 variable instead of concatenating all elements. Used to
                                 get the base for calculating the sensitivities.

    Returns
        ds_concat [xarray.Dataset] : dataset of all Jacobian CH4 at this timestep for
                                     all tracer runs. Has dimensions
                                        - lat
                                        - lon
                                        - lev
                                        - element


    """
    prefix = os.path.expandvars(
        config["OutputPath"] + "/" + config["RunName"] + "/jacobian_runs"
    )
    j_dir = f"{prefix}/{config['RunName']}_{run_id}/OutputDir"
    file_stub = gc_date.strftime("GEOSChem.SpeciesConc.%Y%m%d_0000z.nc4")
    filepath = os.path.join(j_dir, file_stub)
    
    # Construct the list of CH4 vars to request
    keepvars = [f"SpeciesConcVV_CH4_{i:04}" for i in sv_elems]
    if len(keepvars) == 1:
        is_Regional = config["isRegional"]
        is_OH_element = check_is_OH_element(
            sv_elems[0], n_elements, config["OptimizeOH"], is_Regional
        )
        is_BC_element = check_is_BC_element(
            sv_elems[0],
            n_elements,
            config["OptimizeOH"],
            config["OptimizeBCs"],
            is_OH_element,
            is_Regional,
        )
        if is_OH_element or is_BC_element:
            keepvars = ["SpeciesConcVV_CH4"]

    if baserun:
        keepvars = ["SpeciesConcVV_CH4"]

    # It would fail if open all variables with chunks with GCHP,
    # as ncontact is duplicate for GCHP output dimensions
    with xr.open_dataset(filepath, decode_cf=False) as tmp:
        other_vars = [v for v in tmp.variables if "SpeciesConcVV_CH4" not in v]
    
    # Open only these variables
    with xr.open_dataset(
        filepath,
        chunks="auto",
        drop_variables=other_vars,
    ) as dsmf:
        try:
            dsmf = dsmf.isel(time=gc_date.hour, drop=True)
        except Exception as e:
            print(f"Run id {run_id}. Failed at {gc_date} with error: {e}", flush=True)
            raise

        ds_concat = xr.concat((dsmf[v] for v in keepvars), dim="element").rename("ch4")
        ds_concat = ds_concat.to_dataset(name="ch4").assign_attrs(dsmf.attrs)

    if not baserun:
        ds_concat = ds_concat.assign_coords({"element": sv_elems})
    else:
        ds_concat = ds_concat.assign_coords({"element": [0]})

    return ds_concat

def get_gridcell_list(lons, lats):
    """
    Create a 2d array of dictionaries, with each dictionary representing a GC gridcell.
    Dictionaries also initialize the fields necessary to store for tropomi data
    (eg. methane, time, p_sat, etc.)

    Arguments
        lons     [float[]]      : list of gc longitudes for region of interest
        lats     [float[]]      : list of gc latitudes for region of interest

    Returns
        gridcells [dict[][]]     : 2D array of dicts representing a gridcell
    """
    # create array of dictionaries to represent gridcells
    gridcells = []
    for i in range(len(lons)):
        for j in range(len(lats)):
            gridcells.append(
                {
                    "lat": lats[j],
                    "lon": lons[i],
                    "iGC": i,
                    "jGC": j,
                    "methane": [],
                    "p_sat": [],
                    "dry_air_subcolumns": [],
                    "apriori": [],
                    "avkern": [],
                    "time": [],
                    "overlap_area": [],
                    "lat_sat": [],
                    "lon_sat": [],
                    "observation_count": 0,
                    "observation_weights": [],
                }
            )
    gridcells = np.array(gridcells).reshape(len(lons), len(lats))
    return gridcells

def get_gc_lat_lon(gc_cache, start_date):
    """
    get dictionary of lat/lon values for gc gridcells

    Arguments
        gc_cache    [str]   : path to gc data
        start_date  [str]   : start date of the inversion

    Returns
        output      [dict]  : dictionary with the following fields:
                                - lat : list of GC latitudes
                                - lon : list of GC longitudes
    """
    gc_ll = {}
    date = pd.to_datetime(start_date).strftime("%Y%m%d_%H")
    file_species = f"GEOSChem.SpeciesConc.{date}00z.nc4"
    filename = f"{gc_cache}/{file_species}"
    gc_data = xr.open_dataset(filename)
    gc_ll["lon"] = gc_data["lon"].values
    gc_ll["lat"] = gc_data["lat"].values

    gc_data.close()
    return gc_ll


def merge_pressure_grids(p_sat, p_gc):
    """
    Merge TROPOMI & GEOS-Chem vertical pressure grids

    Arguments
        p_sat   [float]    : Pressure edges from TROPOMI (13 edges)     <--- 13-1 = 12 pressure layers
        p_gc    [float]    : Pressure edges from GEOS-Chem (48 edges)   <--- 48-1 = 47 pressure layers

    Returns
        merged  [dict]     : Merged grid dictionary
                                - p_merge       : merged pressure-edge grid
                                - data_type     : for each pressure edge in the merged grid, is it from GEOS-Chem or TROPOMI?
                                - edge_index    : indexes of pressure edges
                                - first_gc_edge : index of first GEOS-Chem pressure edge in the merged grid
    """

    # Combine p_sat and p_gc into merged vertical pressure grid
    p_merge = np.zeros(len(p_sat) + len(p_gc))
    p_merge.fill(np.nan)
    data_type = np.zeros(len(p_sat) + len(p_gc), dtype=int)
    data_type.fill(-99)
    edge_index = []
    i = 0
    j = 0
    k = 0
    while (i < len(p_sat)) or (j < len(p_gc)):
        if i == len(p_sat):
            p_merge[k] = p_gc[j]
            data_type[k] = 2  # geos-chem edge
            j = j + 1
            k = k + 1
            continue
        if j == len(p_gc):
            p_merge[k] = p_sat[i]
            data_type[k] = 1  # tropomi edge
            edge_index.append(k)
            i = i + 1
            k = k + 1
            continue
        if p_sat[i] >= p_gc[j]:
            p_merge[k] = p_sat[i]
            data_type[k] = 1  # tropomi edge
            edge_index.append(k)
            i = i + 1
            k = k + 1
        else:
            p_merge[k] = p_gc[j]
            data_type[k] = 2  # geos-chem edge
            j = j + 1
            k = k + 1

    # Find the first GEOS-Chem pressure edge
    first_gc_edge = -99
    for i in range(len(p_sat) + len(p_gc) - 1):
        if data_type[i] == 2:
            first_gc_edge = i
            break

    # Save data to dictionary
    merged = {}
    merged["p_merge"] = p_merge
    merged["data_type"] = data_type
    merged["edge_index"] = edge_index
    merged["first_gc_edge"] = first_gc_edge

    return merged


def remap(gc_CH4, data_type, p_merge, edge_index, first_gc_edge):
    """
    Remap GEOS-Chem methane to the TROPOMI vertical grid.

    Arguments
        gc_CH4        [float]   : Methane from GEOS-Chem
        p_merge       [float]   : Merged TROPOMI + GEOS-Chem pressure levels, from merge_pressure_grids()
        data_type     [int]     : Labels for pressure edges of merged grid. 1=TROPOMI, 2=GEOS-Chem, from merge_pressure_grids()
        edge_index    [int]     : Indexes of pressure edges, from merge_pressure_grids()
        first_gc_edge [int]     : Index of first GEOS-Chem pressure edge in merged grid, from merge_pressure_grids()

    Returns
        sat_CH4       [float]   : GEOS-Chem methane in TROPOMI pressure coordinates
    """

    # Define CH4 concentrations in the layers of the merged pressure grid
    CH4 = np.zeros(
        len(p_merge) - 1,
    )
    CH4.fill(np.nan)
    k = 0
    for i in range(first_gc_edge, len(p_merge) - 1):
        CH4[i] = gc_CH4[k]
        if data_type[i + 1] == 2:
            k = k + 1
    if first_gc_edge > 0:
        CH4[:first_gc_edge] = CH4[first_gc_edge]

    # Calculate the pressure-weighted mean methane for each TROPOMI layer
    delta_p = p_merge[:-1] - p_merge[1:]
    sat_CH4 = np.zeros(12)
    sat_CH4.fill(np.nan)
    for i in range(len(edge_index) - 1):
        start = edge_index[i]
        end = edge_index[i + 1]
        sum_conc_times_pressure_weights = sum(CH4[start:end] * delta_p[start:end])
        sum_pressure_weights = sum(delta_p[start:end])
        sat_CH4[i] = (
            sum_conc_times_pressure_weights / sum_pressure_weights
        )  # pressure-weighted average

    return sat_CH4


def remap_sensitivities(sensi_lonlat, data_type, p_merge, edge_index, first_gc_edge):
    """
    Remap GEOS-Chem sensitivity data (from perturbation simulations) to the TROPOMI vertical grid.

    Arguments
        sensi_lonlat  [float]   : Sensitivity data from GEOS-Chem perturbation runs, for a specific lon/lat; has dims (lev, grid_element)
        p_merge       [float]   : Combined TROPOMI + GEOS-Chem pressure levels, from merge_pressure_grids()
        data_type     [int]     : Labels for pressure edges of merged grid. 1=TROPOMI, 2=GEOS-Chem, from merge_pressure_grids()
        edge_index    [int]     : Indexes of pressure edges, from merge_pressure_grids()
        first_gc_edge [int]     : Index of first GEOS-Chem pressure edge in merged grid, from merge_pressure_grids()

    Returns
        sat_deltaCH4  [float]   : GEOS-Chem methane sensitivities in TROPOMI pressure coordinates
    """

    # Define DeltaCH4 in the layers of the merged pressure grid, for all perturbed state vector elements
    n_elem = sensi_lonlat.shape[1]
    deltaCH4 = np.zeros((len(p_merge) - 1, n_elem))
    deltaCH4.fill(np.nan)
    k = 0
    for i in range(first_gc_edge, len(p_merge) - 1):
        deltaCH4[i, :] = sensi_lonlat[k, :]
        if data_type[i + 1] == 2:
            k = k + 1
    if first_gc_edge > 0:
        deltaCH4[:first_gc_edge, :] = deltaCH4[first_gc_edge, :]

    # Calculate the weighted mean DeltaCH4 for each layer, for all perturbed state vector elements
    delta_p = p_merge[:-1] - p_merge[1:]
    delta_ps = np.transpose(np.tile(delta_p, (n_elem, 1)))
    sat_deltaCH4 = np.zeros((12, n_elem))
    sat_deltaCH4.fill(np.nan)
    for i in range(len(edge_index) - 1):
        start = edge_index[i]
        end = edge_index[i + 1]
        sum_conc_times_pressure_weights = np.sum(
            deltaCH4[start:end, :] * delta_ps[start:end, :], 0
        )
        sum_pressure_weights = np.sum(delta_p[start:end])
        sat_deltaCH4[i, :] = (
            sum_conc_times_pressure_weights / sum_pressure_weights
        )  # pressure-weighted average

    return sat_deltaCH4

def _ensure_descending_pedges(edges):
    # edges: (N, K+1). Flip to descending (surface -> TOA) if needed.
    need_flip = edges[:, 0] < edges[:, 1]
    if np.any(need_flip):
        edges = edges.copy()
        edges[need_flip] = edges[need_flip, ::-1]
    return edges

def remapping_weights(p_sat_edges, p_gc_edges):
    """
    p_sat_edges : (N, S+1)  TROPOMI pressure edges
    p_gc_edges  : (N, G+1)  GEOS-Chem pressure edges
    Returns     : (N, S, G) weights; for each n,s, sum_g W[n,s,g] == 1 (or 0 if no overlap)
    """
    p_sat_edges = _ensure_descending_pedges(np.asarray(p_sat_edges))
    p_gc_edges  = _ensure_descending_pedges(np.asarray(p_gc_edges))

    sat_bot = p_sat_edges[:, :-1]    # (N, S)
    sat_top = p_sat_edges[:,  1:]    # (N, S)
    gc_bot  = p_gc_edges[:, :-1]     # (N, G)
    gc_top  = p_gc_edges[:,  1:]     # (N, G)

    # Overlap thickness for descending pressures: max(0, min(bots) - max(tops))
    overlap = np.maximum(
        0.0,
        np.minimum(sat_bot[:, :, None], gc_bot[:, None, :])
        - np.maximum(sat_top[:, :, None], gc_top[:, None, :])
    )  # (N, S, G)
    
    # Bottom-fill: only layers *entirely below* GC domain
    # A TROPOMI layer is entirely below GC 
    # if its TOP pressure > GC surface bottom edge pressure.
    gc_bottom_edge = p_gc_edges[:, 0]                  # (N,) surface pressure edge of GC
    below_mask = sat_top >= gc_bottom_edge[:, None]    # (N, S) True -> bottom-fill

    if np.any(below_mask):
        sat_thick = (sat_bot - sat_top)                    # (N,S)
        overlap[below_mask, :]  = 0.0
        overlap[below_mask, 0]  = sat_thick[below_mask]

    denom = overlap.sum(axis=2, keepdims=True)  # (N, S, 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        W = np.where(denom > 0, overlap / denom, 0.0)
        
    return W

def nearest_loc(query_location, reference_grid, tolerance=0.5):
    """Find the index of the nearest grid location to a query location, with some tolerance."""

    distances = np.abs(reference_grid - query_location)
    ind = distances.argmin()
    if distances[ind] >= tolerance:
        return np.nan
    else:
        return ind

def sph2cart(pl, degrees=True):
    if degrees:
        pl = np.deg2rad(pl)

    xyz_shape = list(pl.shape)
    xyz_shape[-1] = 3
    xyz = np.zeros(xyz_shape)

    xyz[..., 0] = np.cos(pl[..., 0]) * np.cos(pl[..., 1])
    xyz[..., 1] = np.cos(pl[..., 0]) * np.sin(pl[..., 1])
    xyz[..., 2] = np.sin(pl[..., 0]) # pl[..., 0] is lat

    return xyz


def spherical_angle(v1, v2, v3):
    p = np.cross(v1, v2)
    q = np.cross(v1, v3)
    #d = np.sqrt(np.dot(p, p) * np.dot(q, q))
    pp = np.sum(p*p, axis=-1)
    qq = np.sum(q*q, axis=-1)
    pq = np.sum(p*q, axis=-1)
    d = np.sqrt(pp * qq)
    return np.arccos(pq/d)


def spherical_excess_area(ll, ul, ur, lr, radius=6371000.):
    v1 = sph2cart(ll)
    v2 = sph2cart(lr)
    v3 = sph2cart(ul)
    a1 = spherical_angle(v1, v2, v3)

    v1 = sph2cart(lr)
    v2 = sph2cart(ur)
    v3 = sph2cart(ll)
    a2 = spherical_angle(v1, v2, v3)

    v1 = sph2cart(ur)
    v2 = sph2cart(ul)
    v3 = sph2cart(lr)
    a3 = spherical_angle(v1, v2, v3)

    v1 = sph2cart(ul)
    v2 = sph2cart(ur)
    v3 = sph2cart(ll)
    a4 = spherical_angle(v1, v2, v3)

    area = (a1 + a2 + a3 + a4 - 2.*np.pi) * radius * radius
    return area

def create_SCRIP_grid(TROPOMI, sat_ind, save_pth):
    grid_rank = 1
    grid_corners = 4
    grid_center_lat = TROPOMI["latitude"]
    grid_center_lon = TROPOMI["longitude"]
    grid_corner_lon = TROPOMI["longitude_bounds"]
    grid_corner_lat = TROPOMI["latitude_bounds"]
    
    grid_imask = np.zeros_like(grid_center_lat, dtype=int)
    
    grid_imask[sat_ind] = 1

    grid_size = len(grid_center_lat.flatten())
    grid_imask = grid_imask.flatten()
    
    # calculate grid_area
    ll = np.stack([grid_corner_lat[...,0], grid_corner_lon[...,0]], axis=-1)
    lr = np.stack([grid_corner_lat[...,1], grid_corner_lon[...,1]], axis=-1)
    ul = np.stack([grid_corner_lat[...,2], grid_corner_lon[...,2]], axis=-1)
    ur = np.stack([grid_corner_lat[...,3], grid_corner_lon[...,3]], axis=-1)
    grid_area = spherical_excess_area(ll, ul, ur, lr)
    grid_area = grid_area.flatten()

    # flatten coordinates
    grid_center_lat = grid_center_lat.reshape((grid_size,))
    grid_center_lon = grid_center_lon.reshape((grid_size,))
    grid_corner_lon = grid_corner_lon.reshape((grid_size,4))
    grid_corner_lat = grid_corner_lat.reshape((grid_size,4))
    
    # Identify invalid cells
    valid_corners = ~np.isnan(grid_corner_lat) & ~np.isnan(grid_corner_lon)

    num_valid_corners = np.sum(valid_corners, axis=1)
    invalid_cells = num_valid_corners < 4

    # Also check for duplicate corners (degenerate cells)
    for i in range(grid_size):
        corners = np.stack([grid_corner_lat[i], grid_corner_lon[i]], axis=-1)
        if len(np.unique(corners, axis=0)) < 4:
            invalid_cells[i] = True

    # Set grid_imask=0 for all invalid cells
    grid_imask[invalid_cells] = 0
    
    # make dataset
    SCRIP_ds = xr.Dataset(
        data_vars={
            "grid_dims": xr.DataArray(
                [grid_size], dims=("grid_rank",)
            ),
            "grid_center_lat": xr.DataArray(
                grid_center_lat, dims=("grid_size",),
                attrs={"units": "degrees"}
            ),
            "grid_center_lon": xr.DataArray(
                grid_center_lon, dims=("grid_size",),
                attrs={"units": "degrees"}
            ),
            "grid_corner_lat": xr.DataArray(
                grid_corner_lat, dims=("grid_size", "grid_corners"),
                attrs={"units": "degrees", "_FillValue": -9999.}
            ),
            "grid_corner_lon": xr.DataArray(
                grid_corner_lon, dims=("grid_size", "grid_corners"),
                attrs={"units": "degrees", "_FillValue": -9999.}
            ),
            "grid_imask": xr.DataArray(
                grid_imask, dims=("grid_size",),
                attrs={"_FillValue": -9999}
            ),
            "grid_area": xr.DataArray(
                grid_area, dims=("grid_size",),
                attrs={"units": "m^2", "long_name": "grid box area"}
            ),
        },
        coords={
            "grid_size": np.arange(grid_size),
            "grid_corners": np.arange(grid_corners),
            "grid_rank": np.arange(grid_rank),
        }
    )
    
    if save_pth is not None:
        print("Saving file {}".format(save_pth))
        SCRIP_ds.to_netcdf(
            save_pth,
            encoding={
                v: {"zlib": True, "complevel": 1} for v in SCRIP_ds.data_vars
            },
        )

def create_ESMF_regridding_weights(TROPOMI, filename, sat_ind, CSgridDir, gridspec_path):
    """Generate the regridding weights from TROPOMI grids to GCHP grids

    Args:
        TROPOMI (dict): TROPOMI dataset
        filename (str): file path to TROPOMI dataset
        sat_ind (array): filtered satellite indices
        CSgridDir (str): directory path to CS_grids
        gridspec_path (str): file path to simulation gridspec file
    """
    
    # create SCRIP grid file
    shortname = re.split(r"\/", filename)[-1]
    date = re.split(r"\.", shortname)[0]

    SCRIP_grid_fpath = f"{date}_SCRIP_grid.nc"
    regrid_weight_fpath = f"{date}_regrid_weights.nc"

    # check if SCRIP grid file exists
    if not os.path.exists(os.path.join(CSgridDir, SCRIP_grid_fpath)):
        create_SCRIP_grid(TROPOMI, sat_ind, os.path.join(CSgridDir, SCRIP_grid_fpath))

    # check if regridding weights file exists
    if not os.path.exists(os.path.join(CSgridDir, regrid_weight_fpath)):
        # remove stale .esmf.nc if it exists
        tmp_esmf = os.path.join(CSgridDir, ".esmf.nc")
        if os.path.exists(tmp_esmf):
            os.remove(tmp_esmf)
        os.environ["ESMF_TMP"] = CSgridDir
        os.environ["TMPDIR"] = CSgridDir
        print(f"Running ESMF_RegridWeightGen for {date}...")
        subprocess.run([
            "mpirun", "-n", "1",
            "ESMF_RegridWeightGen",
            "-s", SCRIP_grid_fpath,
            "-d", gridspec_path,
            "-m", "conserve",
            "--ignore_unmapped",
            "-w", regrid_weight_fpath
        ], check=True, cwd=CSgridDir, stdout=subprocess.DEVNULL)
    return regrid_weight_fpath

def get_overlap_area_CSgrid(TROPOMI, filename, sat_ind, CSgridDir, 
                            gridspec_path, GC_shape):
    # create regridding weights first
    regrid_weight_fpath = create_ESMF_regridding_weights(TROPOMI, filename, sat_ind, 
                                                         CSgridDir, gridspec_path)
    with xr.open_dataset(os.path.join(CSgridDir, regrid_weight_fpath)) as regrid_weight_ds:
        src_ind = regrid_weight_ds['col'].values - 1
        dst_ind = regrid_weight_ds['row'].values - 1
        regrid_weights = regrid_weight_ds['S'].values
        dst_area = regrid_weight_ds['area_b'].values
    Sat_shape = TROPOMI["longitude"].shape
    n_dst = np.product(GC_shape)
    n_obs = np.product(Sat_shape)
    regrid_weights_csr = csr_matrix((regrid_weights, (dst_ind, src_ind)), 
                                    shape=(n_dst, n_obs))
    
    # --- Vectorized obs processing ---
    # flatten indices to valid obs
    sat_ind_flat = np.ravel_multi_index(sat_ind, Sat_shape)
    sat_mask = np.zeros(np.product(Sat_shape), dtype=bool)
    sat_mask[sat_ind_flat] = True
    
    # Subset sparse matrix to valid obs
    W = regrid_weights_csr[:, sat_mask]
    # The regridding weights are defined as: w_ij = f_ij * A_si / A_dj
    # where w_ij is the regridding weights from source grid cell i to destination cell j
    # f_ij is the fraction of source grid cell i contributing to destination grid cell j
    # based on grid cell area overlap
    # A_si is the grid cell area for source grid i
    # A_dj is the grid cell area for destination grid j
    # Thus, overlap_area = w_ij * A_dj
    overlap_area = W.multiply(dst_area[:, None]).tocsr()  # (n_dst × n_valid_obs)
    
    return overlap_area
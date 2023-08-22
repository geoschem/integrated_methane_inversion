import numpy as np
import xarray as xr
import pandas as pd

# common utilities for using different operators

def read_geoschem(date, gc_cache, build_jacobian=False, sensi_cache=None):
    """
    Read GEOS-Chem data and save important variables to dictionary.

    Arguments
        date           [str]   : Date of interest
        gc_cache       [str]   : Path to GEOS-Chem output data
        build_jacobian [log]   : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        sensi_cache    [str]   : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data

    Returns
        dat            [dict]  : Dictionary of important variables from GEOS-Chem:
                                    - CH4
                                    - Latitude
                                    - Longitude
                                    - PEDGE
    """

    # Assemble file paths to GEOS-Chem output collections for input data
    file_species = f"GEOSChem.SpeciesConc.{date}00z.nc4"
    file_pedge = f"GEOSChem.LevelEdgeDiags.{date}00z.nc4"

    # Read lat, lon, CH4 from the SpeciecConc collection
    filename = f"{gc_cache}/{file_species}"
    gc_data = xr.open_dataset(filename)
    LON = gc_data["lon"].values
    LAT = gc_data["lat"].values
    CH4 = gc_data["SpeciesConcVV_CH4"].values[0, :, :, :]
    CH4 = CH4 * 1e9  # Convert to ppb
    CH4 = np.einsum("lij->jil", CH4)
    gc_data.close()

    # Read PEDGE from the LevelEdgeDiags collection
    filename = f"{gc_cache}/{file_pedge}"
    gc_data = xr.open_dataset(filename)
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
        sensitivity_data = xr.open_dataset(f"{sensi_cache}/sensi_{date}.nc")
        sensitivities = sensitivity_data["Sensitivities"].values
        # Reshape so the data have dimensions (lon, lat, lev, grid_element)
        sensitivities = np.einsum("klji->ijlk", sensitivities)
        sensitivity_data.close()
        dat["Sensitivities"] = sensitivities

    return dat

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


def read_all_geoschem(all_strdate, gc_cache, build_jacobian=False, sensi_cache=None):
    """
    Call readgeoschem() for multiple dates in a loop.

    Arguments
        all_strdate    [list, str] : Multiple date strings
        gc_cache       [str]       : Path to GEOS-Chem output data
        build_jacobian [log]       : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        sensi_cache    [str]       : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data

    Returns
        dat            [dict]      : Dictionary of dictionaries. Each sub-dictionary is returned by read_geoschem()
    """

    dat = {}
    for strdate in all_strdate:
        dat[strdate] = read_geoschem(strdate, gc_cache, build_jacobian, sensi_cache)

    return dat


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


def nearest_loc(query_location, reference_grid, tolerance=0.5):
    """Find the index of the nearest grid location to a query location, with some tolerance."""

    distances = np.abs(reference_grid - query_location)
    ind = distances.argmin()
    if distances[ind] >= tolerance:
        return np.nan
    else:
        return ind

import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime
import gc
import pygeohash as pgh
from shapely.geometry import Polygon
from src.inversion_scripts.utils import (
    read_and_filter_satellite,
    mixing_ratio_conv_factor,
    get_strdate,
    check_is_OH_element,
    check_is_BC_element,
)
from src.inversion_scripts.operators.operator_utilities import (
    get_gc_lat_lon,
    read_all_geoschem,
    merge_pressure_grids,
    remap,
    remap_sensitivities,
    remapping_weights,
    get_gridcell_list,
    nearest_loc,
    get_overlap_area_CSgrid,
)
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="xarray")

def apply_average_satellite_operator(
    filename,
    species,
    satellite_product,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    build_jacobian,
    period_i,
    config,
    use_water_obs=False
):
    """
    Apply the averaging satellite operator to map GEOS-Chem data to satellite observation space.

    Arguments
        filename          [str]        : satellite netcdf data file to read
        species           [str]        : The species (CH4 or CO2) to use
        satellite_product [str]        : "BlendedTROPOMI", "TROPOMI", or "Other", specifying the data used in the inversion.
        n_elements        [int]        : Number of state vector elements
        gc_startdate      [datetime64] : First day of inversion period, for GEOS-Chem and satellite
        gc_enddate        [datetime64] : Last day of inversion period, for GEOS-Chem and satellite
        xlim              [float]      : Longitude bounds for simulation domain
        ylim              [float]      : Latitude bounds for simulation domain
        gc_cache          [str]        : Path to GEOS-Chem output data
        build_jacobian    [log]        : Are we trying to map GEOS-Chem sensitivities to satellite observation space?
        period_i       [int]        : kalman filter period
        config         [dict]       : dict of the config file
        use_water_obs  [bool]       : if True, use observations over water

    Returns
        output            [dict]       : Dictionary with:
                                        - obs_GC : GEOS-Chem and satellite data
                                        - satellite gas
                                        - GEOS-Chem gas
                                        - satellite lat, lon
                                        - satellite lat index, lon index
                                          If build_jacobian=True, also include:
                                            - K      : Jacobian matrix
    """

    # Read satellite data
    satellite, sat_ind = read_and_filter_satellite(
        filename, satellite_product, gc_startdate, gc_enddate,
        xlim, ylim, use_water_obs)

    # Number of satellite observations
    n_obs = len(sat_ind[0])
    if n_obs == 0:
        print(f"No satellite observations found in {filename}. Skipping.")
        return None
    print("Found", n_obs, "satellite observations.")

    # Define time threshold (hour 00 after the inversion period)
    date_after_inversion = str(gc_enddate + np.timedelta64(1, "D"))[:10].replace(
        "-", ""
    )
    time_threshold = f"{date_after_inversion}_00"

    # map satellite obs into gridcells and average the observations
    # into each gridcell. Only returns gridcells containing observations
    if config["UseGCHP"]:
        if config['STRETCH_GRID']:
            sf_formatted = f"{config['STRETCH_FACTOR']:.2f}".replace(".", "d")
            target_geohash = pgh.encode(config['TARGET_LAT'], config['TARGET_LON'])
            gridspec_path = f"c{config['CS_RES']}_s{sf_formatted}_t{target_geohash}_gridspec.nc"
        else:
            gridspec_path = f"c{config['CS_RES']}_gridspec.nc"
        GC_shape = (6, config['CS_RES'], config['CS_RES'])
        CSgridDir = f"{os.path.expandvars(config['OutputPath']) }/{config['RunName']}/CS_grids"

        obs_mapped_to_gc = average_tropomi_observations_to_CSgrid(
            TROPOMI, filename, sat_ind, time_threshold, CSgridDir, gridspec_path, GC_shape
        )
    else:
        # get the lat/lons of gc gridcells
        gc_lat_lon = get_gc_lat_lon(gc_cache, gc_startdate)
        obs_mapped_to_gc = average_satellite_observations(
            satellite, species, gc_lat_lon, sat_ind, time_threshold
        )
        GC_shape = (len(gc_lat_lon['lat']), len(gc_lat_lon['lon']))
    n_gridcells = len(obs_mapped_to_gc)

    if build_jacobian:
        # Initialize Jacobian K
        jacobian_K = np.empty([n_gridcells, n_elements], dtype=np.float32)
        jacobian_K.fill(np.nan)

        pertf = os.path.expandvars(
            f'{config["OutputPath"]}/{config["RunName"]}/'
            f"archive_perturbation_sfs/pert_sf_{period_i}.npz"
        )

        emis_perturbations_dict = np.load(pertf, mmap_mode='r')
        emis_perturbations = emis_perturbations_dict["effective_pert_sf"]

        # Calculate sensitivities and save in K matrix
        # determine which elements are for emis,
        # BCs, and OH
        oh_indices = []
        bc_indices = []
        emis_indices = []

        for e in range(n_elements):
            i_elem = e + 1
            # booleans for whether this element is a
            # BC element or OH element
            is_OH_element = check_is_OH_element(
                i_elem, n_elements, config["OptimizeOH"], config["isRegional"]
            )

            is_BC_element = check_is_BC_element(
                i_elem,
                n_elements,
                config["OptimizeOH"],
                config["OptimizeBCs"],
                is_OH_element,
                config["isRegional"],
            )

            if is_OH_element:
                oh_indices.append(e)
            elif is_BC_element:
                bc_indices.append(e)
            else:
                emis_indices.append(e)

    # Initialize array with n_gridcells rows and 5 columns. Columns are
    # satellite species, GEOSChem species, longitude, latitude, observation counts
    obs_GC = np.empty([n_gridcells, 5], dtype=np.float32)
    obs_GC.fill(np.nan)

    if config['UseGCHP']:
        GC_index = np.ravel_multi_index((obs_mapped_to_gc["nfi"],
                                         obs_mapped_to_gc["Ydimi"],
                                         obs_mapped_to_gc["Xdimi"]), GC_shape)
    else:
        GC_index = np.ravel_multi_index((obs_mapped_to_gc["jGC"], #lat
                                         obs_mapped_to_gc["iGC"] # lon
                                         ), GC_shape)

    all_strdate = [gridcell["time"] for gridcell in obs_mapped_to_gc]
    all_strdate = list(set(all_strdate))

    # Read GEOS-Chem data for simulated truth in OSSE simulation
    if config["EnableOSSE"]:
        osse_gc_cache = "./data_geoschem_osse"

        # check if the osse_gc_cache exists
        assert os.path.exists(osse_gc_cache), (
            f"OSSE GEOS-Chem cache directory {osse_gc_cache} does not exist. "
            "Please run the OSSE simulation first."
        )

    for strdate in all_strdate:
        gridcell_dict = obs_mapped_to_gc[obs_mapped_to_gc["time"] == strdate]
        sel_idx = np.where(obs_mapped_to_gc["time"] == strdate)[0]
        if build_jacobian:
            virtual_satellite_pert, virtual_satellite_base, virtual_satellite = get_virtual_satellite(
                strdate, gc_cache, gridcell_dict, n_elements, config, build_jacobian
            )
        else:
            virtual_satellite = get_virtual_satellite(
                strdate, gc_cache, gridcell_dict, n_elements, config, build_jacobian
            )
        if config["EnableOSSE"]:
            synthetic_virtual_satellite = get_virtual_satellite(
                strdate, osse_gc_cache, gridcell_dict, n_elements, config, False
            ) * 1e9  # convert to ppb

        # Save actual and virtual satellite data
        if config["EnableOSSE"]:
            # Synthetic observations if using OSSE, add random noise later
            obs_GC[sel_idx, 0] = synthetic_virtual_satellite
        else:
            # Actual satellite species column observation
            obs_GC[sel_idx, 0] = gridcell_dict[species]

        obs_GC[sel_idx, 1] = virtual_satellite * 1e9  # Virtual satellite column observation and convert to ppb
        obs_GC[sel_idx, 2] = gridcell_dict["lon_sat"]  # satellite longitude
        obs_GC[sel_idx, 3] = gridcell_dict["lat_sat"]  # satellite latitude
        obs_GC[sel_idx, 4] = gridcell_dict["observation_count"]  # observation counts

        if build_jacobian:
            pert_jacobian_xspecies = virtual_satellite_pert # (n_superobs, n_element)
            emis_base_xspecies = virtual_satellite_base # emis_base and BC_base is "RunName_0001" and "SpeciesConcVV_species"
            oh_base_xspecies = virtual_satellite # OH base is "RunName_0000"

            # get perturbations and calculate sensitivities
            perturbations = np.ones((len(gridcell_dict), n_elements), dtype=np.float32)

            # fill pert base array with values
            # array contains 1 entry for each state vector element
            # fill array with nans
            base_xspecies = np.full((len(gridcell_dict), n_elements), np.nan, dtype=np.float32)
            # fill emission elements with the base value
            base_xspecies[:,emis_indices] = np.repeat(emis_base_xspecies,
                                                  np.asarray(emis_indices).size, axis=1)

            # emissions perturbations
            perturbations[:,emis_indices] = np.repeat(emis_perturbations[None,:],
                                                      len(gridcell_dict), axis=0)

            # OH perturbations
            if config["OptimizeOH"]:
                # fill OH elements with the OH base value
                base_xspecies[:,oh_indices] = np.repeat(oh_base_xspecies[:,None],
                                                    np.asarray(oh_indices).size, axis=1)
                # update perturbations array to include OH perturbations
                perturbations[:,oh_indices] = float(config["PerturbValueOH"]) - 1.0

            # BC perturbations
            if config["OptimizeBCs"]:
                # fill BC elements with the base value, which is same as emis value
                base_xspecies[:,bc_indices] = np.repeat(emis_base_xspecies,
                                                    np.asarray(bc_indices).size, axis=1)

                # compute BC perturbation for jacobian construction
                perturbations[:,bc_indices] = config["PerturbValueBCs"]

            # calculate sensitivities
            jacobian_K[sel_idx,:] = ((pert_jacobian_xspecies - base_xspecies) / perturbations).astype(np.float32)

    # add random noise to synthetic observations if using OSSE
    if config["EnableOSSE"]:
        noise = np.random.normal(
            loc=0.0,
            scale=float(config["ObsErrorOSSE"]),
            size=obs_GC[:,0].shape,
        )
        obs_GC[:,0] += noise

    # Output
    output = {}

    # Always return the coincident satellite and GEOS-Chem data
    output["obs_GC"] = obs_GC
    output["GC_index"] = GC_index

    # Optionally return the Jacobian
    if build_jacobian:
        output["K"] = jacobian_K

    return output


def apply_satellite_operator(
    filename,
    species,
    satellite_product,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    period_i,
    config,
    use_water_obs=False,
):
    """
    Apply the satellite operator to map GEOS-Chem species data to satellite observation space.

    Arguments
        filename           [str]        : Satellite netcdf data file to read
        species            [str]        : The species (CH4 or CO2) to use
        satellite_product  [str]        : "BlendedTROPOMI", "TROPOMI", or "Other", specifying the data used in the inversion.
        n_elements         [int]        : Number of state vector elements
        gc_startdate       [datetime64] : First day of inversion period, for GEOS-Chem and satellite
        gc_enddate         [datetime64] : Last day of inversion period, for GEOS-Chem and satellite
        xlim               [float]      : Longitude bounds for simulation domain
        ylim               [float]      : Latitude bounds for simulation domain
        gc_cache           [str]        : Path to GEOS-Chem output data
        period_i       [int]            : kalman filter period
        config         [dict]          : dict of the config file
        use_water_obs  [bool]          : if True, use observations over water

    Returns
        output             [dict]       : Dictionary with one or two fields:
                                                    - obs_GC : GEOS-Chem and satellite species data
                                                    - satellite species
                                                    - GEOS-Chem species
                                                    - satellite lat, lon
                                                    - satellite lat index, lon index
                                                      If build_jacobian=True, also include:
                                                        - K      : Jacobian matrix
    """

    # Read satellite data
    satellite, sat_ind = read_and_filter_satellite(
        filename, satellite_product, gc_startdate, gc_enddate,
        xlim, ylim, use_water_obs)

    # Number of satellite observations
    n_obs = len(sat_ind[0])
    if n_obs == 0:
        return None
    # print("Found", n_obs, "satellite observations.")

    # Initialize a list to store the dates we want to look at
    all_strdate = []

    # Define time threshold (hour 00 after the inversion period)
    date_after_inversion = str(gc_enddate + np.timedelta64(1, "D"))[:10].replace(
        "-", ""
    )
    time_threshold = f"{date_after_inversion}_00"

    # For each satellite observation
    for k in range(n_obs):
        # Get the date and hour
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index
        time = pd.to_datetime(str(satellite["time"][iSat,jSat]))
        strdate = get_strdate(time, time_threshold)
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_gc = read_all_geoschem(
        all_strdate, gc_cache, config
    )

    # Initialize array with n_obs rows and 6 columns. Columns are satellite
    # mixing ratio, GEOSChem mixing ratio, longitude, latitude, II, JJ
    obs_GC = np.zeros([n_obs, 6], dtype=np.float32)
    obs_GC.fill(np.nan)

    if config['UseGCHP']:
        if config['STRETCH_GRID']:
            sf_formatted = f"{config['STRETCH_FACTOR']:.2f}".replace(".", "d")
            target_geohash = pgh.encode(config['TARGET_LAT'], config['TARGET_LON'])
            gridspec_path = f"c{config['CS_RES']}_s{sf_formatted}_t{target_geohash}_gridspec.nc"
        else:
            gridspec_path = f"c{config['CS_RES']}_gridspec.nc"
        GC_shape = (6, config['CS_RES'], config['CS_RES'])
        CSgridDir = f"{os.path.expandvars(config['OutputPath']) }/{config['RunName']}/CS_grids"

        overlap_area_all = get_overlap_area_CSgrid(TROPOMI, filename, sat_ind, CSgridDir,
                            gridspec_path, GC_shape) # (n_dst, n_valid_obs)

    # For each satellite observation:
    for k in range(n_obs):

        # Get GEOS-Chem data for the date of the observation:
        iSat = sat_ind[0][k]
        jSat = sat_ind[1][k]
        longitude_bounds = satellite["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = satellite["latitude_bounds"][iSat, jSat, :]

        p_sat = satellite["pressures"][iSat, jSat, :]
        dry_air_subcolumns = satellite["dry_air_subcolumns"][iSat, jSat, :]  # mol m-2
        apriori = satellite["profile_apriori"][iSat, jSat, :]  # mol m-2
        avkern = satellite["column_AK"][iSat, jSat, :]
        time = pd.to_datetime(str(satellite["time"][iSat,jSat]))
        strdate = get_strdate(time, time_threshold)
        GEOSCHEM = all_date_gc[strdate]

        if config['UseGCHP']:
            overlap_area_csr = overlap_area_all[:, k]
            gc_coords = overlap_area_csr.nonzero()[0]       # row indices of non-zero entries
            overlap_area = overlap_area_csr.data
        else:
            # Polygon representing satellite pixel
            polygon_satellite = Polygon(np.column_stack((longitude_bounds, latitude_bounds)))

            dlon = np.median(np.diff(GEOSCHEM["lon"]))  # GEOS-Chem lon resolution
            dlat = np.median(np.diff(GEOSCHEM["lat"]))  # GEOS-Chem lon resolution

            # Find GEOS-Chem lats & lons closest to the corners of the satellite pixel
            corners_lon_index = []
            corners_lat_index = []
            for l in range(4):
                iGC = nearest_loc(
                    longitude_bounds[l], GEOSCHEM["lon"], tolerance=max(dlon, 0.5)
                )
                jGC = nearest_loc(
                    latitude_bounds[l], GEOSCHEM["lat"], tolerance=max(dlat, 0.5)
                )
                corners_lon_index.append(iGC)
                corners_lat_index.append(jGC)
            # If the tolerance in nearest_loc() is not satisfied, skip the observation
            if np.nan in corners_lon_index + corners_lat_index:
                continue
            # Get lat/lon indexes and coordinates of GEOS-Chem grid cells closest to the satellite corners
            ij_GC = [(x, y) for x in set(corners_lon_index) for y in set(corners_lat_index)]
            gc_coords = [(GEOSCHEM["lon"][i], GEOSCHEM["lat"][j]) for i, j in ij_GC]

            # Compute the overlapping area between the satellite pixel and GEOS-Chem grid cells it touches
            overlap_area = np.zeros(len(gc_coords), dtype=np.float32)
            # For each GEOS-Chem grid cell that touches the satellite pixel:
            for gridcellIndex in range(len(gc_coords)):
                # Define polygon representing the GEOS-Chem grid cell
                coords = gc_coords[gridcellIndex]
                geoschem_corners_lon = [
                    coords[0] - dlon / 2,
                    coords[0] + dlon / 2,
                    coords[0] + dlon / 2,
                    coords[0] - dlon / 2,
                ]
                geoschem_corners_lat = [
                    coords[1] - dlat / 2,
                    coords[1] - dlat / 2,
                    coords[1] + dlat / 2,
                    coords[1] + dlat / 2,
                ]
                # If this is a global 2.0 x 2.5 grid, extend the eastern-most grid cells to 180 degrees
                if (dlon == 2.5) and (coords[0] == 177.5):
                    for i in [1, 2]:
                        geoschem_corners_lon[i] += dlon / 2

                polygon_geoschem = Polygon(
                    np.column_stack((geoschem_corners_lon, geoschem_corners_lat))
                )

                # Calculate overlapping area as the intersection of the two polygons
                if polygon_geoschem.intersects(polygon_satellite):
                    overlap_area[gridcellIndex] = polygon_satellite.intersection(
                        polygon_geoschem
                    ).area

        # If there is no overlap between GEOS-Chem and satellite, skip to next observation:
        if np.sum(overlap_area) == 0:
            continue

        # =======================================================
        #       Map GEOS-Chem to satellite observation space
        # =======================================================

        # Otherwise, initialize satellite virtual mixing ratios and virtual
        #  sensitivity as zero
        area_weighted_virtual_satellite = 0  # virtual satellite mixing ratio
        area_weighted_virtual_satellite_sensitivity = \
            np.zeros(n_elements, dtype=np.float32) if build_jacobian else None  # virtual satellite sensitivity

        # For each GEOS-Chem grid cell that touches the satellite pixel:
        for gridcellIndex in range(len(gc_coords)):
            if overlap_area[gridcellIndex] == 0:
                continue
            if config['UseGCHP']:
                # Get GEOS-Chem lat/lon indices for the cell
                f, j, x = np.unravel_index(gc_coords[gridcellIndex], GC_shape)

                # Get GEOS-Chem pressure edges for the cell
                p_gc = GEOSCHEM["PEDGE"][f, j, x, :]

                # Get GEOS-Chem mixing ratios for the cell
                gc_species = GEOSCHEM[species][f, j, x, :]
            else:
                # Get GEOS-Chem lat/lon indices for the cell
                iGC, jGC = ij_GC[gridcellIndex]

                # Get GEOS-Chem pressure edges for the cell
                p_gc = GEOSCHEM["PEDGE"][iGC, jGC, :]

                # Get GEOS-Chem mixing ratios for the cell
                gc_species = GEOSCHEM[species][iGC, jGC, :]

            # Get merged GEOS-Chem/satellite pressure grid for the cell
            merged = merge_pressure_grids(p_sat, p_gc)

            # Remap GEOS-Chem mixing ratios to satellite pressure levels
            sat_species = remap(
                gc_species,
                merged["data_type"],
                merged["p_merge"],
                merged["edge_index"],
                merged["first_gc_edge"],
            )  # ppb

            # Convert volumetric mixing ratio to mol m-2
            sat_species_molm2 = sat_species * 1/mixing_ratio_conv_factor(species) * dry_air_subcolumns  # mol m-2

            # Derive the column-averaged mixing ratio that satellite would
            # see over this ground cell
            virtual_satellite_gridcellIndex = apply_averaging_kernel(
                apriori, avkern, sat_species_molm2, dry_air_subcolumns, species
            ) # Volumetric mixing ratio

            # Weight by overlapping area (to be divided out later) and add to sum
            area_weighted_virtual_satellite += (
                overlap_area[gridcellIndex] * virtual_satellite_gridcellIndex
            )  # ppb m2

            # If building Jacobian matrix from GEOS-Chem perturbation simulation sensitivity data:
            if build_jacobian:
                if config['UseGCHP']:
                    sensi_lonlat = GEOSCHEM["jacobian_ch4"][f,j,x, :, :]
                else:
                    # Get GEOS-Chem perturbation sensitivities at this lat/lon, for all vertical levels and state vector elements
                    sensi_lonlat = GEOSCHEM["jacobian_ch4"][iGC, jGC, :, :]

                # Map the sensitivities to the satellite pressure levels
                sat_deltaCH4 = remap_sensitivities(
                    sensi_lonlat,
                    merged["data_type"],
                    merged["p_merge"],
                    merged["edge_index"],
                    merged["first_gc_edge"],
                )  # mixing ratio, unitless

                # Derive the change in column-averaged species that the satellite would see over this ground cell
                # in shape of (n_elements,)
                satellite_sensitivity_gridcellIndex = np.sum(avkern[:, None] * sat_deltaCH4 * \
                    dry_air_subcolumns[:, None], axis=0) / np.sum(dry_air_subcolumns)  # mixing ratio, unitless

                # Weight by overlapping area (to be divided out later) and add to sum
                area_weighted_virtual_satellite_sensitivity += (
                    overlap_area[gridcellIndex] * satellite_sensitivity_gridcellIndex
                )  # m2

        # Compute virtual satellite observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_satellite = area_weighted_virtual_satellite / sum(overlap_area)

        # Save actual and virtual satellite data
        obs_GC[k, 0] = satellite[species][
            iSat, jSat
        ]  # Actual satellite mixing ratio column observation
        obs_GC[k, 1] = virtual_satellite  # Virtual satellite mixing ratio column observation
        obs_GC[k, 2] = satellite["longitude"][iSat, jSat]  # satellite longitude
        obs_GC[k, 3] = satellite["latitude"][iSat, jSat]  # satellite latitude
        obs_GC[k, 4] = iSat  # satellite index of longitude
        obs_GC[k, 5] = jSat  # satellite index of latitude

    # Output
    output = {}

    # Always return the coincident satellite and GEOS-Chem data
    output["obs_GC"] = obs_GC

    return output


def average_satellite_observations(
        satellite, species, gc_lat_lon, sat_ind, time_threshold
    ):
    """
    Map satellite observations into appropriate gc gridcells. Then average all
    observations within a gridcell for processing. Use area weighting if
    observation overlaps multiple gridcells.

    Arguments
        satellite      [dict]   : Dict of satellite data
        species        [str]    : Name of species analyzed (CO2 or CH4)
        gc_lat_lon     [list]   : list of dictionaries containing  gc gridcell info
        sat_ind        [int]    : index list of satellite data that passes filters

    Returns
        numpy.ndarray
            Structured array of grid-cell-averaged values with fields:

            - iGC, jGC : int
                GEOS-Chem longitude and latitude indices
            - lat, lon : float
                Grid cell center coordinates
            - lat_sat, lon_sat : float
                Weighted-average satellite footprint center
            - species : float
                Weighted-average species column
            - time : str
                Averaged observation time (string from `get_strdate`)
            - p_sat : float[n_lev]
                Weighted-average vertical pressure profile
            - dry_air_subcolumns : float[n_lev]
                Weighted-average dry-air subcolumn profile
            - apriori : float[n_lev]
                Weighted-average a priori species profile
            - avkern : float[n_lev]
                Weighted-average averaging kernel
            - observation_count : float
                Effective number of contributing observations (fractional if
                split across multiple grid cells)

    """
    n_obs = len(sat_ind[0])
    # print("Found", n_obs, "satellite observations.")
    gc_lats = gc_lat_lon["lat"]
    gc_lons = gc_lat_lon["lon"]
    dlon = np.median(np.diff(gc_lat_lon["lon"]))  # GEOS-Chem lon resolution
    dlat = np.median(np.diff(gc_lat_lon["lat"]))  # GEOS-Chem lon resolution
    gridcell_dicts = get_gridcell_list(gc_lons, gc_lats, species)

    for k in range(n_obs):
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index

        # Find GEOS-Chem lats & lons closest to the corners of the satellite pixel
        longitude_bounds = satellite["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = satellite["latitude_bounds"][iSat, jSat, :]
        corners_lon_index = []
        corners_lat_index = []

        for l in range(4):
            iGC = nearest_loc(longitude_bounds[l], gc_lons, tolerance=max(dlon, 0.5))
            jGC = nearest_loc(latitude_bounds[l], gc_lats, tolerance=max(dlat, 0.5))
            corners_lon_index.append(iGC)
            corners_lat_index.append(jGC)

        # If the tolerance in nearest_loc() is not satisfied, skip the observation
        if np.nan in corners_lon_index + corners_lat_index:
            continue

        # Get lat/lon indexes and coordinates of GEOS-Chem grid cells closest to the satellite corners
        ij_GC = [(x, y) for x in set(corners_lon_index) for y in set(corners_lat_index)]
        gc_coords = [(gc_lons[i], gc_lats[j]) for i, j in ij_GC]

        # Compute the overlapping area between the satellite pixel and GEOS-Chem grid cells it touches
        overlap_area = np.zeros(len(gc_coords))

        # Polygon representing satellite pixel
        polygon_satellite = Polygon(np.column_stack((longitude_bounds, latitude_bounds)))
        for gridcellIndex in range(len(gc_coords)):
            # Define polygon representing the GEOS-Chem grid cell
            coords = gc_coords[gridcellIndex]
            geoschem_corners_lon = [
                coords[0] - dlon / 2,
                coords[0] + dlon / 2,
                coords[0] + dlon / 2,
                coords[0] - dlon / 2,
            ]
            geoschem_corners_lat = [
                coords[1] - dlat / 2,
                coords[1] - dlat / 2,
                coords[1] + dlat / 2,
                coords[1] + dlat / 2,
            ]
            polygon_geoschem = Polygon(
                np.column_stack((geoschem_corners_lon, geoschem_corners_lat))
            )

            # Calculate overlapping area as the intersection of the two polygons
            if polygon_geoschem.intersects(polygon_satellite):
                overlap_area[gridcellIndex] = polygon_satellite.intersection(
                    polygon_geoschem
                ).area
        # If there is no overlap between GEOS-Chem and satellite, skip to next observation:
        total_overlap_area = sum(overlap_area)
        if total_overlap_area == 0:
            continue

        # iterate through any gridcells with observation overlap
        # weight each observation if observation extent overlaps with multiple
        # gridcells
        for index, overlap in enumerate(overlap_area):
            if not overlap == 0:
                # get the matching dictionary for the gridcell with the overlap
                gridcell_dict = gridcell_dicts[ij_GC[index][0]][ij_GC[index][1]]
                gridcell_dict["lat_sat"].append(satellite["latitude"][iSat, jSat])
                gridcell_dict["lon_sat"].append(satellite["longitude"][iSat, jSat])
                gridcell_dict["overlap_area"].append(overlap)
                gridcell_dict["p_sat"].append(satellite["pressures"][iSat, jSat, :])
                gridcell_dict["dry_air_subcolumns"].append(
                    satellite["dry_air_subcolumns"][iSat, jSat, :]
                )
                gridcell_dict["apriori"].append(
                    satellite["profile_apriori"][iSat, jSat, :]
                )
                gridcell_dict["avkern"].append(satellite["column_AK"][iSat, jSat, :])
                gridcell_dict[
                    "time"
                ].append(  # convert times to epoch time to make taking the mean easier
                    int(pd.to_datetime(str(satellite["time"][iSat,jSat])).strftime("%s"))
                )
                gridcell_dict[species].append(
                    satellite[species][iSat, jSat]
                )  # Actual satellite mixing ratio column observation
                # record weights for averaging later
                gridcell_dict["observation_weights"].append(
                    overlap / total_overlap_area
                )
                # increment the observation count based on overlap area
                gridcell_dict["observation_count"] += overlap / total_overlap_area

    # filter out gridcells without any observations
    gridcell_dicts = [
        item for item in gridcell_dicts.flatten() if item["observation_count"] > 0
    ]
    # weighted average observation values for each gridcell
    for gridcell_dict in gridcell_dicts:
        gridcell_dict["lat_sat"] = np.average(
            gridcell_dict["lat_sat"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["lon_sat"] = np.average(
            gridcell_dict["lon_sat"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["overlap_area"] = np.average(
            gridcell_dict["overlap_area"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict[species] = np.average(
            gridcell_dict[species],
            weights=gridcell_dict["observation_weights"],
        )
        # take mean of epoch times and then convert gc filename time string
        time = pd.to_datetime(
            datetime.datetime.fromtimestamp(int(np.mean(gridcell_dict["time"])))
        )
        gridcell_dict["time"] = get_strdate(time, time_threshold)
        # for multi-dimensional arrays, we only take the average across the 0 axis
        gridcell_dict["p_sat"] = np.average(
            gridcell_dict["p_sat"],
            axis=0,
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["dry_air_subcolumns"] = np.average(
            gridcell_dict["dry_air_subcolumns"],
            axis=0,
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["apriori"] = np.average(
            gridcell_dict["apriori"],
            axis=0,
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["avkern"] = np.average(
            gridcell_dict["avkern"],
            axis=0,
            weights=gridcell_dict["observation_weights"],
        )

    if not gridcell_dicts:
        # nothing to return
        return np.zeros(0, dtype=[
            ("iGC","i4"), ("jGC","i4"),
            ("lat_sat","f4"), ("lon_sat","f4"),
            (species,"f4"), ("time","U13"),
            ("p_sat","f4",(0,)),
            ("dry_air_subcolumns","f4",(0,)),
            ("apriori","f4",(0,)),
            ("avkern","f4",(0,)),
            ("observation_count","f4"),
            ("lat","f4"), ("lon","f4"),
        ])

    # infer vertical sizes from the first item
    n_lev_p       = len(gridcell_dicts[0]["p_sat"])
    n_lev_dryair  = len(gridcell_dicts[0]["dry_air_subcolumns"])
    n_lev_apriori = len(gridcell_dicts[0]["apriori"])
    n_lev_avkern  = len(gridcell_dicts[0]["avkern"])

    dtype_latlon = [
        ("iGC","i4"), ("jGC","i4"),
        ("lat_sat","f4"), ("lon_sat","f4"),
        (species,"f4"), ("time","U13"),
        ("p_sat","f4",(n_lev_p,)),
        ("dry_air_subcolumns","f4",(n_lev_dryair,)),
        ("apriori","f4",(n_lev_apriori,)),
        ("avkern","f4",(n_lev_avkern,)),
        ("observation_count","f4"),
        ("lat","f4"), ("lon","f4"),
    ]

    arr = np.zeros(len(gridcell_dicts), dtype=dtype_latlon)

    for idx, cell in enumerate(gridcell_dicts):
        arr["iGC"][idx]   = cell["iGC"]
        arr["jGC"][idx]   = cell["jGC"]
        arr["lat_sat"][idx] = np.float32(cell["lat_sat"])
        arr["lon_sat"][idx] = np.float32(cell["lon_sat"])
        arr[species][idx] = np.float32(cell[species])
        arr["time"][idx]    = cell["time"]  # already a short string from get_strdate
        arr["p_sat"][idx]   = np.asarray(cell["p_sat"], dtype=np.float32)
        arr["dry_air_subcolumns"][idx] = np.asarray(cell["dry_air_subcolumns"], dtype=np.float32)
        arr["apriori"][idx] = np.asarray(cell["apriori"], dtype=np.float32)
        arr["avkern"][idx]  = np.asarray(cell["avkern"], dtype=np.float32)
        arr["observation_count"][idx] = np.float32(cell["observation_count"])
        arr["lat"][idx] = np.float32(cell["lat"])
        arr["lon"][idx] = np.float32(cell["lon"])

    return arr

def average_satellite_observations_to_CSgrid(
        satellite, species, filename, sat_ind, time_threshold,
        CSgridDir, gridspec_path, GC_shape):
    """
    Map satellite observations into appropriate GCHP gridcells. Then average all
    observations within a gridcell for processing.
    Use area weighting (conservative regridding weights x destination area)
    if observation overlaps multiple gridcells.

    Parameters
    ----------
    satellite : dict
        Dictionary of satellite satellite data.
    filename : str
        Path to the satellite dataset file.
    sat_ind : int
        Indices of valid satellite observations (after filtering).
    time_threshold : str or float
        Threshold used to bin/average times into representative strings.
    CSgridDir : str
        Path to cubed-sphere grid definitions.
    gridspec_path : str
        Path to the GCHP gridspec file.
    GC_shape : tuple
        Shape of the GEOS-Chem cubed-sphere grid.

    Returns
    -------
    output_dicts : np.ndarray of dict-like records
        Array of per-gridcell averaged observation data with fields:
            - nfi, Ydimi, Xdimi : cubed-sphere indices
            - lat_sat, lon_sat  : averaged satellite coordinates
            - species           : averaged species column
            - time              : averaged time string
            - p_sat             : averaged pressure grid
            - dry_air_subcolumns: averaged dry air columns
            - apriori           : averaged prior species profiles
            - avkern            : averaged averaging kernels
            - observation_count : number of obs contributing to the average
    """
    Sat_shape = satellite["longitude"].shape
    # flatten indices to valid obs
    sat_ind_flat = np.ravel_multi_index(sat_ind, Sat_shape)
    sat_mask = np.zeros(np.product(Sat_shape), dtype=bool)
    sat_mask[sat_ind_flat] = True

    overlap_area = get_overlap_area_CSgrid(satellite, filename, sat_ind, CSgridDir,
                            gridspec_path, GC_shape) # (n_dst, n_valid_obs)

    # Observation-centric weights
    # observation weights w_ij = overlap_ij / sum_j(overlap_ij)
    # where w_ij is the observation weight for destination grid cell j from the observation cell i
    # overlap_ij is the overlap area between destination cell j and the observation cell i
    # total overlap area for each observation is the sum of the overlapping area
    # between observation cell i and all its overlapping destination grid cell j
    total_overlap_area = np.asarray(overlap_area.sum(axis=0)).ravel() # (n_valid_obs)
    inv_total_overlap_area = np.zeros_like(total_overlap_area, dtype=np.float32)
    valid_obs = total_overlap_area > 0
    inv_total_overlap_area[valid_obs] = 1.0 / total_overlap_area[valid_obs]

    # obs_weight: each obs contributes fully across overlapping cells
    # element-wise multiplication:
    obs_weights = overlap_area.multiply(inv_total_overlap_area).tocsr() # (n_dst × n_valid_obs)

    sum_obs_weights = np.asarray(obs_weights.sum(axis=1)).ravel() # (n_dst)
    GC_mask = sum_obs_weights > 0
    GC_indices = np.where(GC_mask)[0]  # only valid destination cells
    n_valid_GC = len(GC_indices)

    # ---------- Scalars ----------
    sat_lon_flat = satellite["longitude"].ravel()[sat_mask]
    sat_lat_flat = satellite["latitude"].ravel()[sat_mask]
    sat_species_flat = satellite[species].ravel()[sat_mask]
    sat_time_flat = pd.to_datetime(satellite["time"].ravel()[sat_mask])

    # For each destination grid cell j which overlap with at least one observation cell i,
    # weighted_obs_j = sum_i(obs_weights_ij * obs_i) / sum_i(obs_weights_ij)
    sat_lon_avg = (obs_weights[GC_indices, :] @ sat_lon_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    sat_lat_avg = (obs_weights[GC_indices, :] @ sat_lat_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    sat_species_avg = (obs_weights[GC_indices, :] @ sat_species_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    sat_time_avg = (obs_weights[GC_indices, :] @ sat_time_flat.astype(np.int64)) / sum_obs_weights[GC_indices]

    sat_time_avg = pd.to_datetime(sat_time_avg)  # convert back to datetime
    sat_time_str = np.array([get_strdate(t, time_threshold) for t in sat_time_avg])

    # --- Multi-dimensional arrays ---
    def weighted_avg_profiles(var_flat):
        # var_flat: (n_valid_obs, n_lev)
        # obs_weights: (n_dst, n_valid_obs)
        # sum_obs_weights: (n_dst,)
        # result: (n_valid_GC,)
        result = obs_weights[GC_indices, :] @ var_flat
        result /= sum_obs_weights[GC_indices, None]
        return result.astype(np.float32)

    p_sat_flat = species["pressures"].reshape(-1, species["pressures"].shape[-1])[sat_mask, :]
    dryair_flat = species["dry_air_subcolumns"].reshape(-1, species["dry_air_subcolumns"].shape[-1])[sat_mask, :]
    apriori_flat = species["profile_apriori"].reshape(-1, species["profile_apriori"].shape[-1])[sat_mask, :]
    avkern_flat = species["column_AK"].reshape(-1, species["column_AK"].shape[-1])[sat_mask, :]

    p_sat_avg = weighted_avg_profiles(p_sat_flat)
    dryair_avg = weighted_avg_profiles(dryair_flat)
    apriori_avg = weighted_avg_profiles(apriori_flat)
    avkern_avg = weighted_avg_profiles(avkern_flat)

    # --- Fill dictionaries ---
    f_idx, j_idx, x_idx = np.unravel_index(GC_indices, GC_shape)
    n_lev_p = p_sat_avg.shape[1]
    n_lev_dryair = dryair_avg.shape[1]
    n_lev_apriori = apriori_avg.shape[1]
    n_lev_avkern = avkern_avg.shape[1]

    dtype = [
        ("nfi", "i4"), ("Ydimi", "i4"), ("Xdimi", "i4"),
        ("lat_sat", "f4"), ("lon_sat", "f4"), (species, "f4"),
        ("time", "U13"), ("p_sat", "f4", (n_lev_p,)),
        ("dry_air_subcolumns", "f4", (n_lev_dryair,)),
        ("apriori", "f4", (n_lev_apriori,)), ("avkern", "f4", (n_lev_avkern,)),
        ("observation_count", "f4")
    ]

    output_dicts = np.zeros(n_valid_GC, dtype=dtype)
    output_dicts["nfi"] = f_idx
    output_dicts["Ydimi"] = j_idx
    output_dicts["Xdimi"] = x_idx
    output_dicts["lat_sat"] = sat_lat_avg
    output_dicts["lon_sat"] = sat_lon_avg
    output_dicts[species] = sat_species_avg
    output_dicts["time"] = sat_time_str
    output_dicts["p_sat"] = p_sat_avg
    output_dicts["dry_air_subcolumns"] = dryair_avg
    output_dicts["apriori"] = apriori_avg
    output_dicts["avkern"] = avkern_avg
    output_dicts["observation_count"] = sum_obs_weights[GC_indices]

    return output_dicts

def get_virtual_satellite(date, gc_cache, gridcell_dict, n_elements, config, build_jacobian=False):
    """
    Generate virtual satellite species observations from GEOS-Chem.

    Extracts species and pressure from GEOS-Chem, remaps to satellite layers,
    and applies averaging kernels. Optionally computes Jacobian using
    perturbation runs.

    Parameters
    ----------
    date : str
        Date of interest ("YYYYMMDD_HH").
    gc_cache : str
        Path to GEOS-Chem output files.
    gridcell_dict : dict
        Gridcell info with obs indices and satellite data.
    n_elements : int
        Number of state vector elements.
    config : dict
        Inversion configuration options.
    build_jacobian : bool, optional
        Whether to compute sensitivities (default False).

    Returns
    -------
    If build_jacobian=False:
        ndarray (N,) of virtual satellite columns.
    If build_jacobian=True:
        (perturbation columns, base columns, final columns).
    """

    UseGCHP = config['UseGCHP']

    # Assemble file paths to GEOS-Chem output collections for input data
    file_species = f"GEOSChem.SpeciesConc.{date}00z.nc4"
    file_pedge = f"GEOSChem.StateMetLevEdge.{date}00z.nc4"

    # Read lat, lon, species from the SpeciecConc collection
    filename = f"{gc_cache}/{file_species}"

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        with xr.open_dataset(filename, chunks='auto') as gc_data_all:
            if gc_data_all.sizes.get("time", 0) == 0:
                print(f"ERROR: {filename}: empty time dimension", flush=True)
            if UseGCHP:
                nfi   = xr.DataArray(gridcell_dict["nfi"],   dims="obs")
                Ydimi = xr.DataArray(gridcell_dict["Ydimi"], dims="obs")
                Xdimi = xr.DataArray(gridcell_dict["Xdimi"], dims="obs")

                gc_data = gc_data_all.isel(
                    time=0).squeeze().isel(
                    nf=nfi,
                    Ydim=Ydimi,
                    Xdim=Xdimi,
                    drop=True
                )
            else:
                jGC   = xr.DataArray(gridcell_dict["jGC"],   dims="obs")
                iGC   = xr.DataArray(gridcell_dict["iGC"],   dims="obs")
                gc_data = gc_data_all.isel(
                    time=0).squeeze().isel(
                    lat=jGC,
                    lon=iGC,
                    drop=True
                )
            species = gc_data["SpeciesConcVV_{config['Species']}"].transpose("obs","lev").values

    # Read PEDGE from the StateMetLevEdge collection
    filename = f"{gc_cache}/{file_pedge}"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        with xr.open_dataset(filename, chunks='auto') as gc_data_all:
            if gc_data_all.sizes.get("time", 0) == 0:
                print(f"ERROR: {filename}: empty time dimension", flush=True)
            if UseGCHP:
                gc_data = gc_data_all.isel(
                    time=0).squeeze().isel(
                    nf=nfi,
                    Ydim=Ydimi,
                    Xdim=Xdimi,
                    drop=True
                )
            else:
                gc_data = gc_data_all.isel(
                    time=0).squeeze().isel(
                    lat=jGC,
                    lon=iGC,
                    drop=True
                )
            lev_dim = "lev" if "lev" in gc_data["Met_PEDGE"].dims else "ilev"
            PEDGE = gc_data["Met_PEDGE"].transpose("obs", lev_dim).values

    n_superobs = len(gridcell_dict)
    virtual_satellite = np.empty([n_superobs, ], dtype=np.float32)
    virtual_satellite.fill(np.nan)

    p_sat = gridcell_dict["p_sat"]
    dry_air_subcolumns = gridcell_dict["dry_air_subcolumns"]  # mol m-2
    apriori = gridcell_dict["apriori"]  # mol m-2
    avkern = gridcell_dict["avkern"]

    # (N, S, G)  sums to 1 along G, where
    # N is the number of super observations
    # S is the number of satellite pressure edges
    # G is the number of GEOS-Chem pressure edges
    vertical_weights = remapping_weights(p_sat, PEDGE)
    sat_species = np.einsum("nsg,ng->ns", vertical_weights, species)         # (N, S)
    sat_species_molm2 = sat_species * dry_air_subcolumns                     # (N, S)
    virtual_satellite = (
        np.sum(apriori + avkern * (sat_species_molm2 - apriori), axis=1) / \
        np.sum(dry_air_subcolumns, axis=1)
    ).astype(np.float32)                      # (N,), unitless mixing ratio

    # If need to construct Jacobian, read sensitivity data from GEOS-Chem perturbation simulations
    if build_jacobian:
        emis_elements = n_elements
        if config['OptimizeOH']:
            emis_elements -= 2 if config['isRegional'] else 1
        if config['OptimizeBCs']:
            emis_elements -= 4
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
        for e in range(n_elements):
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
        virtual_satellite_pert = [
            get_virtual_satellite_pert(gc_date, k, gridcell_dict, config, v, n_elements, vertical_weights)
            for k, v in pert_simulations_dict.items()
        ]

        virtual_satellite_pert = np.concatenate(virtual_satellite_pert, axis=1)

        virtual_satellite_base = get_virtual_satellite_pert(
            gc_date, "0001", gridcell_dict, config, [0], n_elements, vertical_weights, baserun=True
        )

    if build_jacobian:
        return virtual_satellite_pert, virtual_satellite_base, virtual_satellite
    else:
        return virtual_satellite

def get_virtual_satellite_pert(gc_date, run_id, gridcell_dict, config, sv_elems, n_elements, vertical_weights, baserun=False):
    """
    Compute virtual satellite species columns from a single GEOS-Chem Jacobian run.

    Extracts species tracer(s) from the specified perturbation (or base) simulation,
    remaps them to satellite pressure layers using vertical weights, and applies
    averaging kernels to generate virtual satellite observations.

    Parameters
    ----------
    gc_date : pd.Datetime
        Date and time of the simulation output.
    run_id : str
        ID of the Jacobian GEOS-Chem run (e.g., "0001").
    gridcell_dict : dict
        Gridcell info containing observation indices and satellite data.
    config : dict
        Inversion configuration options.
    sv_elems : list
        State vector elements included in this simulation.
    n_elements : int
        Total number of state vector elements.
    vertical_weights : np.ndarray
        Weights to remap GEOS-Chem levels to satellite layers.
    baserun : bool, optional
        If True, only process the base species (default False).

    Returns
    -------
    virtual_satellite_all : np.ndarray, shape (N, n_tracers)
        Virtual satellite columns for all super-observations (N) and
        tracer elements in this run, in unitless mixing ratio.
    """
    prefix = os.path.expandvars(
        config["OutputPath"] + "/" + config["RunName"] + "/jacobian_runs"
    )
    j_dir = f"{prefix}/{config['RunName']}_{run_id}/OutputDir"
    file_stub = gc_date.strftime("GEOSChem.SpeciesConc.%Y%m%d_0000z.nc4")
    filepath = os.path.join(j_dir, file_stub)

    # Construct the list of species vars to request
    keepvars = [f"SpeciesConcVV_{config['Species']}_{i:04}" for i in sv_elems]
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
            keepvars = [f"SpeciesConcVV_{config['Species']}"]

    if baserun:
        keepvars = ["SpeciesConcVV_{config['Species']}"]

    # It would fail if open all variables with chunks with GCHP,
    # as ncontact is duplicate for GCHP output dimensions
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        with xr.open_dataset(filepath, decode_cf=False) as tmp:
            other_vars = [v for v in tmp.variables if "SpeciesConcVV_{config['Species']}" not in v]

    # Open only these variables
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        with xr.open_dataset(
            filepath,
            drop_variables=other_vars,
            chunks="auto"
        ) as dsmf_all:
            try:
                if dsmf_all.sizes.get("time", 0) == 0:
                    print(f"ERROR: {filepath}: empty time dimension", flush=True)
                if config['UseGCHP']:
                    nfi   = xr.DataArray(gridcell_dict["nfi"],   dims="obs")
                    Ydimi = xr.DataArray(gridcell_dict["Ydimi"], dims="obs")
                    Xdimi = xr.DataArray(gridcell_dict["Xdimi"], dims="obs")
                    dsmf = dsmf_all.isel(
                        time=gc_date.hour).squeeze().isel(
                        nf=nfi,
                        Ydim=Ydimi,
                        Xdim=Xdimi,
                        drop=True
                    )
                else:
                    jGC   = xr.DataArray(gridcell_dict["jGC"],   dims="obs")
                    iGC   = xr.DataArray(gridcell_dict["iGC"],   dims="obs")
                    dsmf = dsmf_all.isel(
                        time=gc_date.hour).squeeze().isel(
                        lat=jGC,
                        lon=iGC,
                        drop=True
                    )
            except Exception as e:
                print(f"Run id {run_id}. Failed at {gc_date} with error: {e}", flush=True)
                raise

            # ---- Batch read all species and standardize to (elem, obs, lev)
            da_list = []
            for v in keepvars:
                da = dsmf[v]
                # Ensure shape is (obs, lev)
                arr = da.transpose("obs","lev").values
                da_list.append(arr)
            species_all = np.stack(da_list, axis=0)   # (elem, obs, lev)

            dry_air_subcolumns = gridcell_dict["dry_air_subcolumns"]  # (N, S)
            apriori = gridcell_dict["apriori"]                        # (N, S)
            avkern = gridcell_dict["avkern"]                          # (N, S)
            denom = np.sum(dry_air_subcolumns, axis=1)                # (N,)

            # ---- Remap levels to satellite layers in batch:
            # vertical_weights:   (N, S, G)
            # species_all:        (E, N, G)
            # -> sat_species_all: (E, N, S)
            sat_species_all = np.einsum("nsg,eng->ens", vertical_weights, species_all)

            # Convert to column units and apply AKs, batched over E
            sat_species_molm2_all = sat_species_all * dry_air_subcolumns[None, :, :]  # (E, N, S)
            numer_all = np.sum(apriori[None, :, :] +
                               avkern[None, :, :] * (sat_species_molm2_all - apriori[None, :, :]),
                               axis=2)  # (E, N)
            virtual_satellite_all = (numer_all / denom[None, :]).T.astype(np.float32)  # (N, E)

    return virtual_satellite_all # unitless mixing ratio

def apply_averaging_kernel(
        apriori,
        avkern,
        sat_species_molm2,
        dry_air_subcolumns,
        species
):
    # Derive the column-averaged mixing ratio that the satellite would see
    # over this ground cell using eq. 46 from TROPOMI Methane ATBD,
    # Hasekamp et al. 2019
    virtual_satellite = (
        sum(apriori + avkern * (sat_species_molm2 - apriori))
        / sum(dry_air_subcolumns)
        * mixing_ratio_conv_factor(species)
    )  # volumetric mixing ratio
    return virtual_satellite

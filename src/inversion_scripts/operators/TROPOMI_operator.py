import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime
from shapely.geometry import Polygon
from src.inversion_scripts.utils import (
    filter_tropomi,
    filter_blended,
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
)

def apply_average_tropomi_operator(
    filename,
    BlendedTROPOMI,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    build_jacobian,
    period_i,
    config,
    use_water_obs=False,
):
    """
    Apply the averaging tropomi operator to map GEOS-Chem methane data to TROPOMI observation space.

    Arguments
        filename       [str]        : TROPOMI netcdf data file to read
        BlendedTROPOMI [bool]       : if True, use blended TROPOMI+GOSAT data
        n_elements     [int]        : Number of state vector elements
        gc_startdate   [datetime64] : First day of inversion period, for GEOS-Chem and TROPOMI
        gc_enddate     [datetime64] : Last day of inversion period, for GEOS-Chem and TROPOMI
        xlim           [float]      : Longitude bounds for simulation domain
        ylim           [float]      : Latitude bounds for simulation domain
        gc_cache       [str]        : Path to GEOS-Chem output data
        build_jacobian [log]        : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        period_i       [int]        : kalman filter period
        config         [dict]       : dict of the config file
        use_water_obs  [bool]       : if True, use observations over water

    Returns
        output         [dict]       : Dictionary with:
                                        - obs_GC : GEOS-Chem and TROPOMI methane data
                                        - TROPOMI methane
                                        - GEOS-Chem methane
                                        - TROPOMI lat, lon
                                        - TROPOMI lat index, lon index
                                          If build_jacobian=True, also include:
                                            - K      : Jacobian matrix
    """

    # Read TROPOMI data
    assert isinstance(BlendedTROPOMI, bool), "BlendedTROPOMI is not a bool"
    if BlendedTROPOMI:
        TROPOMI = read_blended(filename)
    else:
        TROPOMI = read_tropomi(filename)
    if TROPOMI == None:
        print(f"Skipping {filename} due to file processing issue.")
        return TROPOMI

    if BlendedTROPOMI:
        # Only going to consider blended data within lat/lon/time bounds and wihtout problematic coastal pixels
        sat_ind = filter_blended(
            TROPOMI, xlim, ylim, gc_startdate, gc_enddate, use_water_obs
        )
    else:
        # Only going to consider TROPOMI data within lat/lon/time bounds and with QA > 0.5
        sat_ind = filter_tropomi(
            TROPOMI, xlim, ylim, gc_startdate, gc_enddate, use_water_obs
        )

    # Number of TROPOMI observations
    n_obs = len(sat_ind[0])
    if n_obs == 0:
        print(f"No TROPOMI observations found in {filename}. Skipping.")
        return None
    print("Found", n_obs, "TROPOMI observations.")
    # get the lat/lons of gc gridcells
    gc_lat_lon = get_gc_lat_lon(gc_cache, gc_startdate)
    
    # Define time threshold (hour 00 after the inversion period)
    date_after_inversion = str(gc_enddate + np.timedelta64(1, "D"))[:10].replace(
        "-", ""
    )
    time_threshold = f"{date_after_inversion}_00"

    # map tropomi obs into gridcells and average the observations
    # into each gridcell. Only returns gridcells containing observations
    obs_mapped_to_gc = average_tropomi_observations(
        TROPOMI, gc_lat_lon, sat_ind, time_threshold
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

        emis_perturbations_dict = np.load(pertf)
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
    # TROPOMI CH4, GEOSChem CH4, longitude, latitude, observation counts
    obs_GC = np.empty([n_gridcells, 5], dtype=np.float32)
    obs_GC.fill(np.nan)
    
    all_strdate = [gridcell["time"] for gridcell in obs_mapped_to_gc]
    all_strdate = list(set(all_strdate))
    
    for strdate in all_strdate:
        gridcell_dict = obs_mapped_to_gc[obs_mapped_to_gc["time"] == strdate]
        sel_idx = np.where(obs_mapped_to_gc["time"] == strdate)[0]
        if build_jacobian:
            virtual_tropomi_pert, virtual_tropomi_base, virtual_tropomi = get_virtual_tropomi(
                strdate, gc_cache, gridcell_dict, n_elements, config, build_jacobian
            )
        else:
            virtual_tropomi = get_virtual_tropomi(
                strdate, gc_cache, gridcell_dict, n_elements, config, build_jacobian
            )
        # Save actual and virtual TROPOMI data
        obs_GC[sel_idx, 0] = gridcell_dict[
            "methane"
        ]  # Actual TROPOMI methane column observation
        obs_GC[sel_idx, 1] = virtual_tropomi * 1e9  # Virtual TROPOMI methane column observation and convert to ppb
        obs_GC[sel_idx, 2] = gridcell_dict["lon_sat"]  # TROPOMI longitude
        obs_GC[sel_idx, 3] = gridcell_dict["lat_sat"]  # TROPOMI latitude
        obs_GC[sel_idx, 4] = gridcell_dict["observation_count"]  # observation counts
        
        if build_jacobian:
            pert_jacobian_xch4 = virtual_tropomi_pert # (n_superobs, n_element)
            emis_base_xch4 = virtual_tropomi_base # emis_base and BC_base is "RunName_0001" and "SpeciesConcVV_CH4"
            oh_base_xch4 = virtual_tropomi # OH base is "RunName_0000"
            
            # get perturbations and calculate sensitivities
            perturbations = np.ones((len(gridcell_dict), n_elements), dtype=np.float32)

            # fill pert base array with values
            # array contains 1 entry for each state vector element
            # fill array with nans
            base_xch4 = np.full((len(gridcell_dict), n_elements), np.nan, dtype=np.float32)
            # fill emission elements with the base value
            base_xch4[:,emis_indices] = np.repeat(emis_base_xch4[:, None], 
                                                  np.asarray(emis_indices).size, axis=1)

            # emissions perturbations
            perturbations[:,emis_indices] = np.repeat(emis_perturbations[None,:],
                                                      len(gridcell_dict), axis=0)

            # OH perturbations
            if config["OptimizeOH"]:
                # fill OH elements with the OH base value
                base_xch4[:,oh_indices] = np.repeat(oh_base_xch4[:,None],
                                                    np.asarray(oh_indices).size, axis=1)
                # update perturbations array to include OH perturbations
                perturbations[:,oh_indices] = float(config["PerturbValueOH"]) - 1.0

            # BC perturbations
            if config["OptimizeBCs"]:
                # fill BC elements with the base value, which is same as emis value
                base_xch4[:,bc_indices] = np.repeat(emis_base_xch4[:,None],
                                                    np.asarray(bc_indices).size, axis=1)

                # compute BC perturbation for jacobian construction
                perturbations[:,bc_indices] = config["PerturbValueBCs"]

            # calculate sensitivities
            jacobian_K[sel_idx,:] = ((pert_jacobian_xch4 - base_xch4) / perturbations).astype(np.float32)
            
    # Output
    output = {}

    # Always return the coincident TROPOMI and GEOS-Chem data
    output["obs_GC"] = obs_GC

    # Optionally return the Jacobian
    if build_jacobian:
        output["K"] = jacobian_K

    return output


def apply_tropomi_operator(
    filename,
    BlendedTROPOMI,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    build_jacobian,
    period_i,
    config,
    use_water_obs=False,
):
    """
    Apply the tropomi operator to map GEOS-Chem methane data to TROPOMI observation space.

    Arguments
        filename       [str]        : TROPOMI netcdf data file to read
        BlendedTROPOMI [bool]       : if True, use blended TROPOMI+GOSAT data
        n_elements     [int]        : Number of state vector elements
        gc_startdate   [datetime64] : First day of inversion period, for GEOS-Chem and TROPOMI
        gc_enddate     [datetime64] : Last day of inversion period, for GEOS-Chem and TROPOMI
        xlim           [float]      : Longitude bounds for simulation domain
        ylim           [float]      : Latitude bounds for simulation domain
        gc_cache       [str]        : Path to GEOS-Chem output data
        build_jacobian [log]        : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        period_i       [int]        : kalman filter period
        config         [dict]       : dict of the config file
        use_water_obs  [bool]       : if True, use observations over water

    Returns
        output         [dict]       : Dictionary with one or two fields:
                                                        - obs_GC : GEOS-Chem and TROPOMI methane data
                                                    - TROPOMI methane
                                                    - GEOS-Chem methane
                                                    - TROPOMI lat, lon
                                                    - TROPOMI lat index, lon index
                                                      If build_jacobian=True, also include:
                                                        - K      : Jacobian matrix
    """

    # Read TROPOMI data
    assert isinstance(BlendedTROPOMI, bool), "BlendedTROPOMI is not a bool"
    if BlendedTROPOMI:
        TROPOMI = read_blended(filename)
    else:
        TROPOMI = read_tropomi(filename)
    if TROPOMI == None:
        print(f"Skipping {filename} due to file processing issue.")
        return TROPOMI

    if BlendedTROPOMI:
        # Only going to consider blended data within lat/lon/time bounds and wihtout problematic coastal pixels
        sat_ind = filter_blended(
            TROPOMI, xlim, ylim, gc_startdate, gc_enddate, use_water_obs
        )
    else:
        # Only going to consider TROPOMI data within lat/lon/time bounds and with QA > 0.5
        sat_ind = filter_tropomi(
            TROPOMI, xlim, ylim, gc_startdate, gc_enddate, use_water_obs
        )

    # Number of TROPOMI observations
    n_obs = len(sat_ind[0])
    if n_obs == 0:
        return None
    # print("Found", n_obs, "TROPOMI observations.")

    # If need to build Jacobian from GEOS-Chem perturbation simulation sensitivity data:
    if build_jacobian:
        # Initialize Jacobian K
        jacobian_K = np.zeros([n_obs, n_elements], dtype=np.float32)
        jacobian_K.fill(np.nan)

    # Initialize a list to store the dates we want to look at
    all_strdate = []
    
    # Define time threshold (hour 00 after the inversion period)
    date_after_inversion = str(gc_enddate + np.timedelta64(1, "D"))[:10].replace(
        "-", ""
    )
    time_threshold = f"{date_after_inversion}_00"

    # For each TROPOMI observation
    for k in range(n_obs):
        # Get the date and hour
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index
        time = pd.to_datetime(str(TROPOMI["time"][iSat, jSat]))
        strdate = get_strdate(time, time_threshold)
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_gc = read_all_geoschem(
        all_strdate, gc_cache, n_elements, config, build_jacobian
    )
    
    # Initialize array with n_obs rows and 6 columns. Columns are TROPOMI CH4, GEOSChem CH4, longitude, latitude, II, JJ
    obs_GC = np.zeros([n_obs, 6], dtype=np.float32)
    obs_GC.fill(np.nan)

    # For each TROPOMI observation:
    for k in range(n_obs):
        # Get GEOS-Chem data for the date of the observation:
        iSat = sat_ind[0][k]
        jSat = sat_ind[1][k]
        p_sat = TROPOMI["pressures"][iSat, jSat, :]
        dry_air_subcolumns = TROPOMI["dry_air_subcolumns"][iSat, jSat, :]  # mol m-2
        apriori = TROPOMI["methane_profile_apriori"][iSat, jSat, :]  # mol m-2
        avkern = TROPOMI["column_AK"][iSat, jSat, :]
        time = pd.to_datetime(str(TROPOMI["time"][iSat, jSat]))
        strdate = get_strdate(time, time_threshold)
        GEOSCHEM = all_date_gc[strdate]

        dlon = np.median(np.diff(GEOSCHEM["lon"]))  # GEOS-Chem lon resolution
        dlat = np.median(np.diff(GEOSCHEM["lat"]))  # GEOS-Chem lon resolution

        # Find GEOS-Chem lats & lons closest to the corners of the TROPOMI pixel
        longitude_bounds = TROPOMI["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = TROPOMI["latitude_bounds"][iSat, jSat, :]
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
        # Get lat/lon indexes and coordinates of GEOS-Chem grid cells closest to the TROPOMI corners
        ij_GC = [(x, y) for x in set(corners_lon_index) for y in set(corners_lat_index)]
        gc_coords = [(GEOSCHEM["lon"][i], GEOSCHEM["lat"][j]) for i, j in ij_GC]

        # Compute the overlapping area between the TROPOMI pixel and GEOS-Chem grid cells it touches
        overlap_area = np.zeros(len(gc_coords))
        # Polygon representing TROPOMI pixel
        polygon_tropomi = Polygon(np.column_stack((longitude_bounds, latitude_bounds)))
        # For each GEOS-Chem grid cell that touches the TROPOMI pixel:
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
            if (dlon == 2.5) & (coords[0] == 177.5):
                for i in [1, 2]:
                    geoschem_corners_lon[i] += dlon / 2
            polygon_geoschem = Polygon(
                np.column_stack((geoschem_corners_lon, geoschem_corners_lat))
            )
            
            # Calculate overlapping area as the intersection of the two polygons
            if polygon_geoschem.intersects(polygon_tropomi):
                overlap_area[gridcellIndex] = polygon_tropomi.intersection(
                    polygon_geoschem
                ).area

        # If there is no overlap between GEOS-Chem and TROPOMI, skip to next observation:
        if np.sum(overlap_area) == 0:
            continue

        # =======================================================
        #       Map GEOS-Chem to TROPOMI observation space
        # =======================================================

        # Otherwise, initialize tropomi virtual xch4 and virtual sensitivity as zero
        area_weighted_virtual_tropomi = 0  # virtual tropomi xch4
        area_weighted_virtual_tropomi_sensitivity = \
            np.zeros(n_elements, dtype=np.float32) if build_jacobian else None  # virtual tropomi sensitivity

        # For each GEOS-Chem grid cell that touches the TROPOMI pixel:
        for gridcellIndex in range(len(gc_coords)):
            if overlap_area[gridcellIndex] == 0:
                continue
            # Get GEOS-Chem lat/lon indices for the cell
            iGC, jGC = ij_GC[gridcellIndex]

            # Get GEOS-Chem pressure edges for the cell
            p_gc = GEOSCHEM["PEDGE"][iGC, jGC, :]

            # Get GEOS-Chem methane for the cell
            gc_CH4 = GEOSCHEM["CH4"][iGC, jGC, :]

            # Get merged GEOS-Chem/TROPOMI pressure grid for the cell
            merged = merge_pressure_grids(p_sat, p_gc)

            # Remap GEOS-Chem methane to TROPOMI pressure levels
            sat_CH4 = remap(
                gc_CH4,
                merged["data_type"],
                merged["p_merge"],
                merged["edge_index"],
                merged["first_gc_edge"],
            )  # ppb

            # Convert ppb to mol m-2
            sat_CH4_molm2 = sat_CH4 * 1e-9 * dry_air_subcolumns  # mol m-2

            # Derive the column-averaged XCH4 that TROPOMI would see over this ground cell
            # using eq. 46 from TROPOMI Methane ATBD, Hasekamp et al. 2019
            virtual_tropomi_gridcellIndex = (
                sum(apriori + avkern * (sat_CH4_molm2 - apriori))
                / sum(dry_air_subcolumns)
                * 1e9
            )  # ppb

            # Weight by overlapping area (to be divided out later) and add to sum
            area_weighted_virtual_tropomi += (
                overlap_area[gridcellIndex] * virtual_tropomi_gridcellIndex
            )  # ppb m2

            # If building Jacobian matrix from GEOS-Chem perturbation simulation sensitivity data:
            if build_jacobian:
                
                # Get GEOS-Chem perturbation sensitivities at this lat/lon, for all vertical levels and state vector elements
                sensi_lonlat = GEOSCHEM["jacobian_ch4"][iGC, jGC, :, :]

                # Map the sensitivities to TROPOMI pressure levels
                sat_deltaCH4 = remap_sensitivities(
                    sensi_lonlat,
                    merged["data_type"],
                    merged["p_merge"],
                    merged["edge_index"],
                    merged["first_gc_edge"],
                )  # mixing ratio, unitless

                # Derive the change in column-averaged XCH4 that TROPOMI would see over this ground cell
                # in shape of (n_elements,)
                tropomi_sensitivity_gridcellIndex = np.sum(avkern[:, None] * sat_deltaCH4 * \
                    dry_air_subcolumns[:, None], axis=0) / np.sum(dry_air_subcolumns)  # mixing ratio, unitless

                # Weight by overlapping area (to be divided out later) and add to sum
                area_weighted_virtual_tropomi_sensitivity += (
                    overlap_area[gridcellIndex] * tropomi_sensitivity_gridcellIndex
                )  # m2

        # Compute virtual TROPOMI observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_tropomi = area_weighted_virtual_tropomi / sum(overlap_area)
        # For global inversions, area of overlap should equal area of TROPOMI pixel
        # This is because the GEOS-Chem grid is continuous
        if dlon > 2.0:
            assert (
                abs(sum(overlap_area) - polygon_tropomi.area) / polygon_tropomi.area
                < 0.01
            ), f"ERROR: overlap area ({sum(overlap_area)}) /= satellite pixel area ({polygon_tropomi.area})"
        # Save actual and virtual TROPOMI data
        obs_GC[k, 0] = TROPOMI["methane"][
            iSat, jSat
        ]  # Actual TROPOMI methane column observation
        obs_GC[k, 1] = virtual_tropomi  # Virtual TROPOMI methane column observation
        obs_GC[k, 2] = TROPOMI["longitude"][iSat, jSat]  # TROPOMI longitude
        obs_GC[k, 3] = TROPOMI["latitude"][iSat, jSat]  # TROPOMI latitude
        obs_GC[k, 4] = iSat  # TROPOMI index of longitude
        obs_GC[k, 5] = jSat  # TROPOMI index of latitude

        if build_jacobian:
            # Compute TROPOMI sensitivity as weighted mean by overlapping area
            # i.e., need to divide out area [m2] from the previous step
            jacobian_K[k, :] = area_weighted_virtual_tropomi_sensitivity / sum(
                overlap_area
            )

    # Output
    output = {}

    # Always return the coincident TROPOMI and GEOS-Chem data
    output["obs_GC"] = obs_GC

    # Optionally return the Jacobian
    if build_jacobian:
        output["K"] = jacobian_K

    return output


def read_tropomi(filename):
    """
    Read TROPOMI data and save important variables to dictionary.

    Arguments
        filename [str]  : TROPOMI netcdf data file to read

    Returns
        dat      [dict] : Dictionary of important variables from TROPOMI:
                            - CH4
                            - Latitude
                            - Longitude
                            - QA value
                            - UTC time
                            - Time (utc time reshaped for orbit)
                            - Averaging kernel
                            - SWIR albedo
                            - NIR albedo
                            - Blended albedo
                            - CH4 prior profile
                            - Dry air subcolumns
                            - Latitude bounds
                            - Longitude bounds
                            - Vertical pressure profile
    """

    # Initialize dictionary for TROPOMI data
    dat = {}

    # Catch read errors in any of the variables
    try:
        # Store methane, QA, lat, lon, and time
        with xr.open_dataset(filename, group="PRODUCT") as tropomi_data:
            dat["methane"] = tropomi_data["methane_mixing_ratio_bias_corrected"].values[
                0, :, :
            ]
            dat["qa_value"] = tropomi_data["qa_value"].values[0, :, :]
            dat["longitude"] = tropomi_data["longitude"].values[0, :, :]
            dat["latitude"] = tropomi_data["latitude"].values[0, :, :]

            utc_str = tropomi_data["time_utc"].values[0, :]
            utc_str = np.array([d.replace("Z", "") for d in utc_str]).astype(
                "datetime64[ns]"
            )
            dat["time"] = np.repeat(
                utc_str[:, np.newaxis], dat["methane"].shape[1], axis=1
            )

        # Store column averaging kernel, SWIR and NIR surface albedo
        with xr.open_dataset(
            filename, group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS"
        ) as tropomi_data:
            dat["column_AK"] = tropomi_data["column_averaging_kernel"].values[
                0, :, :, ::-1
            ]
            dat["swir_albedo"] = tropomi_data["surface_albedo_SWIR"].values[0, :, :]
            dat["nir_albedo"] = tropomi_data["surface_albedo_NIR"].values[0, :, :]
            dat["blended_albedo"] = 2.4 * dat["nir_albedo"] - 1.13 * dat["swir_albedo"]

        # Store methane prior profile, dry air subcolumns
        with xr.open_dataset(
            filename, group="PRODUCT/SUPPORT_DATA/INPUT_DATA"
        ) as tropomi_data:
            dat["methane_profile_apriori"] = tropomi_data[
                "methane_profile_apriori"
            ].values[
                0, :, :, ::-1
            ]  # mol m-2
            dat["dry_air_subcolumns"] = tropomi_data["dry_air_subcolumns"].values[
                0, :, :, ::-1
            ]  # mol m-2

            # Surface classification values (ubyte type)
            sc_raw = tropomi_data["surface_classification"].values[0, :, :].astype("uint8")
            dat["surface_classification"] = (sc_raw & 3).astype(int)
            dat["surface_classification_249"] = (sc_raw & 249).astype(int)

            # Also get pressure interval and surface pressure for use below
            pressure_interval = (
                tropomi_data["pressure_interval"].values[0, :, :] / 100
            )  # Pa -> hPa
            surface_pressure = (
                tropomi_data["surface_pressure"].values[0, :, :] / 100
            )  # Pa -> hPa

        # Store latitude and longitude bounds for pixels
        with xr.open_dataset(
            filename, group="PRODUCT/SUPPORT_DATA/GEOLOCATIONS"
        ) as tropomi_data:
            dat["longitude_bounds"] = tropomi_data["longitude_bounds"].values[
                0, :, :, :
            ]
            dat["latitude_bounds"] = tropomi_data["latitude_bounds"].values[0, :, :, :]

        # Store vertical pressure profile
        n1 = dat["methane"].shape[
            0
        ]  # length of along-track dimension (scanline) of retrieval field
        n2 = dat["methane"].shape[
            1
        ]  # length of across-track dimension (ground_pixel) of retrieval field
        pressures = np.full([n1, n2, 12 + 1], np.nan, dtype=np.float32)
        for i in range(12 + 1):
            pressures[:, :, i] = surface_pressure - i * pressure_interval
        dat["pressures"] = pressures

    # Return an error if any of the variables were not read correctly
    except Exception as e:
        print(f"Error opening {filename}: {e}")
        return None

    return dat


def read_blended(filename):
    """
    Read Blended TROPOMI+GOSAT data and save important variables to dictionary.
    Arguments
        filename [str]  : Blended TROPOMI+GOSAT netcdf data file to read
    Returns
        dat      [dict] : Dictionary of important variables from Blended TROPOMI+GOSAT:
                            - CH4
                            - Latitude
                            - Longitude
                            - Time (utc time reshaped for orbit)
                            - Averaging kernel
                            - SWIR albedo
                            - NIR albedo
                            - Blended albedo
                            - CH4 prior profile
                            - Dry air subcolumns
                            - Latitude bounds
                            - Longitude bounds
                            - Surface classification
                            - Chi-Square for SWIR
                            - Vertical pressure profile
    """
    assert (
        "BLND" in filename
    ), f"BLND not in filename {filename}, but a blended function is being used"

    try:
        # Initialize dictionary for Blended TROPOMI+GOSAT data
        dat = {}

        # Extract data from netCDF file to our dictionary
        with xr.open_dataset(filename) as blended_data:

            dat["methane"] = blended_data["methane_mixing_ratio_blended"].values[:]
            dat["longitude"] = blended_data["longitude"].values[:]
            dat["latitude"] = blended_data["latitude"].values[:]
            dat["column_AK"] = blended_data["column_averaging_kernel"].values[:, ::-1]
            dat["swir_albedo"] = blended_data["surface_albedo_SWIR"][:]
            dat["nir_albedo"] = blended_data["surface_albedo_NIR"].values[:]
            dat["blended_albedo"] = 2.4 * dat["nir_albedo"] - 1.13 * dat["swir_albedo"]
            dat["methane_profile_apriori"] = blended_data[
                "methane_profile_apriori"
            ].values[:, ::-1]
            dat["dry_air_subcolumns"] = blended_data["dry_air_subcolumns"].values[
                :, ::-1
            ]
            dat["longitude_bounds"] = blended_data["longitude_bounds"].values[:]
            dat["latitude_bounds"] = blended_data["latitude_bounds"].values[:]
            
            # Surface classification values (ubyte type)
            sc_raw = blended_data["surface_classification"].values[:].astype("uint8")
            dat["surface_classification"] = (sc_raw & 3).astype(int)
            dat["surface_classification_249"] = (sc_raw & 249).astype(int)
            
            dat["chi_square_SWIR"] = blended_data["chi_square_SWIR"].values[:]

            # Remove "Z" from time so that numpy doesn't throw a warning
            utc_str = blended_data["time_utc"].values[:]
            dat["time"] = np.array([d.replace("Z", "") for d in utc_str]).astype(
                "datetime64[ns]"
            )

            # Need to calculate the pressure for the 13 TROPOMI levels (12 layer edges)
            pressure_interval = (
                blended_data["pressure_interval"].values[:] / 100
            )  # Pa -> hPa
            surface_pressure = (
                blended_data["surface_pressure"].values[:] / 100
            )  # Pa -> hPa
            n = len(dat["methane"])
            pressures = np.full([n, 12 + 1], np.nan, dtype=np.float32)
            for i in range(12 + 1):
                pressures[:, i] = surface_pressure - i * pressure_interval
            dat["pressures"] = pressures

        # Add an axis here to mimic the (scanline, groundpixel) format of operational TROPOMI data
        # This is so the blended data will be compatible with the TROPOMI operators
        for key in dat.keys():
            dat[key] = np.expand_dims(dat[key], axis=0)

    except Exception as e:
        print(f"Error opening {filename}: {e}")
        return None

    return dat


def average_tropomi_observations(TROPOMI, gc_lat_lon, sat_ind, time_threshold):
    """
    Map TROPOMI observations into appropriate gc gridcells. Then average all
    observations within a gridcell for processing. Use area weighting if
    observation overlaps multiple gridcells.

    Arguments
        TROPOMI        [dict]   : Dict of tropomi data
        gc_lat_lon     [list]   : list of dictionaries containing  gc gridcell info
        sat_ind        [int]    : index list of Tropomi data that passes filters

    Returns
        output         [dict[]]   : flat list of dictionaries the following fields:
                                    - lat                 : gridcell latitude
                                    - lon                 : gridcell longitude
                                    - iGC                 : longitude index value
                                    - jGC                 : latitude index value
                                    - lat_sat             : averaged tropomi latitude
                                    - lon_sat             : averaged tropomi longitude
                                    - overlap_area        : averaged overlap area with gridcell
                                    - p_sat               : averaged pressure for sat
                                    - dry_air_subcolumns  : averaged
                                    - apriori             : averaged
                                    - avkern              : averaged average kernel
                                    - time                : averaged time
                                    - methane             : averaged methane
                                    - observation_count   : number of observations averaged in cell
                                    - observation_weights : area weights for the observation

    """
    n_obs = len(sat_ind[0])
    # print("Found", n_obs, "TROPOMI observations.")
    gc_lats = gc_lat_lon["lat"]
    gc_lons = gc_lat_lon["lon"]
    dlon = np.median(np.diff(gc_lat_lon["lon"]))  # GEOS-Chem lon resolution
    dlat = np.median(np.diff(gc_lat_lon["lat"]))  # GEOS-Chem lon resolution
    gridcell_dicts = get_gridcell_list(gc_lons, gc_lats)

    for k in range(n_obs):
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index

        # Find GEOS-Chem lats & lons closest to the corners of the TROPOMI pixel
        longitude_bounds = TROPOMI["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = TROPOMI["latitude_bounds"][iSat, jSat, :]
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

        # Get lat/lon indexes and coordinates of GEOS-Chem grid cells closest to the TROPOMI corners
        ij_GC = [(x, y) for x in set(corners_lon_index) for y in set(corners_lat_index)]
        gc_coords = [(gc_lons[i], gc_lats[j]) for i, j in ij_GC]

        # Compute the overlapping area between the TROPOMI pixel and GEOS-Chem grid cells it touches
        overlap_area = np.zeros(len(gc_coords))

        # Polygon representing TROPOMI pixel
        polygon_tropomi = Polygon(np.column_stack((longitude_bounds, latitude_bounds)))
        
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
            if polygon_geoschem.intersects(polygon_tropomi):
                overlap_area[gridcellIndex] = polygon_tropomi.intersection(
                    polygon_geoschem
                ).area

        # If there is no overlap between GEOS-Chem and TROPOMI, skip to next observation:
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
                gridcell_dict["lat_sat"].append(TROPOMI["latitude"][iSat, jSat])
                gridcell_dict["lon_sat"].append(TROPOMI["longitude"][iSat, jSat])
                gridcell_dict["overlap_area"].append(overlap)
                gridcell_dict["p_sat"].append(TROPOMI["pressures"][iSat, jSat, :])
                gridcell_dict["dry_air_subcolumns"].append(
                    TROPOMI["dry_air_subcolumns"][iSat, jSat, :]
                )
                gridcell_dict["apriori"].append(
                    TROPOMI["methane_profile_apriori"][iSat, jSat, :]
                )
                gridcell_dict["avkern"].append(TROPOMI["column_AK"][iSat, jSat, :])
                gridcell_dict[
                    "time"
                ].append(  # convert times to epoch time to make taking the mean easier
                    int(pd.to_datetime(str(TROPOMI["time"][iSat, jSat])).strftime("%s"))
                )
                gridcell_dict["methane"].append(
                    TROPOMI["methane"][iSat, jSat]
                )  # Actual TROPOMI methane column observation
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
        gridcell_dict["methane"] = np.average(
            gridcell_dict["methane"],
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
    return gridcell_dicts

def get_virtual_tropomi(date, gc_cache, gridcell_dict, n_elements, config, build_jacobian=False):
    """
    Read GEOS-Chem data and save important variables to dictionary.

    Arguments
        date           [str]   : Date of interest, format "YYYYMMDD_HH"
        gc_cache       [str]   : Path to GEOS-Chem output data
        gridcell_dict  [dict]  : Dictionary of a specific gridcell with sampled satellite observations
        build_jacobian [log]   : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?

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

    with xr.open_dataset(filename) as gc_data_all:
        jGC   = xr.DataArray(gridcell_dict["jGC"],   dims="obs")
        iGC   = xr.DataArray(gridcell_dict["iGC"],   dims="obs")
        gc_data = gc_data_all.isel(
            time=0,
            lat=jGC,
            lon=iGC,
            drop=True
        )
        CH4_da = gc_data["SpeciesConcVV_CH4"]
        if CH4_da.dims[1] == "lev":
            CH4 = CH4_da.values
        else:
            CH4 = np.transpose(CH4_da.values, (1, 0))   # (lev, obs) -> (obs, lev)  

    # Read PEDGE from the LevelEdgeDiags collection
    filename = f"{gc_cache}/{file_pedge}"
    with xr.open_dataset(filename) as gc_data_all:
        gc_data = gc_data_all.isel(
            time=0,
            lat=jGC,
            lon=iGC,
            drop=True
        )
        PEDGE_da = gc_data["Met_PEDGE"]
        if PEDGE_da.dims[1] == "lev":
            PEDGE = PEDGE_da.values
        else:
            PEDGE = np.transpose(PEDGE_da.values, (1, 0))   # (lev, obs) -> (obs, lev)

    n_superobs = len(gridcell_dict)
    virtual_tropomi = np.empty([n_superobs, ], dtype=np.float32)
    virtual_tropomi.fill(np.nan)
    
    p_sat = gridcell_dict["p_sat"]
    dry_air_subcolumns = gridcell_dict["dry_air_subcolumns"]  # mol m-2
    apriori = gridcell_dict["apriori"]  # mol m-2
    avkern = gridcell_dict["avkern"]

    # (N, S, G)  sums to 1 along G, where 
    # N is the number of super observations
    # S is the number of TROPOMI pressure edges
    # G is the number of GEOS-Chem pressure edges
    vertical_weights = remapping_weights(p_sat, PEDGE)
    sat_CH4 = np.einsum("nsg,ng->ns", vertical_weights, CH4)         # (N, S)
    sat_CH4_molm2 = sat_CH4 * dry_air_subcolumns                     # (N, S)
    virtual_tropomi = (
        np.sum(apriori + avkern * (sat_CH4_molm2 - apriori), axis=1) / \
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
        virtual_tropomi_pert = [
            get_virtual_tropomi_pert(gc_date, k, gridcell_dict, config, v, n_elements, vertical_weights)
            for k, v in pert_simulations_dict.items()
        ]

        virtual_tropomi_pert = np.concatenate(virtual_tropomi_pert, axis=1)
        
        virtual_tropomi_base = get_virtual_tropomi_pert(
            gc_date, "0001", gridcell_dict, config, [0], n_elements, vertical_weights, baserun=True
        ).squeeze() # make it 1D
    
    if build_jacobian:
        return virtual_tropomi_pert, virtual_tropomi_base, virtual_tropomi
    else:
        return virtual_tropomi

def get_virtual_tropomi_pert(gc_date, run_id, gridcell_dict, config, sv_elems, n_elements, vertical_weights, baserun=False):
    """
    Concatenate CH4 tracers from all jacobian GEOS-Chem simulations.
    Tracers are assigned a new dimension: "element"

    Arguments
        gc_date    [pd.Datetime] : date object, specifies Ymd_h
        run_id     [str]         : ID for Jacobian GEOS-Chem run, e.g. "0001"
        gridcell_dict  [dict]    : Dictionary of a specific gridcell with sampled satellite observations
        config     [dict]        : dictionary of IMI config file
        sv_elems   [list]        : list of state vector element tracers in this simulations
        n_elements [int]         : number of state vector elements in this inversion
        vertical_weights [np.array]: vertical weights to remap GC to TROPOMI
        baserun    [bool]        : If True, only the base variable in the simulation will
                                 be opened, and the function will just return this one
                                 variable instead of concatenating all elements. Used to
                                 get the base for calculating the sensitivities.

    Returns
        virtual_tropomi_all [np.array] : array of all virtual TROPOMI at this timestep for
                                     all tracer runs. Has dimensions
                                        - n_superobs
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

    with xr.open_dataset(filepath, decode_cf=False) as tmp:
        other_vars = [v for v in tmp.variables if "SpeciesConcVV_CH4" not in v]
    
    # Open only these variables
    with xr.open_dataset(
        filepath,
        drop_variables=other_vars,
    ) as dsmf_all:
        try:
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
        
        n_superobs = len(gridcell_dict)
        n_ele = len(keepvars)
        virtual_tropomi_all = np.empty([n_superobs, n_ele], dtype=np.float32)
        virtual_tropomi_all.fill(np.nan)
        
        for elei in range(n_ele):
            CH4_da = dsmf[keepvars[elei]]
            if CH4_da.dims[1] == "lev":
                CH4 = CH4_da.values
            else:
                CH4 = np.transpose(CH4_da.values, (1, 0))   # (lev, obs) -> (obs, lev)
            
            dry_air_subcolumns = gridcell_dict["dry_air_subcolumns"]  # mol m-2
            apriori = gridcell_dict["apriori"]  # mol m-2
            avkern = gridcell_dict["avkern"]
            
            sat_CH4 = np.einsum("nsg,ng->ns", vertical_weights, CH4)         # (N, S)
            sat_CH4_molm2 = sat_CH4 * dry_air_subcolumns                     # (N, S)
            virtual_tropomi_all[:,elei] = (
                np.sum(apriori + avkern * (sat_CH4_molm2 - apriori), axis=1) / \
                np.sum(dry_air_subcolumns, axis=1)
            ).astype(np.float32)                      # (N,), unitless mixing ratio
                
    return virtual_tropomi_all
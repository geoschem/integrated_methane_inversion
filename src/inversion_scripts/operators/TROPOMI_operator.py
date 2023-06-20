import numpy as np
import xarray as xr
import pandas as pd
import datetime
from shapely.geometry import Polygon
from utils import filter_tropomi
from operators.operator_utilities import (
    get_gc_lat_lon,
    read_all_geoschem,
    merge_pressure_grids,
    remap,
    remap_sensitivities,
    get_gridcell_list,
    nearest_loc,
)


def apply_average_tropomi_operator(
    filename,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    build_jacobian,
    sensi_cache,
):
    """
    Apply the averaging tropomi operator to map GEOS-Chem methane data to TROPOMI observation space.

    Arguments
        filename       [str]        : TROPOMI netcdf data file to read
        n_elements     [int]        : Number of state vector elements
        gc_startdate   [datetime64] : First day of inversion period, for GEOS-Chem and TROPOMI
        gc_enddate     [datetime64] : Last day of inversion period, for GEOS-Chem and TROPOMI
        xlim           [float]      : Longitude bounds for simulation domain
        ylim           [float]      : Latitude bounds for simulation domain
        gc_cache       [str]        : Path to GEOS-Chem output data
        build_jacobian [log]        : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        sensi_cache    [str]        : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data

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
    TROPOMI = read_tropomi(filename)
    if TROPOMI == None:
        print(f"Skipping {filename} due to file processing issue.")
        return TROPOMI

    # We're only going to consider data within lat/lon/time bounds, with QA > 0.5, and with safe surface albedo values
    sat_ind = filter_tropomi(TROPOMI, xlim, ylim, gc_startdate, gc_enddate)

    # get the lat/lons of gc gridcells
    gc_lat_lon = get_gc_lat_lon(gc_cache, gc_startdate)

    # map tropomi obs into gridcells and average the observations
    # into each gridcell. Only returns gridcells containing observations
    obs_mapped_to_gc = average_tropomi_observations(TROPOMI, gc_lat_lon, sat_ind)
    n_gridcells = len(obs_mapped_to_gc)

    if build_jacobian:
        # Initialize Jacobian K
        jacobian_K = np.zeros([n_gridcells, n_elements], dtype=np.float32)
        jacobian_K.fill(np.nan)

    # create list to store the dates/hour of each gridcell
    all_strdate = [gridcell["time"] for gridcell in obs_mapped_to_gc]
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_gc = read_all_geoschem(all_strdate, gc_cache, build_jacobian, sensi_cache)

    # Initialize array with n_gridcells rows and 5 columns. Columns are TROPOMI CH4, GEOSChem CH4, longitude, latitude, observation counts
    obs_GC = np.zeros([n_gridcells, 5], dtype=np.float32)
    obs_GC.fill(np.nan)

    # For each gridcell dict with tropomi obs:
    for i, gridcell_dict in enumerate(obs_mapped_to_gc):

        # Get GEOS-Chem data for the date of the observation:
        p_sat = gridcell_dict["p_sat"]
        dry_air_subcolumns = gridcell_dict["dry_air_subcolumns"]  # mol m-2
        apriori = gridcell_dict["apriori"]  # mol m-2
        avkern = gridcell_dict["avkern"]
        strdate = gridcell_dict["time"]
        GEOSCHEM = all_date_gc[strdate]

        # Get GEOS-Chem pressure edges for the cell
        p_gc = GEOSCHEM["PEDGE"][gridcell_dict["iGC"], gridcell_dict["jGC"], :]
        # Get GEOS-Chem methane for the cell
        gc_CH4 = GEOSCHEM["CH4"][gridcell_dict["iGC"], gridcell_dict["jGC"], :]
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
        virtual_tropomi = (
            sum(apriori + avkern * (sat_CH4_molm2 - apriori))
            / sum(dry_air_subcolumns)
            * 1e9
        )  # ppb

        # If building Jacobian matrix from GEOS-Chem perturbation simulation sensitivity data:
        if build_jacobian:
            # Get GEOS-Chem perturbation sensitivities at this lat/lon, for all vertical levels and state vector elements
            sensi_lonlat = GEOSCHEM["Sensitivities"][
                gridcell_dict["iGC"], gridcell_dict["jGC"], :, :
            ]
            # Map the sensitivities to TROPOMI pressure levels
            sat_deltaCH4 = remap_sensitivities(
                sensi_lonlat,
                merged["data_type"],
                merged["p_merge"],
                merged["edge_index"],
                merged["first_gc_edge"],
            )  # mixing ratio, unitless
            # Tile the TROPOMI averaging kernel
            avkern_tiled = np.transpose(np.tile(avkern, (n_elements, 1)))
            # Tile the TROPOMI dry air subcolumns
            dry_air_subcolumns_tiled = np.transpose(
                np.tile(dry_air_subcolumns, (n_elements, 1))
            )  # mol m-2
            # Derive the change in column-averaged XCH4 that TROPOMI would see over this ground cell
            jacobian_K[i, :] = np.sum(
                avkern_tiled * sat_deltaCH4 * dry_air_subcolumns_tiled, 0
            ) / sum(
                dry_air_subcolumns
            )  # mixing ratio, unitless

        # Save actual and virtual TROPOMI data
        obs_GC[i, 0] = gridcell_dict[
            "methane"
        ]  # Actual TROPOMI methane column observation
        obs_GC[i, 1] = virtual_tropomi  # Virtual TROPOMI methane column observation
        obs_GC[i, 2] = gridcell_dict["lon_sat"]  # TROPOMI longitude
        obs_GC[i, 3] = gridcell_dict["lat_sat"]  # TROPOMI latitude
        obs_GC[i, 4] = gridcell_dict["observation_count"]  # observation counts

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
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    build_jacobian,
    sensi_cache,
):
    """
    Apply the tropomi operator to map GEOS-Chem methane data to TROPOMI observation space.

    Arguments
        filename       [str]        : TROPOMI netcdf data file to read
        n_elements     [int]        : Number of state vector elements
        gc_startdate   [datetime64] : First day of inversion period, for GEOS-Chem and TROPOMI
        gc_enddate     [datetime64] : Last day of inversion period, for GEOS-Chem and TROPOMI
        xlim           [float]      : Longitude bounds for simulation domain
        ylim           [float]      : Latitude bounds for simulation domain
        gc_cache       [str]        : Path to GEOS-Chem output data
        build_jacobian [log]        : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        sensi_cache    [str]        : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data

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
    TROPOMI = read_tropomi(filename)
    if TROPOMI == None:
        print(f"Skipping {filename} due to file processing issue.")
        return TROPOMI

    # We're only going to consider data within lat/lon/time bounds, with QA > 0.5, and with safe surface albedo values
    sat_ind = filter_tropomi(TROPOMI, xlim, ylim, gc_startdate, gc_enddate)

    # Number of TROPOMI observations
    n_obs = len(sat_ind[0])
    print("Found", n_obs, "TROPOMI observations.")

    # If need to build Jacobian from GEOS-Chem perturbation simulation sensitivity data:
    if build_jacobian:
        # Initialize Jacobian K
        jacobian_K = np.zeros([n_obs, n_elements], dtype=np.float32)
        jacobian_K.fill(np.nan)

    # Initialize a list to store the dates we want to look at
    all_strdate = []

    # For each TROPOMI observation
    for k in range(n_obs):
        # Get the date and hour
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index
        time = pd.to_datetime(str(TROPOMI["utctime"][iSat]))
        strdate = time.round("60min").strftime("%Y%m%d_%H")
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_gc = read_all_geoschem(all_strdate, gc_cache, build_jacobian, sensi_cache)

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
        time = pd.to_datetime(str(TROPOMI["utctime"][iSat]))
        strdate = time.round("60min").strftime("%Y%m%d_%H")
        GEOSCHEM = all_date_gc[strdate]
        dlon = np.median(np.diff(GEOSCHEM["lon"])) # GEOS-Chem lon resolution
        dlat = np.median(np.diff(GEOSCHEM["lat"])) # GEOS-Chem lon resolution

        # Find GEOS-Chem lats & lons closest to the corners of the TROPOMI pixel
        longitude_bounds = TROPOMI["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = TROPOMI["latitude_bounds"][iSat, jSat, :]
        corners_lon_index = []
        corners_lat_index = []
        for l in range(4):
            iGC = nearest_loc(longitude_bounds[l], GEOSCHEM["lon"], tolerance=max(dlon,0.5))
            jGC = nearest_loc(latitude_bounds[l], GEOSCHEM["lat"], tolerance=max(dlat,0.5))
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
            polygon_geoschem = Polygon(
                np.column_stack((geoschem_corners_lon, geoschem_corners_lat))
            )
            # Calculate overlapping area as the intersection of the two polygons
            if polygon_geoschem.intersects(polygon_tropomi):
                overlap_area[gridcellIndex] = polygon_tropomi.intersection(
                    polygon_geoschem
                ).area

        # If there is no overlap between GEOS-Chem and TROPOMI, skip to next observation:
        if sum(overlap_area) == 0:
            continue

        # =======================================================
        #       Map GEOS-Chem to TROPOMI observation space
        # =======================================================

        # Otherwise, initialize tropomi virtual xch4 and virtual sensitivity as zero
        area_weighted_virtual_tropomi = 0  # virtual tropomi xch4
        area_weighted_virtual_tropomi_sensitivity = 0  # virtual tropomi sensitivity

        # For each GEOS-Chem grid cell that touches the TROPOMI pixel:
        for gridcellIndex in range(len(gc_coords)):

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
                sensi_lonlat = GEOSCHEM["Sensitivities"][iGC, jGC, :, :]

                # Map the sensitivities to TROPOMI pressure levels
                sat_deltaCH4 = remap_sensitivities(
                    sensi_lonlat,
                    merged["data_type"],
                    merged["p_merge"],
                    merged["edge_index"],
                    merged["first_gc_edge"],
                )  # mixing ratio, unitless

                # Tile the TROPOMI averaging kernel
                avkern_tiled = np.transpose(np.tile(avkern, (n_elements, 1)))

                # Tile the TROPOMI dry air subcolumns
                dry_air_subcolumns_tiled = np.transpose(
                    np.tile(dry_air_subcolumns, (n_elements, 1))
                )  # mol m-2

                # Derive the change in column-averaged XCH4 that TROPOMI would see over this ground cell
                tropomi_sensitivity_gridcellIndex = np.sum(
                    avkern_tiled * sat_deltaCH4 * dry_air_subcolumns_tiled, 0
                ) / sum(
                    dry_air_subcolumns
                )  # mixing ratio, unitless

                # Weight by overlapping area (to be divided out later) and add to sum
                area_weighted_virtual_tropomi_sensitivity += (
                    overlap_area[gridcellIndex] * tropomi_sensitivity_gridcellIndex
                )  # m2

        # Compute virtual TROPOMI observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_tropomi = area_weighted_virtual_tropomi / sum(overlap_area)

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
        # Store methane, QA, lat, lon
        tropomi_data = xr.open_dataset(filename, group="PRODUCT")
        dat["methane"] = tropomi_data["methane_mixing_ratio_bias_corrected"].values[0, :, :]
        dat["qa_value"] = tropomi_data["qa_value"].values[0, :, :]
        dat["longitude"] = tropomi_data["longitude"].values[0, :, :]
        dat["latitude"] = tropomi_data["latitude"].values[0, :, :]

        # Store UTC time
        reference_time = tropomi_data["time"].values
        delta_time = tropomi_data["delta_time"][0].values
        strdate = []
        if delta_time.dtype == "<m8[ns]":
            strdate = reference_time + delta_time
        elif delta_time.dtype == "<M8[ns]":
            strdate = delta_time
        else:
            print(delta_time.dtype)
            pass
        dat["utctime"] = strdate

        # Store time for whole orbit
        times = np.zeros(shape=dat["longitude"].shape, dtype="datetime64[ns]")
        for k in range(dat["longitude"].shape[0]):
            times[k, :] = strdate[k]
        dat["time"] = times
        tropomi_data.close()

        # Store column averaging kernel, SWIR and NIR surface albedo
        tropomi_data = xr.open_dataset(
            filename, group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS"
        )
        dat["column_AK"] = tropomi_data["column_averaging_kernel"].values[0, :, :, ::-1]
        dat["swir_albedo"] = tropomi_data["surface_albedo_SWIR"].values[0, :, :]
        dat["nir_albedo"] = tropomi_data["surface_albedo_NIR"].values[0, :, :]
        dat["blended_albedo"] = 2.4 * dat["nir_albedo"] - 1.13 * dat["swir_albedo"]
        tropomi_data.close()

        # Store methane prior profile, dry air subcolumns
        tropomi_data = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/INPUT_DATA")
        dat["methane_profile_apriori"] = tropomi_data["methane_profile_apriori"].values[
            0, :, :, ::-1
        ]  # mol m-2
        dat["dry_air_subcolumns"] = tropomi_data["dry_air_subcolumns"].values[
            0, :, :, ::-1
        ]  # mol m-2

        # Also get pressure interval and surface pressure for use below
        pressure_interval = (
            tropomi_data["pressure_interval"].values[0, :, :] / 100
        )  # Pa -> hPa
        surface_pressure = (
            tropomi_data["surface_pressure"].values[0, :, :] / 100
        )  # Pa -> hPa
        tropomi_data.close()

        # Store latitude and longitude bounds for pixels
        tropomi_data = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/GEOLOCATIONS")
        dat["longitude_bounds"] = tropomi_data["longitude_bounds"].values[0, :, :, :]
        dat["latitude_bounds"] = tropomi_data["latitude_bounds"].values[0, :, :, :]
        tropomi_data.close()

        # Store vertical pressure profile
        n1 = dat["methane"].shape[
            0
        ]  # length of along-track dimension (scanline) of retrieval field
        n2 = dat["methane"].shape[
            1
        ]  # length of across-track dimension (ground_pixel) of retrieval field
        pressures = np.zeros([n1, n2, 12 + 1], dtype=np.float32)
        pressures.fill(np.nan)
        for i in range(12 + 1):
            pressures[:, :, i] = surface_pressure - i * pressure_interval
        dat["pressures"] = pressures

    # Return an error if any of the variables were not read correctly
    except Exception as e:
        print(f"Error opening {filename}: {e}")
        return None

    return dat

def average_tropomi_observations(TROPOMI, gc_lat_lon, sat_ind):
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
    print("Found", n_obs, "TROPOMI observations.")
    gc_lats = gc_lat_lon["lat"]
    gc_lons = gc_lat_lon["lon"]
    dlon = np.median(np.diff(gc_lat_lon["lon"])) # GEOS-Chem lon resolution
    dlat = np.median(np.diff(gc_lat_lon["lat"])) # GEOS-Chem lon resolution
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
            iGC = nearest_loc(longitude_bounds[l], gc_lons, tolerance=max(dlon,0.5))
            jGC = nearest_loc(latitude_bounds[l], gc_lats, tolerance=max(dlat,0.5))
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
                    int(pd.to_datetime(str(TROPOMI["utctime"][iSat])).strftime("%s"))
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
        gridcell_dict["time"] = (
            pd.to_datetime(
                datetime.datetime.fromtimestamp(int(np.mean(gridcell_dict["time"])))
            )
            .round("60min")
            .strftime("%Y%m%d_%H")
        )
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

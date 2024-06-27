import numpy as np
import xarray as xr
import pandas as pd
import datetime
from shapely.geometry import Polygon
from src.inversion_scripts.utils import (
    read_and_filter_satellite,
    mixing_ratio_conv_factor,
)
from src.inversion_scripts.operators.operator_utilities import (
    get_gc_lat_lon,
    read_all_geoschem,
    merge_pressure_grids,
    remap,
    remap_sensitivities,
    get_gridcell_list,
    nearest_loc,
)


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
    sensi_cache,
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
        sensi_cache       [str]        : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data

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
        filename, satellite_product, gc_startdate, gc_enddate, xlim, ylim)
    
    # Number of satellite observations
    n_obs = len(sat_ind[0])
    print("Found", n_obs, "satellite observations.")

    # get the lat/lons of gc gridcells
    gc_lat_lon = get_gc_lat_lon(gc_cache, gc_startdate)

    # map satellite obs into gridcells and average the observations
    # into each gridcell. Only returns gridcells containing observations
    obs_mapped_to_gc = average_satellite_observations(satellite, species, gc_lat_lon, sat_ind)
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

    # Initialize array with n_gridcells rows and 5 columns. Columns are 
    # satellite gas, GEOSChem gas, longitude, latitude, observation counts
    obs_GC = np.zeros([n_gridcells, 5], dtype=np.float32)
    obs_GC.fill(np.nan)

    # For each gridcell dict with satellite obs:
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
        # Get GEOS-Chem species for the cell
        gc_species = GEOSCHEM[species][gridcell_dict["iGC"], gridcell_dict["jGC"], :]
        # Get merged GEOS-Chem/satellite pressure grid for the cell
        merged = merge_pressure_grids(p_sat, p_gc)
        # Remap GEOS-Chem species to TROPOMI pressure levels
        sat_species = remap(
            gc_species,
            merged["data_type"],
            merged["p_merge"],
            merged["edge_index"],
            merged["first_gc_edge"],
        )  # volumetric mixing ratio
        # Convert volumetric mixing ratio to mol m-2
        sat_species_molm2 = sat_species * 1/mixing_ratio_conv_factor(species) * dry_air_subcolumns  # mol m-2
        # Derive the column-averaged mixing ratio that the satellite would see 
        # over this ground cell 
        virtual_satellite = apply_averaging_kernel(
            apriori, avkern, sat_species_molm2, dry_air_subcolumns, species)
        # Volumetric mixing ratio

        # If building Jacobian matrix from GEOS-Chem perturbation simulation sensitivity data:
        if build_jacobian:
            # Get GEOS-Chem perturbation sensitivities at this lat/lon, for all vertical levels and state vector elements
            sensi_lonlat = GEOSCHEM["Sensitivities"][
                gridcell_dict["iGC"], gridcell_dict["jGC"], :, :
            ]
            # Map the sensitivities to satellite pressure levels
            sat_deltaspecies = remap_sensitivities(
                sensi_lonlat,
                merged["data_type"],
                merged["p_merge"],
                merged["edge_index"],
                merged["first_gc_edge"],
            )  # mixing ratio, unitless
            # Tile the satellite averaging kernel
            avkern_tiled = np.transpose(np.tile(avkern, (n_elements, 1)))
            # Tile the satellite dry air subcolumns
            dry_air_subcolumns_tiled = np.transpose(
                np.tile(dry_air_subcolumns, (n_elements, 1))
            )  # mol m-2
            # Derive the change in column-averaged mixing ratios that TROPOMI would 
            # see over this ground cell
            jacobian_K[i, :] = np.sum(
                avkern_tiled * sat_deltaspecies * dry_air_subcolumns_tiled, 0
            ) / sum(
                dry_air_subcolumns
            )  # mixing ratio, unitless

        # Save actual and virtual satellite data
        obs_GC[i, 0] = gridcell_dict[
            "species"
        ]  # Actual satellite species column observation
        obs_GC[i, 1] = virtual_satellite  # Virtual satellite column observation
        obs_GC[i, 2] = gridcell_dict["lon_sat"]  # satellite longitude
        obs_GC[i, 3] = gridcell_dict["lat_sat"]  # satellite latitude
        obs_GC[i, 4] = gridcell_dict["observation_count"]  # observation counts

    # Output
    output = {}

    # Always return the coincident satellite and GEOS-Chem data
    output["obs_GC"] = obs_GC

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
    build_jacobian,
    sensi_cache,
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
        build_jacobian     [log]        : Are we trying to map GEOS-Chem sensitivities to satellite observation space?
        sensi_cache        [str]        : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data

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
        filename, satellite_product, gc_startdate, gc_enddate, xlim, ylim)

    # Number of satellite observations
    n_obs = len(sat_ind[0])
    # print("Found", n_obs, "satellite observations.")

    # If need to build Jacobian from GEOS-Chem perturbation simulation sensitivity data:
    if build_jacobian:
        # Initialize Jacobian K
        jacobian_K = np.zeros([n_obs, n_elements], dtype=np.float32)
        jacobian_K.fill(np.nan)

    # Initialize a list to store the dates we want to look at
    all_strdate = []

    # For each satellite observation
    for k in range(n_obs):
        # Get the date and hour
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index
        time = pd.to_datetime(str(satellite["time"][iSat,jSat]))
        strdate = time.round("60min").strftime("%Y%m%d_%H")
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_gc = read_all_geoschem(all_strdate, gc_cache, build_jacobian, sensi_cache)

    # Initialize array with n_obs rows and 6 columns. Columns are satellite 
    # mixing ratio, GEOSChem mixing ratio, longitude, latitude, II, JJ
    obs_GC = np.zeros([n_obs, 6], dtype=np.float32)
    obs_GC.fill(np.nan)

    # For each satellite observation:
    for k in range(n_obs):

        # Get GEOS-Chem data for the date of the observation:
        iSat = sat_ind[0][k]
        jSat = sat_ind[1][k]
        p_sat = satellite["pressures"][iSat, jSat, :]
        dry_air_subcolumns = satellite["dry_air_subcolumns"][iSat, jSat, :]  # mol m-2
        apriori = satellite["profile_apriori"][iSat, jSat, :]  # mol m-2
        avkern = satellite["column_AK"][iSat, jSat, :]
        time = pd.to_datetime(str(satellite["time"][iSat,jSat]))
        strdate = time.round("60min").strftime("%Y%m%d_%H")
        GEOSCHEM = all_date_gc[strdate]
        dlon = np.median(np.diff(GEOSCHEM["lon"])) # GEOS-Chem lon resolution
        dlat = np.median(np.diff(GEOSCHEM["lat"])) # GEOS-Chem lon resolution

        # Find GEOS-Chem lats & lons closest to the corners of the satellite pixel
        longitude_bounds = satellite["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = satellite["latitude_bounds"][iSat, jSat, :]
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
        # Get lat/lon indexes and coordinates of GEOS-Chem grid cells closest to the satellite corners
        ij_GC = [(x, y) for x in set(corners_lon_index) for y in set(corners_lat_index)]
        gc_coords = [(GEOSCHEM["lon"][i], GEOSCHEM["lat"][j]) for i, j in ij_GC]

        # Compute the overlapping area between the satellite pixel and GEOS-Chem grid cells it touches
        overlap_area = np.zeros(len(gc_coords))
        # Polygon representing satellite pixel
        polygon_satellite = Polygon(np.column_stack((longitude_bounds, latitude_bounds)))
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
            if (dlon == 2.5) & (coords[0] == 177.5):
                for i in [1,2]:
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
        if sum(overlap_area) == 0:
            continue

        # =======================================================
        #       Map GEOS-Chem to satellite observation space
        # =======================================================

        # Otherwise, initialize satellite virtual mixing ratios and virtual 
        #  sensitivity as zero
        area_weighted_virtual_satellite = 0  # virtual satellite mixing ratio
        area_weighted_virtual_satellite_sensitivity = 0  # virtual satellite sensitivity

        # For each GEOS-Chem grid cell that touches the satellite pixel:
        for gridcellIndex in range(len(gc_coords)):

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

                # Get GEOS-Chem perturbation sensitivities at this lat/lon, for all vertical levels and state vector elements
                sensi_lonlat = GEOSCHEM["Sensitivities"][iGC, jGC, :, :]

                # Map the sensitivities to satellite pressure levels
                sat_deltaspecies = remap_sensitivities(
                    sensi_lonlat,
                    merged["data_type"],
                    merged["p_merge"],
                    merged["edge_index"],
                    merged["first_gc_edge"],
                )  # mixing ratio, unitless

                # Tile the satellite averaging kernel
                avkern_tiled = np.transpose(np.tile(avkern, (n_elements, 1)))

                # Tile the satellite dry air subcolumns
                dry_air_subcolumns_tiled = np.transpose(
                    np.tile(dry_air_subcolumns, (n_elements, 1))
                )  # mol m-2

                # Derive the change in column-averaged mixing ratio that the 
                # satellite would see over this ground cell
                satellite_sensitivity_gridcellIndex = np.sum(
                    avkern_tiled * sat_deltaspecies * dry_air_subcolumns_tiled, 0
                ) / sum(
                    dry_air_subcolumns
                )  # mixing ratio, unitless

                # Weight by overlapping area (to be divided out later) and add to sum
                area_weighted_virtual_satellite_sensitivity += (
                    overlap_area[gridcellIndex] * satellite_sensitivity_gridcellIndex
                )  # m2

        # Compute virtual satellite observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_satellite = area_weighted_virtual_satellite / sum(overlap_area)

        # For global inversions, area of overlap should equal area of satellite pixel
        # This is because the GEOS-Chem grid is continuous
        if dlon > 2.0:
            assert abs(sum(overlap_area)-polygon_satellite.area)/polygon_satellite.area < 0.01, f"ERROR: overlap area ({sum(overlap_area)}) /= satellite pixel area ({polygon_satellite.area})"

        # Save actual and virtual satellite data
        obs_GC[k, 0] = satellite[species][
            iSat, jSat
        ]  # Actual satellite mixing ratio column observation
        obs_GC[k, 1] = virtual_satellite  # Virtual satellite mixing ratio column observation
        obs_GC[k, 2] = satellite["longitude"][iSat, jSat]  # satellite longitude
        obs_GC[k, 3] = satellite["latitude"][iSat, jSat]  # satellite latitude
        obs_GC[k, 4] = iSat  # satellite index of longitude
        obs_GC[k, 5] = jSat  # satellite index of latitude

        if build_jacobian:
            # Compute satellite sensitivity as weighted mean by overlapping area
            # i.e., need to divide out area [m2] from the previous step
            jacobian_K[k, :] = area_weighted_virtual_satellite_sensitivity / sum(
                overlap_area
            )

    # Output
    output = {}

    # Always return the coincident satellite and GEOS-Chem data
    output["obs_GC"] = obs_GC

    # Optionally return the Jacobian
    if build_jacobian:
        output["K"] = jacobian_K

    return output


def average_satellite_observations(satellite, species, gc_lat_lon, sat_ind):
    """
    Map TROPOMI observations into appropriate gc gridcells. Then average all
    observations within a gridcell for processing. Use area weighting if
    observation overlaps multiple gridcells.

    Arguments
        satellite      [dict]   : Dict of satellite data
        species        [str]    : Name of species analyzed (CO2 or CH4)
        gc_lat_lon     [list]   : list of dictionaries containing  gc gridcell info
        sat_ind        [int]    : index list of satellite data that passes filters

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
                                    - $species            : averaged species
                                    - observation_count   : number of observations averaged in cell
                                    - observation_weights : area weights for the observation

    """
    n_obs = len(sat_ind[0])
    # print("Found", n_obs, "satellite observations.")
    gc_lats = gc_lat_lon["lat"]
    gc_lons = gc_lat_lon["lon"]
    dlon = np.median(np.diff(gc_lat_lon["lon"])) # GEOS-Chem lon resolution
    dlat = np.median(np.diff(gc_lat_lon["lat"])) # GEOS-Chem lon resolution
    gridcell_dicts = get_gridcell_list(gc_lons, gc_lats)

    for k in range(n_obs):
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index

        # Find GEOS-Chem lats & lons closest to the corners of the satellite pixel
        longitude_bounds = satellite["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = satellite["latitude_bounds"][iSat, jSat, :]
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
            gridcell_dict["lat_sat"], weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["lon_sat"] = np.average(
            gridcell_dict["lon_sat"], weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["overlap_area"] = np.average(
            gridcell_dict["overlap_area"], weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict[species] = np.average(
            gridcell_dict[species], weights=gridcell_dict["observation_weights"],
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

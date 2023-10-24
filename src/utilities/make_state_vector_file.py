import numpy as np
import xarray as xr
from sklearn.cluster import KMeans
import yaml


def get_grid_bounds(land_cover_pth):
    """
    Get the lat/lon bounds of the grid window for the inversion.
    The land cover file path specifies the window.
    """

    land_cover = xr.load_dataset(land_cover_pth)
    minLat_allowed = np.min(land_cover.lat.values)
    maxLat_allowed = np.max(land_cover.lat.values)
    minLon_allowed = np.min(land_cover.lon.values)
    maxLon_allowed = np.max(land_cover.lon.values)

    return minLat_allowed, maxLat_allowed, minLon_allowed, maxLon_allowed


def check_grid_compatibility(lat_min, lat_max, lon_min, lon_max, land_cover_pth):
    """
    Check whether input lat/lon bounds are compatible with (contained within) the regional grid window.
    The land cover file path specifies the window.
    """

    (
        minLat_allowed,
        maxLat_allowed,
        minLon_allowed,
        maxLon_allowed,
    ) = get_grid_bounds(land_cover_pth)

    if (
        lat_min < minLat_allowed
        or lat_max > maxLat_allowed
        or lon_min < minLon_allowed
        or lon_max > maxLon_allowed
    ):
        compatible = False
    else:
        compatible = True

    return compatible


def make_state_vector_file(
    config_path,
    land_cover_pth,
    hemco_diag_pth,
    save_pth,
):
    """
    Generates the state vector file for an analytical inversion.

    Arguments
        config_path    [str]   : Path to configuration file
        land_cover_pth [str]   : Path to land cover file
        hemco_diag_pth [str]   : Path to initial HEMCO diagnostics file
        save_pth       [str]   : Where to save the state vector file

    Returns
        ds_statevector []     : xarray dataset containing state vector field formatted for HEMCO

    Notes
        - Land cover file looks like 'GEOSFP.20200101.CN.025x03125.NA.nc' (or 0.5-deg equivalent)
        - HEMCO diags file needs to be global, is used to include offshore emissions in state vector
        - Land cover file and HEMCO diags file need to have the same grid resolution
    """

    # Get config
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    lat_min = config["LatMin"]
    lat_max = config["LatMax"]
    lon_min = config["LonMin"]
    lon_max = config["LonMax"]
    is_regional=config["isRegional"]
    buffer_deg = config["BufferDeg"]
    land_threshold = config["LandThreshold"]
    emis_threshold = config["OffshoreEmisThreshold"]
    k_buffer_clust = config["nBufferClusters"]

    # Load land cover data and HEMCO diagnostics
    lc_orig = xr.load_dataset(land_cover_pth)
    hd = xr.load_dataset(hemco_diag_pth)

    # Require hemco diags on same global grid as land cover map
    hd["lon"] = hd["lon"] - 0.03125  # initially offset by 0.03125 degrees

    # Select / group fields together
    lc = (lc_orig["FRLAKE"] + lc_orig["FRLAND"] + lc_orig["FRLANDIC"]).drop("time").squeeze()
    hd = (hd["EmisCH4_Oil"] + hd["EmisCH4_Gas"]).drop("time").squeeze()

    # Check compatibility of region of interest
    if is_regional:
       compatible = check_grid_compatibility(
           lat_min, lat_max, lon_min, lon_max, land_cover_pth
       )
       if not compatible:
           raise ValueError(
               "Region of interest not contained within selected RegionID; see config.yml)."
           )

    # For global inversions exclude Antarctica and limit max_lat to avoid
    # HEMCO error when reading netCDF file
    if not is_regional:
        lat_max = 88.0
        lat_min = -60.0

    # Define bounds of inversion domain
    (
        minLat_allowed,
        maxLat_allowed,
        minLon_allowed,
        maxLon_allowed,
    ) = get_grid_bounds(land_cover_pth)
    lon_min_inv_domain = np.max([lon_min - buffer_deg, minLon_allowed])
    lon_max_inv_domain = np.min([lon_max + buffer_deg, maxLon_allowed])
    lat_min_inv_domain = np.max([lat_min - buffer_deg, minLat_allowed])
    lat_max_inv_domain = np.min([lat_max + buffer_deg, maxLat_allowed])

    # Subset inversion domain for land cover and hemco diagnostics fields
    lc = lc.isel(lon=lc.lon >= lon_min_inv_domain, lat=lc.lat >= lat_min_inv_domain)
    lc = lc.isel(lon=lc.lon <= lon_max_inv_domain, lat=lc.lat <= lat_max_inv_domain)
    hd = hd.isel(lon=hd.lon >= lon_min_inv_domain, lat=hd.lat >= lat_min_inv_domain)
    hd = hd.isel(lon=hd.lon <= lon_max_inv_domain, lat=hd.lat <= lat_max_inv_domain)

    # Initialize state vector from land cover, replacing all values with NaN (to be filled later)
    statevector = lc.where(lc == -9999.0)

    # Set pixels in buffer areas to 0
    if is_regional:
        statevector[:, (statevector.lon < lon_min) | (statevector.lon > lon_max)] = 0
        statevector[(statevector.lat < lat_min) | (statevector.lat > lat_max), :] = 0

    # Also set pixels over water to 0, unless there are offshore emissions
    if land_threshold:
        # Where there is neither land nor emissions, replace with -9999
        land = lc.where((lc > land_threshold) | (hd > emis_threshold))
        statevector.values[land.isnull().values] = -9999

    # Fill in the remaining NaNs with state vector element values
    statevector.values[statevector.isnull().values] = np.arange(
        1, statevector.isnull().sum() + 1
    )[::-1]

    # Assign buffer pixels (the remaining 0's) to state vector
    # -------------------------------------------------------------------------
    if is_regional:
        buffer_area = statevector.values == 0

        # Get image coordinates of all pixels in buffer area
        irows = np.arange(buffer_area.shape[0])
        icols = np.arange(buffer_area.shape[1])
        irows = np.transpose(np.tile(irows, (len(icols), 1)))
        icols = np.tile(icols, (len(irows), 1))
        irows_good = irows[buffer_area > 0]
        icols_good = icols[buffer_area > 0]
        coords = [[icols_good[j], irows_good[j]] for j in range(len(irows_good))]

        # K-means
        X = np.array(coords)
        kmeans = KMeans(n_clusters=k_buffer_clust, random_state=0).fit(X)

        # Assign pixels to state vector
        # -------------------------------------------------------------------------
        highres_statevector_max = np.nanmax(statevector.values)
        n_rows = statevector.shape[0]
        n_cols = statevector.shape[1]
        for r in range(n_rows):
            for c in range(n_cols):
                if statevector[r, c].values == 0:
                    statevector[r, c] = (
                        kmeans.predict([[c, r]])[0] + 1 + highres_statevector_max
                    )
        # -------------------------------------------------------------------------

    # Make dataset
    da_statevector = statevector.copy()
    ds_statevector = da_statevector.to_dataset(name="StateVector")

    # Add attribute metadata
    ds_statevector.lat.attrs["units"] = "degrees_north"
    ds_statevector.lat.attrs["long_name"] = "Latitude"
    ds_statevector.lon.attrs["units"] = "degrees_east"
    ds_statevector.lon.attrs["long_name"] = "Longitude"
    ds_statevector.StateVector.attrs["units"] = "none"
    ds_statevector.StateVector.attrs["missing_value"] = -9999
    ds_statevector.StateVector.attrs["_FillValue"] = -9999

    # Save
    if save_pth is not None:
        print("Saving file {}".format(save_pth))
        ds_statevector.to_netcdf(
            save_pth,
            encoding={
                v: {"zlib": True, "complevel": 9} for v in ds_statevector.data_vars
            },
        )

    if not is_regional:

        # Postprocessing: Open the state vector file
        state_vector_filepath = './StateVector.nc'
        state_vector = xr.load_dataset(state_vector_filepath)

        # make Greenland/large ice-covered region coordinates NaN
        lc_orig = lc_orig.isel(lon=lc_orig.lon >= lon_min_inv_domain, lat=lc_orig.lat >= lat_min_inv_domain)
        lc_orig = lc_orig.isel(lon=lc_orig.lon <= lon_max_inv_domain, lat=lc_orig.lat <= lat_max_inv_domain)
        mask = lc_orig["FRLANDIC"].values < 0.5
        state_vector["StateVector"].values[~mask[0]] = np.nan

        ds_statevector = state_vector["StateVector"].to_dataset()

        # Add attribute metadata
        ds_statevector.lat.attrs["units"] = "degrees_north"
        ds_statevector.lat.attrs["long_name"] = "Latitude"
        ds_statevector.lon.attrs["units"] = "degrees_east"
        ds_statevector.lon.attrs["long_name"] = "Longitude"
        ds_statevector.StateVector.attrs["units"] = "none"
        ds_statevector.StateVector.attrs["missing_value"] = -9999
        ds_statevector.StateVector.attrs["_FillValue"] = -9999

        # Save
        if save_pth is not None:
            print("Saving file {}".format(save_pth))
            ds_statevector.to_netcdf(
                save_pth,
                encoding={
                    v: {"zlib": True, "complevel": 9} for v in ds_statevector.data_vars
                },
            )

    return ds_statevector


if __name__ == "__main__":
    import sys

    config_path = sys.argv[1]
    land_cover_pth = sys.argv[2]
    hemco_diag_pth = sys.argv[3]
    save_pth = sys.argv[4]

    make_state_vector_file(
        config_path,
        land_cover_pth,
        hemco_diag_pth,
        save_pth,
    )

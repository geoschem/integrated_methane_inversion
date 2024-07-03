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

def cluster_buffer_elements(data, num_clusters, offset):
    """
    Description:
        Cluster all 0 valued elements in dataarray into desired num clusters
    arguments:
        data       [][]dataarray : xarrray sensitivity data
        num_clusters         int : number of labels to assign data to
        offset              bool : offset labels by this integer value
    Returns:       [][]dataarray : labeled data
    """
    
    # if using 0 clusters, return data with buffer pixels
    # infilled with -9999
    if num_clusters == 0:
        return xr.where(data == 0, -9999, data)
    
    # Get the latitude and longitude coordinates as separate arrays
    latitudes = data.coords["lat"].values
    longitudes = data.coords["lon"].values
    
    data_copy = data.copy()

    # Get the sensitivity values as a 1D array
    Z = data.values.flatten()
    # labels shape for later
    # labels = np.zeros(Z.shape)
    valid_indices = np.where(Z == 0)[0] 

    # Flatten the latitude and longitude arrays into a 2D grid
    # only keeping valid indices
    X, Y = np.meshgrid(longitudes, latitudes)
    X = X.flatten()[valid_indices]
    Y = Y.flatten()[valid_indices]
    
    # Stack the X, Y, and Z arrays to create a (n_samples, n_features) array
    features = np.column_stack((X, Y))

    # Cluster the features using KMeans
    # Mini-Batch k-means is much faster, but with less accuracy
    kmeans = KMeans(n_clusters=num_clusters, random_state=0)

    cluster_labels = kmeans.fit_predict(features)

    # fill labels on corresponding valid indices of label array
    # adjust labels by offset + 1
    Z[valid_indices] = cluster_labels + offset + 1

    # reconstruct 2D grid
    # cluster_labels = Z.reshape(data.shape)
    data_copy.values = Z.reshape(data.shape)
    return data_copy


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
    buffer_min_lat = 0
    buffer_min_lon = 0
    
    # set minimum buffer degrees based on resolution
    if config["isRegional"]  == "true":
        if config["Res"] == "4.0x5.0":
            deg_lat, deg_lon = 4.0, 5.0 
        elif config["Res"] == "2.0x2.5":
            deg_lat, deg_lon = 2.0, 2.5
        elif config["Res"] == "0.5x0.625":
            deg_lat, deg_lon = 0.5, 0.625
        elif config["Res"] == "0.25x0.3125":
            deg_lat, deg_lon = 0.25, 0.3125
        buffer_min_lat = deg_lat * 3
        buffer_min_lon = deg_lon * 3

    # set the buffer degrees to the maximum to ensure 
    # that the buffer is at least the minimum (3 extra grid cells on each side)
    buffer_deg_lat = max(buffer_deg, buffer_min_lat)
    buffer_deg_lon = max(buffer_deg, buffer_min_lon)
        

    # Load land cover data and HEMCO diagnostics
    lc = xr.load_dataset(land_cover_pth)
    hd = xr.load_dataset(hemco_diag_pth)

    # Require hemco diags on same global grid as land cover map
    # TODO remove this offset once the HEMCO standalone files 
    # are regenerated with recent bugfix that corrects the offset
    if config["Res"] == "0.25x0.3125":
        hd["lon"] = hd["lon"] - 0.03125
    elif config["Res"] == "0.5x0.625":
        hd["lon"] = hd["lon"] - 0.0625

    # Select / group fields together
    lc = (lc["FRLAKE"] + lc["FRLAND"] + lc["FRLANDIC"]).drop("time").squeeze()
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
    lon_min_inv_domain = np.max([lon_min - buffer_deg_lon, minLon_allowed])
    lon_max_inv_domain = np.min([lon_max + buffer_deg_lon, maxLon_allowed])
    lat_min_inv_domain = np.max([lat_min - buffer_deg_lat, minLat_allowed])
    lat_max_inv_domain = np.min([lat_max + buffer_deg_lat, maxLat_allowed])

    # Subset inversion domain for land cover and hemco diagnostics fields
    lc = lc.isel(lon=lc.lon >= lon_min_inv_domain, lat=lc.lat >= lat_min_inv_domain)
    lc = lc.isel(lon=lc.lon <= lon_max_inv_domain, lat=lc.lat <= lat_max_inv_domain)
    hd = hd.isel(lon=hd.lon >= lon_min_inv_domain, lat=hd.lat >= lat_min_inv_domain)
    hd = hd.isel(lon=hd.lon <= lon_max_inv_domain, lat=hd.lat <= lat_max_inv_domain)

    # Initialize state vector from land cover, replacing all values with NaN (to be filled later)
    if is_regional:
        statevector = lc.where(lc == -9999.0)
    else:
        # if global, use hemco file 
        statevector = hd.where(hd == -9999.0)

    # Set pixels in buffer areas to 0
    if is_regional:
        statevector[:, (statevector.lon < lon_min) | (statevector.lon > lon_max)] = 0
        statevector[(statevector.lat < lat_min) | (statevector.lat > lat_max), :] = 0

    # Also set pixels over water to 0, unless there are offshore emissions
    if land_threshold > 0:
        # Where there is neither land nor emissions, replace with 0
        if is_regional:
            land = lc.where((lc > land_threshold) | (hd > emis_threshold))
        else:
            # handle half-width polar grid boxes for global,
            # global files are same shape but different lat
            # at poles in that case
            if (
                np.not_equal(hd.lat.values, lc.lat.values).any() &
                np.equal(hd.lat.shape, lc.lat.shape).all()
            ):
                land = lc.where((lc.values > land_threshold) | (hd.values > emis_threshold))
            else:
                land = lc.where((lc > land_threshold) | (hd > emis_threshold))

        statevector.values[land.isnull().values] = -9999

    # Fill in the remaining NaNs with state vector element values
    statevector.values[statevector.isnull().values] = np.arange(
        1, statevector.isnull().sum() + 1
    )[::-1]

    # Assign buffer pixels (the remaining 0's) to state vector
    # -------------------------------------------------------------------------
    if is_regional:
        statevector = cluster_buffer_elements(statevector, k_buffer_clust, statevector.max().item())

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

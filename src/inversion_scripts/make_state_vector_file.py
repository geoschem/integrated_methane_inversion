import numpy as np
import xarray as xr
from sklearn.cluster import KMeans


def get_nested_grid_bounds(land_cover_pth):
    '''
    Get the lat/lon bounds of the nested grid window (NA, EU, AS) for the inversion.
    The land cover file path specifies the window.
    '''

    if 'NA' in land_cover_pth:
        if 'GEOSFP' in land_cover_pth:
            minLat_allowed = 9.75; maxLat_allowed = 60
            minLon_allowed = -130; maxLon_allowed = -60
        elif 'MERRA2' in land_cover_pth:
            minLat_allowed = 10; maxLat_allowed = 70
            minLon_allowed = -140; maxLon_allowed = -40
    elif 'EU' in land_cover_pth:
        if 'GEOSFP' in land_cover_pth:
            minLat_allowed = 32.75; maxLat_allowed = 61.25
            minLon_allowed = -15; maxLon_allowed = 40
        elif 'MERRA2' in land_cover_pth:
            minLat_allowed = 30; maxLat_allowed = 70
            minLon_allowed = -30; maxLon_allowed = 50
    elif 'AS' in land_cover_pth:
        if 'GEOSFP' in land_cover_pth:
            minLat_allowed = -11; maxLat_allowed = 55
            minLon_allowed = 60; maxLon_allowed = 150
        elif 'MERRA2' in land_cover_pth:
            minLat_allowed = -11; maxLat_allowed = 55
            minLon_allowed = 60; maxLon_allowed = 150
    else:
        raise ValueError("Land cover file path contains no NA/EU/AS token.")

    return minLat_allowed, maxLat_allowed, minLon_allowed, maxLon_allowed


def check_nested_grid_compatibility(lat_min, lat_max, lon_min, lon_max, land_cover_pth):
    '''
    Check whether input lat/lon bounds are compatible with (contained within) the nested grid window (NA, EU, or AS).
    The land cover file path specifies the window.
    '''

    minLat_allowed, maxLat_allowed, minLon_allowed, maxLon_allowed = get_nested_grid_bounds(land_cover_pth)

    if lat_min < minLat_allowed or lat_max > maxLat_allowed or lon_min < minLon_allowed or lon_max > maxLon_allowed:
        compatible = False
    else:
        compatible = True

    return compatible


def make_state_vector_file(land_cover_pth, save_pth, lat_min, lat_max, lon_min, lon_max, buffer_deg=5, land_threshold=0.25, k_buffer_clust=8):
    '''
    Generates the state vector file for an analytical inversion.
    
    Arguments
        land_cover_pth [str]   : Path to land cover file
        save_pth       [str]   : Where to save the state vector file
        lat_min        [float] : Minimum latitude 
        lat_max        [float] : Maximum latitude 
        lon_min        [float] : Minimum longitude 
        lon_max        [float] : Maximum longitude
        buffer_deg     [float] : Width of k-means buffer area in degrees
        land_threshold [float] : Minimum land fraction to include pixel as a state vector element
        k_buffer_clust [int]   : Number of buffer clusters for k-means

    Returns
        ds_statevector []     : xarray dataset containing state vector field formatted for HEMCO

    Notes
        - Land cover file default can be something like 'GEOSFP.20200101.CN.025x03125.NA.nc'
    '''

    # Load land cover data
    lc = xr.open_dataset(land_cover_pth)

    # Group fields together
    lc = (lc['FRLAKE'] + lc['FRLAND'] + lc['FRLANDIC']).drop('time').squeeze()
    
    # Check compatibility of region of interest with nesting window (NA, EU, AS)
    compatible = check_nested_grid_compatibility(lat_min, lat_max, lon_min, lon_max, land_cover_pth)
    if not compatible:
        raise ValueError("Region of interest not contained within selected nesting REGION (NA, EU, AS; see config.yml).")

    # Define bounds of inversion domain
    minLat_allowed, maxLat_allowed, minLon_allowed, maxLon_allowed = get_nested_grid_bounds(land_cover_pth)
    lon_min_inv_domain = np.max([lon_min-buffer_deg, minLon_allowed])
    lon_max_inv_domain = np.min([lon_max+buffer_deg, maxLon_allowed])
    lat_min_inv_domain = np.max([lat_min-buffer_deg, minLat_allowed])
    lat_max_inv_domain = np.min([lat_max+buffer_deg, maxLat_allowed])

    # Subset inversion domain using land cover file
    lc = lc.isel(lon=lc.lon>=lon_min_inv_domain, lat=lc.lat>=lat_min_inv_domain)
    lc = lc.isel(lon=lc.lon<=lon_max_inv_domain, lat=lc.lat<=lat_max_inv_domain)

    # Replace all values with NaN (to be filled later)
    statevector = lc.where(lc == -9999.)

    # Set pixels in buffer areas to 0
    statevector[:, (statevector.lon < lon_min) | (statevector.lon > lon_max)] = 0
    statevector[(statevector.lat < lat_min) | (statevector.lat > lat_max), :] = 0

    # Also set pixels over water to 0
    if land_threshold:
        # Where there is no land, replace with 0
        land = lc.where(lc > land_threshold)
        statevector.values[land.isnull()] = 0

    # Fill in the remaining NaNs with state vector element values
    statevector.values[statevector.isnull()] = np.arange(1, statevector.isnull().sum()+1)[::-1]

    # Now set pixels over water to missing_value = -9999
    if land_threshold:
        # First, where there is no land, replace with NaN
        statevector = statevector.where(lc > land_threshold)
        # Fill with missing_value = -9999
        statevector.values[statevector.isnull()] = -9999

    # Assign buffer pixels (the remaining 0's) to state vector
    # -------------------------------------------------------------------------
    buffer_area = (statevector.values == 0)

    # Get image coordinates of all pixels in buffer area
    irows = np.arange(buffer_area.shape[0])
    icols = np.arange(buffer_area.shape[1])
    irows = np.transpose(np.tile(irows,(len(icols),1)))
    icols = np.tile(icols,(len(irows),1))
    irows_good = irows[buffer_area > 0]
    icols_good = icols[buffer_area > 0]
    coords = [[icols_good[j], irows_good[j]] for j in range(len(irows_good))]
    
    # K-means
    X = np.array(coords)
    kmeans = KMeans(n_clusters=k_buffer_clust, random_state=0).fit(X)
    
    # Assign pixels to state vector
    highres_statevector_max = np.nanmax(statevector.values)
    n_rows = statevector.shape[0]
    n_cols = statevector.shape[1]
    for r in range(n_rows):
        for c in range(n_cols):
            if statevector[r,c].values == 0:
                statevector[r,c] = kmeans.predict([[c,r]])[0] + 1 + highres_statevector_max
    # -------------------------------------------------------------------------

    # Make dataset
    da_statevector = statevector.copy()
    ds_statevector = da_statevector.to_dataset(name='StateVector')

    # Add attribute metadata
    ds_statevector.lat.attrs['units'] = 'degrees_north'
    ds_statevector.lat.attrs['long_name'] = 'Latitude'
    ds_statevector.lon.attrs['units'] = 'degrees_east'
    ds_statevector.lon.attrs['long_name'] = 'Longitude'
    ds_statevector.StateVector.attrs['units'] = 'none'
    ds_statevector.StateVector.attrs['missing_value'] = -9999
    ds_statevector.StateVector.attrs['_FillValue'] = -9999

    # Save
    if save_pth:
        print("Saving file {}".format(save_pth))
        ds_statevector.to_netcdf(save_pth)
    
    return ds_statevector


if __name__ == '__main__':
    import sys

    land_cover_pth = sys.argv[1]
    save_pth = sys.argv[2]
    lat_min = float(sys.argv[3])
    lat_max = float(sys.argv[4])
    lon_min = float(sys.argv[5])
    lon_max = float(sys.argv[6])
    buffer_deg = float(sys.argv[7])
    land_threshold = float(sys.argv[8])
    k_buffer_clust = int(sys.argv[9])
    
    make_state_vector_file(land_cover_pth, save_pth, lat_min, lat_max, lon_min, lon_max, buffer_deg, land_threshold, k_buffer_clust)

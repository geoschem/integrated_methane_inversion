import datetime
import xarray as xr
import os
import numpy as np

def make_gridded_posterior(posterior_SF_path, clusters_path, save_path):
    '''
    Lu's inversion code outputs the posterior scaling factors as a vector (.nc file)
    HEMCO wants the scaling factors as a gridded product, by latitude/longitude
    This script uses the posterior vector file and the clusters file to generate a gridded
    version of the posterior scaling factors

    Arguments
       posterior_SF_path [str] : path to the posterior scaling factors from an inversion
       clusters_path     [str] : path to the cluster file, from which we will take coords information
       save_path         [str] : path where the gridded posterior should be saved

    '''

    # Load clusters and scaling factor data
    clust = xr.load_dataset(clusters_path)
    scfac = xr.load_dataset(posterior_SF_path)

    # Map the scaling factors to the cluster grid
    nlat = len(clust['lat'])
    nlon = len(clust['lon'])
    scfac_arr = np.zeros(clust['Clusters'].shape)
    for ilat in range(nlat):
        for ilon in range(nlon):
            cluster_id = int(clust['Clusters'].values[ilat,ilon])
            scfac_arr[ilat,ilon] = scfac['xhat'].values[cluster_id-1]

    # Convert to data array
    lat = clust['lat'].values
    lon = clust['lon'].values
    scfac_arr = xr.DataArray(scfac_arr, [("lat", list(lat)), ("lon", list(lon))], attrs={'units': "none"})
    
    # Create dataset
    ds_scfac = xr.Dataset({"SF_Nonwetland": (["lat", "lon"], scfac_arr)},
                          coords={"lon": ("lon", lon), "lat": ("lat", lat)})

    # Add attribute metadata
    ds_scfac.lat.attrs['units'] = 'degrees_north'
    ds_scfac.lat.attrs['long_name'] = 'Latitude'
    ds_scfac.lon.attrs['units'] = 'degrees_east'
    ds_scfac.lon.attrs['long_name'] = 'Longitude'
    ds_scfac.SF_Nonwetland.attrs['units'] = '1'

    # Create netcdf
    ds_scfac.to_netcdf(save_path)



if __name__ == '__main__':
    import sys

    posterior_SF_path = sys.argv[1]
    clusters_path = sys.argv[2]
    save_path = sys.argv[3]

    make_gridded_posterior(posterior_SF_path, clusters_path, save_path) 

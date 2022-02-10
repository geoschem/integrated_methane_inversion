import datetime
import xarray as xr
import os
import numpy as np


def do_gridding(vector, statevector):
    '''
    Project input vector onto the inversion grid using information from the state vector file.
    Input vector should be a numpy array of scale factors (SF), diagonal elements of the posterior 
    error covariance matrix (S_post), or diagonal elements of the averaging kernel matrix (A).
    '''
    
    # Map the input vector (e.g., scale factors) to the state vector grid
    nlat = len(statevector['lat'])
    nlon = len(statevector['lon'])
    target_array = np.empty(statevector['StateVector'].shape)
    target_array[:] = np.nan
    for ilat in range(nlat):
        for ilon in range(nlon):
            element_id = statevector['StateVector'].values[ilat,ilon]
            if ~np.isnan(element_id):
                target_array[ilat,ilon] = vector[int(element_id)-1]

    # Convert to data array
    lat = statevector['lat'].values
    lon = statevector['lon'].values
    target_array = xr.DataArray(target_array, 
                                [("lat", list(lat)), ("lon", list(lon))], 
                                attrs={'units': "none"})

    return target_array


def make_gridded_posterior(posterior_SF_path, state_vector_path, save_path):
    '''
    The IMI code outputs the inversion results as vectors and matrices (in a .nc file).
    We (and HEMCO, for scale factors) want the results as a gridded product, by latitude/longitude.
    This script uses the inversion results file and the state vector file to generate a gridded
    version of the posterior scale factors, posterior errors, and averaging kernel sensitivities.

    Arguments
       posterior_SF_path [str] : path to the posterior scale factors from an inversion
       state_vector_path [str] : path to the state vector file, from which we will take coords information
       save_path         [str] : path where the gridded posterior should be saved

    '''

    # Load state vector and inversion results data
    statevector = xr.load_dataset(state_vector_path)
    inv_results = xr.load_dataset(posterior_SF_path)

    # Get the scale factors and the diagonals of the S_post and A matrices
    SF = inv_results['xhat'].values
    S_post = np.diagonal(inv_results['S_post'].values)
    A = np.diagonal(inv_results['A'].values)

    # Do gridding
    gridded_SF = do_gridding(SF, statevector)
    gridded_S_post = do_gridding(S_post, statevector)
    gridded_A = do_gridding(A, statevector)

    # Create dataset
    lat = gridded_SF['lat'].values
    lon = gridded_SF['lon'].values
    ds = xr.Dataset({'ScaleFactor': (["lat", "lon"], gridded_SF),
                     'S_post': (["lat", "lon"], gridded_S_post),
                     'A': (["lat", "lon"], gridded_A)},
                    coords={"lon": ("lon", lon), "lat": ("lat", lat)})

    # Add attribute metadata
    ds.lat.attrs['units'] = 'degrees_north'
    ds.lat.attrs['long_name'] = 'Latitude'
    ds.lon.attrs['units'] = 'degrees_east'
    ds.lon.attrs['long_name'] = 'Longitude'
    ds.ScaleFactor.attrs['units'] = '1'
    ds.S_post.attrs['units'] = '1'
    ds.A.attrs['units'] = '1'
    
    # Create netcdf
    ds.to_netcdf(save_path)

    print(f'Saved gridded file to {save_path}')


if __name__ == '__main__':
    import sys

    posterior_SF_path = sys.argv[1]
    state_vector_path = sys.argv[2]
    save_path = sys.argv[3]

    make_gridded_posterior(posterior_SF_path, state_vector_path, save_path) 

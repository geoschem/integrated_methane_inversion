import numpy as np
import xarray as xr


def do_gridding(vector, statevector):
    """
    Project input vector onto the inversion grid using information from the state vector file.
    Input vector should be a numpy array of scale factors (SF), diagonal elements of the posterior
    error covariance matrix (S_post), or diagonal elements of the averaging kernel matrix (A).
    """

    # Map the input vector (e.g., scale factors) to the state vector grid
    sv_index = np.nan_to_num(statevector.StateVector.values, nan=0).astype(int)
    outarr = vector[sv_index - 1]
    outarr = np.where(
        np.isnan(statevector.StateVector.values)[...,None],
        np.nan,
        outarr
    )

    # to dataarray    
    if statevector.StateVector.dims == ('time', 'lat', 'lon'):
        target_array = xr.DataArray(
            outarr,
            dims = ('time', 'lat', 'lon', 'ensemble'),
            coords = {'time': statevector.time.values, 'lat': statevector.lat.values, 'lon': statevector.lon.values},
            attrs={"units": "1"}
        )
    elif statevector.StateVector.dims == ('time', 'nf', 'Ydim', 'Xdim'):
        target_array = xr.DataArray(
            outarr,
            dims = ('time', 'nf', 'Ydim', 'Xdim', 'ensemble'),
            coords=dict(time=(['time'], statevector.time.values),
                        lats=(['nf', 'Ydim', 'Xdim'], statevector.lats.values),
                        lons=(['nf', 'Ydim', 'Xdim'], statevector.lons.values)),
            attrs={"units": "1"}
        )
    else:
        print('StateVector is not in GCClassic (lat, lon) dimension, nor in GCHP (nf, Ydim, Xdim) dimension')

    return target_array


def make_gridded_posterior(posterior_SF_path, state_vector_path, save_path):
    """
    The IMI code outputs the inversion results as vectors and matrices (in a .nc file).
    We (and HEMCO, for scale factors) want the results as a gridded product, by latitude/longitude.
    This script uses the inversion results file and the state vector file to generate a gridded
    version of the posterior scale factors, posterior errors, and averaging kernel sensitivities.

    Arguments
       posterior_SF_path [str] : path to the posterior scale factors from an inversion
       state_vector_path [str] : path to the state vector file, from which we will take coords information
       save_path         [str] : path where the gridded posterior should be saved

    """

    # Load state vector and inversion results data
    statevector = xr.load_dataset(state_vector_path)
    inv_results = xr.load_dataset(posterior_SF_path)

    target_data_prefixes = ["xhat", "S_post", "A"]
    data_vars = list(inv_results.data_vars)

    # get all vars that match prefixes
    target_data_vars = [
        var
        for var in data_vars
        if any([var.startswith(prefix) for prefix in target_data_prefixes])
    ]

    # Do the gridding for each variable and store in a dictionary
    data_dict = {}
    for var in target_data_vars:
        attrs = inv_results[var].attrs
        attrs["units"] = "1"
        if var.startswith("A") or var.startswith("S_post"):
            # get the diagonals of the S_post and A matrices
            gridded_data = do_gridding(np.diagonal(inv_results[var].values).transpose(), statevector)
            if statevector.StateVector.dims == ('time', 'lat', 'lon'):
                data_dict[var] = (["time", "lat", "lon", "ensemble"], gridded_data.data, attrs)
            elif statevector.StateVector.dims == ('time', 'nf', 'Ydim', 'Xdim'):
                data_dict[var] = (["time", "nf", "Ydim", "Xdim", "ensemble"], gridded_data.data, attrs)
        elif var.startswith("xhat"):
            # get the scale factors
            # fill nan in SF with 1 to prevent GEOS-Chem error
            gridded_data = do_gridding(inv_results[var].values, statevector).fillna(1)
            # change key to ScaleFactor to match HEMCO expectations
            new_SF_key = f"ScaleFactor{var[len('xhat'):]}"
            if statevector.StateVector.dims == ('time', 'lat', 'lon'):
                data_dict[new_SF_key] = (["lat", "lon", "ensemble"], gridded_data.data, attrs)
            elif statevector.StateVector.dims == ('time', 'nf', 'Ydim', 'Xdim'):
                data_dict[new_SF_key] = (["time", "nf", "Ydim", "Xdim", "ensemble"], gridded_data.data, attrs)

    # Create dataset
    if statevector.StateVector.dims == ('time', 'lat', 'lon'):
        time = statevector["time"].values
        lat = statevector["lat"].values
        lon = statevector["lon"].values

        ds = xr.Dataset(
            data_dict,
            coords={"time": ("time", time), "lon": ("lon", lon), "lat": ("lat", lat)},
        )

        # Add attribute metadata for coordinates
        ds.lat.attrs["units"] = "degrees_north"
        ds.lat.attrs["long_name"] = "Latitude"
        ds.lon.attrs["units"] = "degrees_east"
        ds.lon.attrs["long_name"] = "Longitude"
    elif statevector.StateVector.dims == ('time', 'nf', 'Ydim', 'Xdim'):
        ds = xr.Dataset(
            data_dict,
            coords=dict(time=(['time'], statevector.time.values),
                        lats=(['nf', 'Ydim', 'Xdim'], statevector.lats.values),
                        lons=(['nf', 'Ydim', 'Xdim'], statevector.lons.values)),
        )

        # Add attribute metadata for coordinates
        ds.lats.attrs["units"] = "degrees_north"
        ds.lats.attrs["long_name"] = "Latitude"
        ds.lons.attrs["units"] = "degrees_east"
        ds.lons.attrs["long_name"] = "Longitude"
    ds.attrs = inv_results.attrs

    # Create netcdf for ensemble results
    ds.to_netcdf(
        save_path.replace('.nc', '_ensemble.nc'),
        encoding={v: {"zlib": True, "complevel": 1} for v in ds.data_vars}
    )

    # Create netcdf for default results
    ds.isel(ensemble = ds.attrs['default_member_index']).to_netcdf(
        save_path,
        encoding={v: {"zlib": True, "complevel": 1} for v in ds.data_vars}
    )

    print(f"Saved gridded file to {save_path}")


if __name__ == "__main__":
    import sys

    posterior_SF_path = (sys.argv[1]).replace('.nc', '_ensemble.nc')
    state_vector_path = sys.argv[2]
    save_path = sys.argv[3]

    make_gridded_posterior(posterior_SF_path, state_vector_path, save_path)

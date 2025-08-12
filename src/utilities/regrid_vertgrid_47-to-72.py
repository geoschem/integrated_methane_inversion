# regrid_vertgrid_47-to-72.py

import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

def vertical_interp_47_to_72(data_in, p_edge_47, p_edge_72):
    """
    Interpolates 3D GEOS-Chem data (lev47, lat, lon) from 47 vertical levels to 72 levels.

    Parameters:
    -----------
    data_in : np.ndarray
        Input array with shape (47, lat, lon)
    p_edge_47 : np.ndarray
        Pressure at edges of 47 levels (length 48)
    p_edge_72 : np.ndarray
        Pressure at edges of 72 levels (length 73)

    Returns:
    --------
    data_interp : np.ndarray
        Output array with shape (72, lat, lon)
    """

    # Compute pressure at midpoints
    p_mid_47 = 0.5 * (p_edge_47[:-1] + p_edge_47[1:])
    p_mid_72 = 0.5 * (p_edge_72[:-1] + p_edge_72[1:])

    lev47, nlat, nlon = data_in.shape
    data_reshaped = data_in.reshape(lev47, -1)  # (47, nlat*nlon)

    # Interpolate along pressure dimension
    interp_func = interp1d(p_mid_47, data_reshaped, axis=0,
                           bounds_error=False, fill_value='extrapolate')
    data_interp_reshaped = interp_func(p_mid_72)

    return data_interp_reshaped.reshape(72, nlat, nlon)

def regrid_tropomi_bc_vert(input_path, output_path):
    var = 'SpeciesBC_CH4'
    ds = xr.open_dataset(input_path).isel(time=0).squeeze()

    # Ap unit is Pa, Bp is unitless
    Ap_47edge = np.array([0.000000e+00, 4.804826e-00, 6.593752e+02, 1.313480e+03, 
                1.961311e+03, 2.609201e+03, 3.257081e+03, 3.898201e+03, 
                4.533901e+03, 5.169611e+03, 5.805321e+03, 6.436264e+03, 
                7.062198e+03, 7.883422e+03, 8.909992e+03, 9.936521e+03, 
                1.091817e+04, 1.189586e+04, 1.286959e+04, 1.429100e+04, 
                1.562600e+04, 1.696090e+04, 1.816190e+04, 1.930970e+04, 
                2.032590e+04, 2.121500e+04, 2.187760e+04, 2.238980e+04, 
                2.243630e+04, 2.168650e+04, 2.011920e+04, 1.769300e+04, 
                1.503930e+04, 1.278370e+04, 1.086630e+04, 9.236572e+03, 
                7.851231e+03, 5.638791e+03, 4.017541e+03, 2.836781e+03, 
                1.979160e+03, 9.292942e+02, 4.076571e+02, 1.650790e+02, 
                6.167791e+01, 2.113490e+01, 6.600001e+00, 1.000000e+00])
    Bp_47edge = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 
                9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01, 
                8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01, 
                7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 
                6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01, 
                4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01, 
                2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 
                6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])
    Ap_72edge = np.array([0.000000e+00, 4.804826e+00, 6.593752e+02, 1.313480e+03, 
                1.961311e+03, 2.609201e+03, 3.257081e+03, 3.898201e+03, 
                4.533901e+03, 5.169611e+03, 5.805321e+03, 6.436264e+03, 
                7.062198e+03, 7.883422e+03, 8.909992e+03, 9.936521e+03, 
                1.091817e+04, 1.189586e+04, 1.286959e+04, 1.429100e+04, 
                1.562600e+04, 1.696090e+04, 1.816190e+04, 1.930970e+04, 
                2.032590e+04, 2.121500e+04, 2.187760e+04, 2.238980e+04, 
                2.243630e+04, 2.168650e+04, 2.011920e+04, 1.769300e+04, 
                1.503930e+04, 1.278370e+04, 1.086630e+04, 9.236572e+03, 
                7.851231e+03, 6.660341e+03, 5.638791e+03, 4.764391e+03, 
                4.017541e+03, 3.381001e+03, 2.836781e+03, 2.373041e+03, 
                1.979160e+03, 1.645710e+03, 1.364340e+03, 1.127690e+03, 
                9.292942e+02, 7.619842e+02, 6.216801e+02, 5.046801e+02, 
                4.076571e+02, 3.276431e+02, 2.620211e+02, 2.084970e+02, 
                1.650790e+02, 1.300510e+02, 1.019440e+02, 7.951341e+01, 
                6.167791e+01, 4.758061e+01, 3.650411e+01, 2.785261e+01, 
                2.113490e+01, 1.594950e+01, 1.197030e+01, 8.934502e+00, 
                6.600001e+00, 4.758501e+00, 3.270000e+00, 2.000000e+00, 
                1.000000e+00])
    Bp_72edge = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 
                9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01, 
                8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01, 
                7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 
                6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01, 
                4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01, 
                2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 
                6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 
                0.000000e+00])

    PS = 1e5
    p_edge_47 = Ap_47edge + (Bp_47edge * PS)
    p_edge_72 = Ap_72edge + (Bp_72edge * PS)

    data_in = ds[var].values
    lat = ds['lat'].values
    lon = ds['lon'].values

    data_interp = vertical_interp_47_to_72(data_in, p_edge_47, p_edge_72)

    data_out = xr.DataArray(
        data_interp[None, ...],
        dims=['time', 'lev', 'lat', 'lon'],
        coords=[[0.], np.flip(np.arange(1,73)), lat, lon],
        attrs=ds[var].attrs
    )

    dsout = xr.Dataset({var: data_out})
    dsout['lev'].attrs = dict(axis='Z', long_name="GEOS-Chem level", positive="up", units="level")
    dsout['lat'].attrs = ds['lat'].attrs
    dsout['lon'].attrs = ds['lon'].attrs
    dsout['time'].attrs = {
        "long_name": "Time",
        "axis": "T",
        "units": f"days since {np.datetime_as_string(ds['time'].values, unit='D')}"
    }

    dsout.to_netcdf(output_path, engine='netcdf4', format='NETCDF4', unlimited_dims='time',
                    encoding={var: {"zlib": True, "complevel": 1}})

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python regrid_vertgrid_47-to-72.py input.nc4 output.nc4")
        sys.exit(1)
    regrid_tropomi_bc_vert(sys.argv[1], sys.argv[2])
import numpy as np
import xarray as xr
import datetime
from joblib import Parallel, delayed
from utils import zero_pad_num_hour


def zero_pad_num(n):
    nstr = str(n)
    if len(nstr) == 1:
        nstr = '000'+nstr
    if len(nstr) == 2:
        nstr = '00'+nstr
    if len(nstr) == 3:
        nstr = '0'+nstr
    return nstr


def calc_sensi(nelements, perturbation, startday, endday, run_dirs_pth, run_name, sensi_save_pth):
    '''
    Loops over output data from GEOS-Chem perturbation simulations to compute sensitivities 
    for the Jacobian matrix.

    Arguments
        nelements      [int]   : Number of state vector elements
        perturbation   [float] : Size of perturbation (e.g., 0.5)
        startday       [str]   : First day of inversion period; formatted YYYYMMDD
        endday         [str]   : Last day of inversion period; formatted YYYYMMDD
        run_dirs_pth   [str]   : Path to directory containing GC Jacobian run directories
        run_name       [str]   : Simulation run name; e.g. 'CH4_Jacobian'
        sensi_save_pth [str]   : Path to save the sensitivity data

    Resulting 'Sensi' files look like:

        <xarray.Dataset>
        Dimensions:  (grid: 1207, lat: 105, lev: 47, lon: 87)
        Coordinates:
        * lon      (lon) float64 -107.8 -107.5 -107.2 -106.9 ... -81.56 -81.25 -80.94
        * lat      (lat) float64 10.0 10.25 10.5 10.75 11.0 ... 35.25 35.5 35.75 36.0
        * lev      (lev) int32 1 2 3 4 5 6 7 8 9 10 ... 38 39 40 41 42 43 44 45 46 47
        * grid     (grid) int32 1 2 3 4 5 6 7 8 ... 1201 1202 1203 1204 1205 1206 1207
        Data variables:
            Sensi    (grid, lev, lat, lon) float32 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0

    Pseudocode summary:
    
        for each day:
            load the base run SpeciesConc file
            nlon = count the number of longitudes
            nlat = count the number of latitudes
            nlev = count the number of vertical levels
            for each hour:
                base = extract the base run data for the hour
                Sensi = np.empty((nelements, nlev, nlat, nlon))
                Sensi.fill(np.nan)
                for each state vector element:
                    load the SpeciesConc .nc file for the element and day
                    pert = extract the data for the hour
                    sens = pert - base
                    Sensi[element,:,:,:] = sens
                save Sensi as netcdf with appropriate coordinate variables
    '''

    # Make date range
    days = []
    dt = datetime.datetime.strptime(startday, '%Y%m%d')
    dt_max = datetime.datetime.strptime(endday, '%Y%m%d')
    while dt < dt_max:
        dt_str = str(dt)[0:10].replace('-','')
        days.append(dt_str)
        delta = datetime.timedelta(days=1)
        dt += delta

    # Loop over model data to get sensitivities
    hours = range(24)
    elements = range(nelements)

    # For each day
    for d in days:
        # Load the base run SpeciesConc file
        base_data = xr.load_dataset(f'{run_dirs_pth}/{run_name}_0000/OutputDir/GEOSChem.SpeciesConc.{d}_0000z.nc4')
        # Count nlat, nlon, nlev
        nlon = len(base_data['lon']) # 52
        nlat = len(base_data['lat']) # 61
        nlev = len(base_data['lev']) # 47
        # For each hour
        def process(h):
            # Get the base run data for the hour
            base = base_data['SpeciesConc_CH4'][h,:,:,:]
            # Initialize sensitivities array
            Sensi = np.empty((nelements, nlev, nlat, nlon))
            Sensi.fill(np.nan)
            # For each state vector element
            for e in elements:
                # State vector elements are numbered 1..nelements
                elem = zero_pad_num(e+1)
                # Load the SpeciesConc file for the current element and day
                pert_data = xr.load_dataset(f'{run_dirs_pth}/{run_name}_{elem}/OutputDir/GEOSChem.SpeciesConc.{d}_0000z.nc4')
                # Get the data for the current hour
                pert = pert_data['SpeciesConc_CH4'][h,:,:,:]
                # Compute and store the sensitivities
                sens = (pert.values - base.values)/perturbation
                Sensi[e,:,:,:] = sens
            # Save Sensi as netcdf with appropriate coordinate variables
            Sensi = xr.DataArray(Sensi, 
                                 coords=(np.arange(1,nelements+1), np.arange(1,nlev+1), base.lat, base.lon), 
                                 dims=['element','lev','lat','lon'],
                                 name='Sensitivities')
            Sensi = Sensi.to_dataset()
            Sensi.to_netcdf(f'{sensi_save_pth}/Sensi_{d}_{zero_pad_num_hour(h)}.nc')

        results = Parallel(n_jobs=-1) (delayed(process)(hour) for hour in hours)
    print(f'Saved GEOS-Chem sensitivity files to {sensi_save_pth}')





if __name__ == '__main__':
    import sys

    nelements = int(sys.argv[1])
    perturbation = float(sys.argv[2])
    startday = sys.argv[3]
    endday = sys.argv[4]
    run_dirs_pth = sys.argv[5]
    run_name = sys.argv[6]
    sensi_save_pth = sys.argv[7]

    calc_sensi(nelements, perturbation, startday, endday, run_dirs_pth, run_name, sensi_save_pth)

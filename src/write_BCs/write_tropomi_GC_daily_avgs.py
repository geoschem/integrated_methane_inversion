import glob
import os
import re
import datetime
import numpy as np
from joblib import Parallel, delayed
import xarray as xr
from netCDF4 import Dataset
import sys

import yaml
with open("config_write_BCs.yml", "r") as f:
    config = yaml.safe_load(f)

sys.path.insert(1, "../inversion_scripts/operators/")
sys.path.insert(1, "../inversion_scripts/")
from operator_utilities import nearest_loc
from TROPOMI_operator import apply_tropomi_operator
from utils import save_obj, load_obj

def get_TROPOMI_times(filename):
    
    """
    Function that parses the TROPOMI filenames to get the start and end times.
    Example input (str): S5P_RPRO_L2__CH4____20220725T152751_20220725T170921_24775_03_020400_20230201T100624.nc
    Example output (tuple): (np.datetime64('2022-07-25T15:27:51'), np.datetime64('2022-07-25T17:09:21'))
    """

    file_times = re.search(r'(\d{8}T\d{6})_(\d{8}T\d{6})', filename)
    assert file_times is not None, "check TROPOMI filename - wasn't able to find start and end times in the filename"
    start_TROPOMI_time = np.datetime64(datetime.datetime.strptime(file_times.group(1), "%Y%m%dT%H%M%S"))
    end_TROPOMI_time = np.datetime64(datetime.datetime.strptime(file_times.group(2), "%Y%m%dT%H%M%S"))
    
    return start_TROPOMI_time, end_TROPOMI_time

def apply_tropomi_operator_to_one_tropomi_file(filename):
    
    """
    Run apply_tropomi_operator from src/inversion_scripts/operators/TROPOMI_operator.py for a single TROPOMI file (then saves it to a pkl file)
    Example input (str): S5P_RPRO_L2__CH4____20220725T152751_20220725T170921_24775_03_020400_20230201T100624.nc
    Example output: write the file config["workdir"]/step1/S5P_RPRO_L2__CH4____20220725T152751_20220725T170921_24775_03_020400_20230201T100624_GCtoTROPOMI.pkl
    """
    
    result = apply_tropomi_operator(
        filename = filename,
        n_elements = False, # Not relevant
        gc_startdate = start_time_of_interest,
        gc_enddate = end_time_of_interest,
        xlim = [-180, 180],
        ylim = [-90, 90],
        gc_cache = os.path.join(config["workdir"], "runGCC1402", "OutputDir"),
        build_jacobian = False, # Not relevant
        sensi_cache = False) # Not relevant
    
    save_obj(result, os.path.join(config["workdir"], "step1", os.path.basename(filename).replace(".nc","_GCtoTROPOMI.pkl")))

if __name__ == "__main__":

    # From config file, get the start and end times that we will be writing boundary conditions for (+20 days on the end because of our temporal smoothing)
    start_time_of_interest = np.datetime64(datetime.datetime.strptime(config["startdate"], "%Y%m%d"))
    end_time_of_interest = np.datetime64(datetime.datetime.strptime(config["enddate"], "%Y%m%d")) + np.timedelta64(20, 'D')

    # List of all TROPOMI files that interesct our time period of interest
    TROPOMI_files = sorted([file for file in glob.glob(os.path.join(config["tropomi_cache"], "*.nc"))
                            if (start_time_of_interest <= get_TROPOMI_times(file)[0] <= end_time_of_interest)
                            and (start_time_of_interest <= get_TROPOMI_times(file)[1] <= end_time_of_interest)])

    # Make sure some of the TROPOMI files are after the last day you want BCs for
    assert len([file for file in TROPOMI_files
                if (start_time_of_interest <= get_TROPOMI_times(file)[1] <= end_time_of_interest - np.timedelta64(20, 'D'))]) > 20, \
        "There doesn't seem to be many files past the last day you want BCs. This is a problem because of the +/- 15 day temporal smoothing."


    # Run the function across as many cores as you have
    Parallel(n_jobs=-1)(delayed(apply_tropomi_operator_to_one_tropomi_file)(filename) for filename in TROPOMI_files)

    # Read any of the GEOS-Chem files to get the lat/lon grid
    with xr.open_dataset(glob.glob(os.path.join(config["workdir"], "runGCC1402", "OutputDir", "GEOSChem.SpeciesConc*.nc4"))[0]) as data:
        LON = data["lon"].values
        LAT = data["lat"].values
        
    # List of all days in our time range of interest
    alldates = np.arange(start_time_of_interest, end_time_of_interest + np.timedelta64(1, 'D'), dtype='datetime64[D]')
    alldates = [day.astype(datetime.datetime).strftime("%Y%m%d") for day in alldates]

    # Initialize arrays for regridding
    daily_TROPOMI = np.zeros((len(LON), len(LAT), len(alldates)))
    daily_GC = np.zeros((len(LON), len(LAT), len(alldates)))
    daily_count = np.zeros((len(LON), len(LAT), len(alldates)))

    # List of files we wrote above to get data to regrid from
    files_from_step1 = glob.glob(os.path.join(config["workdir"], "step1", "*.pkl"))
    files_from_step1.sort()

    # Perform regridding
    for file in files_from_step1:
        obs_GC = load_obj(file)["obs_GC"]
        NN = obs_GC.shape[0]
        if NN == 0:
            continue
        
        # For each TROPOMI obsevation, assign it to a GEOS-Chem grid cell
        for iNN in range(NN):
            
            # Which day are we on (this is not perfect right now because orbits can cross from one day to the next...
            # but it is the best we can do right now without changing apply_tropomi_operator)
            file_times = re.search(r'(\d{8}T\d{6})_(\d{8}T\d{6})', file)
            assert file_times is not None, "check TROPOMI filename - wasn't able to find start and end times in the filename"
            date = datetime.datetime.strptime(file_times.group(1), "%Y%m%dT%H%M%S").strftime("%Y%m%d")
            time_ind = alldates.index(date)

            c_TROPOMI, c_GC, lon0, lat0 = obs_GC[iNN, :4]
            ii = nearest_loc(lon0, LON, tolerance=5)
            jj = nearest_loc(lat0, LAT, tolerance=4)
            daily_TROPOMI[ii, jj, time_ind] += c_TROPOMI
            daily_GC[ii, jj, time_ind] += c_GC
            daily_count[ii, jj, time_ind] += 1

    # Normalize by how many observations got assigned to a grid cell to finish the regridding
    daily_count[daily_count == 0] = np.nan
    daily_TROPOMI = daily_TROPOMI / daily_count
    daily_GC = daily_GC / daily_count

    # Change order of dimensions
    regrid_CH4 = np.einsum("ijl->lji", daily_TROPOMI) # (lon, lat, time) -> (time, lat, lon)
    regrid_GC = np.einsum("ijl->lji", daily_GC) # (lon, lat, time) -> (time, lat, lon)

    # Write the netCDF file
    outputname = os.path.join(config["workdir"], "step2", "Daily_CH4.nc")
    with Dataset(outputname, "w", format="NETCDF4_CLASSIC") as dataset:
        
        lat = dataset.createDimension("lat", len(LAT))
        lon = dataset.createDimension("lon", len(LON))
        time = dataset.createDimension("time", len(alldates))
        
        latitudes = dataset.createVariable("lat", "f8", ("lat",))
        longitudes = dataset.createVariable("lon", "f8", ("lon",))
        dates = dataset.createVariable("date", "i", ("time",))
        
        nc_CH4 = dataset.createVariable("CH4", "f8", ("time", "lat", "lon"))
        nc_GC = dataset.createVariable("GC", "f8", ("time", "lat", "lon"))
        
        latitudes[:] = LAT
        longitudes[:] = LON
        dates[:] = alldates
        nc_CH4[:, :, :] = regrid_CH4
        nc_GC[:, :, :] = regrid_GC
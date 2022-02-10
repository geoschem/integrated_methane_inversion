import datetime
import xarray as xr
import os
from joblib import Parallel, delayed


def fill_missing_hour(run_name, run_dirs_pth, prev_run_pth, start_day):
    '''
    This script addresses the fact that output files for the first day of a 
    GEOS-Chem simulation do not include data for the first hour of the day; they
    go from 1-23h, instead of 0-23h. The solution is to combine the final output
    file of a previous simulation (the spinup simulation), which contains data 
    for hour 0 of the last day, with the first output file of the more recent 
    simulation (e.g., a sensitivity simulation). This needs to be done for both 
    SpeciesConc and LevelEdgeDiags, and it needs to be done for every run 
    directory: i.e., for every perturbed-state-vector-element simulation, and
    also for the reference simulation. 
    
    Example:
    A GEOS-Chem perturbation simulation for one month, from 20180501 to 20180601. 
    The SpeciesConc and LevelEdgeDiags files (i.e., the output files) for the
    first day, 20180501, are missing data for hour 0. We merge data from the
    final output files of the spinup simulation, which contain data for only hour
    0 of 20180501, into the latest output files for 20180501. Now the data for
    the month are complete. 

    Arguments
        run_name         [str] : Name for this set of runs
        run_dirs_pth     [str] : Path to the parent folder of the GEOS-Chem run
                                 directory or directories where we want to fill
                                 missing data
        prev_run_pth     [str] : Path to the spinup or posterior run directory
        start_day        [str] : First day of simulation, for which the daily
                                 output is missing data; e.g., "20180501"
    '''

    # List run directories
    contents = os.listdir(run_dirs_pth)
    rundirs = [r for r in contents if run_name in r]

    # Process them
    def process(r):
        # Load hour zero from end of spinup run or previous posterior simulation
        prev_file_SC = f'{prev_run_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4'
        prev_file_LE = f'{prev_run_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4'
        prev_data_SC = xr.load_dataset(prev_file_SC)
        prev_data_LE = xr.load_dataset(prev_file_LE)
        
        # Load output SpeciesConc and LevelEdgeDiags file
        output_file_SC = f'{run_dirs_pth}/{r}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0005z.nc4'
        output_data_SC = xr.load_dataset(output_file_SC)
        if '0000' in r:
            output_file_LE = f'{run_dirs_pth}/{r}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0005z.nc4'
            output_data_LE = xr.load_dataset(output_file_LE)
        
        # Merge output and copied datasets and replace original files that were missing the first hour
        merged_data_SC = xr.merge([output_data_SC, prev_data_SC])
        final_file_SC = f'{run_dirs_pth}/{r}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4'
        merged_data_SC.to_netcdf(final_file_SC)
        if '0000' in r:
            merged_data_LE = xr.merge([output_data_LE, prev_data_LE])
            final_file_LE = f'{run_dirs_pth}/{r}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4'
            merged_data_LE.to_netcdf(final_file_LE)

    results = Parallel(n_jobs=-1)(delayed(process)(run) for run in rundirs)


def fill_missing_hour_posterior(run_dirs_pth, prev_run_pth, start_day):
    # Load hour zero from end of spinup run or previous posterior simulation
    prev_file_SC = f'{prev_run_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4'
    prev_file_LE = f'{prev_run_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4'
    prev_data_SC = xr.load_dataset(prev_file_SC)
    prev_data_LE = xr.load_dataset(prev_file_LE)
        
    # Load output SpeciesConc
    output_file_SC = f'{run_dirs_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0005z.nc4'
    output_data_SC = xr.load_dataset(output_file_SC)
    output_file_LE = f'{run_dirs_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0005z.nc4'
    output_data_LE = xr.load_dataset(output_file_LE)

    # Merge output and copied datasets and replace original files that were missing the first hour
    merged_data_SC = xr.merge([output_data_SC, prev_data_SC])
    final_file_SC = f'{run_dirs_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4'
    merged_data_SC.to_netcdf(final_file_SC)
    merged_data_LE = xr.merge([output_data_LE, prev_data_LE])
    final_file_LE = f'{run_dirs_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4'
    merged_data_LE.to_netcdf(final_file_LE)


if __name__ == '__main__':
    import sys

    run_name = sys.argv[1]
    run_dirs_pth = sys.argv[2]
    prev_run_pth = sys.argv[3]
    start_day = sys.argv[4]
    
    if 'posterior' in run_dirs_pth:
        fill_missing_hour_posterior(run_dirs_pth, prev_run_pth, start_day)
    else:
        fill_missing_hour(run_name, run_dirs_pth, prev_run_pth, start_day)

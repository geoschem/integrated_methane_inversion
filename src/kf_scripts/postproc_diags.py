import datetime
import xarray as xr
import os

def fill_missing_hour(run_dirs_pth, postdir, start_day, first_sim_switch):
    '''
    The purpose of this script is to address the problem that output files for the first day of a GEOS-Chem simulation
    do not include data for the first hour of the day; they go from 1-23h, instead of 0-23h. The solution is to make a
    copy of the final output from each posterior simulation, which contains data for hour 0 of the last day, and then 
    combine the copied file with the first output file of the most recent simulation. This needs to be done for both 
    SpeciesConc and LevelEdgeDiags, and it needs to be done for every run directory: i.e., for every perturbed cluster 
    simulation, and also for every posterior simulation.

    Example. I run GEOS-Chem perturbation simulations for week 1, from 20180501 to 20180508. The SpeciesConc and 
    LevelEdgeDiags files (i.e., the output files) for the first day, 20180501, are missing data for hour 0. I copy data 
    from the final output files of my spinup simulation, which contain data for only hour 0 of 20180501, into my latest 
    output files for 20180501. Now my data for week 1 are complete, so I do the inversion for week 1, and run the 
    posterior simulation for week 1. The final SpeciesConc and LevelEdgeDiags files from the posterior simulation are 
    for 20180508, with data only for hour 0. I copy these files and give them new names "... .Copy. ...nc4". I then run 
    perturbation simulations for week 2, from 20180508 to 20180515. Now the output files for 20180508 are missing data
    for hour 0. I get this data from the 20180508 output file from my last posterior simulation. 

    Arguments
        run_dirs_pth     [str] : Path to the parent folder of the GEOS-Chem run directory or directories where we want to fill missing data
        postdir          [str] : Path to the posterior run directory
        start_day        [str] : First day of simulation, for which the daily output is missing data; e.g., "20180501"
        first_sim_switch [str] : Is this the very first simulation? If so, need to fetch data from spinup run; "True" or "False"
    '''

    # List run directories
    contents = os.listdir(run_dirs_pth)
    rundirs = [r for r in contents if 'CH4_' in r]

    # Process them
    for r in rundirs:
        # Load copied SpeciesConc and LevelEdgeDiags files from posterior simulation
        copied_file_SC = f'{postdir}/OutputDir/GEOSChem.SpeciesConc.Copy.{start_day}_0000z.nc4'
        copied_file_LE = f'{postdir}/OutputDir/GEOSChem.LevelEdgeDiags.Copy.{start_day}_0000z.nc4'
        if first_sim_switch == 'True':
            # If this is the very first simulation, use data from spinup run instead of posterior simulation
            copied_file_SC = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/spinup_permian/CH4_spinup/run_dirs/CH4_spinup_0000/OutputDir/GEOSChem.SpeciesConc.'+start_day+'_0000z.nc4'
            copied_file_LE = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/spinup_permian/CH4_spinup/run_dirs/CH4_spinup_0000/OutputDir/GEOSChem.LevelEdgeDiags.'+start_day+'_0000z.nc4'
        copied_data_SC = xr.load_dataset(copied_file_SC)
        copied_data_LE = xr.load_dataset(copied_file_LE)
        # Load output SpeciesConc and LevelEdgeDiags file
        output_file_SC = f'{run_dirs_pth}/{r}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4'
        output_data_SC = xr.load_dataset(output_file_SC)
        output_file_LE = f'{run_dirs_pth}/{r}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4'
        output_data_LE = xr.load_dataset(output_file_LE)
        # Merge output and copied datasets
        merged_data_SC = xr.merge([output_data_SC, copied_data_SC])
        merged_data_LE = xr.merge([output_data_LE, copied_data_LE])
        # Replace original files that were missing the first hour 
        merged_data_SC.to_netcdf(output_file_SC)
        merged_data_LE.to_netcdf(output_file_LE)


if __name__ == '__main__':
    import sys

    run_dirs_pth = sys.argv[1]
    postdir = sys.argv[2]
    start_day = sys.argv[3]
    first_sim_switch = sys.argv[4]

    fill_missing_hour(run_dirs_pth, postdir, start_day, first_sim_switch) 

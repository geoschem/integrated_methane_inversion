import xarray as xr
import os
from joblib import Parallel, delayed

def search_file(file_path, search_string):
    """
    Search for a string in a file and return True if there is a match.
    
    Args:
        file_path (str): Path to the file.
        search_string (str): String to search for in the file.
    
    Returns:
        bool: True if the string is found, False otherwise.
    """
    with open(file_path, 'r') as file:
        for line in file:
            if search_string in line:
                return True
    return False

def fill_missing_hour(run_name, run_dirs_pth, prev_run_pth, start_day, res):
    """
    This script addresses the fact that output files for the first day of a
    GEOS-Chem simulation do not include data for the first hour of the day; they
    go from 1-23h, instead of 0-23h. The solution is to combine the final output
    file of a previous simulation (e.g., the spinup simulation), which contains
    data for hour 0 of the last day, with the first output file of the more recent
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
        res              [str] : Resolution string; e.g., "0.25x0.3125"
    """

    # Get list of run directories
    contents = os.listdir(run_dirs_pth)
    rundirs = [
        r
        for r in contents
        if run_name in r and os.path.isdir(os.path.join(run_dirs_pth, r))
    ]

    # Determine timestamp in filename
    if "0.25x0.3125" in res:
        timestamp = "0005"
    elif "0.5x0.625" in res:
        timestamp = "0005"
    elif "2.0x2.5" in res:
        timestamp = "0010"
    elif "4.0x5.0" in res:
        timestamp = "0010"
    
    # Process them
    def process(r):
        # keep attributes of data variable when arithmetic operations applied
        xr.set_options(keep_attrs=True)
        
        # Load hour zero from end of spinup run or previous posterior simulation
        prev_file_SC = (
            f"{prev_run_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4"
        )
        prev_file_LE = (
            f"{prev_run_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4"
        )
        prev_data_SC = xr.load_dataset(prev_file_SC)
        prev_data_LE = xr.load_dataset(prev_file_LE)

        # Rename SpeciesConcVV_CH4 for the current state vector element
        num=r[-4:]
        
        # For some perturbation simulation we need to scale down the CH4 concentration to 1 ppb
        # Check if this is one of those simulations by checking if HEMCO_Config.rc 
        # reads a 1ppb restart or BC file
        scale_to_1ppb = search_file(f"{run_dirs_pth}/{r}/HEMCO_Config.rc", "jacobian_1ppb_ics_bcs")
        if scale_to_1ppb:
            prev_data_SC["SpeciesConcVV_CH4"] *= 0.0
            prev_data_SC["SpeciesConcVV_CH4"] += 1e-9
            
        prev_data_SC = prev_data_SC.rename({'SpeciesConcVV_CH4':'SpeciesConcVV_CH4_'+num})
        
        # Load output SpeciesConc and LevelEdgeDiags file
        output_file_SC = (
            f"{run_dirs_pth}/{r}/OutputDir/GEOSChem.SpeciesConc.{start_day}_{timestamp}z.nc4"
        )
        output_data_SC = xr.load_dataset(output_file_SC)
        if "0000" in r or "background" in r:
            output_file_LE = f"{run_dirs_pth}/{r}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_{timestamp}z.nc4"
            output_data_LE = xr.load_dataset(output_file_LE)

        # Merge output and copied datasets and replace original files that were missing the first hour
        merged_data_SC = xr.merge([output_data_SC, prev_data_SC])
        final_file_SC = (
            f"{run_dirs_pth}/{r}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4"
        )
        merged_data_SC.to_netcdf(
            final_file_SC,
            encoding={
                v: {"zlib": True, "complevel": 9} for v in merged_data_SC.data_vars
            },
        )
        if "0000" in r or "background" in r:
            merged_data_LE = xr.merge([output_data_LE, prev_data_LE])
            final_file_LE = f"{run_dirs_pth}/{r}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4"
            merged_data_LE.to_netcdf(
                final_file_LE,
                encoding={
                    v: {"zlib": True, "complevel": 9} for v in merged_data_LE.data_vars
                },
            )

    results = Parallel(n_jobs=-1)(delayed(process)(run) for run in rundirs)


def fill_missing_hour_posterior(run_dirs_pth, prev_run_pth, start_day, res):

    # Determine timestamp in filename
    if "0.25x0.3125" in res:
        timestamp = "0005"
    elif "0.5x0.625" in res:
        timestamp = "0005"
    elif "2.0x2.5" in res:
        timestamp = "0010"
    elif "4.0x5.0" in res:
        timestamp = "0010"
    
    # Load hour zero from end of spinup run or previous posterior simulation
    prev_file_SC = (
        f"{prev_run_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4"
    )
    prev_file_LE = (
        f"{prev_run_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4"
    )
    prev_data_SC = xr.load_dataset(prev_file_SC)
    prev_data_LE = xr.load_dataset(prev_file_LE)

    # Load output SpeciesConc
    output_file_SC = (
        f"{run_dirs_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_{timestamp}z.nc4"
    )
    output_data_SC = xr.load_dataset(output_file_SC)
    output_file_LE = (
        f"{run_dirs_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_{timestamp}z.nc4"
    )
    output_data_LE = xr.load_dataset(output_file_LE)

    # Merge output and copied datasets and replace original files that were missing the first hour
    merged_data_SC = xr.merge([output_data_SC, prev_data_SC])
    final_file_SC = (
        f"{run_dirs_pth}/OutputDir/GEOSChem.SpeciesConc.{start_day}_0000z.nc4"
    )
    merged_data_SC.to_netcdf(
        final_file_SC,
        encoding={v: {"zlib": True, "complevel": 9} for v in merged_data_SC.data_vars},
    )
    merged_data_LE = xr.merge([output_data_LE, prev_data_LE])
    final_file_LE = (
        f"{run_dirs_pth}/OutputDir/GEOSChem.LevelEdgeDiags.{start_day}_0000z.nc4"
    )
    merged_data_LE.to_netcdf(
        final_file_LE,
        encoding={v: {"zlib": True, "complevel": 9} for v in merged_data_LE.data_vars},
    )


if __name__ == "__main__":
    import sys

    run_name = sys.argv[1]
    run_dirs_pth = sys.argv[2]
    prev_run_pth = sys.argv[3]
    start_day = sys.argv[4]
    res = sys.argv[5]

    # Check if this is a posterior run, background run, or prior run
    accepted_rundirs = ("posterior_run", f"{run_name}_0000", f"{run_name}_background")
    
    if run_dirs_pth.endswith((accepted_rundirs)):
        fill_missing_hour_posterior(run_dirs_pth, prev_run_pth, start_day, res)
    else:
        fill_missing_hour(run_name, run_dirs_pth, prev_run_pth, start_day, res)

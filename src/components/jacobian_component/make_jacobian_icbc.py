import sys
import xarray as xr

def make_jacobian_icbc(original_file_path, new_file_path, is_BC):
    """
    This function takes a restart or boundary condition file and 
    sets the CH4 concentration to 1 ppb for use in the Jacobian
    simulations.
    Arguments
        original_file_path [str]  : original restart/bc file path
        new_file_path      [str]  : new restart/bc file path
        is_BC              [bool] : True if boundary condition file, 
                                    False if restart file
    Returns
        months_list   [list] : prefix dirs of all months
        days_list     [list] : acceptable file path prefixes
    """
    
    # determine which key to change    
    if is_BC:
        key = "SpeciesBC_CH4"
    else:
        key = "SpeciesRst_CH4"
        
    # read in the original restart/bc file
    orig = xr.load_dataset(original_file_path)
    new_restart = orig.deep_copy()
    
    # set all values to 1 ppb
    new_restart[key] *= 0.0
    new_restart[key] += 1e-9
    
    # write to new file path
    new_restart.to_netcdf(new_file_path)


if __name__ == "__main__":
    original_file_path = sys.argv[1]
    new_file_path = sys.argv[2]
    is_BC = sys.argv[3].lower() == "true"

    make_jacobian_icbc(original_file_path, new_file_path, is_BC)
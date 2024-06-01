import os
import sys
import glob
import xarray as xr

def check_path_and_get_file(path, pattern="*"):
    """
    Check if the path is a file or directory and return the first file
    Arguments
        path      [str]  : path to check
        pattern   [str]  : pattern to match files in the directory
    """
    if os.path.isdir(path):
        print(f"The path '{path}' is a directory.")
        # Find all files matching the pattern in the directory
        matching_files = glob.glob(os.path.join(path, pattern))
        if matching_files:
            # Return the first matching file
            first_file = matching_files[0]
            print(f"First file matching pattern '{pattern}': {first_file}")
            return first_file
        else:
            raise FileNotFoundError(f"No files matching pattern '{pattern}' found in the directory.")
    elif os.path.isfile(path):
        print(f"The path '{path}' is a file.")
        return path
    else:
        print(f"The path '{path}' is neither a file nor a directory.")
        raise FileNotFoundError(f"The path '{path}' is neither a file nor a directory.")

def make_jacobian_icbc(original_file_path, new_file_path):
    """
    This function takes a restart or boundary condition file and 
    sets the CH4 concentration to 1 ppb for use in the Jacobian
    simulations.
    Arguments
        original_file_path [str]  : original restart/bc file path
        new_file_path      [str]  : new restart/bc file path
    Returns
        months_list   [list] : prefix dirs of all months
        days_list     [list] : acceptable file path prefixes
    """
    # read in the original restart/bc file
    orig = xr.load_dataset(original_file_path)
    new_restart = orig.deep_copy()
    
    # determine which data variable to change
    data_vars = list(orig.data_vars)
    if "SpeciesBC_CH4" in data_vars:
        key = "SpeciesBC_CH4"
    elif "SpeciesRst_CH4" in data_vars:
        key = "SpeciesRst_CH4"
    else:
        raise ValueError("No recognized CH4 species found in the file.") 
    
    # set all values to 1 ppb
    new_restart[key] *= 0.0
    new_restart[key] += 1e-9
    
    # write to new file path
    new_restart.to_netcdf(new_file_path)


if __name__ == "__main__":
    original_file_path = sys.argv[1]
    new_file_path = sys.argv[2]
    
    # default to getting the first file in the directory
    # or the file itself if it is a file
    file_path = check_path_and_get_file(original_file_path)

    make_jacobian_icbc(file_path, new_file_path)
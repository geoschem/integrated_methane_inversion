import os
import sys
import glob
import xarray as xr
from src.inversion_scripts.utils import (
    mixing_ratio_conv_factor,
)

def check_path_and_get_file(path, pattern="*"):
    """
    Check if the path is a file or directory and return the first file
    Arguments
        path      [str]  : path to check
        pattern   [str]  : pattern to match files in the directory
    """
    if os.path.isdir(path):
        # Find all files matching the pattern in the directory
        matching_files = glob.glob(os.path.join(path, pattern))
        if matching_files:
            # Return the first matching file
            first_file = matching_files[0]
            return first_file
        else:
            raise FileNotFoundError(f"No files matching pattern '{pattern}' found in the directory.")
    elif os.path.isfile(path):
        return path
    else:
        raise FileNotFoundError(f"The path '{path}' is neither a file nor a directory.")

def make_jacobian_icbc(original_file_path, new_file_path, file_date, species):
    """
    This function takes a restart or boundary condition file and 
    sets the species concentration to 1 mixing ratio unit for use in the 
    Jacobian simulations.
    Arguments
        original_file_path [str]  : original restart/bc file path
        new_file_path      [str]  : new restart/bc file path
    """
    # keep attributes of data variable when arithmetic operations applied
    xr.set_options(keep_attrs=True)
    
    # read in the original restart/bc file
    orig = xr.load_dataset(original_file_path)
    new_restart = orig.copy()
    
    # determine which data variable to change
    data_vars = list(orig.data_vars)
    if f"SpeciesBC_{species}" in data_vars:
        key = f"SpeciesBC_{species}"
        file_prefix = "GEOSChem.BoundaryConditions.lowbg."
    elif f"SpeciesRst_{species}" in data_vars:
        key = f"SpeciesRst_{species}"
        file_prefix = f"GEOSChem.Restart.lowbg."
    else:
        raise ValueError(f"No recognized {species} species found in the file.") 
    
    # set all values to 1 mixing ratio unit, depending on the species
    new_restart[key] *= 0.0
    new_restart[key] += 1/mixing_ratio_conv_factor(species)
        
    write_path = os.path.join(new_file_path, f"{file_prefix}{file_date}_0000z.nc4")
    
    # write to new file path
    new_restart.to_netcdf(write_path)


if __name__ == "__main__":
    original_file_path = sys.argv[1]
    new_file_path = sys.argv[2]
    file_date = sys.argv[3]
    
    # default to getting the first file in the directory
    # or the file itself if it is a file
    file_path = check_path_and_get_file(original_file_path)

    make_jacobian_icbc(file_path, new_file_path, file_date)
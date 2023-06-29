import xarray as xr
import numpy as np
import os
import sys


def make_unit_sf(state_vector_path, save_directory):
    """
    Make an xarray Dataset containing unit emission scale factors.

    Arguments
        state_vector_path  [str] : path to state vector file
        save_directory     [str] : directory in which to save the unit scale factors
    """

    # Make Dataset
    state_vector = xr.load_dataset(state_vector_path)
    scale_factor = state_vector.rename({"StateVector": "ScaleFactor"})
    scale_factor["ScaleFactor"][:] = np.ones(scale_factor["ScaleFactor"].values.shape)
    scale_factor["ScaleFactor"].attrs["units"] = "1"

    # Save to netcdf
    save_path = os.path.join(save_directory, "unit_sf.nc")
    scale_factor.to_netcdf(
        save_path,
        encoding={v: {"zlib": True, "complevel": 9} for v in scale_factor.data_vars},
    )


if __name__ == "__main__":
    # Inputs
    state_vector_path = sys.argv[1]
    save_directory = sys.argv[2]

    # Run the script
    make_unit_sf(state_vector_path, save_directory)

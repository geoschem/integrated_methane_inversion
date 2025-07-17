import xarray as xr
import numpy as np
import os
import sys
import yaml


def make_random_sf(state_vector_path, save_directory, config):
    """
    Make an xarray Dataset containing random emission scale factors.

    Arguments
        state_vector_path  [str] : path to state vector file
        save_directory     [str] : directory in which to save the unit scale factors
    """
    # Make Dataset
    state_vector = xr.load_dataset(state_vector_path)
    state_vector_labels = state_vector["StateVector"]
    mask = state_vector_labels <= last_ROI_element
    scale_factor = state_vector.rename({"StateVector": "ScaleFactor"})
    scale_factor["ScaleFactor"].attrs["units"] = "1"
    last_ROI_element = int(np.nanmax(state_vector.values) - config["nBufferClusters"])
    
    # create random scale factors using settings defined in the config
    np.random.seed(0)
    if config["LognormalErrors"]:
        scale_factor["ScaleFactor"].values = np.random.lognormal(
            mean=1.0,
            sigma=config["EmisRandomPerturbation"],
        )
    else:  # gaussian errors
        scale_factor["ScaleFactor"].values = np.random.uniform(
            1 - config["EmisRandomPerturbation"],
            1 + config["EmisRandomPerturbation"],
            size=scale_factor["ScaleFactor"].values.shape,
        )
        
    # set non-ROI elements to 1.0
    scale_factor["ScaleFactor"].values[~mask] = 1.0

    # Save to netcdf
    save_path = os.path.join(save_directory, "ScaleFactors.nc")
    scale_factor.to_netcdf(
        save_path,
        encoding={v: {"zlib": True, "complevel": 1} for v in scale_factor.data_vars},
    )


if __name__ == "__main__":
    # Inputs
    state_vector_path = sys.argv[1]
    save_directory = sys.argv[2]
    config_file = sys.argv[3]

    # Load configuration
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Run the script
    make_random_sf(state_vector_path, save_directory, config)

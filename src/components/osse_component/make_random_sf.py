import xarray as xr
import numpy as np
import os
import sys
import yaml


def make_random_sf(state_vector_path, config):
    """
    Make an xarray Dataset containing random emission scale factors.

    Arguments
        state_vector_path  [str] : path to state vector file
        config             [str] : config dictionary
    """
    # Make Dataset
    state_vector = xr.load_dataset(state_vector_path)
    state_vector_labels = state_vector["StateVector"]

    last_ROI_element = int(np.nanmax(state_vector_labels.values) - config["nBufferClusters"])
    mask = state_vector_labels <= last_ROI_element
    scale_factor = state_vector.rename({"StateVector": "ScaleFactor"})
    scale_factor["ScaleFactor"].attrs["units"] = "1"
    
    # create random scale factors using settings defined in the config
    np.random.seed(0)
    if config["LognormalErrors"]:
        scale_factor["ScaleFactor"].values = np.random.lognormal(
            mean=1.0,
            sigma=np.log(config["EmisPerturbationOSSE"]),
            size=scale_factor["ScaleFactor"].values.shape,
            ) # Note: set config["EmisPerturbationOSSE"] as GSD value
    else:  
        # gaussian errors
        # Warning: can have negative values
        scale_factor["ScaleFactor"].values = np.random.normal(
            loc=1.0, 
            scale=config["EmisPerturbationOSSE"],
            size=scale_factor["ScaleFactor"].values.shape)
        
    # set non-ROI elements to 1.0
    scale_factor["ScaleFactor"].values[~mask] = 1.0

    # Save to netcdf
    scale_factor.to_netcdf(
        "ScaleFactors.nc",
        encoding={v: {"zlib": True, "complevel": 1} for v in scale_factor.data_vars},
    )

if __name__ == "__main__":
    # Inputs
    state_vector_path = sys.argv[1]
    config_file = sys.argv[2]

    # Load configuration
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Run the script
    if config["CreateAutomaticScaleFactorFileOSSE"]:
        make_random_sf(state_vector_path, config)


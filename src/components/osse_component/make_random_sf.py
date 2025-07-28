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
            sigma=np.log(config["EmisRandomPerturbation"]),
            size=scale_factor["ScaleFactor"].values.shape,
        ) # config["EmisRandomPerturbation"] is GSD value
    else:  # gaussian errors
        # scale_factor["ScaleFactor"].values = np.random.uniform(
        #     1 - config["EmisRandomPerturbation"],
        #     1 + config["EmisRandomPerturbation"],
        #     size=scale_factor["ScaleFactor"].values.shape,
        # )
        scale_factor["ScaleFactor"].values = np.random.normal(
            loc=1.0, 
            scale=config["EmisRandomPerturbation"], 
            size=scale_factor["ScaleFactor"].values.shape) # Warning: has negative values
        
    # set non-ROI elements to 1.0
    scale_factor["ScaleFactor"].values[~mask] = 1.0

    # Save to netcdf
    scale_factor.to_netcdf(
        "ScaleFactors.nc",
        encoding={v: {"zlib": True, "complevel": 1} for v in scale_factor.data_vars},
    )


def make_jacobian_pert_sf(state_vector_path, config):
    """
    Make an array containing perturbation factors for the true emissions.

    Arguments
        state_vector_path  [str] : path to state vector file
        config             [str] : config dictionary
    """

    # StateVector Dataset
    state_vector = xr.load_dataset(state_vector_path)
    # ScaleFactor Dataset
    scale_factor = xr.load_dataset("ScaleFactors.nc")

    # Make Jacobian scale factor file
    prior_sf = xr.Dataset({
        "ScaleFactor": scale_factor["ScaleFactor"],
        "StateVector": state_vector["StateVector"]
    }).to_dataframe().reset_index()
    prior_sf = prior_sf[prior_sf["StateVector"] > 0].sort_values(by="StateVector")
    prior_sf = prior_sf.groupby("StateVector").first()["ScaleFactor"].to_numpy()

    jacobian_pert_sf = 1. / prior_sf
    np.save(f"jacobian_scale_factors.npy", jacobian_pert_sf)  

if __name__ == "__main__":
    # Inputs
    state_vector_path = sys.argv[1]
    config_file = sys.argv[2]

    # Load configuration
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Run the script
    if config["CreateAutomaticScaleFactorFile"]:
        make_random_sf(state_vector_path, config)

    make_jacobian_pert_sf(state_vector_path, config)


"""
Description: Calculate the perturbations for each state vector element and update the perturbation files.
Usage: python make_perturbation_sf.py <config.yml> <period_number> <perturb_value>
Note: period_number is the kalman filter period number for which to calculate the perturbations.
      If not using a kalman filter, set period_number to 1.
"""

import os
import sys
import glob
import yaml
import xarray as xr
import numpy as np
import pandas as pd
from src.inversion_scripts.utils import get_mean_emissions, get_period_mean_emissions


def update_jacobian_perturbation_files(jacobian_dir, state_vector_labels, flat_sf):
    """
    Update the perturbation files in the jacobian run directories with the new scale factors.
    """
    contents = os.listdir(jacobian_dir)
    jacobian_rundirs = [
        os.path.join(jacobian_dir, r)
        for r in contents
        if os.path.isdir(os.path.join(jacobian_dir, r))
    ]

    # loop through each jacobian run directory and update the perturbation file
    # with the new scale factors
    for jacobian_rundir in jacobian_rundirs:
        dir_suffix = jacobian_rundir.split("/")[-1].split("_")[-1]
        perturbation_files = glob.glob(f"{jacobian_rundir}/Perturbations_*.txt")
        for perturbation_file in perturbation_files:

            # check if perturbation file exists
            if os.path.exists(perturbation_file):
                with open(perturbation_file, "r") as file:
                    lines = file.readlines()

                # infer element number based on the file name
                sv_label = perturbation_file[-8:-4]
                sv_element = int(sv_label)
                sv_idx = sv_element - 1

                # make sure we only apply scale factors to emission elements
                if sv_element <= int(state_vector_labels.max().item()):
                    # add the right amount of padding
                    padding = "".ljust(4 - len(str(sv_element)))

                    # construct new perturbation line
                    new_pert_line = (
                        f"ELEM_{sv_label}  {sv_element}  {padding}{flat_sf[sv_idx]}"
                    )

                    # search through perturbations file for element
                    # and replace with new perturbation line
                    for idx, line in enumerate(lines):
                        if line.startswith(f"ELEM_{sv_label}"):
                            lines[idx] = new_pert_line + "\n"
                            break

                    # write out new perturbation file
                    with open(perturbation_file, "w") as file:
                        file.writelines(lines)


def calculate_sfs(state_vector, emis_prior, target_emission=1e-8, prior_sf=None):
    """
    Calculate the scale factors to perturb each state vector
    element by based on the target_emission. Return a dictionary containing
    two flat numpy arrays of perturbation scale factors indexed by state vector
    element. The first array contains the perturbation scale factors to apply
    to the state vector elements in jacobian simulations (based on the
    original prior emissions). The second array contains the effective scale
    factors (based on the nudged prior emissions for kalman mode) used in the
    inversion to calculate the sensitivity of observations to the perturbation.
    For a standalone inversion, the effective scale factors == jacobian scale factors.
    """
    # create a sf dataset with the same structure as the state vector
    sf = state_vector.copy()
    sf = sf.rename({"StateVector": "ScaleFactor"})

    # Calculate scale factors such that applying them to the original emissions
    # will result in a target_emission kg/m2/s2 emission.
    sf["ScaleFactor"] = target_emission / emis_prior["EmisCH4_Total_ExclSoilAbs"]

    # Extract state vector labels
    state_vector_labels = state_vector["StateVector"]

    # Add state vector labels to sf Dataset
    sf = sf.assign(StateVector=state_vector_labels)

    # Use groupby to find the median scale factor for each state vector label
    max_sf_da = sf["ScaleFactor"].groupby(sf["StateVector"]).median()

    # get flat, numpy array by converting to dataframe and sorting based on
    # state vector element
    max_sf_df = max_sf_da.to_dataframe().reset_index()
    max_sf_df = max_sf_df[max_sf_df["StateVector"] > 0].sort_values(by="StateVector")
    jacobian_pert_sf = max_sf_df["ScaleFactor"].values

    # Replace any values greater than the threshold to avoid issues
    # with reaching infinity
    # Replace any NaN values to 1.0
    max_sf_threshold = 15000000.0
    jacobian_pert_sf[jacobian_pert_sf > max_sf_threshold] = max_sf_threshold
    jacobian_pert_sf = np.nan_to_num(jacobian_pert_sf, nan=1.0)

    # If we are using a kalman filter and have nudged prior emissions,
    # calculate the effective scale factors based on the nudged prior emissions
    # these will be used later in the inversion to calculate the sensitivity of
    # observations to the perturbation
    if prior_sf is not None:
        prior_sf = prior_sf.assign(StateVector=state_vector_labels)
        prior_sf_da = prior_sf["ScaleFactor"].groupby(prior_sf["StateVector"]).median()
        prior_sf_df = prior_sf_da.to_dataframe().reset_index()
        prior_sf_df = prior_sf_df[prior_sf_df["StateVector"] > 0].sort_values(
            by="StateVector"
        )
        flat_prior_sf = prior_sf_df["ScaleFactor"].values
    else:
        flat_prior_sf = np.ones(len(jacobian_pert_sf))

    # calculate the effective scale factors
    effective_pert_sf = jacobian_pert_sf / flat_prior_sf

    # return dictionary of scale factor arrays
    perturbation_dict = {
        "effective_pert_sf": effective_pert_sf,
        "jacobian_pert_sf": jacobian_pert_sf,
        "target_emission": target_emission,
    }

    return perturbation_dict


def make_perturbation_sf(config, period_number, perturb_value=1e-8):
    """
    Calculate the perturbations for each state vector element and update the perturbation files.
    Write out an archive of the flat scale factors for later use in sensitivity calculations.
    Args:
        config            [dict] : dictionary of configuration settings
        period_number      [int] : the period number for which to calculate the perturbations
        perturb_value    [float] : the target emission value to perturb each state vector element by
    """
    # make base directory if not already present
    base_directory = os.path.expandvars(
        os.path.join(config["OutputPath"], config["RunName"])
    )

    # get start and end dates
    start_date = str(config["StartDate"])
    end_date = str(config["EndDate"])

    # jacobian rundir path
    jacobian_dir = os.path.join(base_directory, "jacobian_runs")

    # find the hemco emissions file for the period
    prior_cache = os.path.join(base_directory, "hemco_prior_emis/OutputDir")
    if config["KalmanMode"]:
        hemco_emis = get_period_mean_emissions(
            prior_cache, period_number, os.path.join(base_directory, "periods.csv")
        )
        prior_sf = xr.load_dataset(
            os.path.join(
                base_directory, f"archive_sf/prior_sf_period{period_number}.nc"
            )
        )
    else:
        hemco_emis = get_mean_emissions(start_date, end_date, prior_cache)
        prior_sf = None

    # load the state vector dataset
    state_vector = xr.load_dataset(os.path.join(base_directory, "StateVector.nc"))

    # calculate the scale factors to perturb each state vector element by
    perturbation_dict = calculate_sfs(state_vector, hemco_emis, perturb_value, prior_sf)

    # update jacobian perturbation files with new scale factors
    # before we run the jacobian simulations
    update_jacobian_perturbation_files(
        jacobian_dir, state_vector["StateVector"], perturbation_dict["jacobian_pert_sf"]
    )

    # archive npz file of scale factor dictionary for later calculation of sensitivity
    archive_dir = os.path.join(base_directory, "archive_perturbation_sfs")
    os.makedirs(archive_dir, exist_ok=True)
    np.savez(
        os.path.join(archive_dir, f"pert_sf_{period_number}.npz"), **perturbation_dict
    )


if __name__ == "__main__":
    config_path = sys.argv[1]
    period_number = int(sys.argv[2])
    perturb_value = float(sys.argv[3])

    # use appropriate units for perturbation value
    perturb_value = perturb_value * 1e-8

    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    make_perturbation_sf(config, period_number, perturb_value)

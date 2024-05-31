"""
Description: Calculate the perturbations for each state vector element and update the perturbation files.
Usage: python make_perturbation_sf.py <config.yml> <period_number>
Note: period_number is the kalman filter period number for which to calculate the perturbations.
      If not using a kalman filter, set period_number to 1.
"""
import os
import sys
import yaml
import xarray as xr
import numpy as np


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
        perturbation_file = os.path.join(
            jacobian_rundir, f"Perturbations_{dir_suffix}.txt"
        )
        for i in range(state_vector_labels.max().item()):
            # element to adjust pertubation for
            sv_element = int(i + 1)
            sv_label = str(sv_element).zfill(4)

            # construct new perturbation line
            new_pert_line = f"ELEM_{sv_label}  {sv_element}     {flat_sf[i]}"

            # search through perturbations file for element
            # and replace with new perturbation line
            with open(perturbation_file, "r") as file:
                lines = file.readlines()
            for idx, line in enumerate(lines):
                if line.startswith(f"ELEM_{sv_label}"):
                    lines[idx] = new_pert_line + "\n"
                    break

            # write out new perturbation file
            with open(perturbation_file, "w") as file:
                file.write(lines)


def calculate_sfs(base_directory, period_number, state_vector, target_emission=10e-8):
    """
    Calculate the scale factors to perturb each state vector 
    element by based on the target_emission. Return a flat 
    numpy array of the scale factors indexed by state vector 
    element.
    """
    # TODO: set target_emission in configfile
    # load the prior emissions dataset
    # TODO: check if emis_prior already has the soil sink removed
    emis_prior = xr.open_dataset(
        os.path.join(base_directory, f"prior_sf_period{period_number}.nc")
    )

    # create a sf dataset with the same structure as the state vector
    sf = state_vector.copy()
    sf = sf.rename({"StateVector": "ScaleFactor"})

    # Calculate scale factors such that applying them to the original emissions
    # will result in a target_emission kg/m2/s2 emission.
    sf["ScaleFactor"] = target_emission / emis_prior["EmisCH4_Total"]

    # Extract state vector labels
    state_vector_labels = state_vector["StateVector"]

    # Add state vector labels to sf Dataset
    sf = sf.assign(StateVector=state_vector_labels)

    # Use groupby to find the maximum scale factor for each state vector label
    max_sf_da = sf["ScaleFactor"].groupby(sf["StateVector"]).max()

    # get flat, numpy array by converting to dataframe and sorting based on
    # state vector element
    max_sf_df = max_sf_da.to_dataframe().reset_index()
    max_sf_df = max_sf_df[max_sf_df["StateVector"] > 0].sort_values(by="StateVector")
    flat_sf = max_sf_df["ScaleFactor"].values

    return flat_sf


def make_perturbation_sf(config, period_number):
    """
    Calculate the perturbations for each state vector element and update the perturbation files.
    Write out an archive of the flat scale factors for later use in sensitivity calculations.
    """
    # make base directory if not already present
    base_directory = os.path.expandvars(os.path.join(config["OutputPath"], config["RunName"]))

    # load the state vector dataset
    state_vector = xr.load_dataset(os.path.join(base_directory, "StateVector.nc"))

    # calculate the scale factors to perturb each state vector element by
    flat_sf = calculate_sfs(base_directory, period_number, state_vector)

    # jacobian rundir path
    jacobian_dir = os.path.join(base_directory, "jacobian_runs")

    # update jacobian perturbation files with new scale factors 
    # before we run the jacobian simulations
    update_jacobian_perturbation_files(jacobian_dir, state_vector["StateVector"], flat_sf)

    ########################################
    # TODO: cleanup the following code
    # # Calculate scale factors such that applying them to the original emissions
    # # will result in a 20 kg/m2/s2 emission.
    # sf["ScaleFactor"] = 20 / emis_original["EmisCH4_Total"]

    # # loop through the range of state vector labels and set the scale factors
    # state_vector_labels = state_vector["StateVector"].values
    # sf_vals = sf["ScaleFactor"].values
    # flat_sf = []
    # for i in range(state_vector_labels.max().item()):
    #     # find the maximum scale factor for each state vector label
    #     idx = np.where(state_vector_labels == float(i + 1))

    #     # set all scale factors for that state vector label to the maximum
    #     max_sf = sf_vals[idx].max()
    #     sf_vals[idx] = max_sf
    #     # keep track of the maximum scale factor for each state vector label
    #     flat_sf.append(max_sf)

    # # update dataset with new scale factors
    # sf["ScaleFactor"].values = sf_vals

    # archive npy file of flat scale factors for later calculation of sensitivity
    archive_dir = os.path.join(base_directory, "archive_perturbation_sfs")
    os.makedirs(archive_dir, exist_ok=True)
    np.save(
        os.path.join(archive_dir, f"flat_pert_sf_{period_number}.npy"),
        np.array(flat_sf),
    )

    # archive gridded scale factor dataset for use in jacobian simulations
    # archive_path = os.path.join(
    #     archive_dir, f"gridded_pert_sf_{period_number}.nc"
    # )
    # sf.to_netcdf(
    #     archive_path,
    #     encoding={v: {"zlib": True, "complevel": 9} for v in sf.data_vars},
    # )


if __name__ == "__main__":
    config_path = sys.argv[1]
    period_number = sys.argv[2]

    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    make_perturbation_sf(config, period_number)

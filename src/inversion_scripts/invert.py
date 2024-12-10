#!/usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import yaml
import numpy as np
import xarray as xr
from itertools import product
from netCDF4 import Dataset
from src.inversion_scripts.utils import load_obj, calculate_superobservation_error


def ensure_float_list(variable):
    """Make sure the variable is a list of floats."""
    if isinstance(variable, list):
        # Convert each item in the list to a float
        return [float(item) for item in variable]
    elif isinstance(variable, (str, float, int)):
        # Wrap the variable in a list and convert it to float
        return [float(variable)]
    else:
        raise TypeError("Variable must be a string, float, int, or list.")


def do_inversion(
    n_elements,
    jacobian_dir,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    prior_err=0.5,
    obs_err=15,
    gamma=0.25,
    res="0.25x0.3125",
    jacobian_sf=None,
    prior_err_bc=0.0,
    prior_err_oh=0.0,
    is_Regional=True,
    verbose=False,
):
    """
    After running jacobian.py, use this script to perform the inversion and save out results.

    Arguments
        n_elements   [int]   : Number of state vector elements
        jacobian_dir [str]   : Directory where the data from jacobian.py are stored
        lon_min      [float] : Minimum longitude
        lon_max      [float] : Maximum longitude
        lat_min      [float] : Minimum latitude
        lat_max      [float] : Maximum latitude
        prior_err    [float] : Prior error standard deviation (default 0.5)
        obs_err      [float] : Observational error standard deviation (default 15 ppb)
        gamma        [float] : Regularization parameter (default 0.25)
        res          [str]   : Resolution string from config.yml (default '0.25x0.3125')
        jacobian_sf  [str]   : Path to Jacobian scale factors file if using precomputed K
        prior_err_bc [float] : Prior error standard deviation (default 0.0)
        prior_err_oh [float] : Prior error standard deviation (default 0.0)
        is_Regional  [bool]  : Is this a regional simulation?

    Returns
        xhat         [float] : Posterior scaling factors
        delta_optimized        [float] : Change from prior     [xhat = 1 + delta_optimized]
        KTinvSoK     [float] : K^T*inv(S_o)*K        [part of inversion equation]
        KTinvSoyKxA  [float] : K^T*inv(S_o)*(y-K*xA) [part of inversion equation]
        S_post       [float] : Posterior error covariance matrix
        A            [float] : Averaging kernel matrix

    """
    # boolean for whether we are optimizing boundary conditions
    bc_optimization = prior_err_bc > 0.0
    oh_optimization = prior_err_oh > 0.0

    # Need to ignore data in the GEOS-Chem 3 3 3 3 buffer zone
    # Shave off one or two degrees of latitude/longitude from each side of the domain
    # ~1 degree if 0.25x0.3125 resolution, ~2 degrees if 0.5x0.6125 resolution
    # This assumes 0.25x0.3125 and 0.5x0.625 simulations are always regional
    if "0.25x0.3125" in res:
        degx = 4 * 0.3125
        degy = 4 * 0.25
    elif "0.5x0.625" in res:
        degx = 4 * 0.625
        degy = 4 * 0.5
    else:
        degx = 0
        degy = 0

    xlim = [lon_min + degx, lon_max - degx]
    ylim = [lat_min + degy, lat_max - degy]

    # Read output data from jacobian.py (virtual & true TROPOMI columns, Jacobian matrix)
    files = glob.glob(f"{jacobian_dir}/*.pkl")
    files.sort()

    # ==========================================================================================
    # Now we will assemble two different expressions needed for the analytical inversion.
    #
    # These expressions are from eq. (5) and (6) in Zhang et al. (2018) ACP:
    # "Monitoring global OH concentrations using satellite observations of atmospheric methane".
    #
    # Specifically, we are going to solve:
    #   xhat = xA + G*(y-K*xA)
    #        = xA + inv(gamma * K^T*inv(S_o)*K + inv(S_a)) * gamma * K^T*inv(S_o) * (y-K*xA)
    #                          (--------------)                     (-----------------------)
    #                            Expression 1                             Expression 2
    #
    # Expression 1 = "KTinvSoK"
    # Expression 2 = "KTinvSoyKxA"
    #
    # In the code below this becomes
    #   xhat = xA + inv(gamma*KTinvSoK + inv(S_a)) * gamma*KTinvSoyKxA
    #        = xA + delta_optimized
    #        = 1  + delta_optimized      [since xA=1 when optimizing scale factors]
    #
    # We build KTinvSoK and KTinvSoyKxA "piece by piece", loading one jacobian .pkl file at a
    # time. This is so that we don't need to assemble or invert the full Jacobian matrix, which
    # can be very large.
    # ==========================================================================================

    # Initialize two expressions from the inversion equation
    KTinvSoK = np.zeros(
        [n_elements, n_elements], dtype=float
    )  # expression 1: K^T * inv(S_o) * K
    KTinvSoyKxA = np.zeros(
        [n_elements], dtype=float
    )  # expression 2: K^T * inv(S_o) * (y-K*xA)

    # Initialize
    # For each .pkl file generated by jacobian.py:
    for fi in files:
        if verbose:
            print(fi)

        # Load TROPOMI/GEOS-Chem and Jacobian matrix data from the .pkl file
        dat = load_obj(fi)

        # Skip if there aren't any TROPOMI observations on this day
        if dat["obs_GC"].shape[0] == 0:
            continue

        # Otherwise, grab the TROPOMI/GEOS-Chem data
        obs_GC = dat["obs_GC"]

        # Only consider data within the new latitude and longitude bounds
        ind = np.where(
            (obs_GC[:, 2] >= xlim[0])
            & (obs_GC[:, 2] <= xlim[1])
            & (obs_GC[:, 3] >= ylim[0])
            & (obs_GC[:, 3] <= ylim[1])
        )[0]

        # Skip if no data in bounds
        if len(ind) == 0:
            continue

        # TROPOMI and GEOS-Chem data within bounds
        obs_GC = obs_GC[ind, :]

        # weight obs_err based on the observation count to prevent overfitting
        # Note: weighting function defined by Zichong Chen for his
        # middle east inversions. May need to be tuned based on region.
        # From Chen et al. 2023:
        # "Satellite quantification of methane emissions and oil/gas methane
        # intensities from individual countries in the Middle East and North
        # Africa: implications for climate action"
        s_superO_1 = calculate_superobservation_error(obs_err, 1)
        s_superO_p = np.array(
            [
                calculate_superobservation_error(obs_err, p) if p >= 1 else s_superO_1
                for p in obs_GC[:, 4]
            ]
        )
        # Define observational errors (diagonal entries of S_o matrix)
        obs_error = np.power(obs_err, 2)
        gP = s_superO_p**2 / s_superO_1**2
        # scale error variance by gP
        obs_error = gP * obs_error

        # check to make sure obs_err isn't negative, set 1 as default value
        obs_error = [obs if obs > 0 else 1 for obs in obs_error]

        # Jacobian entries for observations within bounds [ppb]
        if jacobian_sf is None:
            K = 1e9 * dat["K"][ind, :]
        else:
            # Get Jacobian from reference inversion
            fi_ref = fi.replace("data_converted", "data_converted_reference")
            dat_ref = load_obj(fi_ref)
            K = 1e9 * dat_ref["K"][ind, :]

        # Number of observations
        if verbose:
            print("Sum of Jacobian entries:", np.sum(K))

        # Apply scaling matrix if using precomputed Jacobian
        if jacobian_sf is not None:
            scale_factors = np.load(jacobian_sf)
            if bc_optimization:
                # add (unit) scale factors for BCs
                # as the last 4 elements of the scaling matrix
                scale_factors = np.append(scale_factors, np.ones(4))
            reps = K.shape[0]
            scaling_matrix = np.tile(scale_factors, (reps, 1))
            if oh_optimization:
                K[:, :-2] *= scaling_matrix
            else:
                K *= scaling_matrix

        # Measurement-model mismatch: TROPOMI columns minus GEOS-Chem virtual TROPOMI columns
        # This is (y - F(xA)), i.e., (y - (K*xA + c)) or (y - K*xA) in shorthand
        delta_y = obs_GC[:, 0] - obs_GC[:, 1]  # [ppb]

        # If there are any nans in the data, abort
        if (
            np.any(np.isnan(delta_y))
            or np.any(np.isnan(K))
            or np.any(np.isnan(obs_error))
        ):
            print("missing values", fi)
            break

        # Define KTinvSo = K^T * inv(S_o)
        KT = K.transpose()
        KTinvSo = np.zeros(KT.shape, dtype=float)
        for k in range(KT.shape[1]):
            KTinvSo[:, k] = KT[:, k] / obs_error[k]

        # Parts of inversion equation
        partial_KTinvSoK = KTinvSo @ K  # expression 1: K^T * inv(S_o) * K
        partial_KTinvSoyKxA = (
            KTinvSo @ delta_y
        )  # expression 2: K^T * inv(S_o) * (y-K*xA)

        # Add partial expressions to sums
        KTinvSoK += partial_KTinvSoK
        KTinvSoyKxA += partial_KTinvSoyKxA

    # Inverse of prior error covariance matrix, inv(S_a)
    Sa_diag = np.zeros(n_elements)
    Sa_diag.fill(prior_err**2)

    # Number of elements to apply scale factor to
    scale_factor_idx = n_elements

    # If optimizing OH, adjust for it in the inversion
    if oh_optimization:
        # Add prior error for OH as the last element(s) of the diagonal
        # Following Masakkers et al. (2019, ACP) weight the OH term by the
        # ratio of the number of elements (n_OH_elements/n_emission_elements)
        if is_Regional:
            OH_weight = 1 / (n_elements - 1)
            Sa_diag[-1:] = OH_weight * prior_err_oh**2
            scale_factor_idx -= 1
        else:
            OH_weight = 2 / (n_elements - 2)
            Sa_diag[-2:] = OH_weight * prior_err_oh**2
            scale_factor_idx -= 2

    # If optimizing boundary conditions, adjust for it in the inversion
    if bc_optimization:
        scale_factor_idx -= 4

        # add prior error for BCs as the last 4 elements of the diagonal
        if oh_optimization:
            if is_Regional:
                Sa_diag[-5:-1] = prior_err_bc**2
            else:
                Sa_diag[-6:-1] = prior_err_bc**2
        else:
            Sa_diag[-4:] = prior_err_bc**2

    inv_Sa = np.diag(1 / Sa_diag)  # Inverse of prior error covariance matrix

    # Solve for posterior scale factors xhat
    delta_optimized = np.linalg.inv(gamma * KTinvSoK + inv_Sa) @ (gamma * KTinvSoyKxA)

    # Update scale factors by 1 to match what GEOS-Chem expects
    # xhat = 1 + delta_optimized
    # Notes:
    #  - If optimizing BCs, the last 4 elements are in concentration space,
    #    so we do not need to add 1
    #  - If optimizing OH, the last element also needs to be updated by 1
    xhat = delta_optimized.copy()
    xhat[:scale_factor_idx] += 1
    if oh_optimization:
        if is_Regional:
            xhat[-1] += 1
            print(f"xhat[OH] = {xhat[-1]}")
        else:
            xhat[-2:] += 1
            print(f"xhat[OH] = {xhat[-2:]}")

    # Posterior error covariance matrix
    S_post = np.linalg.inv(gamma * KTinvSoK + inv_Sa)

    # Averaging kernel matrix
    A = np.identity(n_elements) - S_post @ inv_Sa

    # Calculate J_A, where delta_optimized = xhat - xA
    # J_A = (xhat - xA)^T * inv_Sa * (xhat - xA)
    delta_optimizedT = delta_optimized.transpose()
    J_A = delta_optimizedT @ inv_Sa @ delta_optimized
    J_A_normalized = J_A / n_elements

    # Print some statistics
    print(
        f"hyperparameters: (prior_err: {prior_err}, obs_err: {obs_err}, gamma: {gamma}, "
        + f"prior_err_bc: {prior_err_bc}, prior_err_oh: {prior_err_oh})"
    )
    print(
        f"Normalized J_A: {J_A_normalized}"
    )  # ideal gamma is where this is close to 1
    print(
        "Min:",
        xhat[:scale_factor_idx].min(),
        "Mean:",
        xhat[:scale_factor_idx].mean(),
        "Max",
        xhat[:scale_factor_idx].max(),
    )

    return xhat, delta_optimized, KTinvSoK, KTinvSoyKxA, S_post, A, J_A_normalized


def do_inversion_ensemble(
    n_elements,
    jacobian_dir,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    prior_errs,
    obs_errs,
    gammas,
    res,
    jacobian_sf,
    prior_errs_bc,
    prior_errs_oh,
    is_Regional,
):
    """
    Run series of inversions with hyperparameter vectors and save out the results.
    """
    hyperparam_ensemble = list(
        product(prior_errs, obs_errs, gammas, prior_errs_bc, prior_errs_oh)
    )

    results_dict = {
        "KTinvSoK": [],
        "KTinvSoyKxA": [],
        "ratio": [],
        "xhat": [],
        "S_post": [],
        "A": [],
        "J_A_normalized": [],
        "labels": [],
        "hyperparameters": [],
    }
    for member in hyperparam_ensemble:
        prior_err, obs_err, gamma, prior_err_bc, prior_err_oh = member
        params = {
            "prior_err": prior_err, 
            "obs_err": obs_err,   
            "gamma": gamma,     
            "prior_err_bc": prior_err_bc,
            "prior_err_oh": prior_err_oh,
        }
        # create a label for the output file variables
        label = f"_prior_err_{prior_err}_obs_err_{obs_err}_gamma_{gamma}"
        if prior_err_bc > 0:
            label = label + f"_prior_err_bc_{prior_err_bc}"
        if prior_err_oh > 0:
            label = label + f"_prior_err_oh_{prior_err_oh}"
        xhat, delta_optimized, KTinvSoK, KTinvSoyKxA, S_post, A, J_A_normalized = (
            do_inversion(
                n_elements,
                jacobian_dir,
                lon_min,
                lon_max,
                lat_min,
                lat_max,
                prior_err,
                obs_err,
                gamma,
                res,
                jacobian_sf,
                prior_err_bc,
                prior_err_oh,
                is_Regional,
                verbose=False,
            )
        )
        results_dict["KTinvSoK"].append(KTinvSoK)
        results_dict["KTinvSoyKxA"].append(KTinvSoyKxA)
        results_dict["ratio"].append(delta_optimized)
        results_dict["xhat"].append(xhat)
        results_dict["S_post"].append(S_post)
        results_dict["A"].append(A)
        results_dict["J_A_normalized"].append(J_A_normalized)
        results_dict["labels"].append(label)
        results_dict["hyperparameters"].append(params)

    # Find the ensemble member that is closest to 1 following Lu et al. (2021)
    idx_best_Ja = np.argmin(np.abs(np.array(results_dict["J_A_normalized"]) - 1))
    print(
        f"Best J_A: {results_dict['J_A_normalized'][idx_best_Ja]} with"
        + f" (prior_err, obs_err, gamma, prior_err_bc, prior_err_oh) = {hyperparam_ensemble[idx_best_Ja]}"
    )

    # Create an xarray.Dataset
    dataset = xr.Dataset()
    for key, values in results_dict.items():
        if key == "labels":
            continue  # Skip labels; used only for naming variables
        for idx, (label, value) in enumerate(zip(results_dict["labels"], values)):
            if isinstance(value, np.ndarray) and value.ndim == 2:
                dims = ("nvar1", "nvar2")
            elif isinstance(value, np.ndarray):
                dims = ("nvar1",)
            else:
                dims = ()  # Scalar

            dataset[f"{key}{label}"] = (dims, value)

            # Use the best J_A as the default data variable values
            if idx == idx_best_Ja:
                dataset[key] = (dims, value)

    # reorder the variables, so the default vars are at the top
    best_vars = list(results_dict.keys())
    best_vars.remove("labels")
    best_vars.remove("hyperparameters")
    ensemble_vars = [item for item in list(dataset.data_vars) if item not in best_vars]

    new_dataset = xr.Dataset(
        {var: dataset[var] for var in best_vars + ensemble_vars},
        coords=dataset.coords,
    )

    dataset = new_dataset

    return dataset


if __name__ == "__main__":
    import sys

    config_path = sys.argv[1]
    n_elements = int(sys.argv[2])
    jacobian_dir = sys.argv[3]
    output_path = sys.argv[4]
    lon_min = float(sys.argv[5])
    lon_max = float(sys.argv[6])
    lat_min = float(sys.argv[7])
    lat_max = float(sys.argv[8])
    res = sys.argv[9]
    jacobian_sf = sys.argv[10]

    # read in config file
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # set parameters based on config file
    is_Regional = ensure_float_list(config["isRegional"])
    prior_err = ensure_float_list(config["PriorError"])
    obs_err = ensure_float_list(config["ObsError"])
    gamma = ensure_float_list(config["Gamma"])

    # 0.0 if not optimizing BCs or OH
    prior_err_BC = config["PriorErrorBCs"] if config["OptimizeBCs"] else 0.0
    prior_err_OH = config["PriorErrorOH"] if config["OptimizeOH"] else 0.0
    prior_err_BC = ensure_float_list(prior_err_BC)
    prior_err_OH = ensure_float_list(prior_err_OH)

    # Reformat Jacobian scale factor input
    if jacobian_sf == "None":
        jacobian_sf = None

    # Run the inversion code
    out_ds = do_inversion_ensemble(
        n_elements,
        jacobian_dir,
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        prior_err,
        obs_err,
        gamma,
        res,
        jacobian_sf,
        prior_err_BC,
        prior_err_OH,
        is_Regional,
    )

    # Save the results of the ensemble inversion
    out_ds.to_netcdf(output_path)
    print(f"Saved results to {output_path}")

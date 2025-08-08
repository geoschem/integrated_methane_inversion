# Description: Script to perform inversion using lognormal errors
# Usage: python lognormal_invert.py <path_to_config_file> <path_to_state_vector_file> <jacobian_sf>
# Inputs:
#       path_to_config_file: path to yaml config file
#       path_to_state_vector_file: path to state vector netcdf file
#       jacobian_sf: (optional) path to numpy array of scale factors for jacobian

import sys
import yaml
from itertools import product
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from scipy.sparse import spdiags
from src.inversion_scripts.utils import ensure_float_list
from src.inversion_scripts.make_gridded_posterior import make_gridded_posterior


def lognormal_invert(config, state_vector_filepath, jacobian_sf):
    """
    Description:
        Run inversion using lognormal errors following method from eqn 2 of
        Chen et al., 2022 https://doi.org/10.5194/acp-22-10809-2022
        Outputs inversion results to netcdf files.
    Arguments:
        config                [Dict]   : dictionary of config variables
        state_vector_filepath [String] : path to state vector netcdf file
        jacobian_sf           [String] : path to numpy array of scale factors
    """
    results_save_path = f"inversion_result_ln.nc"
    # dictionary to store inversion results
    results_dict = {
        "xhat": [],
        "lnxn": [],
        "S_post": [],
        "A": [],
        "DOFS": [],
        "Ja_normalized": [],
        "prior_err": [],
        "obs_err": [],
        "gamma": [],
        "prior_err_bc": [],
        "prior_err_oh": [],
        "prior_err_buffer": [],
    }

    state_vector = xr.load_dataset(state_vector_filepath)
    state_vector_labels = state_vector["StateVector"]
    lats, lons = state_vector_labels.lat, state_vector_labels.lon

    # used to determine convergence of xn 5e-3 is .5%
    convergence_threshold = 5e-3

    # Load in the observation and background data
    ds = np.load("obs_ch4_tropomi.npz")
    y = np.array(ds["obs_tropomi"])
    ds = np.load("gc_ch4_bkgd.npz")
    ybkg = np.array(ds["gc_ch4_bkgd"])

    # We only solve using lognormal errors for state vector elements
    # within the domain of interest, not the buffer elements, the
    # BC elements, or OH optimization. So, to do this we split K into
    # two matrices, one for the lognormal elements, and one for the
    # normal elements.
    optimize_bcs = config["OptimizeBCs"]
    optimize_oh = config["OptimizeOH"]
    is_regional = config["isRegional"]
    if optimize_oh:
        if is_regional:
            OH_element_num = 1
        else:
            OH_element_num = 2
    else:
        OH_element_num = 0
    BC_element_num = 4 if optimize_bcs else 0
    num_sv_elems = (
        int(state_vector_labels.max().item()) + BC_element_num + OH_element_num
    )
    num_buffer_elems = int(config["nBufferClusters"])
    num_normal_elems = num_buffer_elems + BC_element_num + OH_element_num
    ds = np.load("full_jacobian_K.npz")
    K_temp = np.array(ds["K"]) * 1e9

    # Apply scaling matrix if using precomputed Jacobian
    if jacobian_sf is not None:
        scale_factors = np.load(jacobian_sf)
        # apply unit scaling for BC elements and OH elements if using
        if optimize_bcs or optimize_oh:
            scale_factors = np.append(
                scale_factors, np.ones(BC_element_num + OH_element_num)
            )
        reps = K_temp.shape[0]
        scaling_matrix = np.tile(scale_factors, (reps, 1))
        K_temp *= scaling_matrix

    # Define Sa, gamma, So, and Sa_bc values to iterate through
    prior_errors = ensure_float_list(config["PriorError"])
    sa_buffer_elems = ensure_float_list(config["PriorErrorBufferElements"])
    sa_bc_vals = ensure_float_list(config["PriorErrorBCs"]) if optimize_bcs else [0.0]
    sa_oh_vals = ensure_float_list(config["PriorErrorOH"]) if optimize_oh else [0.0]
    gamma_vals = ensure_float_list(config["Gamma"])
    obs_err_keys = [
        f"so_{obs_err}" for obs_err in ensure_float_list(config["ObsError"])
    ]

    so_dict = np.load("so_super.npz")

    # Calculate the difference between tropomi and the background
    # simulation, which has no emissions
    y_ybkg_diff = y - ybkg
    ybkg, y, y_ybkg_diff = (
        np.swapaxes(ybkg, 0, 1),
        np.swapaxes(y, 0, 1),
        np.swapaxes(y_ybkg_diff, 0, 1),
    )

    # fixed kappa of 10 following Chen et al., 2022 https://doi.org/10.5194/acp-22-10809-2022
    kappa = 10

    # iterate through different combination of gamma, lnsa, and sa_bc
    # TODO: parallelize this once we allow vectorization of these values
    combinations = list(
        product(
            gamma_vals,
            prior_errors,
            obs_err_keys,
            sa_bc_vals,
            sa_buffer_elems,
            sa_oh_vals,
        )
    )
    for gamma, sa, so_key, sa_bc, sa_buffer, sa_oh in combinations:
        # params dict to store hyperparameters
        params = {
            "prior_err": sa,
            "obs_err": float(so_key.removeprefix("so_")),
            "gamma": gamma,
            "prior_err_bc": sa_bc,
            "prior_err_oh": sa_oh,
            "prior_err_buffer": sa_buffer,
        }
        # The levenberg-marquardt method assumes that the prior emissions is
        # the median prior emissions, but typically priors are the mean emission.
        # To account for this we convert xa to a median. This can be done by
        # scaling the lognormal part of K by 1/exp((lnsa**2)/2).
        # Here, we calculate this scaling factor
        prior_scale = 1 / np.exp((np.log(float(sa)) ** 2) / 2)

        # split K based on whether we are solving for lognormal or normal elements
        # K_ROI is the matrix for the lognormal elements (the region of interest)
        # the lognormal part of K gets scaled by prior_scale to convert to median
        K_ROI = K_temp[:, :-num_normal_elems]
        K_normal = K_temp[:, -num_normal_elems:]
        K_ROI = prior_scale * K_ROI
        K_full = np.concatenate((K_ROI, K_normal), axis=1)

        m, n = np.shape(K_ROI)

        # Create base xa and lnxa matrices
        # Note: the resulting xa vector has lognormal elements until the
        # final Buffer, BCs, and OH elements
        xa = np.ones((n, 1)) * 1.0
        lnxa = np.log(xa)

        # Create normal elements for buffer, BCs, and OH
        # BC elements are relative to 0 because they are in concentration space
        # Other elements are in scale factor space where 1 is the prior
        xa_normal_buffer = np.ones((num_buffer_elems, 1)) * 1.0
        xa_normal_BCs = np.ones((BC_element_num, 1)) * 0.0
        xa_normal_OH = np.ones((OH_element_num, 1)) * 1.0
        xa_normal = np.concatenate(
            (xa_normal_buffer, xa_normal_BCs, xa_normal_OH), axis=0
        )

        # concatenate normal elements to xa and lnxa
        xa = np.concatenate((xa, xa_normal), axis=0)
        lnxa = np.concatenate((lnxa, xa_normal), axis=0)

        # get the So matrix
        so = so_dict[so_key]

        # Create inverted So matrix
        Soinv = spdiags(1 / so, 0, m, m)

        lnsa_val = np.log(sa)

        # Create lnSa matrix
        # lnsa = lnsa_val**2 * np.ones((n, 1))
        lnsa = lnsa_val**2 * np.ones((n, 1))

        # For the buffer elems, BCs, and OH elements
        # we apply a different Sa value
        # In the most basic we only generate Sa for buffer elements
        base_sa_normal = sa_buffer**2 * np.ones(
            (num_normal_elems - (BC_element_num + OH_element_num), 1)
        )

        # conditionally add BC and OH elements
        if optimize_bcs:
            bc_errors = sa_bc**2 * np.ones((BC_element_num, 1))
            base_sa_normal = np.concatenate((base_sa_normal, bc_errors), axis=0)

        if optimize_oh:
            oh_errors = sa_oh**2 * np.ones((OH_element_num, 1))
            # weight the OH term(s) following Maasakkers et al. (2019)
            oh_weight = OH_element_num / (num_normal_elems - OH_element_num)
            oh_errors_constraint = (oh_weight * sa_oh**2) * np.ones((OH_element_num, 1))
            sa_normal = np.concatenate(
                (base_sa_normal, oh_errors), axis=0
            )  # unweighted Sa vector
            sa_normal_constraint = np.concatenate(
                (base_sa_normal, oh_errors_constraint), axis=0
            )  # weighted Sa vector
        else:
            sa_normal = base_sa_normal.copy()
            sa_normal_constraint = base_sa_normal.copy()

        # concatenate lognormal prior errors with normal prior errors
        lnsa_arr = np.concatenate((lnsa, sa_normal), axis=0)
        lnsa_arr_constraint = np.concatenate((lnsa, sa_normal_constraint), axis=0)

        # Create two separate lnSa matrices
        lnsa = np.zeros((n + num_normal_elems, n + num_normal_elems))
        lnsa_constraint = lnsa.copy()
        np.fill_diagonal(lnsa, lnsa_arr)  # unweighted
        np.fill_diagonal(lnsa_constraint, lnsa_arr_constraint)  # weighted
        invlnsa = np.linalg.inv(lnsa)
        invlnsa_constraint = np.linalg.inv(lnsa_constraint)

        # we start with lnxa using the prior values (scale factors of ln(1))
        lnxn = lnxa

        # start with arbitrary value for xn_iteration_pct_diff above .05%
        xn_iteration_pct_diff = 1

        # Iterate for calculation of ln(xn) until convergence threshold is met (5e-3)
        # We decompose eqn 2 from chen et al into 4 terms
        # term 1: gamma*K'.T@inv(So)@K'
        # term 2: inv((1+kappa)*inv(ln(sa)))
        # term 3: gamma*K'.T@inv(So)@(y_ybkg_diff - K@xn)
        # term 4: -inv(ln(sa))@(ln(xn) - ln(xa))
        # We can then solve for xn iteratively by doing:
        # ln(xn) = x(n-1) + inv(term1+term2)@(term3 + term4)
        # where x(n-1) is the previous iteration of xn until convergence
        print("Status: Iterating to calculate ln(xn)")
        while xn_iteration_pct_diff >= convergence_threshold:

            # we need to transform lnxn to xn to calculate K_prime
            xn = np.concatenate(
                (np.exp(lnxn[:-num_normal_elems]), lnxn[-num_normal_elems:]),
                axis=0,
            )
            # K_prime is the updated jacobian using the new xn from the previous iteration
            K_prime = np.concatenate(
                (K_ROI * xn[:-num_normal_elems].T, K_normal), axis=1
            )

            # commonly used term for term1 and term3
            gamma_K_prime_transpose_Soinv = gamma * K_prime.T @ Soinv

            # Compute the next xn_update
            term1 = gamma_K_prime_transpose_Soinv @ K_prime
            term2 = (1 + kappa) * invlnsa_constraint
            inv_term = np.linalg.inv(term1 + term2)

            term3 = gamma_K_prime_transpose_Soinv @ (y_ybkg_diff - K_full @ xn)
            term4 = invlnsa_constraint @ (lnxn - lnxa)

            # put it all together to calculate lnxn_update
            lnxn_update = lnxn + inv_term @ (term3 - term4)

            # Check for convergence
            xn_iteration_pct_diff = max(
                abs(
                    np.exp(lnxn_update[:-num_normal_elems])
                    - np.exp(lnxn[:-num_normal_elems])
                )
                / np.exp(lnxn[:-num_normal_elems])
            )

            lnxn = lnxn_update

        print("Status: Done Iterating")
        xn = np.concatenate(
            (np.exp(lnxn[:-num_normal_elems]), lnxn[-num_normal_elems:]), axis=0
        )

        # Calculate averaging kernel and degrees of freedom for signal
        K_primeT_so = gamma * np.transpose(K_prime) @ Soinv
        # posterior error covariance matrix (uses unweighted Sa)
        lns = np.linalg.inv(K_primeT_so @ K_prime + invlnsa)
        # Averaging kernel (uses unweighted Sa)
        G = lns @ K_primeT_so
        ak = G @ K_prime
        dofs = np.trace(ak)
        print(f"DOFS: {dofs}")

        # Calculate posterior mean xhat
        dlns = np.diag(lns[:-num_normal_elems, :-num_normal_elems])
        xnmean = np.concatenate(
            (
                xn[:-num_normal_elems]
                * np.expand_dims(np.exp(dlns * (0.5)) * prior_scale, axis=1),
                xn[-num_normal_elems:],
            )
        )

        # Calculate Ja diagnostic only for domain of interest (ignoring buffer and BC elements)
        # Ja diagnostic is useful for determining regurlarization parameter (gamma)
        Ja = (
            np.transpose(lnxn[:-num_normal_elems] - lnxa[:-num_normal_elems])
            @ invlnsa[:-num_normal_elems, :-num_normal_elems]
            @ (lnxn[:-num_normal_elems] - lnxa[:-num_normal_elems])
        )

        print(
            f"Diagnostics:\n  (Ja: {Ja}, gamma: {gamma}, "
            + f"sa: {sa}, sa_bc: {sa_bc}, sa_oh: {sa_oh_vals} sa_buffer: {sa_buffer})"
        )

        # Append results to results_dict
        xhat = xnmean
        print(f"xhat = {xhat.sum()}")
        results_dict["xhat"].append(xhat.flatten()),
        results_dict["lnxn"].append(lnxn.flatten()),
        results_dict["S_post"].append(lns),
        results_dict["A"].append(ak),
        results_dict["DOFS"].append(dofs),
        results_dict["Ja_normalized"].append(Ja.item() / num_sv_elems),
        for k, v in params.items():
            results_dict[k].append(v)

    # Define the default data variables as those with normalized Ja closest to 1
    idx_default_Ja = np.argmin(np.abs(np.array(results_dict["Ja_normalized"]) - 1))

    # Filter ensemble members to only members with Ja between 0.5 and 2
    filter_ens_members = True  # set to False to turn off filtering
    include_ens_members = [
        i for i, Ja in enumerate(results_dict["Ja_normalized"]) if 0.5 <= Ja <= 2.0
    ]

    if filter_ens_members and len(include_ens_members) > 0:
        for k in results_dict.keys():
            results_dict[k] = [results_dict[k][i] for i in include_ens_members]
    elif len(include_ens_members) == 0:
        print(
            "Warning: No ensemble members with 0.5 <= J_A/n <= 2.0, "
            + "Returning all members in ensemble. This may lead to suboptimal results."
            + " Consider adding additional ensemble members with different hyperparameters."
        )
    else:
        print(
            "Warning: Returning all members in ensemble without filtering "
            + "Ja/n thresholds [0.5, 2.0]. This may lead to suboptimal results."
            + " Consider adding ensemble filters."
        )

    # Create an xarray dataset to store inversion results
    dataset = xr.Dataset()
    for k, v in results_dict.items():
        v = np.array(v)
        dims = ["ensemble"] + [f"nvar{i}" for i in range(1, v.ndim)]
        dataset[k] = (dims, v)

    # save index number of ens member with J_A/n
    # closes to 1 as the default member
    dataset.attrs = {"default_member_index": idx_default_Ja}

    # ensemble dimension to end
    dataset = dataset.transpose(..., "ensemble")

    # also calculate the mean of the ensemble as the main result
    dataset_mean = dataset.mean(dim="ensemble")

    dataset.to_netcdf(
        results_save_path.replace(".nc", "_ensemble.nc"),
        encoding={v: {"zlib": True, "complevel": 1} for v in dataset.data_vars},
    )

    dataset_mean.to_netcdf(
        results_save_path,
        encoding={v: {"zlib": True, "complevel": 1} for v in dataset_mean.data_vars},
    )

    # make gridded posterior
    make_gridded_posterior(
        results_save_path.replace(".nc", "_ensemble.nc"),
        state_vector_filepath,
        "gridded_posterior_ln.nc",
    )


if __name__ == "__main__":
    config_path = sys.argv[1]
    state_vector_filepath = sys.argv[2]
    jacobian_sf = None if sys.argv[3] == "None" else sys.argv[3]

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    lognormal_invert(config, state_vector_filepath, jacobian_sf)

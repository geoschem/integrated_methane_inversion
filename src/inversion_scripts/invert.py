#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import numpy as np
import xarray as xr
from itertools import product
from pathlib import Path
from collections import defaultdict, deque
from src.inversion_scripts.utils import (
    load_obj,
    calculate_superobservation_error,
    ensure_float_list,
    get_mean_emissions,
    update_prior_error_for_OptimizeSoil,
    map_files_to_reference,
)
from src.utilities.config_utils import load_config


def align_obs_rows_with_reference(obs_GC, obs_GC_ref):
    """Match target observations to reference rows using shared metadata columns (lat, lon, obs_count)."""
    ncols = min(obs_GC.shape[1], obs_GC_ref.shape[1])

    def make_key(row):
        # Ignore the leading observed/model xCH4 columns and match on shared metadata.
        return tuple(np.round(row[2:ncols], decimals=6))

    ref_lookup = defaultdict(deque)
    for idx, row in enumerate(obs_GC_ref):
        ref_lookup[make_key(row)].append(idx)

    obs_indices = []
    ref_indices = []
    for idx, row in enumerate(obs_GC):
        key = make_key(row)
        if ref_lookup[key]:
            obs_indices.append(idx)
            ref_indices.append(ref_lookup[key].popleft())

    return np.asarray(obs_indices, dtype=int), np.asarray(ref_indices, dtype=int)

def get_prior_sigma_vector(
    n_elements,
    prior_err,
    OptimizeSoil=False,
    prior_ds=None,
    StateVectorFile=None,
):
    """Return the prior standard deviation for each state-vector element."""
    sigma = np.full(n_elements, prior_err, dtype=float)
    if OptimizeSoil:
        prior_err_new = update_prior_error_for_OptimizeSoil(
            prior_ds, prior_err, StateVectorFile, n_elements
        )
        sigma[: len(prior_err_new)] = prior_err_new
    return sigma


def get_expected_state_vector_ids(StateVectorFile):
    """Return the sorted state-vector IDs defined in the active state-vector file."""
    state_vector = xr.load_dataset(StateVectorFile)
    state_vector_ids = state_vector["StateVector"].values.reshape(-1)
    state_vector_ids = state_vector_ids[np.isfinite(state_vector_ids)]
    state_vector_ids = state_vector_ids[state_vector_ids > 0]
    state_vector_ids = np.unique(state_vector_ids.astype(np.int32))
    return np.sort(state_vector_ids)


def build_prior_covariance(
    n_elements,
    prior_err,
    OptimizeSoil=False,
    prior_ds=None,
    StateVectorFile=None,
    prebuilt_prior_err_covariance=False,
):
    """
    Build the prior covariance in either full or diagonal form.

    Returns the covariance, the constraint covariance, and a boolean indicating
    whether the covariance should be treated as a full matrix downstream.
    """
    if prebuilt_prior_err_covariance:
        # Load prebuilt covariance matrix with off-diagonal elements
        Sa = np.zeros((n_elements, n_elements), dtype=float)
        covariance_path = Path("prior_norm_error_covariance.npz")
        if not covariance_path.exists():
            raise FileNotFoundError(f"Covariance matrix file not found: {covariance_path}")
        with np.load(covariance_path) as prebuilt:
            Sa_prebuilt = prebuilt["covariance"]
            state_vector_ids_prebuilt = prebuilt["state_vector_ids"]

        expected_state_vector_ids = get_expected_state_vector_ids(StateVectorFile)
        state_vector_ids_prebuilt = np.asarray(state_vector_ids_prebuilt, dtype=np.int32)
        if Sa_prebuilt.shape[0] != Sa_prebuilt.shape[1]:
            raise ValueError(
                f"Prior covariance must be square, got shape {Sa_prebuilt.shape}"
            )
        if Sa_prebuilt.shape[0] != state_vector_ids_prebuilt.size:
            raise ValueError(
                "Prior covariance size does not match the saved state-vector IDs: "
                f"{Sa_prebuilt.shape[0]} vs {state_vector_ids_prebuilt.size}"
            )
        if not np.array_equal(state_vector_ids_prebuilt, expected_state_vector_ids):
            raise ValueError(
                "Saved prior covariance state-vector IDs do not match the active "
                "StateVectorFile."
            )
        if Sa_prebuilt.shape[0] > n_elements:
            raise ValueError(
                "Prior covariance block is larger than the inversion state vector: "
                f"{Sa_prebuilt.shape[0]} > {n_elements}"
            )

        Sa_prebuilt_elems = Sa_prebuilt.shape[0]
        # The prebuilt matrix can define only a leading subset of the full state vector,
        # so we scale and insert it into the top-left block.
        sigma_prebuilt = get_prior_sigma_vector(
            Sa_prebuilt_elems,
            prior_err,
            OptimizeSoil=OptimizeSoil,
            prior_ds=prior_ds,
            StateVectorFile=StateVectorFile,
        )
        # The prebuilt matrix stores only the normalized covariance structure, so we
        # apply sigma_i * sigma_j here. This reduces to a scalar prior_err**2 factor
        # when all sigmas are the same, but supports element-wise prior_err_new values.
        Sa[:Sa_prebuilt_elems, :Sa_prebuilt_elems] = (
            sigma_prebuilt[:, None] * Sa_prebuilt * sigma_prebuilt[None, :]
        )
        return Sa, Sa.copy(), True

    # Otherwise, build only a diagonal covariance matrix
    sigma = get_prior_sigma_vector(
        n_elements,
        prior_err,
        OptimizeSoil=OptimizeSoil,
        prior_ds=prior_ds,
        StateVectorFile=StateVectorFile,
    )
    Sa_diag = sigma**2
    return Sa_diag, Sa_diag.copy(), False


def apply_diagonal_prior(matrix, indices, value):
    """Write a prior variance value onto either a 1D diagonal vector or 2D matrix."""
    if np.isscalar(matrix[0]) or getattr(matrix, "ndim", 1) == 1:
        matrix[indices] = value
    else:
        # Convert slices like [-4:] into explicit diagonal positions for 2D matrices.
        diag_indices = np.arange(matrix.shape[0])[indices]
        matrix[diag_indices, diag_indices] = value


def get_oh_index_slice(n_elements, is_Regional):
    """Return the slice occupied by OH state-vector elements."""
    return slice(-1, None) if is_Regional else slice(-2, None)


def get_bc_index_slice(optimize_oh, is_Regional):
    """Return the slice occupied by boundary-condition state-vector elements."""
    if optimize_oh:
        return slice(-5, -1) if is_Regional else slice(-6, -2)
    return slice(-4, None)


def apply_oh_prior(Sa, Sa_constraint, n_elements, prior_err_oh, is_Regional):
    """Apply OH prior variances to both the unweighted and weighted constraint priors."""
    if is_Regional:
        OH_weight = 1 / (n_elements - 1)
    else:
        OH_weight = 2 / (n_elements - 2)
    oh_slice = get_oh_index_slice(n_elements, is_Regional)
    apply_diagonal_prior(Sa, oh_slice, prior_err_oh**2)
    apply_diagonal_prior(Sa_constraint, oh_slice, OH_weight * prior_err_oh**2)


def apply_bc_prior(Sa, Sa_constraint, prior_err_bc, optimize_oh, is_Regional):
    """Apply BC prior variances to both the unweighted and weighted constraint priors."""
    bc_slice = get_bc_index_slice(optimize_oh, is_Regional)
    apply_diagonal_prior(Sa, bc_slice, prior_err_bc**2)
    apply_diagonal_prior(Sa_constraint, bc_slice, prior_err_bc**2)


def invert_prior_covariance(Sa, Sa_constraint, use_full_prior_covariance):
    """Invert either full prior covariances or diagonal prior-variance vectors."""
    if use_full_prior_covariance:
        return np.linalg.inv(Sa_constraint), np.linalg.inv(Sa)
    return np.diag(1 / Sa_constraint), np.diag(1 / Sa)


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
    OptimizeSoil=False,
    prior_ds=None,
    StateVectorFile=None,
    verbose=False,
    prebuilt_prior_err_covariance=False,
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
        OptimizeSoil [bool]  : Optimize soil sink?
        prior_ds     [xr.Dataset]: prior emission dataset
        StateVectorFile [str]: Path to gridded state vector file

    Returns
        xhat         [float] : Posterior scaling factors
        delta_optimized        [float] : Change from prior     [xhat = 1 + delta_optimized]
        KTinvSoK     [float] : K^T*inv(S_o)*K        [part of inversion equation]
        KTinvSoyKxA  [float] : K^T*inv(S_o)*(y-K*xA) [part of inversion equation]
        S_post       [float] : Posterior error covariance matrix
        A            [float] : Averaging kernel matrix

    """
    # make mapping of target files to reference files if using precomputed Jacobian
    if jacobian_sf is not None:
        reference_dir = jacobian_dir.replace(
            "data_converted", "data_converted_reference"
        )
        K_ref_file_mappings = map_files_to_reference(jacobian_dir, reference_dir)

    # boolean for whether we are optimizing boundary conditions
    optimize_bc = prior_err_bc > 0.0
    optimize_oh = prior_err_oh > 0.0

    # Need to ignore data in the GEOS-Chem 3 3 3 3 buffer zone
    # Shave off one or two degrees of latitude/longitude from each side of the domain
    # ~1 degree if 0.25x0.3125 resolution, ~2 degrees if 0.5x0.6125 resolution
    # This assumes 0.125x0.15625, 0.25x0.3125, and 0.5x0.625 simulations are always regional
    if "0.125x0.15625" in res:
        degx = 4 * 0.15625
        degy = 4 * 0.125
    elif "0.25x0.3125" in res:
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

    # Read output data from jacobian.py (virtual & true satellite columns, Jacobian matrix)
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

        # Load satellite/GEOS-Chem and Jacobian matrix data from the .pkl file
        dat = load_obj(fi)

        # Skip if there aren't any satellite observations on this day
        if dat["obs_GC"].shape[0] == 0:
            continue

        # Otherwise, grab the satellite/GEOS-Chem data
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

        # Satellite and GEOS-Chem data within bounds
        obs_GC = obs_GC[ind, :]

        ref_ind = None
        if jacobian_sf is not None:
            # Precomputed Jacobians come from a reference run, so align observations
            # before indexing K to ensure each retained row maps to the same scene.
            fi_ref = str(K_ref_file_mappings.get(Path(fi)))
            if fi_ref is None:
                print(f"No reference file found for {fi} in {jacobian_dir}")
                continue
            dat_ref = load_obj(fi_ref)
            obs_ind, ref_ind = align_obs_rows_with_reference(obs_GC, dat_ref["obs_GC"])
            if len(ref_ind) == 0:
                print(f"No overlapping reference observations found for {fi_ref}")
                continue
            obs_GC = obs_GC[obs_ind, :]

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
            K = 1e9 * dat_ref["K"][ref_ind, :]

        # Number of observations
        if verbose:
            print("Sum of Jacobian entries:", np.sum(K))

        # Apply scaling matrix if using precomputed Jacobian
        if jacobian_sf is not None:
            scale_factors = np.load(jacobian_sf)
            if optimize_bc:
                # add (unit) scale factors for BCs
                # as the last 4 elements of the scaling matrix
                scale_factors = np.append(scale_factors, np.ones(4))
            reps = K.shape[0]
            scaling_matrix = np.tile(scale_factors, (reps, 1))
            if optimize_oh:
                if is_Regional:
                    K[:, :-1] *= scaling_matrix
                else:
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

    # Build either a full precomputed prior covariance or the original diagonal form.
    Sa, Sa_constraint, use_full_prior_covariance = build_prior_covariance(
        n_elements,
        prior_err,
        OptimizeSoil=OptimizeSoil,
        prior_ds=prior_ds,
        StateVectorFile=StateVectorFile,
        prebuilt_prior_err_covariance=prebuilt_prior_err_covariance,
    )

    # Number of elements to apply scale factor to
    scale_factor_idx = n_elements

    # If optimizing OH, adjust for it in the inversion
    if optimize_oh:
        # Add prior error for OH as the last element(s) of the diagonal
        # Following Masakkers et al. (2019, ACP) weight the OH term by the
        # ratio of the number of elements (n_OH_elements/n_emission_elements)
        # use this weighted constraint matrix to calculate the solution only
        apply_oh_prior(Sa, Sa_constraint, n_elements, prior_err_oh, is_Regional)
        scale_factor_idx -= 1 if is_Regional else 2

    # If optimizing boundary conditions, adjust for it in the inversion
    if optimize_bc:
        scale_factor_idx -= 4
        apply_bc_prior(Sa, Sa_constraint, prior_err_bc, optimize_oh, is_Regional)

    # The inversion uses the weighted constraint prior for the state estimate and the
    # unweighted prior for posterior diagnostics such as S_post and the averaging kernel.
    inv_Sa_constraint, inv_Sa = invert_prior_covariance(
        Sa, Sa_constraint, use_full_prior_covariance
    )

    # Solve for posterior scale factors xhat using the weighted constraint matrix
    delta_optimized = np.linalg.inv(gamma * KTinvSoK + inv_Sa_constraint) @ (
        gamma * KTinvSoyKxA
    )

    # Update scale factors by 1 to match what GEOS-Chem expects
    # xhat = 1 + delta_optimized
    # Notes:
    #  - If optimizing BCs, the last 4 elements are in concentration space,
    #    so we do not need to add 1
    #  - If optimizing OH, the last element also needs to be updated by 1
    xhat = delta_optimized.copy()
    xhat[:scale_factor_idx] += 1
    if optimize_oh:
        if is_Regional:
            xhat[-1] += 1
            print(f"xhat[OH] = {xhat[-1]}")
        else:
            xhat[-2:] += 1
            print(f"xhat[OH] = {xhat[-2:]}")

    # Posterior error covariance matrix (use unweighted Sa)
    S_post = np.linalg.inv(gamma * KTinvSoK + inv_Sa)

    # Averaging kernel matrix (use unweighted Sa)
    A = np.identity(n_elements) - S_post @ inv_Sa

    # Calculate J_A, where delta_optimized = xhat - xA
    # J_A = (xhat - xA)^T * inv_Sa * (xhat - xA)
    delta_optimizedT = delta_optimized.transpose()
    J_A = delta_optimizedT @ inv_Sa @ delta_optimized
    Ja_normalized = J_A / n_elements

    # Print some statistics
    print(
        f"hyperparameters: (prior_err: {prior_err}, obs_err: {obs_err}, gamma: {gamma}, "
        + f"prior_err_bc: {prior_err_bc}, prior_err_oh: {prior_err_oh})"
    )
    print(f"Normalized J_A: {Ja_normalized}")  # ideal gamma is where this is close to 1
    print(
        "Min:",
        xhat[:scale_factor_idx].min(),
        "Mean:",
        xhat[:scale_factor_idx].mean(),
        "Max",
        xhat[:scale_factor_idx].max(),
    )

    return xhat, delta_optimized, KTinvSoK, KTinvSoyKxA, S_post, A, Ja_normalized


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
    OptimizeSoil=False,
    prior_ds=None,
    StateVectorFile=None,
    prebuilt_prior_err_covariance=False,
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
        "Ja_normalized": [],
        "prior_err": [],
        "obs_err": [],
        "gamma": [],
        "prior_err_bc": [],
        "prior_err_oh": [],
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
        xhat, delta_optimized, KTinvSoK, KTinvSoyKxA, S_post, A, Ja_normalized = (
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
                OptimizeSoil,
                prior_ds,
                StateVectorFile,
                verbose=False,
                prebuilt_prior_err_covariance=prebuilt_prior_err_covariance,
            )
        )
        results_dict["KTinvSoK"].append(KTinvSoK)
        results_dict["KTinvSoyKxA"].append(KTinvSoyKxA)
        results_dict["ratio"].append(delta_optimized)
        results_dict["xhat"].append(xhat)
        results_dict["S_post"].append(S_post)
        results_dict["A"].append(A)
        results_dict["Ja_normalized"].append(Ja_normalized)
        for k, v in params.items():
            results_dict[k].append(v)

    # Find the ensemble member that is closest to 1 following Lu et al. (2021)
    idx_default_Ja = np.argmin(np.abs(np.array(results_dict["Ja_normalized"]) - 1))
    print(
        f"J_A/n closest to 1: {results_dict['Ja_normalized'][idx_default_Ja]} with"
        + f" (prior_err, obs_err, gamma, prior_err_bc, prior_err_oh) = {hyperparam_ensemble[idx_default_Ja]}"
    )

    # Filter ensemble members to only members with Ja between 0.5 and 2
    filter_ens_members = True  # set to False to turn off filtering
    include_ens_members = [
        i for i, Ja in enumerate(results_dict["Ja_normalized"]) if 0.5 <= Ja <= 2.0
    ]
    if filter_ens_members and len(include_ens_members) > 0:
        for k in results_dict.keys():
            results_dict[k] = [results_dict[k][i] for i in include_ens_members]
        # map idx_default_Ja to new filtered index
        idx_default_Ja = include_ens_members.index(idx_default_Ja)
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

    # Create an xarray.Dataset
    dataset = xr.Dataset()
    for k, v in results_dict.items():
        v = np.array(v)
        dims = ["ensemble"] + [f"nvar{i}" for i in range(1, v.ndim)]
        dataset[k] = (dims, v)

    # ensemble dimension to end
    dataset = dataset.transpose(..., "ensemble")

    # Specify attributes
    dataset.xhat.attrs["long_name"] = "Posterior scaling factors"
    dataset.xhat.attrs["units"] = "1"
    dataset.S_post.attrs["long_name"] = "Posterior error covariance matrix"
    dataset.S_post.attrs["units"] = "1"
    dataset.A.attrs["long_name"] = "Averaging kernel matrix"
    dataset.A.attrs["units"] = "1"
    dataset.Ja_normalized.attrs["long_name"] = "Normalized cost function Ja/n"
    dataset.Ja_normalized.attrs["units"] = "1"
    dataset.prior_err.attrs["long_name"] = "Prior error (Sa)"
    dataset.prior_err.attrs["units"] = "1"
    dataset.obs_err.attrs["long_name"] = "Observation error (So)"
    dataset.obs_err.attrs["units"] = "ppb"
    dataset.gamma.attrs["long_name"] = "Regularization parameter"
    dataset.gamma.attrs["units"] = "1"
    dataset.prior_err_bc.attrs["long_name"] = "Prior error for BC elements"
    dataset.prior_err_bc.attrs["units"] = "ppb"
    dataset.prior_err_oh.attrs["long_name"] = "Prior error for OH elements"
    dataset.prior_err_oh.attrs["units"] = "1"
    dataset.KTinvSoK.attrs["long_name"] = "K^T * inv(So) * K expression from inversion equation"
    dataset.KTinvSoK.attrs["units"] = "1"
    dataset.KTinvSoyKxA.attrs["long_name"] = "K^T * inv(So) * (y-K*xA) expression from inversion equation"
    dataset.KTinvSoyKxA.attrs["units"] = "1"
    dataset.ratio.attrs["long_name"] = "Change from prior (xhat - xA)"
    dataset.ratio.attrs["units"] = "1"

    # Calculate the mean of the ensemble as the main result
    dataset_mean = dataset.mean(dim="ensemble")

    dataset_mean.xhat.attrs["long_name"] = "Posterior scaling factors"
    dataset_mean.xhat.attrs["units"] = "1"
    dataset_mean.S_post.attrs["long_name"] = "Posterior error covariance matrix"
    dataset_mean.S_post.attrs["units"] = "1"
    dataset_mean.A.attrs["long_name"] = "Averaging kernel matrix"
    dataset_mean.A.attrs["units"] = "1"
    dataset_mean.Ja_normalized.attrs["long_name"] = "Normalized cost function Ja/n"
    dataset_mean.Ja_normalized.attrs["units"] = "1"
    dataset_mean.prior_err.attrs["long_name"] = "Prior error (Sa)"
    dataset_mean.prior_err.attrs["units"] = "1"
    dataset_mean.obs_err.attrs["long_name"] = "Observation error (So)"
    dataset_mean.obs_err.attrs["units"] = "ppb"
    dataset_mean.gamma.attrs["long_name"] = "Regularization parameter"
    dataset_mean.gamma.attrs["units"] = "1"
    dataset_mean.prior_err_bc.attrs["long_name"] = "Prior error for BC elements"
    dataset_mean.prior_err_bc.attrs["units"] = "ppb"
    dataset_mean.prior_err_oh.attrs["long_name"] = "Prior error for OH elements"
    dataset_mean.prior_err_oh.attrs["units"] = "1"
    dataset_mean.KTinvSoK.attrs["long_name"] = "K^T * inv(So) * K expression from inversion equation"
    dataset_mean.KTinvSoK.attrs["units"] = "1"
    dataset_mean.KTinvSoyKxA.attrs["long_name"] = "K^T * inv(So) * (y-K*xA) expression from inversion equation"
    dataset_mean.KTinvSoyKxA.attrs["units"] = "1"
    dataset_mean.ratio.attrs["long_name"] = "Change from prior (xhat - xA)"
    dataset_mean.ratio.attrs["units"] = "1"

    return dataset, dataset_mean


if __name__ == "__main__":
    import sys
    import os

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
    StateVectorFile = sys.argv[11]

    # read in config file
    config = load_config(config_path)

    # set parameters based on config file
    is_Regional = config["isRegional"]
    prior_err = ensure_float_list(config["PriorError"])
    obs_err = ensure_float_list(config["ObsError"])
    gamma = ensure_float_list(config["Gamma"])

    # 0.0 if not optimizing BCs or OH
    prior_err_BC = config["PriorErrorBCs"] if config["OptimizeBCs"] else 0.0
    prior_err_OH = config["PriorErrorOH"] if config["OptimizeOH"] else 0.0
    prior_err_BC = ensure_float_list(prior_err_BC)
    prior_err_OH = ensure_float_list(prior_err_OH)
    prebuilt_prior_err_covariance = config["OffDiagonalPriorCov"]
    
    OptimizeSoil = config["OptimizeSoil"]
    if OptimizeSoil:
        # prior emissions
        prior_cache = f"{os.path.expandvars(config['OutputPath']) }/{config['RunName']}/hemco_prior_emis/OutputDir/"
        start_date = config["StartDate"]
        end_date = config["EndDate"]
        prior_ds = get_mean_emissions(start_date, end_date, prior_cache)
    else:
        prior_ds = None

    # Reformat Jacobian scale factor input
    if jacobian_sf == "None":
        jacobian_sf = None

    # Run the inversion code
    out_ds, out_ds_mean = do_inversion_ensemble(
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
        OptimizeSoil,
        prior_ds,
        StateVectorFile,
        prebuilt_prior_err_covariance,
    )

    # Save the results of the ensemble inversion
    out_ds.to_netcdf(
        output_path.replace(".nc", "_ensemble.nc"),
        encoding={v: {"zlib": True, "complevel": 1} for v in out_ds.data_vars},
    )

    out_ds_mean.to_netcdf(
        output_path,
        encoding={v: {"zlib": True, "complevel": 1} for v in out_ds_mean.data_vars},
    )
    print(f"Saved results to {output_path}")

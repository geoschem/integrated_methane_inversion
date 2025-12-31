#!/usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import yaml
import re
import datetime
import numpy as np
import xarray as xr
from itertools import product
from scipy.sparse import hstack
from src.inversion_scripts.utils import (
    load_obj,
    save_obj,
    calculate_superobservation_error,
    ensure_float_list,
    get_mean_emissions,
    get_period_mean_emissions,
    update_prior_error_for_OptimizeSoil,
)
from src.inversion_scripts.regrid_precomputed_jacobian import(
    get_regrid_weights_jacobian_row,
    get_regrid_weights_jacobian_col,
    regrid_jacobian_row_col,
    median_and_sort_along_statevector,
)

def do_inversion(
    config,
    n_elements,
    jacobian_files,
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
    period_number=1,
    verbose=False,
):
    """
    After running jacobian.py, use this script to perform the inversion and save out results.

    Arguments
        n_elements   [int]   : Number of state vector elements
        jacobian_files [list]  : Jacobian files generated from jacobian.py
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
    # boolean for whether we are optimizing boundary conditions
    optimize_bc = config["OptimizeBCs"]
    optimize_oh = config["OptimizeOH"]

    # Need to ignore data in the GEOS-Chem 3 3 3 3 buffer zone
    # Shave off one or two degrees of latitude/longitude from each side of the domain
    # ~1 degree if 0.25x0.3125 resolution, ~2 degrees if 0.5x0.6125 resolution
    # This assumes 0.25x0.3125 and 0.5x0.625 simulations are always regional
    # initialize shaved off degrees
    degx = 0
    degy = 0
    if config['isRegional']:
        if "0.125x0.15625" in res:
            degx = 4 * 0.15625
            degy = 4 * 0.125
        elif "0.25x0.3125" in res:
            degx = 4 * 0.3125
            degy = 4 * 0.25
        elif "0.5x0.625" in res:
            degx = 4 * 0.625
            degy = 4 * 0.5
    xlim = [lon_min + degx, lon_max - degx]
    ylim = [lat_min + degy, lat_max - degy]

    # get prior emissions, destination area in square radians and all reference directories if using precomputed Jacobian
    if config['PrecomputedJacobian']:
        if (config['MultiPrecomputedJacobian']) and (config['OnlyEmisPrecomputedK']):
            # prior_ds of the destination grid
            RunDirs=f"{os.path.expandvars(config['OutputPath']) }/{config['RunName']}"
            
            # get start and end dates
            start_date = str(config["StartDate"])
            end_date = str(config["EndDate"])
            
            OptimizeSoil = config['OptimizeSoil']
            prior_cache = os.path.join(RunDirs, "hemco_prior_emis/OutputDir")
            if config["KalmanMode"]:
                prior_ds = get_period_mean_emissions(
                    prior_cache, period_number, os.path.join(RunDirs, "periods.csv")
                )
            else:
                prior_ds = get_mean_emissions(start_date, end_date, prior_cache)
            if not OptimizeSoil:
                prior_emis = prior_ds["EmisCH4_Total_ExclSoilAbs"].values
            else:
                prior_emis = prior_ds["EmisCH4_Total"].values
        
            ref_prior_emis_sv = []
            
            CSgridDir=f"{RunDirs}/CS_grids"
            ref_directory_path_all = config['ReferenceRunDir']
            with open(ref_directory_path_all, "r") as f:
                ref_directory_path_list_all = f.read()

            ref_directory_path_list = ref_directory_path_list_all.splitlines()
            ref_parent_dir = os.path.dirname(ref_directory_path_list[0])
            ref_dir_basename = os.path.basename(ref_directory_path_list[0])
            ref_dir_prename = re.sub(r'\d+$', '', ref_dir_basename)
            
            ref_sv_fpath = config['ReferenceStateVectorFile']
            if config['isRegional']:
                ref_sv_fname = os.path.basename(ref_sv_fpath).replace('combined', 'subset')
                ref_sv_fpath = f"{CSgridDir}/{ref_sv_fname}"
            # squeeze out time dimension
            ref_sv_ds = xr.open_dataset(ref_sv_fpath).squeeze("time")
            
            # # only use specific target face(s) for the inversion
            # # Find the index where target_face == 74
            # idx = np.where(ref_sv_ds['target_face'].values == 74)[0]
            # ref_sv_ds = ref_sv_ds.isel(target_face=idx)
            
            target_face = ref_sv_ds['target_face'].values
            num_ref_dir = len(target_face)
            print(f"Number for reference directories containing precomputed Jacobian: {num_ref_dir}")
            # get all reference prior emissions sorted by reference state vector
            for i in range(num_ref_dir):
                # reference directory naming is 1-based
                ref_RunName = f"{ref_dir_prename}{target_face[i]+1:03d}"
                ref_dir = os.path.join(ref_parent_dir, ref_RunName)
                ref_config_path = os.path.join(ref_dir, f"config_{ref_RunName}.yml")
                ref_config = yaml.load(open(ref_config_path), Loader=yaml.FullLoader)
                
                ref_OptimizeSoil = ref_config['OptimizeSoil']
                ref_prior_cache = os.path.join(ref_dir, "hemco_prior_emis/OutputDir")
                if ref_config["KalmanMode"]:
                    ref_prior_ds = get_period_mean_emissions(
                        ref_prior_cache, period_number, os.path.join(ref_dir, "periods.csv")
                    )
                else:
                    ref_prior_ds = get_mean_emissions(start_date, end_date, ref_prior_cache)
                if not ref_OptimizeSoil:
                    ref_prior_emis = ref_prior_ds["EmisCH4_Total_ExclSoilAbs"].values
                else:
                    ref_prior_emis = ref_prior_ds["EmisCH4_Total"].values
                ref_prior_emis = median_and_sort_along_statevector(ref_prior_emis, ref_sv_ds['StateVector'].values[i,...])
                ref_prior_emis_sv.append(ref_prior_emis)
            
    # Read output data from jacobian.py (virtual & true TROPOMI columns, Jacobian matrix)
    
    files = jacobian_files

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

        GC_index = dat['GC_index']
                        
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

        GC_index = GC_index[ind]

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
        if config['PrecomputedJacobian']:
            if (config['MultiPrecomputedJacobian']) and (config['OnlyEmisPrecomputedK']):
                # regrid Jacobian row (super observations) first, 
                # and then regrid Jacobian column (state vector)
                
                # regrid jacobian row
                jacobian_RegridRow = []
                overlap_area_jacobian_ratio_src = []
                overlap_area_src = []
                for ti in range(num_ref_dir):
                    # Sanity check:
                    # reference directory naming is 1-based
                    ref_RunName = f"{ref_dir_prename}{target_face[ti]+1:03d}"
                    ref_dir = os.path.join(ref_parent_dir, ref_RunName)
                    ref_config_path = os.path.join(ref_dir, f"config_{ref_RunName}.yml")
                    ref_config = yaml.load(open(ref_config_path), Loader=yaml.FullLoader)
                    assert (
                        np.isclose(ref_sv_ds['TARGET_LAT'].values[ti], ref_config['TARGET_LAT']) and
                        np.isclose(ref_sv_ds['TARGET_LON'].values[ti], ref_config['TARGET_LON'])
                    ), \
                        f"The TARGET_LAT/LON in the reference directory is not consistent with the reference config \
                        for {ref_RunName} with target_face id of {target_face[ti]+1} \
                        (TARGET_LAT: {ref_sv_ds['TARGET_LAT'].values[ti]:.2f} vs. {ref_config['TARGET_LAT']:.2f}, \
                        TARGET_LON: {ref_sv_ds['TARGET_LON'].values[ti]:.2f} vs. {ref_config['TARGET_LON']:.2f})"
                    
                    jacobian_ref_path = os.path.join(ref_dir, "inversion", "data_converted", os.path.basename(fi))
                    jacobian_ref_dat = load_obj(jacobian_ref_path)
                    
                    # assume the reference Jacobian is global
                    assert not (ref_config['OptimizeBCs'] or ref_config['isRegional']), \
                        "The reference precomputed Jacobian (stretched GCHP with ensemble target faces) \
                            must be global (not BC-optimized or regional)"
                    # regrid jacobian row
                    if (ref_config['OptimizeOH']):
                        # discard Jacobian column(s) for optimizing OH
                        jacobian_ref = jacobian_ref_dat["K"][:,:-2]
                    else:
                        jacobian_ref = jacobian_ref_dat["K"]
                    ref_GC_index = jacobian_ref_dat["GC_index"]
                    # GC_index is already subsetted to the region and 
                    # thus later for K_emis, we do not need to use [ind,:] to subset anymore
                    jacobian_regridding_weights_row = get_regrid_weights_jacobian_row(config, RunDirs, GC_index, ref_config, ref_GC_index)
                    jacobian_RegridRow_temp = jacobian_regridding_weights_row.dot(jacobian_ref).astype('float32')
                    jacobian_RegridRow.append(jacobian_RegridRow_temp)
                    # get inputs needed for regridding jacobian col
                    overlap_area_jacobian_ratio_src_temp, overlap_area_src_temp = get_regrid_weights_jacobian_col(config, RunDirs, ref_config, ref_dir, ref_prior_emis_sv[ti])
                    overlap_area_src.append(overlap_area_src_temp)
                    overlap_area_jacobian_ratio_src.append(overlap_area_jacobian_ratio_src_temp)

                # dense matrix in shape of (n_dst_valid_superobs, n_ref_sv_total)
                jacobian_RegridRow = np.concatenate(jacobian_RegridRow, axis=1)
                # sparse matrix in shape of (n_dst_sv, n_ref_sv_total)
                overlap_area_jacobian_ratio_src = hstack(overlap_area_jacobian_ratio_src, format="csr")
                # sparse matrix in shape of (n_dst_sv, n_ref_sv_total)
                overlap_area_src = hstack(overlap_area_src, format="csr")

                # regrid jacobian column
                # the jacoban is reconciled by the jacobian ratio and then 
                # get the area-weighted mean over reference state vector elements that overlap with the destination state vector
                # regrid_jacobian_row_col = sum ( jacobian_RegridRow * overlap_area_jacobian_ratio_src ) / sum(overlap_area_src)
                K_emis = 1e9 * regrid_jacobian_row_col(
                    jacobian_RegridRow, overlap_area_jacobian_ratio_src, overlap_area_src, RunDirs, config, prior_emis
                )

                if config["OptimizeBCs"] or config["OptimizeOH"]:
                    K_noemis = 1e9 * dat["K_noEmis"][ind, :]
                    K = np.concatenate((K_emis, K_noemis), axis=1)
                else:
                    K = K_emis
                
                if config.get('SaveRegriddedK', False):
                    regridded_K_path = fi.replace('data_converted', 'data_converted_regridded')
                    os.makedirs(os.path.dirname(regridded_K_path), exist_ok=True)
                    if not os.path.exists(regridded_K_path):
                        output = {}
                        output['GC_index'] = GC_index
                        output['K'] = K/1e9 # convert back to unitless mixing ratio
                        output['obs_GC'] = obs_GC
                        print(f"Saved regridded K_emis to {regridded_K_path}")
                        save_obj(output, regridded_K_path)
                        del output
                if config.get('SaveRegriddedRowK', False):
                    regridded_K_path = fi.replace('data_converted', 'data_converted_regriddedRow')
                    os.makedirs(os.path.dirname(regridded_K_path), exist_ok=True)
                    if not os.path.exists(regridded_K_path):
                        output = {}
                        output['GC_index'] = GC_index
                        output['K'] = jacobian_RegridRow # convert back to unitless mixing ratio
                        output['obs_GC'] = obs_GC
                        print(f"Saved regridded K_emis to {regridded_K_path}")
                        save_obj(output, regridded_K_path)
                        del output
            else:
                # Get Jacobian from reference inversion
                fi_ref = fi.replace("data_converted", "data_converted_reference")
                dat_ref = load_obj(fi_ref)
                K = 1e9 * dat_ref["K"][ind, :]
                
                # Apply scaling matrix if using precomputed Jacobian
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
        else:
            K = 1e9 * dat["K"][ind, :]

        # Number of observations
        if verbose:
            print("Sum of Jacobian entries:", np.sum(K))

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
    if OptimizeSoil:
        prior_err_new = update_prior_error_for_OptimizeSoil(prior_ds, prior_err, StateVectorFile, n_elements)
        Sa_diag = prior_err_new**2
    else:
        Sa_diag.fill(prior_err**2)
    Sa_diag_constraint = Sa_diag.copy() # constraint matrix to calculate the solution only

    # Number of elements to apply scale factor to
    scale_factor_idx = n_elements

    # If optimizing OH, adjust for it in the inversion
    if optimize_oh:
        # Add prior error for OH as the last element(s) of the diagonal
        # Following Masakkers et al. (2019, ACP) weight the OH term by the
        # ratio of the number of elements (n_OH_elements/n_emission_elements)
        # use this weighted constraint matrix to calculate the solution only
        if is_Regional:
            OH_weight = 1 / (n_elements - 1)
            Sa_diag_constraint[-1:] = OH_weight * prior_err_oh**2
            Sa_diag[-1:] = prior_err_oh**2
            scale_factor_idx -= 1
        else:
            OH_weight = 2 / (n_elements - 2)
            Sa_diag_constraint[-2:] = OH_weight * prior_err_oh**2 # weighted constraint matrix
            Sa_diag[-2:] = prior_err_oh**2 # unweighted matrix
            scale_factor_idx -= 2

    # If optimizing boundary conditions, adjust for it in the inversion
    if optimize_bc:
        scale_factor_idx -= 4

        # add prior error for BCs as the last 4 elements of the diagonal to both the unweighted (Sa_diag)
        # and weighted constraint (Sa_diag)
        if optimize_oh:
            if is_Regional:
                Sa_diag[-5:-1] = prior_err_bc**2
                Sa_diag_constraint[-5:-1] = prior_err_bc**2
            else:
                Sa_diag[-6:-1] = prior_err_bc**2
                Sa_diag_constraint[-6:-2] = prior_err_bc**2
        else:
            Sa_diag[-4:] = prior_err_bc**2
            Sa_diag_constraint[-4:] = prior_err_bc**2
    
    inv_Sa_constraint = np.diag(
        1 / Sa_diag_constraint
    )  # Inverse of weighted constraint matrix
    inv_Sa = np.diag(1 / Sa_diag)  # Inverse of unweighted prior error covariance matrix

    # Solve for posterior scale factors xhat using the weighted constraint matrix
    delta_optimized = np.linalg.inv(gamma * KTinvSoK + inv_Sa_constraint) @ (
        gamma * KTinvSoyKxA
    )
    # Posterior error covariance matrix (use unweighted Sa)
    S_post = np.linalg.inv(gamma * KTinvSoK + inv_Sa)
    
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
    config,
    n_elements,
    jacobian_files,
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
    period_number=1,
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
                config,
                n_elements,
                jacobian_files,
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
                period_number,
                verbose=False,
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

    # Create an xarray.Dataset
    dataset = xr.Dataset()
    for k, v in results_dict.items():
        v = np.array(v)
        dims = ["ensemble"] + [f"nvar{i}" for i in range(1, v.ndim)]
        dataset[k] = (dims, v)

    # save index number of ens member with J_A/n
    # closest to 1 as the default member
    dataset.attrs = {"default_member_index": idx_default_Ja}

    # ensemble dimension to end
    dataset = dataset.transpose(..., "ensemble")

    dataset_default = dataset.isel(ensemble=idx_default_Ja)

    return dataset, dataset_default


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
    period_number = sys.argv[12]

    # read in config file
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

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

    # make it robust to read output files in the defined time range
    # Get TROPOMI data filenames for the desired date range
    gc_startdate = np.datetime64(datetime.datetime.strptime(str(config['StartDate']), "%Y%m%d"))
    gc_enddate = np.datetime64(datetime.datetime.strptime(str(config['EndDate']), "%Y%m%d"))
    allfiles = glob.glob(f"{jacobian_dir}/*.pkl")
    jacobian_files = []
    for index in range(len(allfiles)):
        filename = allfiles[index]
        shortname = re.split(r"\/", filename)[-1]
        shortname = re.split(r"\.", shortname)[0]
        strdate = re.split(r"\.|_+|T", shortname)[4]
        strdate = datetime.datetime.strptime(strdate, "%Y%m%d")
        if (strdate >= gc_startdate) and (strdate < gc_enddate):
            jacobian_files.append(filename)
    jacobian_files.sort()
    
    # Run the inversion code
    out_ds, out_ds_default = do_inversion_ensemble(
        config,
        n_elements,
        jacobian_files,
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
        period_number,
    )

    # add atributes for stretching GCHP simulation
    if config['STRETCH_GRID']:
        out_ds.attrs['STRETCH_FACTOR'] = np.float32(config['STRETCH_FACTOR'])
        out_ds.attrs['TARGET_LAT'] = np.float32(config['TARGET_LAT'])
        out_ds.attrs['TARGET_LON'] = np.float32(config['TARGET_LON'])
        
        out_ds_default.attrs['STRETCH_FACTOR'] = np.float32(config['STRETCH_FACTOR'])
        out_ds_default.attrs['TARGET_LAT'] = np.float32(config['TARGET_LAT'])
        out_ds_default.attrs['TARGET_LON'] = np.float32(config['TARGET_LON'])
    # Save the results of the ensemble inversion
    out_ds.to_netcdf(
        output_path.replace(".nc", "_ensemble.nc"),
        encoding={v: {"zlib": True, "complevel": 1} for v in out_ds.data_vars},
    )

    out_ds_default.to_netcdf(
        output_path,
        encoding={v: {"zlib": True, "complevel": 1} for v in out_ds.data_vars},
    )
    print(f"Saved results to {output_path}")

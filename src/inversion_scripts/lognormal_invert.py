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
    OH_element_num = 1 if optimize_oh else 0
    BC_element_num = 4 if optimize_bcs else 0
    num_sv_elems = (
        int(state_vector_labels.max().item()) + OH_element_num + BC_element_num
    )
    num_buffer_elems = int(config["nBufferClusters"])
    num_normal_elems = num_buffer_elems + OH_element_num + BC_element_num
    ds = np.load("full_jacobian_K.npz")
    K_temp = np.array(ds["K"]) * 1e9

    # Apply scaling matrix if using precomputed Jacobian
    if jacobian_sf is not None:
        scale_factors = np.load(jacobian_sf)
        # apply unit scaling for BC elements and OH elements if using
        if optimize_bcs or optimize_oh:
            scale_factors = np.append(scale_factors, np.ones(BC_element_num + OH_element_num))
        reps = K_temp.shape[0]
        scaling_matrix = np.tile(scale_factors, (reps, 1))
        K_temp *= scaling_matrix

    # The levenberg-marquardt method assumes that the prior emissions is 
    # the median prior emissions, but typically priors are the mean emission.
    # To account for this we convert xa to a median. This can be done by 
    # scaling the lognormal part of K by 1/exp((lnsa**2)/2).
    # Here, we calculate this scaling factor
    prior_scale = 1/np.exp((np.log(float(config["PriorError"]))**2)/2)
    
    # split K based on whether we are solving for lognormal or normal elements
    # K_ROI is the matrix for the lognormal elements (the region of interest)
    # the lognormal part of K gets scaled by prior_scale to convert to median
    K_ROI = K_temp[:, :-num_normal_elems]
    K_ROI = prior_scale * K_ROI
    K_normal = K_temp[:, -num_normal_elems:]
    K_full = np.concatenate((K_ROI, K_normal), axis=1)

    # get the So matrix
    ds = np.load("so_super.npz")
    so = ds["so"]

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
    m, n = np.shape(K_ROI)

    # Create base xa and lnxa matrices
    # Note: the resulting xa vector has lognormal elements until the final bc elements
    xa = np.ones((n, 1)) * 1.0
    lnxa = np.log(xa)
    xa_normal = np.zeros((num_normal_elems, 1)) * 1.0
    xa = np.concatenate((xa, xa_normal), axis=0)
    lnxa = np.concatenate((lnxa, xa_normal), axis=0)

    # Create inverted So matrix
    Soinv = spdiags(1 / so, 0, m, m)
    
    # Define Sa, gamma, and Sa_bc values to iterate through
    prior_errors = [float(config["PriorError"])]
    sa_buffer_elems = [float(config["PriorErrorBufferElements"])]
    sa_bc_vals = [float(config["PriorErrorBCs"])] if optimize_bcs else [None]
    sa_oh_vals = [float(config["PriorErrorOH"])] if optimize_oh else [None]
    gamma_vals = [float(config["Gamma"])]

    # iterate through different combination of gamma, lnsa, and sa_bc
    # TODO: for now we will only allow one value for each of these
    # TODO: parallelize this once we allow vectorization of these values
    combinations = list(product(gamma_vals, prior_errors, sa_bc_vals, sa_buffer_elems, sa_oh_vals))
    for gamma, sa, sa_bc, sa_buffer, sa_oh in combinations:
        lnsa_val = np.log(sa)
        results_save_path = f"inversion_result_ln.nc"

        # Create lnSa matrix
        # lnsa = lnsa_val**2 * np.ones((n, 1))
        lnsa = lnsa_val**2 * np.ones((n, 1))

        
        # For the buffer elems, BCs, and OH elements
        # we apply a different Sa value
        # In the most basic we only generate Sa for buffer elements
        sa_normal = sa_buffer**2 * np.ones((num_normal_elems - (BC_element_num + OH_element_num), 1))
        
        # conditionally add BC and OH elements
        if optimize_bcs:
            bc_errors = sa_bc**2 * np.ones((BC_element_num, 1))
            sa_normal = np.concatenate((sa_normal, bc_errors), axis=0)   
        if optimize_oh:
            oh_errors = sa_oh**2 * np.ones((OH_element_num, 1))
            sa_normal = np.concatenate((sa_normal, oh_errors), axis=0)

        # concatenate lognormal prior errors with normal prior errors
        lnsa_arr = np.concatenate((lnsa, sa_normal), axis=0)

        # Create lnSa matrix
        lnsa = np.zeros((n + num_normal_elems, n + num_normal_elems))
        np.fill_diagonal(lnsa, lnsa_arr)
        invlnsa = np.linalg.inv(lnsa)

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
            K_prime = np.concatenate((K_ROI * xn[:-num_normal_elems].T, K_normal), axis=1)
            
            # commonly used term for term1 and term3
            gamma_K_prime_transpose_Soinv = gamma * K_prime.T @ Soinv
            
            # Compute the next xn_update
            term1 = gamma_K_prime_transpose_Soinv @ K_prime
            term2 = (1 + kappa) * invlnsa
            inv_term = np.linalg.inv(term1 + term2)

            term3 = gamma_K_prime_transpose_Soinv @ (y_ybkg_diff - K_full @ xn)
            term4 = invlnsa @ (lnxn - lnxa)

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
        lns = np.linalg.inv(K_primeT_so @ K_prime + invlnsa)
        G = lns @ K_primeT_so
        ak = G @ K_prime
        dofs = np.trace(ak)
        print(f"DOFS: {dofs}")

        # Calculate posterior mean xhat
        dlns = np.diag(lns[:-num_normal_elems, :-num_normal_elems])
        xnmean = np.concatenate(
            (
                xn[:-num_normal_elems] * np.expand_dims(np.exp(dlns * (0.5)) * prior_scale, axis=1),
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

        # Create gridded datarrays of S_post, xhat, and A
        xhat = xnmean
        print(f"xhat = {xhat.sum()}")
        xhat_arr = np.zeros((len(lats), len(lons)))
        ak_sensitivities = np.diagonal(ak)
        ak_arr = np.zeros((len(lats), len(lons)))
        lnS_post = np.diagonal(lns)
        lnS_post_arr = np.zeros((len(lats), len(lons)))
        for i in range(np.shape(xhat)[0]):
            idx = np.where(state_vector_labels == float(i + 1))
            xhat_arr[idx] = xhat[i]
            ak_arr[idx] = ak_sensitivities[i]
            lnS_post_arr[idx] = lnS_post[i]

        # lambda function to make a DataArray from numpy array
        make_dataArray = lambda arr: xr.DataArray(
            data=arr,
            dims=["lat", "lon"],
            coords=dict(
                lon=(["lon"], lons.values),
                lat=(["lat"], lats.values),
            ),
        )

        # transform to data arrays
        scale_factors = make_dataArray(xhat_arr)
        ak_arr = make_dataArray(ak_arr)
        lnS_post_arr = make_dataArray(lnS_post_arr)

        # create gridded posterior
        ds = xr.Dataset(
            {
                "ScaleFactor": (["lat", "lon"], scale_factors.data),
                "A": (["lat", "lon"], ak_arr.data),
                "S_post": (["lat", "lon"], lnS_post_arr.data),
            },
            coords={"lon": ("lon", lons.data), "lat": ("lat", lats.data)},
        )

        # Add attribute metadata
        ds.lat.attrs["units"] = "degrees_north"
        ds.lat.attrs["long_name"] = "Latitude"
        ds.lon.attrs["units"] = "degrees_east"
        ds.lon.attrs["long_name"] = "Longitude"
        ds.ScaleFactor.attrs["units"] = "1"

        # save to netcdf file
        ds.to_netcdf("gridded_posterior_ln.nc")

        # Save (ungridded) inversion results
        dataset = Dataset(results_save_path, "w", format="NETCDF4_CLASSIC")
        dataset.createDimension("nvar", num_sv_elems)
        dataset.createDimension("float", 1)
        nc_xn = dataset.createVariable("xhat", np.float32, ("nvar"))
        nc_lnxn = dataset.createVariable("lnxn", np.float32, ("nvar"))
        nc_lnS_post = dataset.createVariable("S_post", np.float32, ("nvar", "nvar"))
        nc_A = dataset.createVariable("A", np.float32, ("nvar", "nvar"))
        nc_dofs = dataset.createVariable("DOFS", np.float32, ("float",))
        nc_Ja = dataset.createVariable("Ja", np.float32, ("float",))
        nc_xn[:] = xnmean.flatten()
        nc_lnxn[:] = lnxn.flatten()
        nc_lnS_post[:, :] = lns
        nc_A[:, :] = ak
        nc_dofs[0] = dofs
        nc_Ja[0] = Ja
        dataset.close()


if __name__ == "__main__":
    config_path = sys.argv[1]
    state_vector_filepath = sys.argv[2]
    jacobian_sf = None if sys.argv[3] == "None" else sys.argv[3]

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    lognormal_invert(config, state_vector_filepath, jacobian_sf)

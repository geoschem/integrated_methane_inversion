# Description: Script to perform inversion using lognormal errors
# Usage: python lognormal_invert.py <path_to_config_file> <path_to_state_vector_file>

# TODO: merge this script with invert.py to avoid redundancy
# This script performs the inversion but using lognormal
# errors instead of normal errors. As an alternative to
# the invert.py script.

import sys
import yaml
from itertools import product
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from scipy.sparse import spdiags


def lognormal_invert(config, state_vector_filepath):
    """
    Description:
        Run inversion using lognormal errors following method from eqn 2 of
        Chen et al., 2022 https://doi.org/10.5194/acp-22-10809-2022
        Outputs inversion results to netcdf files.
    Arguments:
        state_vector_filepath [String] : path to state vector netcdf file
        config                [Dict]   : dictionary of config variables
    """
    state_vector = xr.load_dataset(state_vector_filepath)
    state_vector_labels = state_vector["StateVector"]
    lats, lons = state_vector_labels.lat, state_vector_labels.lon

    # used to determine convergence of xn 5e-3 is .5%
    convergence_threshold = 5e-3

    # Load in the observation and background data
    ds = np.load("obs_ch4_tropomi.npz")
    y = np.asmatrix(ds["obs_tropomi"])
    ds = np.load("gc_ch4_bkgd.npz")
    ybkg = np.asmatrix(ds["gc_ch4_bkgd"])

    # We only solve using lognormal errors for state vector elements
    # within the domain of interest, not the buffer elements, or the
    # BC elements. So, to do this we split K into two matrices, one
    # for the lognormal elements, and one for the normal elements.
    optimize_bcs = config["OptimizeBCs"]
    num_sv_elems = (
        int(state_vector_labels.max().item()) + 4
        if optimize_bcs
        else state_vector_labels.max().item()
    )
    num_buffer_elems = int(config["nBufferClusters"])
    num_normal_elems = num_buffer_elems + 4 if optimize_bcs else num_buffer_elems
    ds = np.load("full_jacobian_K.npz")
    K_lognormal = np.asmatrix(ds["K"][:, :-num_normal_elems]) * 1e9
    K_normal = np.asmatrix(ds["K"][:, -num_normal_elems:]) * 1e9

    # get the So matrix
    ds = np.load("so_super.npz")
    so = ds["so"] ** 2

    # Calculate the difference between tropomi and the background
    # simulation, which has no emissions
    y_ybkg_diff = y - ybkg
    ybkg, y, y_ybkg_diff = (
        np.swapaxes(ybkg, 0, 1),
        np.swapaxes(y, 0, 1),
        np.swapaxes(y_ybkg_diff, 0, 1),
    )
    K_full = np.concatenate((K_lognormal, K_normal), axis=1)

    # fixed kappa of 10 following Chen et al., 2022 https://doi.org/10.5194/acp-22-10809-2022
    kappa = 10
    m, n = np.shape(K_lognormal)

    # Create base xa and lnxa matrices
    # Note: the resulting xa matrix has lognormal elements until the final bc elements
    xa = np.ones((n, 1)) * 1.0
    lnxa = np.asmatrix(np.log(xa))
    xa_bcs = np.zeros((num_normal_elems, 1)) * 1.0
    xa = np.asmatrix(np.concatenate((xa, xa_bcs), axis=0))
    lnxa = np.asmatrix(np.concatenate((lnxa, xa_bcs), axis=0))

    # Create inverted So matrix
    Soinv = spdiags(1 / so, 0, m, m)

    # Define Sa, gamma, and Sa_bc values to iterate through
    prior_errors = [float(config["PriorError"])]
    sa_buffer_elems = [float(config["PriorErrorBufferElements"])]
    sa_bc_vals = [float(config["PriorErrorBCs"])] if optimize_bcs else [None]
    gamma_vals = [float(config["Gamma"])]

    # iterate through different combination of gamma, lnsa, and sa_bc
    # TODO: for now we will only allow one value for each of these
    # TODO: parallelize this once we allow vectorization of these values
    combinations = list(product(gamma_vals, prior_errors, sa_bc_vals, sa_buffer_elems))
    for gamma, sa, sa_bc, sa_buffer in combinations:
        lnsa_val = np.log(sa)
        results_save_path = f"inversion_result_ln.nc"

        invso_over_gamma = 1 / (so / gamma)
        So_over_gamma_inv = spdiags(invso_over_gamma, 0, m, m)

        lnsa = lnsa_val**2 * np.ones((n, 1))

        # For the buffer elems and BCs we apply a different Sa value
        if optimize_bcs:
            bc_errors = sa_bc**2 * np.ones((4, 1))
            buffer_errors = sa_buffer**2 * np.ones((num_normal_elems - 4, 1))
            sa_normal = np.concatenate((buffer_errors, bc_errors), axis=0)
        else:
            sa_normal = sa_buffer**2 * np.ones((num_normal_elems, 1))

        # concatenate lognormal prior errors with normal prior errors
        lnsa_arr = np.concatenate((lnsa, sa_normal), axis=0)

        # Create lnSa matrix
        lnsa = np.zeros((n + num_normal_elems, n + num_normal_elems))
        np.fill_diagonal(lnsa, lnsa_arr)
        invlnsa = np.linalg.inv(lnsa)

        # we start with xa and lnxa using the prior values (scale factors of 1)
        xn = xa
        lnxn = lnxa

        lnk = np.concatenate(
            (
                np.multiply(K_lognormal, np.transpose(xn[:-num_normal_elems])),
                K_normal,
            ),
            axis=1,
        )

        # We decompose eqn 2 into 3 terms
        # term 1: gamma*ln(K).T@inv(So)@ln(K) + inv((1+kappa)*inv(ln(sa)))
        # term 2: gamma*ln(K).T@inv(So)@(y_ybkg_diff - K@xa)
        # term 3: -inv(ln(sa))@(ln(xn) - ln(xa))
        # We can then solve for xn iteratively by doing:
        # xn = x(n-1) + term1@(term2 + term3)
        # where x(n-1) is the previous iteration of xn until convergence
        gamma_lnk_transpose_Soinv = (
            gamma * np.transpose(lnk) @ Soinv
        )  # Used in term1 and term2

        term1 = np.linalg.inv(
            gamma_lnk_transpose_Soinv @ lnk + np.multiply((1 + kappa), invlnsa)
        )

        term2 = gamma_lnk_transpose_Soinv @ (y_ybkg_diff - (K_full @ xn))

        term3 = -1 * invlnsa @ (lnxn - lnxa)

        lnxn_update = lnxn + term1 @ (term2 + term3)

        xn_iteration_pct_diff = max(
            abs(
                np.exp(lnxn_update[:-num_normal_elems])
                - np.exp(lnxn[:-num_normal_elems])
            )
            / np.exp(lnxn[:-num_normal_elems])
        )

        # Iterate for calculation of ln(xn) until convergence threshold is met (5e-3)
        print("Status: Iterating to calculate ln(xn)")
        while xn_iteration_pct_diff >= convergence_threshold:
            lnxn = lnxn_update
            xn = np.concatenate(
                (np.exp(lnxn[:-num_normal_elems]), lnxn[-num_normal_elems:]),
                axis=0,
            )
            lnk = np.concatenate(
                (
                    np.multiply(K_lognormal, np.transpose(xn[:-num_normal_elems])),
                    K_normal,
                ),
                axis=1,
            )
            gamma_lnk_transpose_Soinv = gamma * np.transpose(lnk) @ Soinv  #
            term1 = np.linalg.inv(
                gamma_lnk_transpose_Soinv @ lnk + (1 + kappa) * invlnsa
            )
            term2 = gamma_lnk_transpose_Soinv @ (y_ybkg_diff - (K_full @ xn))
            term3 = -1 * invlnsa @ (lnxn - lnxa)
            lnxn_update = lnxn + term1 @ (term2 + term3)
            xn_iteration_pct_diff = max(  # percent diff between xn and xn_update
                abs(
                    np.exp(lnxn_update[:-num_normal_elems])
                    - np.exp(lnxn[:-num_normal_elems])
                )
                / np.exp(lnxn[:-num_normal_elems])
            )

        print("Status: Done Iterating")
        lnxn = lnxn_update
        xn = np.concatenate(
            (np.exp(lnxn[:-num_normal_elems]), lnxn[-num_normal_elems:]), axis=0
        )
        lnk = np.concatenate(
            (
                np.multiply(K_lognormal, np.transpose(xn[:-num_normal_elems])),
                K_normal,
            ),
            axis=1,
        )
        kso = np.transpose(lnk) @ So_over_gamma_inv
        lns = np.linalg.inv(kso @ lnk + invlnsa)
        G = lns @ kso
        ak = G @ lnk
        dofs = np.trace(ak)

        dlns = np.diag(lns[:-num_normal_elems, :-num_normal_elems])
        xnmean = np.concatenate(
            (
                np.multiply(
                    xn[:-num_normal_elems],
                    np.expand_dims(np.exp(dlns * (0.5)), axis=1),
                ),
                xn[-num_normal_elems:],
            )
        )

        Ja = (
            np.transpose(lnxn[:-num_normal_elems] - lnxa[:-num_normal_elems])
            @ invlnsa[:-num_normal_elems, :-num_normal_elems]
            @ (lnxn[:-num_normal_elems] - lnxa[:-num_normal_elems])
        )

        print(f"Diagnostics:\n  (Ja: {Ja}, gamma: {gamma}, sa: {sa}, sa_bc: {sa_bc}, sa_buffer: {sa_buffer})")

        # Diagnostic for Ja with BCs
        # Ja_with_BCs = np.transpose(lnxn - lnxa) @ invlnsa @ (lnxn - lnxa)
        # print("Ja with BCs", Ja_with_BCs)

        # Create gridded datarrays of S_post, xhat, and A
        xhat = xnmean
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

        # handy lambda function to make a dataArray from numpy array
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

        # save inversion results
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

        # Save out gridded posterior
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

        # Create netcdf
        ds.to_netcdf("gridded_posterior_ln.nc")


if __name__ == "__main__":
    config_path = sys.argv[1]
    state_vector_filepath = sys.argv[2]
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    lognormal_invert(config, state_vector_filepath)
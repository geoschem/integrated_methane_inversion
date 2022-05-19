#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pickle
from utils import load_obj


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
        res          [str]   : '0.25x0.3125' or '0.5x0.625' -- from config.yml

    Returns
        xhat         [float] : Posterior scaling factors
        ratio        [float] : Change from prior     [xhat = 1 + ratio]
        KTinvSoK     [float] : K^T*inv(S_o)*K        [part of inversion equation]
        KTinvSoyKxA  [float] : K^T*inv(S_o)*(y-K*xA) [part of inversion equation]
        S_post       [float] : Posterior error covariance matrix
        A            [float] : Averaging kernel matrix

    """

    # Need to ignore data in the GEOS-Chem 3 3 3 3 buffer zone
    # Shave off one or two degrees of latitude/longitude from each side of the domain
    # ~1 degree if 0.25x0.3125 resolution, ~2 degrees if 0.5x0.6125 resolution
    if "0.25x0.3125" in res:
        degx = 4 * 0.3125
        degy = 4 * 0.25
    elif "0.5x0.625" in res:
        degx = 4 * 0.625
        degy = 4 * 0.5
    else:
        msg = "Bad input for res; must be '0.25x0.3125' or '0.5x0.625' "
        raise ValueError(msg)

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
    #        = xA + ratio
    #        = 1  + ratio      [since xA=1 when optimizing scale factors]
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
    
    # For each .pkl file generated by jacobian.py:
    for fi in files:

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
        # Note: weighting function defined by Zichong for his 
        # middle east inversions. May need to be tuned based on region.
        obsWeight = 1.02 - (obs_GC[:, 4] * .02)
        obs_error = obsWeight * obs_err
        
        # check to make sure obs_err isn't negative, set 1 as default value
        obs_error = [obs if obs > 0 else 1 for obs in obs_error]

        # Jacobian entries for observations within bounds [ppb]
        K = 1e9 * dat["K"][ind, :]

        # Number of observations
        print("Sum of Jacobian entries:", np.sum(K))

        # Define observational errors (diagonal entries of S_o matrix)
        obs_error = np.power(obs_error, 2)
        
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
    inv_Sa = np.diag(1 / Sa_diag)  # Inverse of prior error covariance matrix

    # Solve for posterior scale factors xhat
    ratio = np.linalg.inv(gamma * KTinvSoK + inv_Sa) @ (gamma * KTinvSoyKxA)
    xhat = 1 + ratio

    # Posterior error covariance matrix
    S_post = np.linalg.inv(gamma * KTinvSoK + inv_Sa)

    # Averaging kernel matrix
    A = np.identity(n_elements) - S_post @ inv_Sa

    return xhat, ratio, KTinvSoK, KTinvSoyKxA, S_post, A


if __name__ == "__main__":
    import sys

    n_elements = int(sys.argv[1])
    jacobian_dir = sys.argv[2]
    output_path = sys.argv[3]
    lon_min = float(sys.argv[4])
    lon_max = float(sys.argv[5])
    lat_min = float(sys.argv[6])
    lat_max = float(sys.argv[7])
    prior_err = float(sys.argv[8])
    obs_err = float(sys.argv[9])
    gamma = float(sys.argv[10])
    res = sys.argv[11]

    # Run the inversion code
    out = do_inversion(
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
    )
    xhat = out[0]
    ratio = out[1]
    KTinvSoK = out[2]
    KTinvSoyKxA = out[3]
    S_post = out[4]
    A = out[5]

    # Print some statistics
    print("Min:", xhat.min(), "Mean:", xhat.mean(), "Max", xhat.max())

    # Save results
    dataset = Dataset(output_path, "w", format="NETCDF4_CLASSIC")
    nvar = dataset.createDimension("nvar", n_elements)
    nc_KTinvSoK = dataset.createVariable("KTinvSoK", np.float32, ("nvar", "nvar"))
    nc_KTinvSoyKxA = dataset.createVariable("KTinvSoyKxA", np.float32, ("nvar"))
    nc_ratio = dataset.createVariable("ratio", np.float32, ("nvar"))
    nc_xhat = dataset.createVariable("xhat", np.float32, ("nvar"))
    nc_S_post = dataset.createVariable("S_post", np.float32, ("nvar", "nvar"))
    nc_A = dataset.createVariable("A", np.float32, ("nvar", "nvar"))
    nc_KTinvSoK[:, :] = KTinvSoK
    nc_KTinvSoyKxA[:] = KTinvSoyKxA
    nc_ratio[:] = ratio
    nc_xhat[:] = xhat
    nc_S_post[:, :] = S_post
    nc_A[:, :] = A
    dataset.close()

    print(f"Saved results to {output_path}")

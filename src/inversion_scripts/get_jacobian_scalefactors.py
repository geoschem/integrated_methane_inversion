#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import xarray as xr


def get_jacobian_scalefactors(period_number, base_directory, sf_archive, ref_archive):
    """
    Running sensitivity inversions with pre-constructed Jacobian may require scaling the 
    Jacobian to match updated prior estimates. This is because the Jacobian is defined
    by a relative (50%) perturbation of emissions, so that the sensitivities depend on
    the magnitude of emissions in the prior inventory. 

    This function generates scale factors to perform that scaling.

    Arguments
        period_number  [int]   : Current inversion period, starting from 1
        base_directory [str]   : The base directory for the inversion, where e.g., "preview_sim/" resides
        sf_archive     [str]   : Path to archive of prior/posterior scale factors
        ref_archive    [str]   : Path to archive of scale factors from original inversion
    
    Returns
        sf_K           [float] : Scale factors to be applied to the Jacobian matrix
    """

    # Get target and reference scale factors (prior) for current period
    sf_path = os.path.join(sf_archive, f"prior_sf_period{period_number}.nc")
    sf_path_ref = os.path.join(ref_archive, f"prior_sf_period{period_number}.nc")
    sf = xr.load_dataset(sf_path)["ScaleFactor"]
    sf_ref = xr.load_dataset(sf_path_ref)["ScaleFactor"]

    # Get scale factors to apply to the Jacobian matrix K
    statevector_path = os.path.join(base_directory, "StateVector.nc")
    statevector = xr.load_dataset(statevector_path)["StateVector"]
    n_elements = int(np.nanmax(statevector.data))
    sf_ratio = sf / sf_ref
    sf_K = []
    for e in range(1, n_elements + 1):
        # Get scale factor for state vector element e
        sf_e = np.nanmean(sf_ratio.where(statevector == e).values)
        # Append
        sf_K.append(sf_e)
    sf_K = np.asarray(sf_K)

    # Reset buffer area to 1 # TODO Do we want this feature?
    # config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    # n_buff = config["nBufferClusters"]
    # sf_K[-n_buff:] = 1

    return sf_K


if __name__ == "__main__":
    import sys

    period_number = int(sys.argv[1])
    base_directory = sys.argv[2]
    sf_archive = sys.argv[3]
    ref_archive = sys.argv[4]

    # Get the scale factors
    out = get_jacobian_scalefactors(
        period_number, base_directory, sf_archive, ref_archive
    )

    # Save them
    save_path_1 = os.path.join(
        sf_archive, f"jacobian_scale_factors_period{period_number}.npy"
    )
    save_path_2 = os.path.join(
        base_directory,
        f"kf_inversions/period{period_number}",
        f"jacobian_scale_factors.npy",
    )
    np.save(save_path_1, out)
    np.save(save_path_2, out)

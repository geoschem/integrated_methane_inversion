#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import yaml
import numpy as np
import xarray as xr


def get_jacobian_scalefactors(period_number, inv_directory, ref_directory):
    """
    Running sensitivity inversions with pre-constructed Jacobian may require scaling the
    Jacobian to match updated prior estimates. This is because the Jacobian is defined
    by relative perturbations of emissions, so that the sensitivities depend on
    the magnitude of emissions in the prior inventory.

    This function generates scale factors to perform that scaling.

    Arguments
        period_number  [int]   : Current inversion period, starting from 1
        inv_directory  [str]   : The base directory for the inversion, where e.g., "preview/" resides
        ref_directory  [str]   : The base directory for the reference inversion

    Returns
        sf_K           [float] : Scale factors to be applied to the Jacobian matrix
    """

    # Get target and reference perturbation scale factors applied for jacobian calculation
    # These are the perturbations applied in K = (perturbed_sim - base_sim) / perturbation
    pert_sf_path = os.path.join(
        inv_directory, "archive_perturbation_sfs", f"pert_sf_{period_number}.npz"
    )
    ref_pert_sf_path = os.path.join(
        ref_directory, "archive_perturbation_sfs", f"pert_sf_{period_number}.npz"
    )
    pert_sf_dict = np.load(pert_sf_path)
    ref_pert_sf_dict = np.load(ref_pert_sf_path)
    pert_sf = pert_sf_dict["effective_pert_sf"]
    ref_pert_sf = ref_pert_sf_dict["effective_pert_sf"]

    # Get the ratio of the targetted emissions in the target and reference inversions
    target_emis_ratio = (
        pert_sf_dict["target_emission"] / ref_pert_sf_dict["target_emission"]
    )

    # Note 1: This line assumes the spatial dist of emissions is the same within each state
    # vector element. Be careful switching prior emission inventories (especially if grid
    # cells are clustered). For this to work we need to use the same state vector in both the
    # target and reference inversions.
    # Note 2: This also assumes that the temporal variabily of the swapped emissions is the same.
    # If the temporal variability is different, there will be error associated with scaling
    # the Jacobian.
    sf_K = pert_sf / ref_pert_sf

    # Apply the target_emis_ratio to the scale factors
    sf_K = sf_K * target_emis_ratio

    return np.asarray(sf_K)


if __name__ == "__main__":
    import sys

    period_number = int(sys.argv[1])
    inv_directory = sys.argv[2]
    ref_directory = sys.argv[3]
    kalman_mode = sys.argv[4].lower() == "true"

    # Get the scale factors
    out = get_jacobian_scalefactors(period_number, inv_directory, ref_directory)

    # Save them
    save_path_1 = os.path.join(
        os.path.join(inv_directory, "archive_sf"),
        f"jacobian_scale_factors_period{period_number}.npy",
    )
    path_prefix = f"kf_inversions/period{period_number}" if kalman_mode else "inversion"
    save_path_2 = os.path.join(
        inv_directory,
        path_prefix,
        f"jacobian_scale_factors.npy",
    )
    np.save(save_path_1, out)
    np.save(save_path_2, out)

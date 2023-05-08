#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import xarray as xr


def get_jacobian_scalefactors(period_number, inv_directory, ref_directory):
    """
    Running sensitivity inversions with pre-constructed Jacobian may require scaling the 
    Jacobian to match updated prior estimates. This is because the Jacobian is defined
    by a relative (50%) perturbation of emissions, so that the sensitivities depend on
    the magnitude of emissions in the prior inventory. 

    This function generates scale factors to perform that scaling.

    Arguments
        period_number  [int]   : Current inversion period, starting from 1
        inv_directory [str]   : The base directory for the inversion, where e.g., "preview_sim/" resides
        ref_directory  [str]   : The base directory for the reference inversion
    
    Returns
        sf_K           [float] : Scale factors to be applied to the Jacobian matrix
    """

    # Get target and reference scale factors (prior) for current period
    sf_archive = os.path.join(inv_directory, "archive_sf")
    ref_archive = os.path.join(ref_directory, "archive_sf")
    sf_path = os.path.join(sf_archive, f"prior_sf_period{period_number}.nc")
    sf_path_ref = os.path.join(ref_archive, f"prior_sf_period{period_number}.nc")
    sf = xr.load_dataset(sf_path)["ScaleFactor"]
    sf_ref = xr.load_dataset(sf_path_ref)["ScaleFactor"]

    # Reset buffer area to 1?
    #   TODO Do we want this feature? See similar comment in prepare_sf.py
    #        If we do, need to update prepare_sf.py and then uncomment lines here:
    #        (untested)
    # config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    # n_buff = config["nBufferClusters"]
    # statevector_path = os.path.join(inv_directory, "StateVector.nc")
    # statevector = xr.load_dataset(statevector_path)["StateVector"]
    # n_elements = int(np.nanmax(statevector.data))
    # sf = sf.where(statevector <= n_elements - n_buff) # Replace buffers with nan
    # sf = sf.fillna(1) # Fill nan with 1

    # Get HEMCO diagnostics for current inversion and reference inversion
    #  The HEMCO diags emissions are needed to calculate Jacobian scale
    #  factors for sensitivity inversions that change the prior inventory
    run_name = inv_directory.split("/")[-1]
    ref_run_name = ref_directory.split("/")[-1]
    prior_sim = os.path.join(
        inv_directory, "jacobian_runs", f"{run_name}_0000", "OutputDir"
    )
    ref_prior_sim = os.path.join(
        ref_directory, "jacobian_runs", f"{ref_run_name}_0000", "OutputDir"
    )
    hemco_list = [f for f in os.listdir(prior_sim) if "HEMCO" in f]
    hemco_list.sort()
    ref_hemco_list = [f for f in os.listdir(ref_prior_sim) if "HEMCO" in f]
    ref_hemco_list.sort()
    pth = os.path.join(prior_sim, hemco_list[period_number - 1])
    ref_pth = os.path.join(ref_prior_sim, ref_hemco_list[period_number - 1])
    hemco_emis = xr.load_dataset(pth)["EmisCH4_Total"].isel(time=0, drop=True)
    ref_hemco_emis = xr.load_dataset(ref_pth)["EmisCH4_Total"].isel(time=0, drop=True)

    # Get scale factors to apply to the Jacobian matrix K
    statevector_path = os.path.join(inv_directory, "StateVector.nc")
    statevector = xr.load_dataset(statevector_path)["StateVector"]
    n_elements = int(np.nanmax(statevector.data))
    # TODO - This line assumes the spatial dist of emissions in the buffer elements is
    #        the same in hemco_emis and ref_hemco_emis.... I think.
    #        Should we be summing emissions in each buffer element and then dividing
    #        the sums to get the appropriate scale factor?
    sf_ratio = (sf * hemco_emis) / (sf_ref * ref_hemco_emis)
    sf_K = []
    for e in range(1, n_elements + 1):
        # Get scale factor for state vector element e
        sf_e = np.nanmean(sf_ratio.where(statevector == e).values)
        # Append
        sf_K.append(sf_e)
    sf_K = np.asarray(sf_K)

    return sf_K


if __name__ == "__main__":
    import sys

    period_number = int(sys.argv[1])
    inv_directory = sys.argv[2]
    ref_directory = sys.argv[3]

    # Get the scale factors
    out = get_jacobian_scalefactors(period_number, inv_directory, ref_directory)

    # Save them
    save_path_1 = os.path.join(
        os.path.join(inv_directory, "archive_sf"),
        f"jacobian_scale_factors_period{period_number}.npy",
    )
    save_path_2 = os.path.join(
        inv_directory,
        f"kf_inversions/period{period_number}",
        f"jacobian_scale_factors.npy",
    )
    np.save(save_path_1, out)
    np.save(save_path_2, out)

import xarray as xr
import pandas as pd
import os
import sys
import numpy as np
import yaml
from src.inversion_scripts.utils import sum_total_emissions, get_posterior_emissions
from src.inversion_scripts.utils import get_period_mean_emissions


def prepare_sf(config_path, period_number, base_directory, nudge_factor):
    """
    Function to prepare scale factors for HEMCO emissions.

    This is done at the beginning of an inversion, to generate the appropriate scale factors for the
    prior simulation. In the first period, the scale factors are just unit scale factors (i.e., use
    the emissions from HEMCO directly). In following periods, the scale factors are derived from
    previous inversion results, including nudging to the original prior emission estimates.

    Arguments
        period_number   [int]   : What period are we on? For the first period, period_number = 1
        base_directory  [str]   : The base directory for the inversion, where e.g., "prior_run/" resides
        nudge_factor    [float] : Weight applied to original prior when nudging (default = 0.1)
    """

    # Read config file
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)

    # Fix nudge_factor type
    nudge_factor = float(nudge_factor)

    # Define some useful paths
    unit_sf_path = os.path.join(base_directory, "unit_sf.nc")
    statevector_path = os.path.join(base_directory, "StateVector.nc")
    original_prior_cache = os.path.join(base_directory, "prior_run/OutputDir")
    
    # Get the original emissions for the first inversion period
    periods_csv_path = os.path.join(base_directory, "periods.csv")
    original_emis_ds = get_period_mean_emissions(original_prior_cache, 1, periods_csv_path)

    # Get state vector, grid-cell areas, mask
    statevector = xr.load_dataset(statevector_path)
    areas = original_emis_ds["AREA"]
    state_vector_labels = statevector["StateVector"]
    last_ROI_element = int(
        np.nanmax(state_vector_labels.values) - config["nBufferClusters"]
    )
    mask = state_vector_labels <= last_ROI_element

    # Initialize unit scale factors
    sf = xr.load_dataset(unit_sf_path)
    posterior_scale_ds = sf.copy()

    # If we are past the first inversion period, need to use previous inversion results to construct
    # the initial scale factors for the current period.
    period_number = int(period_number)
    if period_number > 1:
        # For each period up to (but not including) the current one
        for p in range(period_number - 1):
            # Add one since we're counting from period 1, not 0
            p = p + 1

            # Get the original HEMCO emissions for period p
            # Note: we remove soil absorption from the prior for our nudging operations
            # since it is not optimized in the inversion.
            original_emis_ds = get_period_mean_emissions(original_prior_cache, p, periods_csv_path)
            original_emis = original_emis_ds["EmisCH4_Total_ExclSoilAbs"]

            # Get the gridded posterior for period p
            gridded_posterior_filename = (
                "gridded_posterior_ln.nc"
                if config["LognormalErrors"]
                else "gridded_posterior.nc"
            )
            gridded_posterior_path = os.path.join(
                base_directory, f"kf_inversions/period{p}/{gridded_posterior_filename}"
            )
            posterior_p = xr.load_dataset(gridded_posterior_path)

            # Get posterior emissions multiplied up to current period p, and apply nudging
            posterior_scale_ds["ScaleFactor"] = (
                posterior_p["ScaleFactor"] * sf["ScaleFactor"]
            )
            current_posterior_emis = original_emis * posterior_scale_ds["ScaleFactor"]

            nudged_posterior_emis = (
                nudge_factor * original_emis
                + (1 - nudge_factor) * current_posterior_emis
            )  # TODO nudge_factor is currently inverse of what's in the paper, i.e. 0.1 instead of 0.9

            # Sum emissions
            current_total = sum_total_emissions(current_posterior_emis, areas, mask)
            nudged_total = sum_total_emissions(nudged_posterior_emis, areas, mask)

            # Get the final posterior emissions
            lambda_scaler = current_total / nudged_total
            nudged_posterior_roi = nudged_posterior_emis * mask
            nudged_posterior_buf = nudged_posterior_emis * abs(mask - 1)
            scaled_nudged_posterior_emis = (
                nudged_posterior_buf + nudged_posterior_roi * lambda_scaler
            )

            # Get the final posterior scale factors
            sf["ScaleFactor"] = scaled_nudged_posterior_emis / original_emis

            # Reset buffer area to 1
            # Note: resetting the buffer area to 1 seems to prevent issues where
            # emissions can become negative.
            sf["ScaleFactor"] = sf["ScaleFactor"].where(
                state_vector_labels <= last_ROI_element
            )  # Replace buffers with nan
            sf["ScaleFactor"] = sf["ScaleFactor"].fillna(1)  # Fill nan with 1

        print(
            f"Used HEMCO emissions up to week {p} to prepare prior scaling factors for this week."
        )

    # Print the current total emissions in the region of interest
    emis = get_posterior_emissions(original_emis_ds, sf)["EmisCH4_Total"]
    total_emis = sum_total_emissions(emis, areas, mask)
    print(f"Total prior emission = {total_emis} Tg a-1")

    # Ensure good netcdf attributes for HEMCO
    sf.lat.attrs["units"] = "degrees_north"
    sf.lat.attrs["long_name"] = "Latitude"
    sf.lon.attrs["units"] = "degrees_east"
    sf.lon.attrs["long_name"] = "Longitude"
    sf.ScaleFactor.attrs["units"] = "1"

    # Save final scale factors
    save_path = os.path.join(base_directory, "ScaleFactors.nc")
    sf.to_netcdf(
        save_path,
        encoding={v: {"zlib": True, "complevel": 9} for v in sf.data_vars},
    )

    # Archive scale factors
    archive_path = os.path.join(base_directory, "archive_sf")
    sf.to_netcdf(
        os.path.join(archive_path, f"prior_sf_period{period_number}.nc"),
        encoding={v: {"zlib": True, "complevel": 9} for v in sf.data_vars},
    )


if __name__ == "__main__":
    config_path = sys.argv[1]
    period_number = sys.argv[2]
    base_directory = sys.argv[3]
    nudge_factor = sys.argv[4]

    prepare_sf(config_path, period_number, base_directory, nudge_factor)

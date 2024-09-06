import xarray as xr
import numpy as np
import os
import sys
import yaml
from src.inversion_scripts.utils import (
    sum_total_emissions,
    get_posterior_emissions,
    get_period_mean_emissions,
)


def print_posterior_emissions(config_path, period_number, base_directory):
    """
    Simple function to print out the total posterior emissions after an inversion.
    """
    # Read config file
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)

    # Useful paths, directories, file names
    statevector_path = os.path.join(base_directory, "StateVector.nc")
    sf_archive_path = os.path.join(base_directory, "archive_sf")
    post_sf_path = os.path.join(
        sf_archive_path, f"posterior_sf_period{period_number}.nc"
    )
    prior_cache_path = os.path.join(base_directory, "hemco_prior_emis/OutputDir")
    periods_csv_path = os.path.join(base_directory, "periods.csv")
    hemco_diags = get_period_mean_emissions(
        prior_cache_path, period_number, periods_csv_path
    )

    # Get state vector, grid-cell areas, mask
    statevector = xr.load_dataset(statevector_path)
    areas = hemco_diags["AREA"]
    state_vector_labels = statevector["StateVector"]
    last_ROI_element = int(
        np.nanmax(state_vector_labels.values) - config["nBufferClusters"]
    )
    mask = state_vector_labels <= last_ROI_element

    # Emissions
    hemco_emis = hemco_diags
    posterior_sf = xr.load_dataset(post_sf_path)
    posterior_emis_ds = get_posterior_emissions(hemco_emis, posterior_sf)
    if "time" in posterior_emis_ds.dims:
        posterior_emis = posterior_emis_ds["EmisCH4_Total"].isel(time=0, drop=True)
    else:
        posterior_emis = posterior_emis_ds["EmisCH4_Total"].squeeze(drop=True)
    total_emis = sum_total_emissions(posterior_emis, areas, mask)

    # Print
    print(f"Total posterior emission = {total_emis} Tg a-1")


if __name__ == "__main__":

    config_path = sys.argv[1]
    period_number = sys.argv[2]
    base_directory = sys.argv[3]

    print_posterior_emissions(config_path, period_number, base_directory)

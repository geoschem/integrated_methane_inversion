import xarray as xr
import numpy as np
import os
import sys
import yaml

sys.path.append("../inversion_scripts/")
from utils import sum_total_emissions


def print_posterior_emissions(config_path, period_number, base_directory):
    """
    Simple function to print out the total posterior emissions after an inversion.
    """
    # Read config file
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)

    # Useful paths, directories, file names
    statevector_path = os.path.join(base_directory, "StateVector.nc")
    sf_archive_path = os.path.join(base_directory, "archive_sf")
    post_sf_path = os.path.join(sf_archive_path, f"prior_sf_period{period_number}.nc")
    jacobian_dir = os.path.join(base_directory, "jacobian_runs")
    prior_sim = [r for r in os.listdir(jacobian_dir) if "0000" in r][0]
    prior_cache = os.path.join(base_directory, f"jacobian_runs/{prior_sim}/OutputDir")
    hemco_list = [f for f in os.listdir(prior_cache) if "HEMCO" in f]
    hemco_list.sort()
    hemco_diags_path = os.path.join(prior_cache, hemco_list[period_number - 1])
    hemco_diags = xr.load_dataset(hemco_diags_path)

    # Get state vector, grid-cell areas, mask
    statevector = xr.load_dataset(statevector_path)
    areas = xr.load_dataset(hemco_diags)["AREA"]
    state_vector_labels = statevector["StateVector"]
    last_ROI_element = int(
        np.nanmax(state_vector_labels.values) - config["nBufferClusters"]
    )
    mask = state_vector_labels <= last_ROI_element

    # Emissions
    hemco_emis = hemco_diags["EmisCH4_Total"].isel(time=0, drop=True)
    posterior_sf = xr.load_dataset(post_sf_path)["ScaleFactor"]
    # Recall that posterior_sf from archive_sf multiplies sf over the full inversion record
    posterior_emis = posterior_sf * hemco_emis
    total_emis = sum_total_emissions(posterior_emis, areas, mask)

    # Print
    print(f"Total posterior emission = {total_emis} Tg a-1")

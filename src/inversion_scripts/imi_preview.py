#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import os
import datetime
import cartopy.crs as ccrs
import colorcet as cc
from utils import (
    calculate_gridcell_areas,
    sum_total_emissions,
    count_obs_in_mask,
    plot_field,
)
from joblib import Parallel, delayed
from jacobian import read_tropomi
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)


def get_TROPOMI_data(file_path, xlim, ylim, startdate_np64, enddate_np64):
    """
    Returns a dict with the lat, lon, xch4, and albedo_swir observations
    extracted from the given tropomi file. Filters are applied to remove
    unsuitable observations
    Args:
        file_path : string
            path to the tropomi file
        xlim: list
            longitudinal bounds for region of interest
        ylim: list
            latitudinal bounds for region of interest
        startdate_np64: datetime64
            start date for time period of interest
        enddate_np64: datetime64
            end date for time period of interest
    Returns:
         tropomi_data: dict
            dictionary of the extracted values
    """
    # tropomi data dictionary
    tropomi_data = {"lat": [], "lon": [], "xch4": [], "swir_albedo": []}

    # Load the TROPOMI data
    TROPOMI = read_tropomi(file_path)

    # We're only going to consider data within lat/lon/time bounds, with QA > 0.5, and with safe surface albedo values
    sat_ind = np.where(
        (TROPOMI["longitude"] > xlim[0])
        & (TROPOMI["longitude"] < xlim[1])
        & (TROPOMI["latitude"] > ylim[0])
        & (TROPOMI["latitude"] < ylim[1])
        & (TROPOMI["time"] >= startdate_np64)
        & (TROPOMI["time"] <= enddate_np64)
        & (TROPOMI["qa_value"] >= 0.5)
        & (TROPOMI["swir_albedo"] > 0.05)
        & (TROPOMI["blended_albedo"] < 0.85)
    )

    # Loop over observations and archive
    num_obs = len(sat_ind[0])
    for k in range(num_obs):
        lat_idx = sat_ind[0][k]
        lon_idx = sat_ind[1][k]
        tropomi_data["lat"].append(TROPOMI["latitude"][lat_idx, lon_idx])
        tropomi_data["lon"].append(TROPOMI["longitude"][lat_idx, lon_idx])
        tropomi_data["xch4"].append(TROPOMI["methane"][lat_idx, lon_idx])
        tropomi_data["swir_albedo"].append(TROPOMI["swir_albedo"][lat_idx, lon_idx])

    return tropomi_data


def imi_preview(config_path, state_vector_path, preview_dir, tropomi_cache, cpu_count):
    """
    Function to perform preview
    Requires preview simulation to have been run already (to generate HEMCO diags)
    Requires TROPOMI data to have been downloaded already
    """

    # ----------------------------------
    # Setup
    # ----------------------------------

    # Read config file
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)

    # Open the state vector file
    state_vector = xr.open_dataset(state_vector_path)
    state_vector_labels = state_vector["StateVector"]

    # Identify the last element of the region of interest
    last_ROI_element = int(
        np.nanmax(state_vector_labels.values) - config["nBufferClusters"]
    )

    # Define mask for ROI, to be used below
    mask = state_vector_labels <= last_ROI_element

    # ----------------------------------
    # Total prior emissions
    # ----------------------------------

    # Prior emissions
    preview_cache = os.path.join(preview_dir, "OutputDir")
    hemco_diags_file = [
        f for f in os.listdir(preview_cache) if "HEMCO_diagnostics" in f
    ][0]
    prior_pth = os.path.join(preview_cache, hemco_diags_file)

    prior = xr.open_dataset(prior_pth)["EmisCH4_Total"].isel(time=0)

    # Compute total emissions
    if config["Res"] == "0.25x0.3125":
        dlat = 0.25 / 2
        dlon = 0.3125 / 2
    elif config["Res"] == "0.5x0.625":
        dlat = 0.5 / 2
        dlon = 0.625 / 2

    # Compute grid cell areas and sum total emissions in the region of interest
    areas = calculate_gridcell_areas(state_vector, mask, dlat, dlon)
    total_prior_emissions = sum_total_emissions(
        prior, areas, state_vector_labels, last_ROI_element
    )
    outstring1 = (
        f"Total prior emissions in region of interest = {total_prior_emissions} Tg/y"
    )
    print(outstring1)

    # ----------------------------------
    # Observations in region of interest
    # ----------------------------------

    # Paths to tropomi data files
    tropomi_files = [f for f in os.listdir(tropomi_cache) if ".nc" in f]
    tropomi_paths = [os.path.join(tropomi_cache, f) for f in tropomi_files]

    # Latitude/longitude bounds of the inversion domain
    xlim = [float(state_vector.lon.min()), float(state_vector.lon.max())]
    ylim = [float(state_vector.lat.min()), float(state_vector.lat.max())]

    # Start and end dates of the inversion
    startday = str(config["StartDate"])
    endday = str(config["EndDate"])
    start = f"{startday[0:4]}-{startday[4:6]}-{startday[6:8]} 00:00:00"
    end = f"{endday[0:4]}-{endday[4:6]}-{endday[6:8]} 23:59:59"
    startdate_np64 = np.datetime64(
        datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    )
    enddate_np64 = np.datetime64(
        datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")
        - datetime.timedelta(days=1)
    )

    # Only consider tropomi files within date range (in case more are present)
    tropomi_paths = [
        p
        for p in tropomi_paths
        if int(p.split("____")[1][0:8]) >= int(startday)
        and int(p.split("____")[1][0:8]) < int(endday)
    ]
    tropomi_paths.sort()

    # Open tropomi files and filter data
    lat = []
    lon = []
    xch4 = []
    albedo = []

    # read in and filter tropomi observations (uses parallel processing)
    observation_dicts = Parallel(n_jobs=cpu_count)(
        delayed(get_TROPOMI_data)(file_path, xlim, ylim, startdate_np64, enddate_np64)
        for file_path in tropomi_paths
    )
    for dict in observation_dicts:
        lat.extend(dict["lat"])
        lon.extend(dict["lon"])
        xch4.extend(dict["xch4"])
        albedo.extend(dict["swir_albedo"])

    # Assemble in dataframe
    df = pd.DataFrame()
    df["lat"] = lat
    df["lon"] = lon
    df["xch4"] = xch4
    df["swir_albedo"] = albedo

    # Count the number of observations in the region of interest
    num_obs = count_obs_in_mask(mask, df)
    if num_obs < 1:
        sys.exit("Error: No observations found in region of interest")
    outstring2 = f"Found {num_obs} observations in the region of interest"
    print("\n" + outstring2)

    # ----------------------------------
    # Estimate information content
    # ----------------------------------

    # State vector, observations
    n = last_ROI_element  # Number of state vector elements in the ROI
    m = num_obs / n  # Number of observations per state vector element

    # Other parameters
    if config["Res"] == "0.25x0.3125":
        L = 25 * 1000  # Rough length scale of state vector element [m]
    elif config["Res"] == "0.5x0.625":
        L = 50 * 1000  # Rough length scale of state vector element [m]
    U = 5 * (1000 / 3600)  # 5 km/h uniform wind speed in m/s
    p = 101325  # Surface pressure [Pa = kg/m/s2]
    g = 9.8  # Gravity [m/s2]
    Mair = 0.029  # Molar mass of air [kg/mol]
    Mch4 = 0.01604  # Molar mass of methane [kg/mol]
    alpha = 0.4  # Simple parameterization of turbulence

    # Change units of total prior emissions
    total_prior_emissions_kgs = (
        total_prior_emissions * 1e9 / (3600 * 24 * 365)
    )  # kg/s from Tg/y
    total_prior_emissions_kgs_per_element = (
        total_prior_emissions_kgs / L**2 / n
    )  # kg/m2/s from kg/s, per element

    # Error standard deviations with updated units
    sA = config["PriorError"] * total_prior_emissions_kgs_per_element
    sO = config["ObsError"] * 1e-9

    # Averaging kernel sensitivity for each grid element, and dofs
    k = alpha * (Mair * L * g / (Mch4 * U * p))
    a = sA**2 / (sA**2 + (sO / k) ** 2 / m)
    dofs = n * a

    outstring3 = f"k = {np.round(k,5)} kg-1 m2 s"
    outstring4 = f"a = {np.round(a,5)}"
    outstring5 = f"expectedDOFS: {np.round(dofs,5)}"
    print(outstring3)
    print(outstring4)
    print(outstring5)

    # ----------------------------------
    # Estimate dollar cost
    # ----------------------------------

    # Estimate cost by scaling reference cost of $20 for one-month Permian inversion
    # Reference number of state variables = 243
    # Reference number of days = 31
    # Reference cost for EC2 storage = $50 per month
    reference_cost = 20
    reference_num_compute_hours = 10
    hours_in_month = 31 * 24
    reference_storage_cost = 50 * reference_num_compute_hours / hours_in_month
    num_state_variables = np.nanmax(state_vector_labels.values)
    num_days = np.round((enddate_np64 - startdate_np64) / np.timedelta64(1, "D"))
    if config["Res"] == "0.25x0.3125":
        res_factor = 1
    elif config["Res"] == "0.5x0.625":
        res_factor = 0.5
    additional_storage_cost = ((num_days / 31) - 1) * reference_storage_cost
    expected_cost = (
        (reference_cost + additional_storage_cost)
        * (num_state_variables / 243) ** 2
        * (num_days / 31)
        * res_factor
    )

    outstring6 = (
        f"approximate cost = ${np.round(expected_cost,2)} for on-demand instance"
    )
    outstring7 = f"                 = ${np.round(expected_cost/3,2)} for spot instance"
    print(outstring6)
    print(outstring7)

    # ----------------------------------
    # Output
    # ----------------------------------

    # Write preview diagnostics to text file
    outputtextfile = open(os.path.join(preview_dir, "preview_diagnostics.txt"), "w+")
    outputtextfile.write("##" + outstring1 + "\n")
    outputtextfile.write("##" + outstring2 + "\n")
    outputtextfile.write("##" + outstring3 + "\n")
    outputtextfile.write("##" + outstring4 + "\n")
    outputtextfile.write("##" + outstring6 + "\n")
    outputtextfile.write("##" + outstring7 + "\n")
    outputtextfile.write(outstring5)
    outputtextfile.close()

    # Prepare plot data for prior
    prior_kgkm2h = prior * (1000**2) * 60 * 60  # Units kg/km2/h

    # Prepare plot data for observations
    df_means = df.copy(deep=True)
    df_means["lat"] = np.round(df_means["lat"], 1)  # Bin to 0.1x0.1 degrees
    df_means["lon"] = np.round(df_means["lon"], 1)
    df_means = df_means.groupby(["lat", "lon"]).mean()
    ds = df_means.to_xarray()

    # Prepare plot data for observation counts
    df_counts = df.copy(deep=True).drop(["xch4", "swir_albedo"], axis=1)
    df_counts["counts"] = 1
    df_counts["lat"] = np.round(df_counts["lat"], 1)  # Bin to 0.1x0.1 degrees
    df_counts["lon"] = np.round(df_counts["lon"], 1)
    df_counts = df_counts.groupby(["lat", "lon"]).sum()
    ds_counts = df_counts.to_xarray()

    plt.rcParams.update({"font.size": 18})

    # Plot prior emissions
    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        prior_kgkm2h,
        cmap=cc.cm.linear_kryw_5_100_c67_r,
        plot_type="pcolormesh",
        vmin=0,
        vmax=14,
        lon_bounds=None,
        lat_bounds=None,
        levels=21,
        title="Prior emissions",
        cbar_label="Emissions (kg km$^{-2}$ h$^{-1}$)",
        mask=mask,
        only_ROI=False,
    )
    plt.savefig(os.path.join(preview_dir, "preview_prior_emissions.png"), dpi=150)

    # Plot observations
    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        ds["xch4"],
        cmap="Spectral_r",
        plot_type="imshow",
        vmin=1800,
        vmax=1850,
        lon_bounds=None,
        lat_bounds=None,
        title="TROPOMI $X_{CH4}$",
        cbar_label="Column mixing ratio (ppb)",
        mask=mask,
        only_ROI=False,
    )
    plt.savefig(os.path.join(preview_dir, "preview_observations.png"), dpi=150)

    # Plot albedo
    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        ds["swir_albedo"],
        cmap="magma",
        plot_type="imshow",
        vmin=0,
        vmax=0.4,
        lon_bounds=None,
        lat_bounds=None,
        title="SWIR Albedo",
        cbar_label="Albedo",
        mask=mask,
        only_ROI=False,
    )
    plt.savefig(os.path.join(preview_dir, "preview_albedo.png"), dpi=150)

    # Plot observation density
    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        ds_counts["counts"],
        cmap="Blues",
        plot_type="imshow",
        vmin=0,
        vmax=np.nanmax(ds_counts["counts"].values),
        lon_bounds=None,
        lat_bounds=None,
        title="Observation density",
        cbar_label="Number of observations",
        mask=mask,
        only_ROI=False,
    )
    plt.savefig(os.path.join(preview_dir, "preview_observation_density.png"), dpi=150)


if __name__ == "__main__":
    import sys

    config_path = sys.argv[1]
    state_vector_path = sys.argv[2]
    preview_dir = sys.argv[3]
    tropomi_cache = sys.argv[4]
    cpu_count = int(sys.argv[5])

    imi_preview(config_path, state_vector_path, preview_dir, tropomi_cache, cpu_count)

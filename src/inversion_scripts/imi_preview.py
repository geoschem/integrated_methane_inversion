#!/usr/bin/env python
# -*- coding: utf-8 -*-

#SBATCH -N 1

import sys
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml
import os
import datetime
import time
import cartopy.crs as ccrs
import colorcet as cc
from utils import (
    sum_total_emissions,
    count_obs_in_mask,
    plot_field,
    filter_tropomi,
    calculate_area_in_km,
)
from joblib import Parallel, delayed
from operators.TROPOMI_operator import read_tropomi
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

    # Handle unreadable files
    if TROPOMI == None:
        print(f"Skipping {file_path} due to error")
        return TROPOMI

    # We're only going to consider data within lat/lon/time bounds, with QA > 0.5, and with safe surface albedo values
    sat_ind = filter_tropomi(TROPOMI, xlim, ylim, startdate_np64, enddate_np64)

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


def imi_preview(
    inversion_path, config_path, state_vector_path, preview_dir, tropomi_cache
):
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
    # redirect output to log file
    output_file = open(f"{inversion_path}/imi_output.log", "a")
    sys.stdout = output_file
    sys.stderr = output_file

    # Open the state vector file
    state_vector = xr.load_dataset(state_vector_path)
    state_vector_labels = state_vector["StateVector"]

    # Identify the last element of the region of interest
    last_ROI_element = int(
        np.nanmax(state_vector_labels.values) - config["nBufferClusters"]
    )

    # # Define mask for ROI, to be used below
    a, df, num_days, prior, outstrings = estimate_averaging_kernel(
        config, state_vector_path, preview_dir, tropomi_cache, preview=True
    )
    mask = state_vector_labels <= last_ROI_element

    # Count the number of observations in the region of interest
    num_obs = count_obs_in_mask(mask, df)
    if num_obs < 1:
        sys.exit("Error: No observations found in region of interest")
    outstring2 = f"Found {num_obs} observations in the region of interest"
    print("\n" + outstring2)

    # ----------------------------------
    # Estimate dollar cost
    # ----------------------------------

    # Estimate cost by scaling reference cost of $20 for one-month Permian inversion
    # Reference number of state variables = 243
    # Reference number of days = 31
    # Reference cost for EC2 storage = $50 per month
    # Reference area = area of 24-39 N 95-111W
    reference_cost = 20
    reference_num_compute_hours = 10
    reference_area_km = calculate_area_in_km(
        [(-111, 24), (-95, 24), (-95, 39), (-111, 39)]
    )
    hours_in_month = 31 * 24
    reference_storage_cost = 50 * reference_num_compute_hours / hours_in_month
    num_state_variables = np.nanmax(state_vector_labels.values)

    lats = [float(state_vector.lat.min()), float(state_vector.lat.max())]
    lons = [float(state_vector.lon.min()), float(state_vector.lon.max())]
    coords = [
        (lons[0], lats[0]),
        (lons[1], lats[0]),
        (lons[1], lats[1]),
        (lons[0], lats[1]),
    ]
    inversion_area_km = calculate_area_in_km(coords)

    if config["Res"] == "0.25x0.3125":
        res_factor = 1
    elif config["Res"] == "0.5x0.625":
        res_factor = 0.5
    additional_storage_cost = ((num_days / 31) - 1) * reference_storage_cost
    expected_cost = (
        (reference_cost + additional_storage_cost)
        * (num_state_variables / 243)
        * (inversion_area_km / reference_area_km)
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
    outputtextfile.write("##" + outstring2 + "\n")
    outputtextfile.write("##" + outstring6 + "\n")
    outputtextfile.write("##" + outstring7 + "\n")
    outputtextfile.write(outstrings)
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
    plt.savefig(
        os.path.join(preview_dir, "preview_prior_emissions.png"),
        bbox_inches="tight",
        dpi=150,
    )

    # Plot observations
    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        ds["xch4"],
        cmap="Spectral_r",
        plot_type="pcolormesh",
        vmin=1800,
        vmax=1850,
        lon_bounds=None,
        lat_bounds=None,
        title="TROPOMI $X_{CH4}$",
        cbar_label="Column mixing ratio (ppb)",
        mask=mask,
        only_ROI=False,
    )
    plt.savefig(
        os.path.join(preview_dir, "preview_observations.png"),
        bbox_inches="tight",
        dpi=150,
    )

    # Plot albedo
    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        ds["swir_albedo"],
        cmap="magma",
        plot_type="pcolormesh",
        vmin=0,
        vmax=0.4,
        lon_bounds=None,
        lat_bounds=None,
        title="SWIR Albedo",
        cbar_label="Albedo",
        mask=mask,
        only_ROI=False,
    )
    plt.savefig(
        os.path.join(preview_dir, "preview_albedo.png"), bbox_inches="tight", dpi=150
    )

    # Plot observation density
    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        ds_counts["counts"],
        cmap="Blues",
        plot_type="pcolormesh",
        vmin=0,
        vmax=np.nanmax(ds_counts["counts"].values),
        lon_bounds=None,
        lat_bounds=None,
        title="Observation density",
        cbar_label="Number of observations",
        mask=mask,
        only_ROI=False,
    )
    plt.savefig(
        os.path.join(preview_dir, "preview_observation_density.png"),
        bbox_inches="tight",
        dpi=150,
    )

    sensitivities_da = map_sensitivities_to_sv(a, state_vector, last_ROI_element)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        sensitivities_da["Sensitivities"],
        cmap=cc.cm.CET_L19,
        lon_bounds=None,
        lat_bounds=None,
        title="Estimated Averaging kernel sensitivities",
        cbar_label="Sensitivity",
        only_ROI=True,
        state_vector_labels=state_vector_labels,
        last_ROI_element=last_ROI_element,
    )
    plt.savefig(
        os.path.join(preview_dir, "preview_estimated_sensitivities.png"),
        bbox_inches="tight",
        dpi=150,
    )


def map_sensitivities_to_sv(sensitivities, sv, last_ROI_element):
    """
    maps sensitivities onto 2D xarray Datarray for visualization
    """
    s = sv.copy().rename({"StateVector": "Sensitivities"})
    mask = s["Sensitivities"] <= last_ROI_element
    s["Sensitivities"] = s["Sensitivities"].where(mask)
    # map sensitivities onto corresponding xarray DataArray
    for i in range(1, last_ROI_element + 1):
        mask = sv["StateVector"] == i
        s = xr.where(mask, sensitivities[i - 1], s)

    return s


def estimate_averaging_kernel(
    config, state_vector_path, preview_dir, tropomi_cache, preview=False
):
    """
    Estimates the averaging kernel sensitivities using prior emissions
    and the number of observations available in each grid cell
    """

    # ----------------------------------
    # Setup
    # ----------------------------------

    # Open the state vector file
    state_vector = xr.load_dataset(state_vector_path)
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
    prior = xr.load_dataset(prior_pth)["EmisCH4_Total"].isel(time=0)

    # Compute total emissions in the region of interest
    areas = xr.load_dataset(prior_pth)["AREA"]
    total_prior_emissions = sum_total_emissions(prior, areas, mask)
    outstring1 = (
        f"Total prior emissions in region of interest = {total_prior_emissions} Tg/y \n"
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
    observation_dicts = Parallel(n_jobs=-1)(
        delayed(get_TROPOMI_data)(file_path, xlim, ylim, startdate_np64, enddate_np64)
        for file_path in tropomi_paths
    )
    # remove any problematic observation dicts (eg. corrupted data file)
    observation_dicts = list(filter(None, observation_dicts))

    for dict in observation_dicts:
        lat.extend(dict["lat"])
        lon.extend(dict["lon"])
        xch4.extend(dict["xch4"])
        albedo.extend(dict["swir_albedo"])

    # Assemble in dataframe
    df = pd.DataFrame()
    df["lat"] = lat
    df["lon"] = lon
    df["count"] = np.ones(len(lat))
    df["swir_albedo"] = albedo
    df["xch4"] = xch4

    # set resolution specific variables
    if config["Res"] == "0.25x0.3125":
        L_native = 25 * 1000  # Rough length scale of native state vector element [m]
        lat_step = 0.25
        lon_step = 0.3125
    elif config["Res"] == "0.5x0.625":
        lat_step = 0.5
        lon_step = 0.625
        L_native = 50 * 1000  # Rough length scale of native state vector element [m]

    # bin observations into gridcells and map onto statevector
    observation_counts = add_observation_counts(df, state_vector, lat_step, lon_step)

    # parallel processing function
    def process(i):
        mask = state_vector_labels == i
        # prior emissions for each element (in Tg/y)
        emissions_temp = sum_total_emissions(prior, areas, mask)
        # append the calculated length scale of element
        L_temp = L_native * state_vector_labels.where(mask).count().item()
        # append the number of obs in each element
        num_obs_temp = np.nansum(observation_counts["count"].where(mask).values)
        return emissions_temp, L_temp, num_obs_temp

    # in parallel, create lists of emissions, number of observations,
    # and rough length scale for each cluster element in ROI
    result = Parallel(n_jobs=-1)(
        delayed(process)(i) for i in range(1, last_ROI_element + 1)
    )

    # unpack list of tuples into individual lists
    emissions, L, num_obs = [list(item) for item in zip(*result)]

    # ----------------------------------
    # Estimate information content
    # ----------------------------------

    # State vector, observations
    emissions = np.array(emissions)
    m = np.array(num_obs)  # Number of observations per state vector element
    L = np.array(L)

    # Other parameters
    U = 5 * (1000 / 3600)  # 5 km/h uniform wind speed in m/s
    p = 101325  # Surface pressure [Pa = kg/m/s2]
    g = 9.8  # Gravity [m/s2]
    Mair = 0.029  # Molar mass of air [kg/mol]
    Mch4 = 0.01604  # Molar mass of methane [kg/mol]
    alpha = 0.4  # Simple parameterization of turbulence

    # Change units of total prior emissions
    emissions_kgs = emissions * 1e9 / (3600 * 24 * 365)  # kg/s from Tg/y
    emissions_kgs_per_m2 = emissions_kgs / np.power(
        L, 2
    )  # kg/m2/s from kg/s, per element

    time_delta = enddate_np64 - startdate_np64
    num_days = np.round((time_delta) / np.timedelta64(1, "D"))

    # Error standard deviations with updated units
    sA = config["PriorError"] * emissions_kgs_per_m2
    sO = config["ObsError"] * 1e-9

    # Averaging kernel sensitivity for each grid element
    k = alpha * (Mair * L * g / (Mch4 * U * p))
    a = sA**2 / (sA**2 + (sO / k) ** 2 / m)

    outstring3 = f"k = {np.round(k,5)} kg-1 m2 s"
    outstring4 = f"a = {np.round(a,5)} \n"
    outstring5 = f"expectedDOFS: {np.round(sum(a),5)}"
    print(outstring3)
    print(outstring4)
    print(outstring5)

    if preview:
        outstrings = (
            f"##{outstring1}\n" + f"##{outstring3}\n" + f"##{outstring4}\n" + outstring5
        )
        return a, df, num_days, prior, outstrings
    else:
        return a


def add_observation_counts(df, state_vector, lat_step, lon_step):
    """
    Given arbitrary observation coordinates in a pandas df, group
    them by gridcell and return the number of observations mapped
    onto the statevector dataset
    """
    to_lon = lambda x: np.floor(x / lon_step) * lon_step
    to_lat = lambda x: np.floor(x / lat_step) * lat_step

    df = df.rename(columns={"lon": "old_lon", "lat": "old_lat"})

    df["lat"] = to_lat(df.old_lat)
    df["lon"] = to_lon(df.old_lon)
    groups = df.groupby(["lat", "lon"])

    counts_ds = groups.sum().to_xarray().drop_vars(["old_lat", "old_lon"])
    return xr.merge([counts_ds, state_vector])


if __name__ == "__main__":
    inversion_path = sys.argv[1]
    config_path = sys.argv[2]
    state_vector_path = sys.argv[3]
    preview_dir = sys.argv[4]
    tropomi_cache = sys.argv[5]

    imi_preview(
        inversion_path, config_path, state_vector_path, preview_dir, tropomi_cache
    )

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SBATCH -N 1

import os
import sys
import yaml
import time
import warnings
import datetime
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib
import colorcet as cc
import cartopy.crs as ccrs

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from src.inversion_scripts.point_sources import get_point_source_coordinates
from src.inversion_scripts.utils import (
    sum_total_emissions,
    count_obs_in_mask,
    plot_field,
    filter_tropomi,
    filter_blended,
    calculate_area_in_km,
    calculate_superobservation_error,
    get_mean_emissions,
    get_posterior_emissions,
)
from joblib import Parallel, delayed
from src.inversion_scripts.operators.TROPOMI_operator import (
    read_tropomi,
    read_blended,
)

warnings.filterwarnings("ignore", category=FutureWarning)


def get_TROPOMI_data(
    file_path, BlendedTROPOMI, xlim, ylim, startdate_np64, enddate_np64
):
    """
    Returns a dict with the lat, lon, xch4, and albedo_swir observations
    extracted from the given tropomi file. Filters are applied to remove
    unsuitable observations
    Args:
        file_path : string
            path to the tropomi file
        BlendedTROPOMI : bool
            if True, use blended TROPOMI+GOSAT data
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
    assert isinstance(BlendedTROPOMI, bool), "BlendedTROPOMI is not a bool"
    if BlendedTROPOMI:
        TROPOMI = read_blended(file_path)
    else:
        TROPOMI = read_tropomi(file_path)
    if TROPOMI == None:
        print(f"Skipping {file_path} due to error")
        return TROPOMI

    if BlendedTROPOMI:
        # Only going to consider data within lat/lon/time bounds and without problematic coastal pixels
        sat_ind = filter_blended(TROPOMI, xlim, ylim, startdate_np64, enddate_np64)
    else:
        # Only going to consider data within lat/lon/time bounds, with QA > 0.5, and with safe surface albedo values
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
    for key in config.keys():
        if isinstance(config[key],str):
            config[key] = os.path.expandvars(config[key])

    # Open the state vector file
    state_vector = xr.load_dataset(state_vector_path)
    state_vector_labels = state_vector["StateVector"]

    # Identify the last element of the region of interest
    last_ROI_element = int(
        np.nanmax(state_vector_labels.values) - config["nBufferClusters"]
    )

    # # Define mask for ROI, to be used below
    a, df, num_days, prior, outstrings = estimate_averaging_kernel(
        config,
        state_vector_path,
        preview_dir,
        tropomi_cache,
        preview=True,
        kf_index=None,
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
    elif config["Res"] == "2.0x2.5":
        res_factor = 0.125
    elif config["Res"] == "4.0x5.0":
        res_factor = 0.0625
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
        point_sources=get_point_source_coordinates(config),
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
    
    # plot state vector
    num_colors = state_vector_labels.where(mask).max().item()
    sv_cmap = matplotlib.colors.ListedColormap(np.random.rand(int(num_colors),3))
    fig = plt.figure(figsize=(8, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_field(
        ax,
        state_vector_labels,
        cmap=sv_cmap,
        lon_bounds=None,
        lat_bounds=None,
        title="State Vector Elements",
        cbar_label="Element ID",
        only_ROI=True,
        state_vector_labels=state_vector_labels,
        last_ROI_element=last_ROI_element,
    )
    plt.savefig(
        os.path.join(preview_dir, "preview_state_vector.png"),
        bbox_inches="tight",
        dpi=150,
    )
    
    # plot estimated averaging kernel sensitivities
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
    
    # calculate expected DOFS
    expectedDOFS = np.round(sum(a), 5)
    if expectedDOFS < config["DOFSThreshold"]:
        print(
            f"\nExpected DOFS = {expectedDOFS} are less than DOFSThreshold = {config['DOFSThreshold']}. Exiting.\n"
        )
        print(
            "Consider increasing the inversion period, increasing the prior error, or using another prior inventory.\n"
        )
        # if run with sbatch this ensures the exit code is not lost.
        file = open(".error_status_file.txt", "w")
        file.write("Error Status: 1")
        file.close()
        sys.exit(1)


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
    config, state_vector_path, preview_dir, tropomi_cache, preview=False, kf_index=None
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
    # Start and end dates of the inversion
    startday = str(config["StartDate"])
    endday = str(config["EndDate"])

    # Prior emissions
    prior_cache = os.path.join(config["OutputPath"], config["RunName"], "prior_run/OutputDir")
    # adjustments for when performing for dynamic kf clustering
    if kf_index is not None:
        # use different date range for KF inversion if kf_index is not None
        rundir_path = preview_dir.split("preview_run")[0]
        periods = pd.read_csv(f"{rundir_path}periods.csv")
        startday = str(periods.iloc[kf_index - 1]["Starts"])
        endday = str(periods.iloc[kf_index - 1]["Ends"])

        # use the nudged (prior) emissions for generating averaging kernel estimate
        sf = xr.load_dataset(f"{rundir_path}archive_sf/prior_sf_period{kf_index}.nc")
        prior_ds = get_mean_emissions(startday, endday, prior_cache)
        prior_ds = get_posterior_emissions(prior_ds, sf)
    else:
        prior_ds = get_mean_emissions(startday, endday, prior_cache)
        
    prior = prior_ds["EmisCH4_Total"]

    # Compute total emissions in the region of interest
    areas = prior_ds["AREA"]
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

    # Use blended TROPOMI+GOSAT data or operational TROPOMI data?
    BlendedTROPOMI = config["BlendedTROPOMI"]

    # Open tropomi files and filter data
    lat = []
    lon = []
    xch4 = []
    albedo = []

    # Read in and filter tropomi observations (uses parallel processing)
    observation_dicts = Parallel(n_jobs=-1)(
        delayed(get_TROPOMI_data)(
            file_path, BlendedTROPOMI, xlim, ylim, startdate_np64, enddate_np64
        )
        for file_path in tropomi_paths
    )
    # Remove any problematic observation dicts (eg. corrupted data file)
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

    # Set resolution specific variables
    # L_native = Rough length scale of native state vector element [m]
    if config["Res"] == "0.25x0.3125":
        L_native = 25 * 1000
        lat_step = 0.25
        lon_step = 0.3125
    elif config["Res"] == "0.5x0.625":
        L_native = 50 * 1000
        lat_step = 0.5
        lon_step = 0.625
    elif config["Res"] == "2.0x2.5":
        L_native = 200 * 1000
        lat_step = 2.0
        lon_step = 2.5
    elif config["Res"] == "4.0x5.0":
        L_native = 400 * 1000
        lat_step = 4.0
        lon_step = 5.0

    # bin observations into gridcells and map onto statevector
    observation_counts = add_observation_counts(df, state_vector, lat_step, lon_step)

    # parallel processing function
    def process(i):
        mask = state_vector_labels == i
        # prior emissions for each element (in Tg/y)
        emissions_temp = sum_total_emissions(prior, areas, mask)
        # number of native state vector elements in each element
        size_temp = state_vector_labels.where(mask).count().item()
        # append the calculated length scale of element
        L_temp = L_native * size_temp
        # append the number of obs in each element
        num_obs_temp = np.nansum(observation_counts["count"].where(mask).values)
        return emissions_temp, L_temp, size_temp, num_obs_temp

    # in parallel, create lists of emissions, number of observations,
    # and rough length scale for each cluster element in ROI
    result = Parallel(n_jobs=-1)(
        delayed(process)(i) for i in range(1, last_ROI_element + 1)
    )

    # unpack list of tuples into individual lists
    emissions, L, num_native_elements, num_obs = [list(item) for item in zip(*result)]

    if np.sum(num_obs) < 1:
        sys.exit("Error: No observations found in region of interest")
    outstring2 = f"Found {np.sum(num_obs)} observations in the region of interest"

    # ----------------------------------
    # Estimate information content
    # ----------------------------------

    time_delta = enddate_np64 - startdate_np64
    num_days = np.round((time_delta) / np.timedelta64(1, "D"))

    # State vector, observations
    emissions = np.array(emissions)
    m = np.array(num_days)  # Number of observation days
    L = np.array(L)
    num_native_elements = np.array(num_native_elements)

    # If Kalman filter mode, count observations per inversion period
    if config["KalmanMode"]:
        startday_dt = datetime.datetime.strptime(startday, "%Y%m%d")
        endday_dt = datetime.datetime.strptime(endday, "%Y%m%d")
        n_periods = np.floor((endday_dt - startday_dt).days / config["UpdateFreqDays"])
        n_obs_per_period = np.round(num_obs / n_periods)
        outstring2 = f"Found {int(np.sum(n_obs_per_period))} observations in the region of interest per inversion period, for {int(n_periods)} period(s)"
        m = config["UpdateFreqDays"]  # number of days in inversion period

    print("\n" + outstring2)

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

    # Error standard deviations with updated units
    sA = config["PriorError"] * emissions_kgs_per_m2
    sO = config["ObsError"]

    # Calculate superobservation error to use in averaging kernel sensitivity equation
    # from P observations per grid cell = number of observations per grid cell / m days
    # P is number of observations per grid cell (native state vector element)
    # Note: to account for clustering we do num_obs / num_native_elements / num_days
    P = np.array(num_obs) / num_native_elements / num_days
    s_superO_1 = calculate_superobservation_error(
        sO, 1
    )  # for handling cells with 0 observations (avoid divide by 0)

    # list containing superobservation error per state vector element
    s_superO_p = [
        calculate_superobservation_error(sO, element) if element >= 1.0 else s_superO_1
        for element in P
    ]
    s_superO = np.array(s_superO_p) * 1e-9  # convert to ppb

    # Averaging kernel sensitivity for each grid element
    # Note: m is the number of superobservations (just days if native),
    # but if clustered it is days * num_native_elements
    k = alpha * (Mair * L * g / (Mch4 * U * p))
    a = sA**2 / (sA**2 + (s_superO / k) ** 2 / (m * num_native_elements))

    outstring3 = f"k = {np.round(k,5)} kg-1 m2 s"
    outstring4 = f"a = {np.round(a,5)} \n"
    outstring5 = f"expectedDOFS: {np.round(sum(a),5)}"

    if config["KalmanMode"]:
        outstring5 += " per inversion period"

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

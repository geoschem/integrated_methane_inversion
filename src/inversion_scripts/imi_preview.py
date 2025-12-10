#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SBATCH -N 1

import os
import sys
import yaml
import warnings
import datetime
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib
import colorcet as cc
import cartopy.crs as ccrs
from scipy.ndimage import binary_dilation
import gc

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from src.inversion_scripts.point_sources import get_point_source_coordinates
from src.inversion_scripts.utils import (
    sum_total_emissions,
    plot_field,
    plot_field_gchp, # note we need to set vmin and vmax to make it proper for all cubic faces
    filter_tropomi,
    filter_blended,
    calculate_superobservation_error,
    get_mean_emissions,
    get_posterior_emissions,
    sum_and_sort_along_statevector,
)
from src.inversion_scripts.operators.TROPOMI_operator import (
    read_tropomi,
    read_blended,
)
from src.inversion_scripts.classify_TROPOMI_obs_to_CSgrids import (
    latlon_to_cartesian,
    build_kdtree,
    classify_obs_to_cs_grid,
    map_obs_to_CSgrid,
)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

def get_TROPOMI_data(
    file_path, BlendedTROPOMI, xlim, ylim, startdate_np64, enddate_np64, use_water_obs
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
        use_water_obs: bool
            if True, use observations over water
    Returns:
         tropomi_data: dict
            dictionary of the extracted values
    """
    # tropomi data dictionary
    tropomi_data = {"lat": [], "lon": [], "xch4": [], "swir_albedo": [], "time": []}

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
        sat_ind = filter_blended(
            TROPOMI, xlim, ylim, startdate_np64, enddate_np64, use_water_obs
        )
    else:
        # Only going to consider data within lat/lon/time bounds, with QA > 0.5, and with safe surface albedo values
        sat_ind = filter_tropomi(
            TROPOMI, xlim, ylim, startdate_np64, enddate_np64, use_water_obs
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
        tropomi_data["time"].append(TROPOMI["time"][lat_idx, lon_idx])

    return tropomi_data


def imi_preview(
    config_path, state_vector_path, preview_dir, tropomi_cache
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
        if isinstance(config[key], str):
            config[key] = os.path.expandvars(config[key])

    # Open the state vector file and squeeze time dimension
    state_vector = xr.load_dataset(state_vector_path).squeeze()
    state_vector_labels = state_vector["StateVector"]

    # Identify the last element of the region of interest
    last_ROI_element = int(
        np.nanmax(state_vector_labels.values) - config["nBufferClusters"]
    )

    if config['UseGCHP']:
        basedir = os.path.expandvars(
            os.path.join(config["OutputPath"], config["RunName"])
        )
        gridfpath = f'{basedir}/CS_grids/grids.c{config["CS_RES"]}.nc'
        gridds = xr.open_dataset(gridfpath)
        corner_lons = gridds['corner_lons']
        corner_lats = gridds['corner_lats']
        
    # Set latitude/longitude bounds for plots
    if not config['UseGCHP']:
        # Trim 1-2.5 degrees to remove GEOS-Chem buffer zone
        if config["Res"] == "0.25x0.3125":
            degx = 4 * 0.3125
            degy = 4 * 0.25
        elif config["Res"] == "0.5x0.625":
            degx = 4 * 0.625
            degy = 4 * 0.5
        elif config["Res"] == "2.0x2.5":
            degx = 4 * 2.5
            degy = 4 * 2.0

        lon_bounds = [
            np.min(state_vector.lon.values) + degx,
            np.max(state_vector.lon.values) - degx,
        ]
        lat_bounds = [
            np.min(state_vector.lat.values) + degy,
            np.max(state_vector.lat.values) - degy,
        ]
    elif config['STRETCH_GRID']:
        buffer_bounds = 0.
        temp_lons = gridds['corner_lons'].values[5,...].copy()
        temp_lons[temp_lons>180] -= 360
        lon_min = max(temp_lons.min() - buffer_bounds, -180)
        lon_max = min(temp_lons.max() + buffer_bounds, 180)
        lat_min = max(gridds['corner_lats'].values[5,...].min(), -90)
        lat_max = min(gridds['corner_lats'].values[5,...].max(), 90)
        lon_bounds = [lon_min, lon_max]
        lat_bounds = [lat_min, lat_max]
    else:
        lon_bounds = [-180, 180]
        lat_bounds = [-90, 90]
    
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

    # ----------------------------------
    # Estimate dollar cost
    # ----------------------------------

    # Estimate cost by scaling reference cost of $20 for one-month Permian inversion
    # Reference number of state variables = 243
    # Reference number of days = 31
    # Reference cost for EC2 storage = $50 per month
    # Reference area = area of 24-39 N 95-111W
    # Note: calculate_area_in_km in src.inversion_scripts.utils cannot get the surface area correctly when it is nearly global coverage
    #       Thus, here we turn to get the ratio of the number of grid boxes relative to reference grid
    reference_cost = 20
    reference_num_compute_hours = 10
    ref_nbox = ((39 - 24) / 0.25) * ((-95 + 111) / 0.3125)

    hours_in_month = 31 * 24
    reference_storage_cost = 50 * reference_num_compute_hours / hours_in_month
    num_state_variables = np.nanmax(state_vector_labels.values)

    if config['UseGCHP']:
        nbox = 6 * config['CS_RES'] ** 2
    else:
        if config["Res"] == "0.125x0.15625":
            deltalat = 0.125
            deltalon = 0.15625
        if config["Res"] == "0.25x0.3125":
            deltalat = 0.25
            deltalon = 0.3125
        elif config["Res"] == "0.5x0.625":
            deltalat = 0.5
            deltalon = 0.625
        elif config["Res"] == "2.0x2.5":
            deltalat = 2.0
            deltalon = 2.5
        elif config["Res"] == "4.0x5.0":
            deltalat = 4.0
            deltalon = 5.0
        lats = [float(state_vector.lat.min()), float(state_vector.lat.max())]
        lons = [float(state_vector.lon.min()), float(state_vector.lon.max())]
        nbox = (lats[1] - lats[0]) / deltalat * (lons[1] - lons[0]) / deltalon
    nbox_factor = nbox / ref_nbox
    additional_storage_cost = ((num_days / 31) - 1) * reference_storage_cost
    expected_cost = (
        (reference_cost + additional_storage_cost)
        * (num_state_variables / 243)
        * nbox_factor
        * (num_days / 31)
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
    if config['UseGCHP']:
        plot_field_gchp(
            ax,
            corner_lons,
            corner_lats,
            prior_kgkm2h,
            cmap=cc.cm.linear_kryw_5_100_c67_r,
            plot_type="pcolormesh",
            vmin=0,
            vmax=14,
            lon_bounds=lon_bounds,
            lat_bounds=lat_bounds,
            levels=21,
            title="Prior emissions",
            point_sources=get_point_source_coordinates(config),
            cbar_label="Emissions (kg km$^{-2}$ h$^{-1}$)",
            only_ROI=False,
        )
    else:
        plot_field(
            ax,
            prior_kgkm2h,
            cmap=cc.cm.linear_kryw_5_100_c67_r,
            plot_type="pcolormesh",
            vmin=0,
            vmax=14,
            lon_bounds=lon_bounds,
            lat_bounds=lat_bounds,
            levels=21,
            title="Prior emissions",
            point_sources=get_point_source_coordinates(config),
            cbar_label="Emissions (kg km$^{-2}$ h$^{-1}$)",
            mask=mask if config["isRegional"] else None,
            only_ROI=False,
        )
    plt.savefig(
        os.path.join(preview_dir, "preview_prior_emissions.png"),
        bbox_inches="tight",
        dpi=150,
    )

    # simple function to find the dynamic range for colorbar
    dynamic_range = lambda vals: (
        np.round(np.nanmedian(vals) / 25.0) * 25 - 25,
        np.round(np.nanmedian(vals) / 25.0) * 25 + 25,
    )
    # Plot observations
    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    xch4_min, xch4_max = dynamic_range(ds["xch4"].values)
    plot_field(
        ax,
        ds["xch4"],
        cmap="Spectral_r",
        plot_type="pcolormesh",
        vmin=xch4_min,
        vmax=xch4_max,
        lon_bounds=lon_bounds,
        lat_bounds=lat_bounds,
        title="TROPOMI $X_{CH4}$",
        cbar_label="Column mixing ratio (ppb)",
        mask=mask if config["isRegional"] else None,
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
        lon_bounds=lon_bounds,
        lat_bounds=lat_bounds,
        title="SWIR Albedo",
        cbar_label="Albedo",
        mask=mask if config["isRegional"] else None,
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
        lon_bounds=lon_bounds,
        lat_bounds=lat_bounds,
        title="Observation density",
        cbar_label="Number of observations",
        mask=mask if config["isRegional"] else None,
        only_ROI=False,
    )
    plt.savefig(
        os.path.join(preview_dir, "preview_observation_density.png"),
        bbox_inches="tight",
        dpi=150,
    )

    # plot state vector
    num_colors = state_vector_labels.where(mask).max().item()
    sv_cmap = matplotlib.colors.ListedColormap(np.random.rand(int(num_colors), 3))
    fig = plt.figure(figsize=(8, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    if config['UseGCHP']:
        plot_field_gchp(
            ax,
            corner_lons,
            corner_lats,
            state_vector_labels,
            cmap=sv_cmap,
            vmin=1,
            vmax=num_colors,
            lon_bounds=lon_bounds,
            lat_bounds=lat_bounds,
            title="State Vector Elements",
            cbar_label="Element ID",
            only_ROI=True,
            state_vector_labels=state_vector_labels,
            last_ROI_element=last_ROI_element,
        )
    else:
        plot_field(
            ax,
            state_vector_labels,
            cmap=sv_cmap,
            vmin=1,
            vmax=num_colors,
            lon_bounds=lon_bounds,
            lat_bounds=lat_bounds,
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
    sensitivities = map_sensitivities_to_sv(a, state_vector_labels, last_ROI_element)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()})
    if config['UseGCHP']:
        plot_field_gchp(
            ax,
            corner_lons,
            corner_lats,
            sensitivities,
            cmap=cc.cm.CET_L19,
            vmin=0,
            vmax=np.nanpercentile(sensitivities.values, 95),
            lon_bounds=lon_bounds,
            lat_bounds=lat_bounds,
            title="Estimated Averaging kernel sensitivities",
            cbar_label="Sensitivity",
            only_ROI=True,
            state_vector_labels=state_vector_labels,
            last_ROI_element=last_ROI_element,
        )
    else:
        plot_field(
            ax,
            sensitivities,
            cmap=cc.cm.CET_L19,
            vmin=0,
            vmax=np.nanpercentile(sensitivities.values, 95),
            lon_bounds=lon_bounds,
            lat_bounds=lat_bounds,
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


def map_sensitivities_to_sv(sensitivities, state_vector_lables, last_ROI_element):
    """
    Map 1D sensitivities onto a label grid.

    Parameters
    ----------
    sensitivities : array-like, shape (last_ROI_element,)
        Sensitivity value corresponding to ROI label 1..last_ROI_element.
    state_vector_lables : xr.DataArray
        The StateVector array containing integer labels and NaNs.
        Can have shape (lat, lon), (nf, Ydim, Xdim), (time, nf, Ydim, Xdim), etc.
    last_ROI_element : int

    Returns
    -------
    xr.DataArray with the same dims/coords as labels_da
    """
    labels = state_vector_lables.values
    sens = np.asarray(sensitivities)

    # Valid ROI labels: 1..last_ROI_element
    valid = np.isfinite(labels) & (labels <= last_ROI_element)

    # Output array
    out = np.full(labels.shape, np.nan, dtype=sens.dtype)

    # Convert labels → 0-based indices
    idx = labels[valid].astype(int) - 1

    # Fill output
    out[valid] = sens[idx]

    # Wrap back into a DataArray
    return xr.DataArray(
        out,
        coords=state_vector_lables.coords,
        dims=state_vector_lables.dims,
        name="Sensitivities"
    )

def get_sectoral_outputs(prior_ds, areas, mask, preview_dir):
    """
    Get sectoral emissions from the prior dataset
    """
    # Plot sectoral emissions
    sectors = [
        var
        for var in list(prior_ds.keys())
        if "EmisCH4" in var and not ("Total" in var or "Excl" in var)
    ]

    # Calculate total emissions for each sector
    prior_sector_vals = []
    positive_sectors = []
    for sector in sectors:
        prior_val = sum_total_emissions(prior_ds[sector], areas, mask)
        if prior_val > 0:
            prior_sector_vals.append(prior_val)
            positive_sectors.append(sector.replace("EmisCH4_", ""))

    # Combine the lists into tuples and sort them based on prior_sector_vals
    combined = list(zip(positive_sectors, prior_sector_vals))
    combined_sorted = sorted(combined, key=lambda x: x[1])
    positive_sectors, prior_sector_vals = zip(*combined_sorted)

    # Plot bars for prior emissions
    fig = plt.figure(figsize=(10, 5))
    ax = fig.subplots(1, 1)
    bar_height = 0.35
    ind = np.arange(len(positive_sectors))
    bars1 = ax.barh(
        ind,
        prior_sector_vals,
        bar_height,
        color="goldenrod",
        label="Prior Emissions",
    )

    # Add labels and title
    ax.set_xlabel("Emissions ($Tg\ a^{-1}$)")
    ax.set_ylabel("Sector")
    ax.set_title("Sectoral Emissions (Prior Inventory)")
    ax.set_yticks(ind)
    ax.set_yticklabels(positive_sectors)

    plt.savefig(f"{preview_dir}/prior_sectoral_emissions.png", bbox_inches="tight")

    sector_totals = {}

    for item in combined_sorted:
        category = item[0]
        sector_prior = item[1]
        sector_totals[f"{category}Prior"] = sector_prior

    # Save the statistics to a file
    stats_pd = pd.DataFrame(sector_totals, index=[0])
    stats_pd.to_csv(f"{preview_dir}/prior_sectoral_statistics.csv", index=False)

    return

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

    # Open the state vector file and squeeze time dimension
    state_vector = xr.load_dataset(state_vector_path).squeeze()
    state_vector_labels = state_vector["StateVector"]

    # Identify the last element of the region of interest
    last_ROI_element = int(
        np.nanmax(state_vector_labels.values) - config["nBufferClusters"]
    )

    # Whether to use observations over water?
    use_water_obs = config["UseWaterObs"] if "UseWaterObs" in config.keys() else False

    # Define mask for ROI, to be used below
    mask = state_vector_labels <= last_ROI_element

    # ----------------------------------
    # Total prior emissions
    # ----------------------------------
    # Start and end dates of the inversion
    startday = str(config["StartDate"])
    endday = str(config["EndDate"])

    # Prior emissions
    prior_cache = os.path.expandvars(
        os.path.join(config["OutputPath"], config["RunName"], "hemco_prior_emis/OutputDir")
    )

    # adjustments for when performing for dynamic kf clustering
    if kf_index is not None:
        # use different date range for KF inversion if kf_index is not None
        rundir_path = preview_dir.split("preview")[0]
        periods = pd.read_csv(f"{rundir_path}periods.csv")
        startday = str(periods.iloc[kf_index - 1]["Starts"])
        endday = str(periods.iloc[kf_index - 1]["Ends"])

        # use the nudged (prior) emissions for generating averaging kernel estimate
        sf = xr.load_dataset(f"{rundir_path}archive_sf/prior_sf_period{kf_index}.nc")
        prior_ds = get_mean_emissions(startday, endday, prior_cache)
        prior_ds = get_posterior_emissions(prior_ds, sf, config["OptimizeSoil"])
    else:
        prior_ds = get_mean_emissions(startday, endday, prior_cache)

    prior = prior_ds["EmisCH4_Total"]

    # Compute total emissions in the region of interest
    if config['UseGCHP']:
        basedir = os.path.expandvars(
            os.path.join(config["OutputPath"], config["RunName"])
        )
        gridfpath = f'{basedir}/CS_grids/grids.c{config["CS_RES"]}.nc'
        gridds = xr.open_dataset(gridfpath)
        areas = gridds['area']
    else:
        areas = prior_ds["AREA"]
    total_prior_emissions = sum_total_emissions(prior, areas, mask)
    outstring1 = (
        f"Total prior emissions in region of interest = {total_prior_emissions} Tg/y \n"
    )
    print(outstring1)


    # calculate sectoral totals if running preview
    if preview:
        get_sectoral_outputs(prior_ds, areas, mask, preview_dir)

    # ----------------------------------
    # Observations in region of interest
    # ----------------------------------

    # Paths to tropomi data files
    tropomi_files = [f for f in os.listdir(tropomi_cache) if ".nc" in f]
    tropomi_paths = [os.path.join(tropomi_cache, f) for f in tropomi_files]

    if config['UseGCHP']:
        xlim = [-180, 180]
        ylim = [-90, 90]
    else:
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
    trtime = []

    # Read in and filter tropomi observations (uses parallel processing)
    observation_dicts = Parallel(n_jobs=-1)(
        delayed(get_TROPOMI_data)(
            file_path,
            BlendedTROPOMI,
            xlim,
            ylim,
            startdate_np64,
            enddate_np64,
            use_water_obs,
        )
        for file_path in tropomi_paths
    )
    # Remove any problematic observation dicts (eg. corrupted data file)
    observation_dicts = list(filter(None, observation_dicts))

    for obs_dict in observation_dicts:
        lat.extend(obs_dict["lat"])
        lon.extend(obs_dict["lon"])
        xch4.extend(obs_dict["xch4"])
        albedo.extend(obs_dict["swir_albedo"])
        trtime.extend(obs_dict["time"])

    # Assemble in dataframe
    df = pd.DataFrame()
    df["lat"] = lat
    df["lon"] = lon
    df["obs_count"] = np.ones(len(lat))
    df["swir_albedo"] = albedo
    df["xch4"] = xch4
    df["time"] = trtime

    # Set resolution specific variables
    # L_native = Rough length scale of native state vector element [m]
    if config['UseGCHP']:
        df_super = classify_obs_to_cs_grid(df, gridfpath)
        daily_observation_counts = map_obs_to_CSgrid(df_super, gridfpath)
    else:
        # Set resolution specific variables
        # L_native = Rough length scale of native state vector element [m]
        if config["Res"] == "0.125x0.15625":
            lat_step = 0.125
            lon_step = 0.15625
        elif config["Res"] == "0.25x0.3125":
            lat_step = 0.25
            lon_step = 0.3125
        elif config["Res"] == "0.5x0.625":
            lat_step = 0.5
            lon_step = 0.625
        elif config["Res"] == "2.0x2.5":
            lat_step = 2.0
            lon_step = 2.5
        elif config["Res"] == "4.0x5.0":
            lat_step = 4.0
            lon_step = 5.0

        # bin observations into gridcells and map onto statevector
        to_lon = lambda x: np.floor(x / lon_step) * lon_step
        to_lat = lambda x: np.floor(x / lat_step) * lat_step

        df_super = df.rename(columns={"lon": "old_lon", "lat": "old_lat"})

        df_super["lat"] = to_lat(df_super.old_lat)
        df_super["lon"] = to_lon(df_super.old_lon)

        # extract relevant fields and group by lat, lon, date
        df_super = df_super[["lat", "lon", "time", "obs_count"]].copy()
        df_super["date"] = df_super["time"].dt.floor("D")
        grouped = (
            df_super.groupby(["lat", "lon", "date"]).size().reset_index(name="obs_count")
        )

        # convert the grouped DataFrame to an xarray Dataset
        daily_observation_counts = grouped.set_index(["lat", "lon", "date"]).to_xarray()

    # create a daily superobservation count as well
    daily_observation_counts["superobs_count"] = daily_observation_counts["obs_count"]

    # set the nans to 0 if there are no observations. For superobs each day is 1 superob
    daily_observation_counts["superobs_count"].values = np.where(
        np.isnan(np.array(daily_observation_counts["obs_count"].values)), 0, 1
    )
    daily_observation_counts["obs_count"] = daily_observation_counts["obs_count"].fillna(0)

    flux_per_sv, L, num_obs, m_superi = compute_sv_element_stats(
        state_vector_labels=state_vector_labels.values,
        areas=areas.values,
        prior=prior.values,
        daily_observation_counts=daily_observation_counts,
        config=config,
        mask=mask.values,
        last_ROI_element=last_ROI_element,
        sum_and_sort_along_statevector=sum_and_sort_along_statevector,
    )

    if np.sum(num_obs) < 1:
        sys.exit("Error: No observations found in region of interest")
    outstring2 = f"Found {np.sum(num_obs)} observations in the region of interest"
    print("\n" + outstring2)

    # ----------------------------------
    # Estimate information content
    # ----------------------------------

    time_delta = enddate_np64 - startdate_np64
    num_days = np.round((time_delta) / np.timedelta64(1, "D"))

    # If Kalman filter mode, count observations per inversion period
    if config["KalmanMode"]:
        startday_dt = datetime.datetime.strptime(startday, "%Y%m%d")
        endday_dt = datetime.datetime.strptime(endday, "%Y%m%d")
        if not config["MakePeriodsCSV"]:
            rundir_path = preview_dir.split("preview")[0]
            periods = pd.read_csv(f"{rundir_path}periods.csv")
            n_periods = periods.iloc[-1]["period_number"]
        else:
            n_periods = np.floor(
                (endday_dt - startday_dt).days / config["UpdateFreqDays"]
            )
        # average number of successful observation days in each inversion period
        m_superi = m_superi / n_periods
        n_obs_per_period = np.round(num_obs / n_periods)
        outstring2 = f"Found {int(np.sum(n_obs_per_period))} observations in the region of interest per inversion period, for {int(n_periods)} period(s)"
        print("\n" + outstring2)

    # Other parameters
    U = 5 * (1000 / 3600)  # 5 km/h uniform wind speed in m/s
    p = 101325  # Surface pressure [Pa = kg/m/s2]
    g = 9.8  # Gravity [m/s2]
    Mair = 0.029  # Molar mass of air [kg/mol]
    Mch4 = 0.01604  # Molar mass of methane [kg/mol]
    alpha = 0.4  # Simple parameterization of turbulence

    # Use the first element of the error list if multiple values are provided
    sigmaA = config["PriorError"][0] if isinstance(config["PriorError"], list) else config["PriorError"]
    # Error standard deviations with updated units
    sA = sigmaA * flux_per_sv
    sO = config["ObsError"][0] if isinstance(config["ObsError"], list) else config["ObsError"]

    # Calculate superobservation error to use in averaging kernel sensitivity equation
    # from P observations per grid cell = number of observations per grid cell / number of super-observations
    # P is number of observations per grid cell (native state vector element)
    P = np.array(num_obs) / m_superi
    P = np.nan_to_num(P)  # replace nan with 0
    s_superO_1 = calculate_superobservation_error(
        sO, 1
    )  # for handling cells with 0 observations (avoid divide by 0)

    # list containing superobservation error per state vector element
    s_superO_p = [
        calculate_superobservation_error(sO, element) if element >= 1.0 else s_superO_1
        for element in P
    ]
    s_superO = np.array(s_superO_p) * 1e-9  # convert to ppb

    # TODO: add eqn number from Estrada et al. 2024 once published
    # Averaging kernel sensitivity for each grid element
    # Note: m_superi is the number of superobservations,
    # defined as sum of days in each grid cell with >0 successful obs
    # in the state vector element
    # a is set to 0 where m_superi is 0
    m_superi = np.array(m_superi)
    k = alpha * (Mair * L * g / (Mch4 * U * p))
    a = sA**2 / (sA**2 + (s_superO / k) ** 2 / (m_superi))

    # Places with 0 superobs should be 0
    a = np.where(np.equal(m_superi, 0), float(0), a)

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
            f"##{outstring1}\n"
            + f"##{outstring2}\n"
            + f"##{outstring3}\n"
            + f"##{outstring4}\n"
            + outstring5
        )
        return a, df.drop(columns=["time"]), num_days, prior, outstrings
    else:
        return a

def compute_sv_element_stats(
    state_vector_labels,
    areas,
    prior,
    daily_observation_counts,
    config,
    mask,
    last_ROI_element,
    sum_and_sort_along_statevector,
):
    """
    Compute per–state-vector-element quantities:
      - emissions (sum(prior * area), kg/s)
      - L_native (sqrt(mean cell area), in same units as sqrt(areas))
      - num_native_elements (cell counts)
      - num_obs (sum of obs_count over dilated masks)
      - n_success_days (sum of superobs_count over dilated masks)

    Parameters
    ----------
    state_vector_labels : xarray.DataArray or ndarray
        State vector label grid, with NaN for background, integers for elements.
        Shape must match `areas` and the spatial dims of `daily_observation_counts`.
    areas : xarray.DataArray or ndarray
        Grid-cell areas (same shape as state_vector_labels), in m2 or whatever unit.
    prior : xarray.DataArray or ndarray
        Prior flux field in same grid and units so prior*areas is kg/s (or equivalent).
    daily_observation_counts : xarray.Dataset
        Must contain:
          - "obs_count" and "superobs_count"
          - For GCHP: dims (date, nf, Y, X) and vars "lats", "lons"
          - For non-GCHP: dims (lat, lon, date) and coords "lat", "lon"
    config : dict
        Must contain key "UseGCHP" (bool).
    mask : ndarray or xarray.DataArray (bool)
        ROI mask: True where state_vector_labels <= last_ROI_element.
    last_ROI_element : int
        Largest label index in ROI (no buffers).
    sum_and_sort_along_statevector : callable
        Function (val, sv, fill_value=np.nan) -> per-label sums, in ascending label order.

    Returns
    -------
    emissions : np.ndarray, shape (last_ROI_element,)
    L_native : np.ndarray, shape (last_ROI_element,)
    num_native_elements : np.ndarray, shape (last_ROI_element,), int
    num_obs : np.ndarray, shape (last_ROI_element,)
    n_success_days : np.ndarray, shape (last_ROI_element,)
    """

    # ------------------------------------------------------------------
    # 0. Flatten state-vector labels (background = NaN)
    # ------------------------------------------------------------------
    sv = np.asarray(state_vector_labels)
    sv_labels_flat = sv.ravel()
    Ncells = sv_labels_flat.size

    # Sanity check: ROI labels must be exactly 1..last_ROI_element
    unique_labels = np.unique(sv[mask])
    unique_labels = unique_labels[~np.isnan(unique_labels)]
    unique_labels = unique_labels.astype(int)
    assert unique_labels[0] == 1
    assert unique_labels[-1] == last_ROI_element
    assert unique_labels.size == last_ROI_element

    # ------------------------------------------------------------------
    # 1. Per–state-vector-element static quantities
    # ------------------------------------------------------------------
    areas_arr = np.asarray(areas)
    prior_arr = np.asarray(prior)

    # (a) total area per SV element (ROI + buffers, then slice ROI)
    area_per_sv_all = sum_and_sort_along_statevector(
        val=areas_arr,
        sv=sv,
    )
    area_per_sv = area_per_sv_all[:last_ROI_element]

    # (b) number of native grid cells per SV element
    ones_arr = np.ones_like(areas_arr, dtype=float)
    cell_count_per_sv_all = sum_and_sort_along_statevector(
        val=ones_arr,
        sv=sv,
    )
    cell_count_per_sv = cell_count_per_sv_all[:last_ROI_element]

    # (c) native length scale L_native = sqrt(mean cell area)
    mean_area_per_sv = area_per_sv / np.maximum(cell_count_per_sv, 1.0)
    L_native_per_sv = np.sqrt(mean_area_per_sv)

    # (d) emission flux per SV element (kg/s)
    flux_per_sv_all = sum_and_sort_along_statevector(
        val=prior_arr,
        sv=sv,
    )
    flux_per_sv = flux_per_sv_all[:last_ROI_element]

    # ------------------------------------------------------------------
    # 2. Collapse obs + superobs over time for each grid cell
    # ------------------------------------------------------------------
    use_gchp = config['UseGCHP']
    obs_da = daily_observation_counts["obs_count"]
    superobs_da = daily_observation_counts["superobs_count"]

    if use_gchp:
        # obs dims: (date, nf, Y, X)
        obs = obs_da.values
        superobs = superobs_da.values
        T, nf, Ny, Nx = obs.shape
        assert Ncells == nf * Ny * Nx

        obs_flat = obs.reshape(T, Ncells)
        superobs_flat = superobs.reshape(T, Ncells)

        obs_per_cell = np.nansum(obs_flat, axis=0)
        super_per_cell = np.nansum(superobs_flat, axis=0)

        # KDTree coordinates must match flatten order
        CSlats = daily_observation_counts["lats"].values  # (nf, Y, X)
        CSlons = daily_observation_counts["lons"].values
        kdtree, grid_shape = build_kdtree(CSlats, CSlons)
        assert grid_shape == sv.shape

        lat_flat = CSlats.ravel()
        lon_flat = CSlons.ravel()

    else:
        # obs dims: (lat, lon, date)  -> reorder to (T, Ny, Nx)
        obs_raw = obs_da.values
        superobs_raw = superobs_da.values

        obs = np.moveaxis(obs_raw, -1, 0)
        superobs = np.moveaxis(superobs_raw, -1, 0)
        T, Ny, Nx = obs.shape
        assert Ncells == Ny * Nx

        obs_flat = obs.reshape(T, Ncells)
        superobs_flat = superobs.reshape(T, Ncells)

        obs_per_cell = np.nansum(obs_flat, axis=0)
        super_per_cell = np.nansum(superobs_flat, axis=0)

        lats_1d = daily_observation_counts["lat"].values
        lons_1d = daily_observation_counts["lon"].values
        kdtree, grid_shape = build_kdtree(lats_1d, lons_1d)
        assert grid_shape == sv.shape

        lat_grid, lon_grid = np.meshgrid(lats_1d, lons_1d, indexing="ij")
        lat_flat = lat_grid.ravel()
        lon_flat = lon_grid.ravel()

    assert sv_labels_flat.shape[0] == Ncells
    assert obs_per_cell.shape[0] == Ncells

    # ------------------------------------------------------------------
    # 3. KDTree query for ROI-labeled cells only
    # ------------------------------------------------------------------
    # mask: state_vector_labels <= last_ROI_element
    sv_mask_flat = np.asarray(mask).ravel()
    sv_indices_flat = np.where(sv_mask_flat)[0]     # (M,) indices of ROI cells
    sv_labels_nonsorted = sv_labels_flat[sv_indices_flat].astype(int)  # (M,)

    # Safety: ensure only ROI labels appear
    assert sv_labels_nonsorted.min() >= 1
    assert sv_labels_nonsorted.max() <= last_ROI_element

    # Coordinates for ROI cells
    query_cart = latlon_to_cartesian(
        lat_flat[sv_indices_flat],
        lon_flat[sv_indices_flat],
    )

    # KDTree neighbor search
    n_neighbors = 25
    _, neighbor_idxs = kdtree.query(query_cart, k=n_neighbors)
    M, K = neighbor_idxs.shape  # M = #ROI cells, K = n_neighbors

    # ------------------------------------------------------------------
    # 4. UNION semantics: deduplicate (label, neighbor_cell) pairs
    # ------------------------------------------------------------------
    # For each ROI cell, repeat its label K times to align with neighbors
    labels_rep = np.repeat(sv_labels_nonsorted, K)    # (M*K,)
    cells_neighbors = neighbor_idxs.ravel()           # (M*K,)

    # Encode (label, cell) uniquely as: code = label * Ncells + cell
    labels_64 = labels_rep.astype(np.int64, copy=False)
    cells_64 = cells_neighbors.astype(np.int64, copy=False)
    pair_codes = labels_64 * np.int64(Ncells) + cells_64

    # Optimization: if each label appears only once among ROI cells,
    # duplicates are mathematically impossible, so skip np.unique.
    label_counts = np.bincount(
        sv_labels_nonsorted.astype(int),
        minlength=last_ROI_element + 1,
    )
    if label_counts[1:].max() == 1:
        pair_codes_unique = pair_codes
    else:
        pair_codes_unique = np.unique(pair_codes)

    # Decode back into (label, cell)
    labels_unique = (pair_codes_unique // Ncells).astype(int)
    cells_unique = (pair_codes_unique % Ncells).astype(int)

    # Convert 1-based labels to 0-based SV indices
    sv_indices = labels_unique - 1

    # ------------------------------------------------------------------
    # 5. Aggregate obs & superobs per SV element
    # ------------------------------------------------------------------
    num_obs_per_sv = np.bincount(
        sv_indices,
        weights=obs_per_cell[cells_unique],
        minlength=last_ROI_element,
    )
    n_success_days_per_sv = np.bincount(
        sv_indices,
        weights=super_per_cell[cells_unique],
        minlength=last_ROI_element,
    )

    # ------------------------------------------------------------------
    # 6. Return NumPy arrays aligned by SV index (0→label1, ..., N-1→labelN)
    # ------------------------------------------------------------------
    emissions = flux_per_sv                     # (Nsv,)
    L_native = L_native_per_sv                  # (Nsv,)
    num_obs = num_obs_per_sv                    # (Nsv,)
    n_success_days = n_success_days_per_sv      # (Nsv,)

    return emissions, L_native, num_obs, n_success_days

if __name__ == "__main__":
    try:
        config_path = sys.argv[1]
        state_vector_path = sys.argv[2]
        preview_dir = sys.argv[3]
        tropomi_cache = sys.argv[4]

        imi_preview(
            config_path, state_vector_path, preview_dir, tropomi_cache
        )
    except Exception as err:
        with open(".preview_error_status.txt", "w") as file1:
            # Writing data to a file
            file1.write(
                "This file is used to tell the controlling script that the imi_preview failed"
            )
        print(err)
        sys.exit(1)

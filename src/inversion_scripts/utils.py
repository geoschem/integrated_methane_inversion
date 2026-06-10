from __future__ import annotations

import os
import sys
import pickle
from datetime import datetime, timedelta
import numpy as np
import xarray as xr
import cartopy
import cartopy.crs as ccrs
from shapely.geometry.polygon import Polygon
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pyproj import Geod
import pandas as pd
import re
import warnings
from src.inversion_scripts.classify_TROPOMI_obs_to_CSgrids import(
    latlon_to_cartesian,
    build_kdtree,
)
from pathlib import Path


def save_obj(obj, name):
    """Save something with Pickle."""

    with open(name, "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def save_netcdf(ds, save_path, comp_level=1):
    """Save an xarray dataset to netcdf."""
    ds.to_netcdf(
        save_path,
        encoding={v: {"zlib": True, "complevel": comp_level} for v in ds.data_vars},
    )


def load_obj(name):
    """Load something with Pickle."""

    with open(name, "rb") as f:
        return pickle.load(f)


def zero_pad_num_hour(n):
    nstr = str(n)
    if len(nstr) == 1:
        nstr = "0" + nstr
    return nstr


def mixing_ratio_conv_factor(species):
    if species == "CH4":
        return 1e9
    elif species == "CO2":
        return 1e6
    else:
        raise ValueError(
            f"{species} is not recognized. Please add a line to "
            "mixing_ratio_conv_factor in src/inversion_scripts/utils.py"
        )


def species_molar_mass(species):
    if species == "CH4":
        return 0.01604  # Molar mass of methane [kg/mol]
    elif species == "CO2":
        return 0.04401  # Molar mass of CO2 [kg/mol]
    else:
        raise ValueError(
            f"{species} is not recognized. Please add a line to "
            "species_molar_mass in src/inversion_scripts/utils.py"
        )


def sum_total_emissions(emissions, areas, mask):
    """
    Function to sum total emissions across the region of interest.

    Arguments:
        emissions : xarray data array for emissions across inversion domain
        areas     : xarray data array for grid-cell areas across inversion domain
        mask      : xarray data array binary mask for the region of interest

    Returns:
        Total emissions in Tg/y
    """

    s_per_d = 86400
    d_per_y = 365
    tg_per_kg = 1e-9
    emissions_in_kg_per_s = emissions * areas * mask
    total = emissions_in_kg_per_s.sum() * s_per_d * d_per_y * tg_per_kg
    return float(total)


def get_default_sector_groups(ds, species="CH4", include_soil=False):
    """
    Return a simple default sector grouping for emissions variables.

    Wastewater and landfills are grouped as Waste; oil and gas are grouped as
    Oil/Gas. Other sector variables are returned individually.
    """

    prefix = f"Emis{species}_"
    grouped_suffixes = {
        "Wetlands": ["Wetlands"],
        "Livestock": ["Livestock"],
        "Waste": ["Landfills", "Wastewater"],
        "Oil/Gas": ["Oil", "Gas"],
        "Coal": ["Coal"],
        "Reservoirs": ["Reservoirs"],
        "Rice": ["Rice"],
        "Other": ["OtherAnth"],
    }

    emis_vars = [
        var
        for var in ds.data_vars
        if var.startswith(prefix)
        and "Total" not in var
        and "Excl" not in var
        and (include_soil or "Soil" not in var)
    ]

    sector_groups = {}
    used_vars = set()
    for sector, suffixes in grouped_suffixes.items():
        sector_vars = [
            f"{prefix}{suffix}"
            for suffix in suffixes
            if f"{prefix}{suffix}" in emis_vars
        ]
        if sector_vars:
            sector_groups[sector] = sector_vars
            used_vars.update(sector_vars)

    for var in sorted(emis_vars):
        if var not in used_vars:
            sector_groups[var.replace(prefix, "")] = [var]

    return sector_groups


def aggregate_state_vector_field(field, state_vector_labels, last_ROI_element):
    """
    Sum a gridded field over each ROI state-vector element.
    """

    labels = np.asarray(state_vector_labels).ravel()
    values = np.asarray(field).ravel()
    valid = (
        np.isfinite(labels)
        & np.isfinite(values)
        & (labels >= 1)
        & (labels <= last_ROI_element)
    )
    labels = labels[valid].astype(int)
    values = values[valid]
    sums = np.bincount(labels, weights=values, minlength=last_ROI_element + 1)
    return sums[1 : last_ROI_element + 1]


def build_sector_weight_matrix(
    prior_ds,
    state_vector_labels,
    last_ROI_element,
    sector_groups=None,
    areas=None,
    species="CH4",
    include_soil=False,
):
    """
    Build normalized sector weights over ROI state-vector elements.

    Rows are sectors and columns are state-vector elements from 1 through
    ``last_ROI_element``. Each nonzero sector row sums to one, so multiplying
    the matrix by ROI state-vector scale factors gives sector-level scale
    factors.
    """

    if sector_groups is None:
        sector_groups = get_default_sector_groups(
            prior_ds, species=species, include_soil=include_soil
        )
    if areas is None:
        areas = prior_ds["AREA"]

    weights = {}
    for sector, variables in sector_groups.items():
        field = sum((prior_ds[var] for var in variables if var in prior_ds), 0)
        if not hasattr(field, "dims"):
            continue
        sector_emissions = aggregate_state_vector_field(
            field * areas, state_vector_labels, last_ROI_element
        )
        total = np.sum(sector_emissions)
        if total > 0:
            weights[sector] = sector_emissions / total

    if len(weights) == 0:
        raise ValueError("No nonzero sector emissions were found in the ROI.")

    return pd.DataFrame(
        weights,
        index=np.arange(1, last_ROI_element + 1),
    ).T


def reduce_matrix_by_sector(matrix, sector_weights, use_pseudoinverse=False):
    """
    Reduce a state-vector matrix to sector space.

    Use ``use_pseudoinverse=True`` for averaging kernels (W A W*). Use the
    default for covariance-like matrices (W S W.T).
    """

    W = sector_weights.to_numpy(dtype=float)
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("Expected a two-dimensional state-vector matrix.")
    n_roi = W.shape[1]
    if matrix.shape[0] < n_roi or matrix.shape[1] < n_roi:
        raise ValueError(
            "Matrix dimensions are smaller than the number of ROI state-vector "
            "elements in sector_weights."
        )
    matrix = matrix[:n_roi, :n_roi]

    if use_pseudoinverse:
        W_star = W.T @ np.linalg.pinv(W @ W.T)
        reduced = W @ matrix @ W_star
    else:
        reduced = W @ matrix @ W.T

    return pd.DataFrame(
        reduced,
        index=sector_weights.index,
        columns=sector_weights.index,
    )


def covariance_to_correlation(covariance):
    """Convert a covariance matrix to a correlation matrix."""

    covariance = np.asarray(covariance, dtype=float)
    std = np.sqrt(np.diag(covariance))
    with np.errstate(divide="ignore", invalid="ignore"):
        correlation = covariance / np.outer(std, std)
    correlation[~np.isfinite(correlation)] = np.nan
    return correlation


def filter_obs_with_mask(mask, df, UseGCHP=False):
    """
    Select observations lying within a boolean mask
    mask is boolean xarray data array
    df is pandas dataframe with lat, lon, etc.
    """

    # Query lats/lons
    query_lats = df["lat"].values
    query_lons = df["lon"].values

    if UseGCHP:
        lats = mask["lats"].values
        lons = mask["lons"].values
        
        # Loop
        bad_ind = []
        for k in range(len(df)):
            # Find closest reference coordinates to selected lat/lon bounds
            # Build KDTree and polygons
            kdtree, shape = build_kdtree(lats, lons)
            query_cart = latlon_to_cartesian(query_lats[k], query_lons[k])
            _, neighbor_idx = kdtree.query(query_cart, k=1)
            neighbor_idx = neighbor_idx.flatten()
            f, j, i = np.unravel_index(neighbor_idx, shape)
            # If not in mask, save as bad index
            if mask[f,j,i] == 0:
                bad_ind.append(k)

    else:
        reference_lat_grid = mask["lat"].values
        reference_lon_grid = mask["lon"].values
        # Loop
        bad_ind = []
        for k in range(len(df)):
            # Find closest reference coordinates to selected lat/lon bounds
            ref_lat_ind = np.abs(reference_lat_grid - query_lats[k]).argmin()
            ref_lon_ind = np.abs(reference_lon_grid - query_lons[k]).argmin()
            # If not in mask, save as bad index
            if mask[ref_lat_ind, ref_lon_ind] == 0:
                bad_ind.append(k)

    # Drop bad indexes and count remaining entries
    df_filtered = df.copy()
    df_filtered = df_filtered.drop(df_filtered.index[bad_ind])

    return df_filtered


def count_obs_in_mask(mask, df, UseGCHP=False):
    """
    Count the number of observations in a boolean mask
    mask is boolean xarray data array
    df is pandas dataframe with lat, lon, etc.
    """

    df_filtered = filter_obs_with_mask(mask, df, UseGCHP)
    n_obs = len(df_filtered)

    return n_obs


def check_is_OH_element(sv_elem, nelements, opt_OH, is_regional):
    """
    Determine if the current state vector element is the OH element
    """
    return opt_OH and (
        ((not is_regional) and (sv_elem > (nelements - 2)))
        or (is_regional and (sv_elem == nelements))
    )


def check_is_BC_element(sv_elem, nelements, opt_OH, opt_BC, is_OH_element, is_regional):
    """
    Determine if the current state vector element is a boundary condition element
    """
    return (
        not is_OH_element
        and opt_BC
        and (
            (opt_OH and not(is_regional) and (sv_elem > (nelements - 6)))
            or (opt_OH and is_regional and (sv_elem > (nelements - 5)))
            or ((not opt_OH) and (sv_elem > (nelements - 4)))
        )
    )


def plot_hyperparameter_analysis(axs, ens_totals_posterior, params_dict):
    """
    Function to plot the sensitivity of the inversion to the hyperparameters
    """
    ens_totals_std = ens_totals_posterior.std()
    ens_mean_emis = ens_totals_posterior.mean()

    num_axes = len(params_dict.keys())

    # Only proceed if there are axes to plot
    if num_axes > 0:
        # Flatten axs in case of multiple subplots
        axs = axs.flatten()

        for i, key in enumerate(params_dict.keys()):
            ax = axs[i]
            ax.plot(
                params_dict[key],
                ens_totals_posterior,
                marker="o",
                linestyle="",
                label="Ensemble Member",
            )
            ax.axhline(y=ens_mean_emis, linestyle="-", label="Ensemble Mean")
            ax.axhline(
                y=ens_mean_emis - ens_totals_std,
                linestyle="--",
                label="Ensemble Standard Deviation",
            )
            ax.axhline(y=ens_mean_emis + ens_totals_std, linestyle="--")
            ax.set_title(f"Inversion sensitivity to {key}")
            ax.set_ylabel("Total emissions (Tg/yr)")
            ax.set_xlabel(key)
            if i == 0:
                ax.legend()

        plt.tight_layout()
    else:
        print("Not enough ensemble members to plot")


def plot_ensemble(
    ax,
    ens_posterior_totals,
    total_prior_emissions,
    default_emission_totals=None,
    plot_save_path=None,
):
    """
    Function to plot the total emissions from the ensemble members and the prior
    emissions. Optionally, the default emission totals can be plotted as well.

    Arguments
        ax                      : matplotlib axis object
        ens_posterior_totals    : list of total emissions from the ensemble members
        total_prior_emissions   : total emissions from the prior
        default_emission_totals : total emissions from the default member (optional)
        plot_save_path          : path to save the plot (optional)
    """
    # calculate the mean and standard deviation of the ensemble posterior totals
    ens_mean_emis = np.mean(ens_posterior_totals)
    ens_totals_min = np.min(ens_posterior_totals)
    ens_totals_max = np.max(ens_posterior_totals)

    # Plot the prior emissions
    ax.bar(
        0, total_prior_emissions, width=0.5, color="goldenrod", label="Prior", zorder=1
    )
    # Plot the ensemble mean and min/max error bars
    ax.bar(
        1, ens_mean_emis, width=0.5, color="steelblue", label="Ensemble mean", zorder=2
    )
    # plot the min/max error bars
    ax.errorbar(
        1,
        ens_mean_emis,
        yerr=[[ens_mean_emis - ens_totals_min], [ens_totals_max - ens_mean_emis]],
        fmt="none",
        color="k",
        capsize=4,
        zorder=5,
    )
    # Plot the ensemble members
    ax.plot(
        np.ones(len(ens_posterior_totals)),
        ens_posterior_totals,
        marker="o",
        linestyle="",
        alpha=0.5,
        color="darkblue",
        label="Ensemble member",
        zorder=3,
    )
    if default_emission_totals:
        ax.plot(
            1,
            default_emission_totals,
            marker="o",
            linestyle="",
            color="c",
            label="Default member (Ja/n closest to 1)",
            zorder=4,
        )

    # Labeling
    ax.set_xticks([1, 0])
    ax.set_xticklabels(["Posterior", "Prior"])
    ax.set_ylabel("Emissions ($Tg\ a^{-1}$)")
    ax.set_title("Total Emissions")
    ax.legend()
    if plot_save_path:
        plt.savefig(os.path.join(plot_save_path, "total_emis_ensemble.png"))


def plot_field(
    ax,
    field,
    cmap,
    plot_type="pcolormesh",
    lon_bounds=None,
    lat_bounds=None,
    levels=None,
    vmin=None,
    vmax=None,
    title=None,
    point_sources=None,
    cbar_label=None,
    mask=None,
    only_ROI=False,
    state_vector_labels=None,
    last_ROI_element=None,
    is_regional=True,
    save_path=None,
    clean_title=None,
    UseGCHP=False,
):
    """
    Function to plot inversion results.

    Arguments
        ax         : matplotlib axis object
        field      : xarray dataarray
        cmap       : colormap to use, e.g. 'viridis'
        plot_type  : 'pcolormesh' or 'imshow'
        lon_bounds : [lon_min, lon_max]
        lat_bounds : [lat_min, lat_max]
        levels     : number of colormap levels (None for continuous)
        vmin       : colorbar lower bound
        vmax       : colorbar upper bound
        title      : plot title
        point_sources: plot given point sources on map
        cbar_label : colorbar label
        mask       : mask for region of interest, boolean dataarray
        only_ROI   : zero out data outside the region of interest, true or false
        save_path  : path to save the plot
        clean_title: title without special characters
    """
    field = field.squeeze()
    # Select map features
    if is_regional:
        oceans_50m = cartopy.feature.NaturalEarthFeature("physical", "ocean", "50m")
        lakes_50m = cartopy.feature.NaturalEarthFeature("physical", "lakes", "50m")
        states_provinces_50m = cartopy.feature.NaturalEarthFeature(
            "cultural", "admin_1_states_provinces_lines", "50m"
        )
        ax.add_feature(cartopy.feature.BORDERS, facecolor="none")
        ax.add_feature(oceans_50m, facecolor="none", edgecolor="black")
        ax.add_feature(lakes_50m, facecolor="none", edgecolor="black")
        ax.add_feature(states_provinces_50m, facecolor="none", edgecolor="black")
    else:
        ax.coastlines(resolution="110m")

    # Show only ROI values?
    if only_ROI:
        field = field.where((state_vector_labels <= last_ROI_element))

    # Plot
    if plot_type == "pcolormesh":
        field.plot.pcolormesh(
            cmap=cmap,
            levels=levels,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
            cbar_kwargs={"label": cbar_label, "fraction": 0.041, "pad": 0.04},
        )
    elif plot_type == "imshow":
        field.plot.imshow(
            cmap=cmap,
            levels=levels,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
            cbar_kwargs={"label": cbar_label, "fraction": 0.041, "pad": 0.04},
        )
    else:
        raise ValueError('plot_type must be "pcolormesh" or "imshow"')

    # Zoom on ROI?
    if lon_bounds and lat_bounds:
        extent = [lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]]
        ax.set_extent(extent, crs=ccrs.PlateCarree())

    # Show boundary of ROI?
    if (mask is not None) and (not UseGCHP) and (is_regional):
        mask.plot.contour(levels=1, colors="k", linewidths=4, ax=ax)

    # Remove duplicated axis labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0)
    gl.right_labels = False
    gl.top_labels = False

    # Title
    if title:
        ax.set_title(title)

    # Marks any specified high-resolution coordinates on the preview observation density map
    if point_sources:
        for coord in point_sources:
            ax.plot(coord[1], coord[0], marker="x", markeredgecolor="black")
        point = Line2D(
            [0],
            [0],
            label="point source",
            marker="x",
            markersize=10,
            markeredgecolor="black",
            markerfacecolor="k",
            linestyle="",
        )
        ax.legend(handles=[point])

    # Save plot
    if save_path:
        # Ensure the directory exists
        os.makedirs(save_path, exist_ok=True)

        # Replace spaces in the title or clean title with underscores
        if clean_title:
            sanitized_title = clean_title.replace(" ", "_") + ".png"
        else:
            sanitized_title = title.replace(" ", "_") + ".png"

        # Construct the full file path
        full_save_path = os.path.join(save_path, sanitized_title.lower())

        # Save the plot
        plt.savefig(full_save_path, format="png", bbox_inches="tight")
        print(f"Plot saved to {full_save_path}")

def plot_field_gchp(
    ax,
    corner_lons,
    corner_lats,
    field,
    cmap,
    plot_type="pcolormesh",
    lon_bounds=None,
    lat_bounds=None,
    levels=None,
    vmin=None,
    vmax=None,
    title=None,
    point_sources=None,
    cbar_label=None,
    only_ROI=False,
    state_vector_labels=None,
    last_ROI_element=None,
    is_regional=True,
    stretch_grid=False,
    save_path=None,
    clean_title=None,
):
    """
    Function to plot inversion results.

    Arguments
        ax         : matplotlib axis object
        corner_lons: xarray dataarraycorner_lons in GCHP
        corner_lats: xarray dataarray corner_lats in GCHP
        field      : xarray dataarray
        cmap       : colormap to use, e.g. 'viridis'
        plot_type  : 'pcolormesh' or 'imshow'
        lon_bounds : [lon_min, lon_max]
        lat_bounds : [lat_min, lat_max]
        levels     : number of colormap levels (None for continuous)
        vmin       : colorbar lower bound
        vmax       : colorbar upper bound
        title      : plot title
        point_sources: plot given point sources on map
        cbar_label : colorbar label
        mask       : mask for region of interest, boolean dataarray
        only_ROI   : zero out data outside the region of interest, true or false
        save_path  : path to save the plot
        clean_title: title without special characters
    """

    # Select map features
    if (is_regional | stretch_grid):
        oceans_50m = cartopy.feature.NaturalEarthFeature("physical", "ocean", "50m")
        lakes_50m = cartopy.feature.NaturalEarthFeature("physical", "lakes", "50m")
        states_provinces_50m = cartopy.feature.NaturalEarthFeature(
            "cultural", "admin_1_states_provinces_lines", "50m"
        )
        ax.add_feature(cartopy.feature.BORDERS, facecolor="none")
        ax.add_feature(oceans_50m, facecolor="none", edgecolor="black")
        ax.add_feature(lakes_50m, facecolor="none", edgecolor="black")
        ax.add_feature(states_provinces_50m, facecolor="none", edgecolor="black")
    else:
        ax.coastlines(resolution="110m")

    # Show only ROI values?
    if only_ROI:
        field = field.where((state_vector_labels <= last_ROI_element))

    # Plot
    if plot_type == "pcolormesh":
        for face in range(6):
            x = corner_lons.isel(nf=face)
            y = corner_lats.isel(nf=face)
            v = field.squeeze().isel(nf=face)
            mesh = ax.pcolormesh(
                x, y, v, 
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                transform=ccrs.PlateCarree()
            )
        # Add colorbar
        cbar = ax.get_figure().colorbar(
            mesh,
            ax=ax,
            label=cbar_label,
            fraction=0.041,
            pad=0.04
        )
        
    elif plot_type == "imshow":
        for face in range(6):
            x = corner_lons.isel(nf=face)
            y = corner_lats.isel(nf=face)
            v = field.squeeze().isel(nf=face)

            img = ax.imshow(
                v,
                origin="lower",
                extent=[x.min(), x.max(), y.min(), y.max()],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                transform=ccrs.PlateCarree(),
                interpolation="none"  
            )

        # Mimic xarray's automatic colorbar
        cbar = ax.get_figure().colorbar(
            img,
            ax=ax,
            fraction=0.041,
            pad=0.04
        )
        cbar.set_label(cbar_label)
    else:
        raise ValueError('plot_type must be "pcolormesh" or "imshow"')

    # Zoom on ROI?
    if lon_bounds and lat_bounds:
        extent = [lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]]
        ax.set_extent(extent, crs=ccrs.PlateCarree())

    # Remove duplicated axis labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0)
    gl.right_labels = False
    gl.top_labels = False

    # Title
    if title:
        ax.set_title(title)

    # Marks any specified high-resolution coordinates on the preview observation density map
    if point_sources:
        for coord in point_sources:
            ax.plot(coord[1], coord[0], marker="x", markeredgecolor="black")
        point = Line2D(
            [0],
            [0],
            label="point source",
            marker="x",
            markersize=10,
            markeredgecolor="black",
            markerfacecolor="k",
            linestyle="",
        )
        ax.legend(handles=[point])

    # Save plot
    if save_path:
        # Ensure the directory exists
        os.makedirs(save_path, exist_ok=True)

        # Replace spaces in the title or clean title with underscores
        if clean_title:
            sanitized_title = clean_title.replace(" ", "_") + ".png"
        else:
            sanitized_title = title.replace(" ", "_") + ".png"

        # Construct the full file path
        full_save_path = os.path.join(save_path, sanitized_title.lower())

        # Save the plot
        plt.savefig(full_save_path, format="png", bbox_inches="tight")
        print(f"Plot saved to {full_save_path}")

def plot_time_series(
    x_data,
    y_data,
    line_labels,
    title,
    y_label,
    x_label="Date",
    DOFS=None,
    fig_size=(15, 6),
    x_rotation=45,
    y_sci_notation=True,
):
    """
    Function to plot inversion time series results.

    Arguments
        x_data         : x data datetimes to plot
        y_data         : list of y data to plot
        line_labels    : line label string for each y data
        title          : plot title
        y_label        : label for y axis
        x_label        : label for x axis
        DOFS           : DOFs for each interval
        fig_size       : tuple for figure size
        x_rotation     : rotation of x axis labels
        y_sci_notation : whether to use scientific notation for y axis
    """
    assert len(y_data) == len(line_labels)
    plt.clf()
    # Set the figure size
    _, ax1 = plt.subplots(figsize=fig_size)

    # Plot emissions time series
    for i in range(len(y_data)):
        # only use line for moving averages
        if "moving" in line_labels[i].lower():
            ax1.plot(x_data, y_data[i], label=line_labels[i])
        else:
            ax1.plot(
                x_data, y_data[i], linestyle="None", marker="o", label=line_labels[i]
            )

    lines = ax1.get_lines()

    # Plot DOFS time series using red
    if DOFS is not None:
        # Create a twin y-axis
        ax2 = ax1.twinx()
        ax2.plot(x_data, DOFS, linestyle="None", marker="o", color="red", label="DOFS")
        ax2.set_ylabel("DOFS", color="red")
        ax2.tick_params(axis="y", labelcolor="red")
        ax2.set_ylim(0, 1)
        # add DOFS line to legend
        lines = ax1.get_lines() + ax2.get_lines()

    # use a date string for the x axis locations
    plt.gca().xaxis.set_major_locator(mdates.WeekdayLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    # tilt the x axis labels
    plt.xticks(rotation=x_rotation)
    # scientific notation for y axis
    if y_sci_notation:
        plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    plt.title(title)
    plt.legend(lines, [line.get_label() for line in lines])
    plt.show()


def filter_tropomi(tropomi_data, xlim, ylim, startdate, enddate, use_water_obs=False):
    """
    Description:
        Filter out any data that does not meet the following
        criteria: We only consider data within lat/lon/time bounds,
        with QA > 0.5 and that don't cross the antimeridian.
        Also, we filter out pixels south of 60S and (optionally) over water.
    Returns:
        numpy array with satellite indices for filtered tropomi data.
    """
    valid_idx = (
        (tropomi_data["longitude"] > xlim[0])
        & (tropomi_data["longitude"] < xlim[1])
        & (tropomi_data["latitude"] > ylim[0])
        & (tropomi_data["latitude"] < ylim[1])
        & (tropomi_data["time"] >= startdate)
        & (tropomi_data["time"] <= enddate)
        & (tropomi_data["qa_value"] >= 0.5)
        & (tropomi_data["longitude_bounds"].ptp(axis=2) < 100)
        & ~(
            tropomi_data["surface_classification_0xF9"] == 184
        )  # exclude land_snow_or_ice
        & (tropomi_data["latitude"] > -60)
    )

    if use_water_obs:
        return np.where(valid_idx)
    else:
        return np.where(valid_idx & (tropomi_data["surface_classification"] != 1))


def filter_blended(blended_data, xlim, ylim, startdate, enddate, use_water_obs=False):
    """
    Description:
        Filter out any data that does not meet the following
        criteria: We only consider data within lat/lon/time bounds,
        that don't cross the antimeridian, and we filter out all
        coastal pixels (surface classification 3) and inland water
        pixels with a poor fit (surface classifcation 2,
        SWIR chi-2 > 20000) (recommendation from Balasus et al. 2023).
        Also, we filter out pixels south of 60S and (optionally) over water.
    Returns:
        numpy array with satellite indices for filtered tropomi data.
    """
    valid_idx = (
        (blended_data["longitude"] > xlim[0])
        & (blended_data["longitude"] < xlim[1])
        & (blended_data["latitude"] > ylim[0])
        & (blended_data["latitude"] < ylim[1])
        & (blended_data["time"] >= startdate)
        & (blended_data["time"] <= enddate)
        & (blended_data["longitude_bounds"].ptp(axis=2) < 100)
        & ~(
            (blended_data["surface_classification"] == 3)
            | (
                (blended_data["surface_classification"] == 2)
                & (blended_data["chi_square_SWIR"][:] > 20000)
            )
        )
        & ~(
            blended_data["surface_classification_0xF9"] == 184
        )  # exclude land_snow_or_ice
        & (blended_data["latitude"] > -60)
    )

    if use_water_obs:
        return np.where(valid_idx)
    else:
        return np.where(valid_idx & (blended_data["surface_classification"] != 1))


def calculate_area_in_km(coordinate_list):
    """
    Description:
        Calculate area in km of a polygon given a list of coordinates
    Arguments
        coordinate_list  [tuple]: list of lat/lon coordinates.
                         coordinates must be in correct polygon order
    Returns:
        int: area in km of polygon
    """

    polygon = Polygon(coordinate_list)

    geod = Geod(ellps="clrk66")
    poly_area, _ = geod.geometry_area_perimeter(polygon)

    return abs(poly_area) * 1e-6


def calculate_superobservation_error(sO, p):
    """
    Returns the estimated observational error accounting for superobservations.
    Using eqn (5) from Chen et al., 2023, https://doi.org/10.5194/egusphere-2022-1504
    Args:
        sO : float
            observational error specified in config file
        p  : float
            average number of observations contained within each superobservation
    Returns:
         s_super: float
            observational error for superobservations
    """
    # values from Chen et al., 2023, https://doi.org/10.5194/egusphere-2022-1504
    r_retrieval = 0.55
    s_transport = 4.5
    s_super = np.sqrt(sO**2 * (((1 - r_retrieval) / p) + r_retrieval) + s_transport**2)
    return s_super


def read_tropomi(filename):
    """
    Read TROPOMI data and save important variables to dictionary.

    Arguments
        filename [str]  : TROPOMI netcdf data file to read

    Returns
        dat      [dict] : Dictionary of important variables from TROPOMI:
                            - CH4
                            - Latitude
                            - Longitude
                            - QA value
                            - UTC time
                            - Time (utc time reshaped for orbit)
                            - Averaging kernel
                            - SWIR albedo
                            - NIR albedo
                            - Blended albedo
                            - CH4 prior profile
                            - Dry air subcolumns
                            - Latitude bounds
                            - Longitude bounds
                            - Vertical pressure profile
    """

    # Initialize dictionary for TROPOMI data
    dat = {}

    # Catch read errors in any of the variables
    try:
        # Store methane, QA, lat, lon, and time
        with xr.open_dataset(filename, group="PRODUCT") as tropomi_data:
            dat["CH4"] = tropomi_data["methane_mixing_ratio_bias_corrected"].values[
                0, :, :
            ]
            dat["qa_value"] = tropomi_data["qa_value"].values[0, :, :]
            dat["longitude"] = tropomi_data["longitude"].values[0, :, :]
            dat["latitude"] = tropomi_data["latitude"].values[0, :, :]

            utc_str = tropomi_data["time_utc"].values[0, :]
            utc_str = np.array([d.replace("Z", "") for d in utc_str]).astype(
                "datetime64[ns]"
            )
            dat["time"] = np.repeat(utc_str[:, np.newaxis], dat["CH4"].shape[1], axis=1)

        # Store column averaging kernel, SWIR and NIR surface albedo
        with xr.open_dataset(
            filename, group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS"
        ) as tropomi_data:
            dat["column_AK"] = tropomi_data["column_averaging_kernel"].values[
                0, :, :, ::-1
            ]
            dat["swir_albedo"] = tropomi_data["surface_albedo_SWIR"].values[0, :, :]
            dat["nir_albedo"] = tropomi_data["surface_albedo_NIR"].values[0, :, :]
            dat["blended_albedo"] = 2.4 * dat["nir_albedo"] - 1.13 * dat["swir_albedo"]

        # Store methane prior profile, dry air subcolumns
        with xr.open_dataset(
            filename, group="PRODUCT/SUPPORT_DATA/INPUT_DATA"
        ) as tropomi_data:
            dat["profile_apriori"] = tropomi_data["methane_profile_apriori"].values[
                0, :, :, ::-1
            ]  # mol m-2
            dat["dry_air_subcolumns"] = tropomi_data["dry_air_subcolumns"].values[
                0, :, :, ::-1
            ]  # mol m-2

            # Surface classification values of NaN will be filtered out
            # due to a low QA value. We can make sure by setting them to 5
            # and making sure nothing outside of [0,1,2,3] makes it through.
            sc = tropomi_data["surface_classification"].values[0, :, :]
            nan_mask = np.isnan(sc)
            sc_no_nans = np.nan_to_num(sc, nan=0)
            sc = (sc_no_nans.astype("uint8") & 0x03).astype(int)
            sc[nan_mask] = 5
            dat["surface_classification"] = sc
            sc_0xF9 = (sc_no_nans.astype("uint8") & 0xF9).astype(int)
            sc_0xF9[nan_mask] = 5
            dat["surface_classification_0xF9"] = sc_0xF9

            # Also get pressure interval and surface pressure for use below
            pressure_interval = (
                tropomi_data["pressure_interval"].values[0, :, :] / 100
            )  # Pa -> hPa
            surface_pressure = (
                tropomi_data["surface_pressure"].values[0, :, :] / 100
            )  # Pa -> hPa

        # Store latitude and longitude bounds for pixels
        with xr.open_dataset(
            filename, group="PRODUCT/SUPPORT_DATA/GEOLOCATIONS"
        ) as tropomi_data:
            dat["longitude_bounds"] = tropomi_data["longitude_bounds"].values[
                0, :, :, :
            ]
            dat["latitude_bounds"] = tropomi_data["latitude_bounds"].values[0, :, :, :]

        # Store vertical pressure profile
        n1 = dat["CH4"].shape[
            0
        ]  # length of along-track dimension (scanline) of retrieval field
        n2 = dat["CH4"].shape[
            1
        ]  # length of across-track dimension (ground_pixel) of retrieval field
        pressures = np.full([n1, n2, 12 + 1], np.nan, dtype=np.float32)
        for i in range(12 + 1):
            pressures[:, :, i] = surface_pressure - i * pressure_interval
        dat["pressures"] = pressures

    # Return an error if any of the variables were not read correctly
    except Exception as e:
        print(f"Error opening {filename}: {e}")
        return None

    return dat


def read_blended(filename):
    """
    Read Blended TROPOMI+GOSAT data and save important variables to dictionary.
    Arguments
        filename [str]  : Blended TROPOMI+GOSAT netcdf data file to read
    Returns
        dat      [dict] : Dictionary of important variables from Blended TROPOMI+GOSAT:
                            - CH4
                            - Latitude
                            - Longitude
                            - Time (utc time reshaped for orbit)
                            - Averaging kernel
                            - SWIR albedo
                            - NIR albedo
                            - Blended albedo
                            - CH4 prior profile
                            - Dry air subcolumns
                            - Latitude bounds
                            - Longitude bounds
                            - Surface classification
                            - Chi-Square for SWIR
                            - Vertical pressure profile
    """
    assert (
        "BLND" in filename
    ), f"BLND not in filename {filename}, but a blended function is being used"
    try:
        # Initialize dictionary for Blended TROPOMI+GOSAT data
        dat = {}

        # Extract data from netCDF file to our dictionary
        with xr.open_dataset(filename) as blended_data:
            dat["CH4"] = blended_data["methane_mixing_ratio_blended"].values[:]
            dat["longitude"] = blended_data["longitude"].values[:]
            dat["latitude"] = blended_data["latitude"].values[:]
            dat["column_AK"] = blended_data["column_averaging_kernel"].values[:, ::-1]
            dat["swir_albedo"] = blended_data["surface_albedo_SWIR"][:]
            dat["nir_albedo"] = blended_data["surface_albedo_NIR"].values[:]
            dat["blended_albedo"] = 2.4 * dat["nir_albedo"] - 1.13 * dat["swir_albedo"]
            dat["profile_apriori"] = blended_data["methane_profile_apriori"].values[
                :, ::-1
            ]
            dat["dry_air_subcolumns"] = blended_data["dry_air_subcolumns"].values[
                :, ::-1
            ]
            dat["longitude_bounds"] = blended_data["longitude_bounds"].values[:]
            dat["latitude_bounds"] = blended_data["latitude_bounds"].values[:]
            dat["surface_classification"] = (
                blended_data["surface_classification"].values[:].astype("uint8") & 0x03
            ).astype(int)
            dat["surface_classification_0xF9"] = (
                blended_data["surface_classification"].values[:].astype("uint8") & 0xF9
            ).astype(int)
            dat["chi_square_SWIR"] = blended_data["chi_square_SWIR"].values[:]

            # Remove "Z" from time so that numpy doesn't throw a warning
            utc_str = blended_data["time_utc"].values[:]
            dat["time"] = np.array([d.replace("Z", "") for d in utc_str]).astype(
                "datetime64[ns]"
            )

            # Need to calculate the pressure for the 13 TROPOMI levels (12 layer edges)
            pressure_interval = (
                blended_data["pressure_interval"].values[:] / 100
            )  # Pa -> hPa
            surface_pressure = (
                blended_data["surface_pressure"].values[:] / 100
            )  # Pa -> hPa
            n = len(dat["CH4"])
            pressures = np.full([n, 12 + 1], np.nan, dtype=np.float32)
            for i in range(12 + 1):
                pressures[:, i] = surface_pressure - i * pressure_interval
            dat["pressures"] = pressures

        # Add an axis here to mimic the (scanline, groundpixel) format of operational TROPOMI data
        # This is so the blended data will be compatible with the TROPOMI operators
        for key in dat.keys():
            dat[key] = np.expand_dims(dat[key], axis=0)

    except Exception as e:
        print(f"Error opening {filename}: {e}")
        return None

    return dat


def read_and_filter_satellite(
    filename, satellite_str, gc_startdate, gc_enddate, xlim, ylim, use_water_obs
):

    # Read TROPOMI data
    if satellite_str == "BlendedTROPOMI":
        satellite = read_blended(filename)
    elif satellite_str == "TROPOMI":
        satellite = read_tropomi(filename)
    else:
        print("Other data source is not currently supported")
        sys.exit(1)

    # If empty, skip this file
    if satellite == None:
        print(f"Skipping {filename} due to file processing issue.")
        return satellite

    # Filter the data
    if satellite_str == "BlendedTROPOMI":
        # Only going to consider blended data within lat/lon/time bounds and wihtout problematic coastal pixels
        sat_ind = filter_blended(
            satellite, xlim, ylim, gc_startdate, gc_enddate, use_water_obs
        )
    elif satellite_str == "TROPOMI":
        # Only going to consider TROPOMI data within lat/lon/time bounds and with QA > 0.5
        sat_ind = filter_tropomi(
            satellite, xlim, ylim, gc_startdate, gc_enddate, use_water_obs
        )

    else:
        print("Other data source filtering is not currently supported --HON")
        sys.exit(1)

    return satellite, sat_ind


def get_posterior_emissions(prior, scale, species, OptimizeSoil=False):
    """
    Function to calculate the posterior emissions from the prior
    and the scale factors. Properly accounting for no optimization
    of the soil sink.
    Args:
        prior  : xarray dataset
            prior emissions
        scales : xarray dataset or datarray of scale factors
    Returns:
        posterior : xarray dataset
            posterior emissions
    """
    # Make copies to avoid modifying the original data
    prior = prior.copy()
    scale = scale.copy()

    # keep attributes of data even when arithmetic operations applied
    xr.set_options(keep_attrs=True)

    # if xarray datarray
    if isinstance(scale, xr.DataArray):
        scale_factors = scale
    # if xarray dataset
    elif isinstance(scale, xr.Dataset):
        scale_factors = scale["ScaleFactor"]
    else:
        raise ValueError("Scale factors must be an xarray DataArray or Dataset")

    posterior = prior.copy()
    if not OptimizeSoil:
        # we do not optimize soil absorbtion in the inversion. This
        # means that we need to keep the soil sink constant and properly
        # account for it in the posterior emissions calculation.
        # To do this, we:
        # make a copy of the original soil sink
        prior_soil_sink = prior[f"Emis{species}_SoilAbsorb"].copy()

        filtered_keys = [
            key
            for key in prior.keys()
            if "EmisCH4" in key
            and key != "EmisCH4_Total"
            and key != "EmisCH4_SoilAbsorb"
        ]
        # scale the prior emissions for all sectors except soil using the scale factors
        for ds_var in filtered_keys:
            posterior[ds_var] = prior[ds_var] * scale_factors

        # But reset the soil sink to the original value
        posterior[f"Emis{species}_SoilAbsorb"] = prior_soil_sink

        # Add the original soil sink back to the total emissions
        posterior[f"Emis{species}_Total"] = (
            posterior[f"Emis{species}_Total_ExclSoilAbs"]
            + posterior[f"Emis{species}_SoilAbsorb"]
        )
    else:
        filtered_keys = [key for key in prior.keys() if "EmisCH4" in key]
        # scale the prior emissions for all sectors using the scale factors
        for ds_var in filtered_keys:
            posterior[ds_var] = prior[ds_var] * scale_factors

    return posterior

def get_strdate(current_time, date_threshold):
    # round observation time to nearest hour
    strdate = current_time.round("60min").strftime("%Y%m%d_%H")
    # Unless it equals the date threshold (hour 00 after the inversion period)
    if strdate == date_threshold:
        strdate = current_time.floor("60min").strftime("%Y%m%d_%H")

    return strdate


def filter_prior_files(filenames, start_date, end_date):
    """
    Filter a list of HEMCO diagnostic files based on the specified date range.
    """
    # Parse the input dates
    start_date = datetime.strptime(start_date, "%Y%m%d")
    end_date = datetime.strptime(end_date, "%Y%m%d") - timedelta(days=1)

    filtered_files = []
    for file in filenames:
        match1 = re.search(r"\.(\d{8}_\d{4})z", file)
        match2 = re.search(r"\.(\d{12})", file)
        # Extract the date part from the filename
        file_date = None
        if match1:
            file_date = datetime.strptime(match1.group(1), "%Y%m%d_%H%M")
        elif match2:
            file_date = datetime.strptime(match2.group(1), "%Y%m%d%H%M")
        
        # Check if the file date is within the specified range
        if start_date <= file_date <= end_date:
            filtered_files.append(file)

    return filtered_files


def get_mean_emissions(start_date, end_date, prior_cache_path):
    """
    Calculate the mean emissions for the specified date range.
    """
    # find all prior files in the specified date range
    prior_files = [
        f for f in os.listdir(prior_cache_path)
        if "HEMCO_sa_diagnostics" in f or "GEOSChem.Emissions" in f
    ]
    prior_files = filter_prior_files(prior_files, str(start_date), str(end_date))

    # drop anchor with duplicate dimensions of ncontact from GCHP outputs
    hemco_diags = []
    for f in prior_files:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
            ds = xr.load_dataset(os.path.join(prior_cache_path, f))

        if 'anchor' in ds.data_vars:
            ds = ds.drop_vars('anchor')
        hemco_diags.append(ds)

    # concatenate all datasets and aggregate into the mean prior
    # emissions for the specified date range
    prior_ds = xr.concat(hemco_diags, dim="time")
    return prior_ds.mean(dim=["time"])


def get_period_mean_emissions(prior_cache_path, period, periods_csv_path):
    """
    Calculate the mean emissions for the specified kalman period.
    """
    period_df = pd.read_csv(periods_csv_path)
    period_df = period_df[period_df["period_number"].astype(int) == int(period)]
    period_df = period_df.reset_index(drop=True)
    start_date = str(period_df.loc[0, "Starts"])
    end_date = str(period_df.loc[0, "Ends"])
    return get_mean_emissions(start_date, end_date, prior_cache_path)


def ensure_float_list(variable):
    """Make sure the variable is a list of floats."""
    if isinstance(variable, list):
        # Convert each item in the list to a float
        return [float(item) for item in variable]
    elif isinstance(variable, (str, float, int)):
        # Wrap the variable in a list and convert it to float
        return [float(variable)]
    else:
        raise TypeError("Variable must be a string, float, int, or list.")


def compute_min_max_ensemble_sectors(
    prior_ds: xr.Dataset,
    ens_posterior_ds: xr.Dataset,
    ens_scale_ds: xr.Dataset,
    areas: xr.DataArray,
    mask: xr.DataArray,
    num_ensemble_members: int,
) -> tuple[dict[str, float], dict[str, float]]:
    """
    Compute the minimum and maximum posterior emissions for each sector across the ensemble.

    Returns
    -------
    sector_mins : dict
        Minimum posterior emissions for each sector.
    sector_maxes : dict
        Maximum posterior emissions for each sector.
    """
    # Identify sectors (EmisCH4 variables, excluding totals/exclusions)
    sectors = [
        var
        for var in list(ens_posterior_ds.keys())
        if "EmisCH4" in var and not ("Total" in var or "Excl" in var)
    ]

    sector_mins = {sector: float("inf") for sector in sectors}
    sector_maxes = {sector: float("-inf") for sector in sectors}

    for member in range(num_ensemble_members):
        active_ds = get_posterior_emissions(
            prior_ds, ens_scale_ds.isel(ensemble=member)
        )

        for sector in sectors:
            post_val = sum_total_emissions(active_ds[sector], areas, mask)
            # Update sector min and max values
            if np.isfinite(post_val):
                sector_mins[sector] = min(sector_mins[sector], post_val)
                sector_maxes[sector] = max(sector_maxes[sector], post_val)

    # Remove sectors with inf or -inf values from the results
    sector_mins = {
        f"{k.replace('EmisCH4_', '')}_ensMin": v
        for k, v in sector_mins.items()
        if np.isfinite(v)
    }
    sector_maxes = {
        f"{k.replace('EmisCH4_', '')}_ensMax": v
        for k, v in sector_maxes.items()
        if np.isfinite(v)
    }

    return sector_mins, sector_maxes


def export_visualization_outputs(
    *,
    plot_save_path: str,
    start_date,
    end_date,
    total_prior_emissions: float,
    total_posterior_emissions: float,
    DOFS: float,
    prior_bias: float,
    posterior_bias: float,
    prior_std: float,
    posterior_std: float,
    prior_RMSE: float,
    posterior_RMSE: float,
    df_for_count: pd.DataFrame,
    mask: xr.DataArray,
    combined_sorted: list,
    prior_ds: xr.Dataset,
    posterior_ds: xr.Dataset,
    areas: xr.DataArray,
    # Ensemble (optional)
    num_ensemble_members: int = 0,
    ens_totals_posterior: np.ndarray | None = None,
    ens_inv_result: xr.Dataset | None = None,
    ens_scale_ds: xr.Dataset | None = None,
    # Artifacts to save (optional)
    state_vector_labels: xr.DataArray | None = None,
    scale: xr.DataArray | None = None,
    avkern_ROI: xr.DataArray | None = None,
    ds_obs: xr.Dataset | None = None,
    ds_obs_regridded: xr.Dataset | None = None,
) -> pd.DataFrame:
    """
    Export statistics CSVs and NetCDF artifacts for the visualization notebook.

    Returns
    -------
    pandas.DataFrame
        The single-row statistics DataFrame written to `viz_notebook_statistics.csv`.
    """
    # Create output dirs
    os.makedirs(plot_save_path, exist_ok=True)

    # 1) Build the statistics dictionary (mirrors the notebook)
    statistics: dict[str, object] = {
        "StartDate": start_date,
        "EndDate": end_date,
        "PriorEmissions": float(total_prior_emissions),
        "PosteriorEmissions": float(total_posterior_emissions),
        "DOFS": float(DOFS),
        "BiasPrior": float(prior_bias),
        "BiasPosterior": float(posterior_bias),
        "BiasPrior_std": float(prior_std),
        "BiasPosterior_std": float(posterior_std),
        "RMSEPrior": float(prior_RMSE),
        "RMSEPosterior": float(posterior_RMSE),
        "NumSuperObservations": int(count_obs_in_mask(mask, df_for_count)),
    }

    # 2) Sector totals (use provided combined_sorted and append domain totals like the notebook)
    sector_totals: dict[str, float] = {}
    combined_plus = list(combined_sorted)

    if "EmisCH4_Total" in prior_ds and "EmisCH4_Total" in posterior_ds:
        combined_plus.append(
            (
                "Total",
                float(sum_total_emissions(prior_ds["EmisCH4_Total"], areas, mask)),
                float(sum_total_emissions(posterior_ds["EmisCH4_Total"], areas, mask)),
            )
        )

    if (
        "EmisCH4_Total_ExclSoilAbs" in prior_ds
        and "EmisCH4_Total_ExclSoilAbs" in posterior_ds
    ):
        combined_plus.append(
            (
                "TotalExclSoilAbs",
                float(
                    sum_total_emissions(
                        prior_ds["EmisCH4_Total_ExclSoilAbs"], areas, mask
                    )
                ),
                float(
                    sum_total_emissions(
                        posterior_ds["EmisCH4_Total_ExclSoilAbs"], areas, mask
                    )
                ),
            )
        )

    for name, prior_val, post_val in combined_plus:
        sector_totals[f"{name}Prior"] = float(prior_val)
        sector_totals[f"{name}Posterior"] = float(post_val)

    statistics.update(sector_totals)

    # 3) Ensemble stats (optional) + CSV
    if (
        num_ensemble_members
        and num_ensemble_members > 1
        and ens_totals_posterior is not None
        and ens_inv_result is not None
    ):

        statistics["EnsemblePosterior_Mean"] = float(np.mean(ens_totals_posterior))
        statistics["EnsemblePosterior_Std"] = float(np.std(ens_totals_posterior))
        statistics["EnsemblePosterior_Max"] = float(np.max(ens_totals_posterior))
        statistics["EnsemblePosterior_Min"] = float(np.min(ens_totals_posterior))
        statistics["EnsemblePosteriorsArray"] = ",".join(
            [str(x) for x in ens_totals_posterior]
        )

        rows = []
        value_names = [
            "Ja_normalized",
            "prior_err",
            "obs_err",
            "gamma",
            "prior_err_bc",
            "prior_err_oh",
        ]
        for member in ens_inv_result.ensemble.values:
            row = {"EnsembleMember": int(member)}
            # Mirror notebook's indexing assumption (member integer indexes into ens_totals_posterior)
            try:
                row["posterior"] = float(ens_totals_posterior[int(member)])
            except Exception:
                # Fallback to last element if out-of-range / non-contiguous IDs
                row["posterior"] = float(ens_totals_posterior[-1])
            for vn in value_names:
                if vn in ens_inv_result:
                    row[vn] = np.asarray(ens_inv_result.sel(ensemble=member)[vn]).item()
            rows.append(row)

        pd.DataFrame(rows).to_csv(
            os.path.join(plot_save_path, "ensemble_statistics.csv"), index=False
        )

        sector_mins, sector_maxes = compute_min_max_ensemble_sectors(
            prior_ds,
            get_posterior_emissions(prior_ds, ens_scale_ds.isel(ensemble=0)),
            ens_scale_ds,
            areas,
            mask,
            num_ensemble_members,
        )

        statistics.update(sector_mins)
        statistics.update(sector_maxes)

    # 4) Write the main statistics CSV
    stats_pd = pd.DataFrame(statistics, index=[0])
    stats_pd.to_csv(
        os.path.join(plot_save_path, "viz_notebook_statistics.csv"), index=False
    )

    # 5) NetCDF artifacts (optional)
    netcdf_path = os.path.join(plot_save_path, "netCDF")
    os.makedirs(netcdf_path, exist_ok=True)

    if state_vector_labels is not None:
        state_vector_labels.to_netcdf(os.path.join(netcdf_path, "state_vector.nc"))
    if prior_ds is not None:
        prior_ds.to_netcdf(os.path.join(netcdf_path, "prior.nc"))
    if posterior_ds is not None:
        posterior_ds.to_netcdf(os.path.join(netcdf_path, "posterior.nc"))
    if scale is not None:
        scale.to_netcdf(os.path.join(netcdf_path, "scale.nc"))
    if avkern_ROI is not None:
        avkern_ROI.to_netcdf(os.path.join(netcdf_path, "avkern.nc"))
    if ds_obs is not None:
        ds_obs.to_netcdf(os.path.join(netcdf_path, "obs_ds.nc"))
    if ds_obs_regridded is not None:
        ds_obs_regridded.to_netcdf(os.path.join(netcdf_path, "obs_ds_regridded.nc"))

    return stats_pd


def update_prior_error_for_OptimizeSoil(
    prior_ds, org_prior_error, StateVectorFile, n_elements
):
    """
    Update prior error for the case when OptimizeSoil is turned on.

    Args:
        prior_ds (xarray.Dataset): prior emission dataset
        org_prior_error (float): relative prior error
        StateVectorFile (str): Path to gridded state vector file
        n_elements (int): number of state vector elements

    Returns:
        np.ndarray: Updated relative prior error for each state vector element
    """
    prior_soil = prior_ds["EmisCH4_SoilAbsorb"].values
    prior_flux = prior_ds["EmisCH4_Total"].values
    prior_emis = prior_flux - prior_soil

    state_vector = xr.open_dataset(StateVectorFile).squeeze()
    state_vector_labels = state_vector["StateVector"].fillna(-9999).values.astype(int)
    last_ROI_element = np.nanmax(state_vector_labels)

    prior_err = np.zeros(n_elements)

    for i in range(1, last_ROI_element + 1):
        mask = state_vector_labels == i

        # mean emissions & soil sinks for this state vector element
        emisi = np.nanmean(prior_emis[mask])
        soili = np.nanmean(prior_soil[mask])
        fluxi = np.nanmean(prior_flux[mask])

        if abs(fluxi) > 0:
            prior_err[i - 1] = (
                np.sqrt((org_prior_error * emisi) ** 2 + (org_prior_error * soili) ** 2)
                / fluxi
            )
    return prior_err


def get_orbit_num(filename: str | Path) -> int:
    """
    Extract the orbit number from a TROPOMI or blended data filename.
    Assumes the filename format includes the orbit number as the 5th element from the end
    when split by underscores (e.g., S5P_RPRO_L2__CH4____20210704T174744_20210704T192914_19300_03_020400_20221126T123158_GCtoTROPOMI.pkl).

    Args:
        filename: The filename from which to extract the orbit number.

    Returns:
        The extracted orbit number as an integer.
    """

    name = Path(filename).name
    parts = name.split("_")

    try:
        orbit_num = parts[-5]
    except IndexError as exc:
        raise ValueError(f"Unexpected filename format: {filename}") from exc

    if not orbit_num.isdigit():
        raise ValueError(f"Orbit number not found in filename: {filename}")

    return int(orbit_num)


def build_reference_index(ref_dir: str | Path) -> dict[int, Path]:
    """
    Build an index mapping orbit numbers to reference file paths for quick lookup.

    Args:
        ref_dir: Directory containing reference files (e.g., data_converted_reference).

    Returns:
        A dictionary mapping orbit numbers to their corresponding reference file paths.
    """
    index = {}
    for path in Path(ref_dir).glob("*.pkl"):
        orbit = get_orbit_num(path)
        if orbit in index:
            raise ValueError(f"Duplicate reference files found for orbit {orbit}")
        index[orbit] = path
    return index


def map_files_to_reference(
    target_dir: str | Path, ref_dir: str | Path
) -> dict[Path, Path]:
    """
    Map target files to their corresponding reference files based on orbit numbers.
    Assumes that both target and reference files have the orbit number in the same position in their filenames.

    Args:
        target_dir: Directory containing target files (e.g., data_converted).
        ref_dir: Directory containing reference files (e.g., data_converted_reference).

    Returns:
        A dictionary mapping each target file path to its corresponding reference file path.
    """

    ref_index = build_reference_index(ref_dir)
    mapping = {}

    for target_file in Path(target_dir).glob("*.pkl"):
        try:
            orbit = get_orbit_num(target_file)
            mapping[target_file] = ref_index[orbit]
        except KeyError:
            mapping[target_file] = None
            print(f"No reference file found for orbit {orbit} in {target_file}")
        except ValueError as e:
            print(e)

    return mapping

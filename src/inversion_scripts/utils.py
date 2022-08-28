import numpy as np
import xarray as xr
from functools import partial
import pyproj
from shapely.geometry.polygon import Polygon
import shapely.ops as ops
import cartopy
import cartopy.crs as ccrs
import pickle


def save_obj(obj, name):
    """Save something with Pickle."""

    with open(name, "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    """Load something with Pickle."""

    with open(name, "rb") as f:
        return pickle.load(f)


def zero_pad_num_hour(n):
    nstr = str(n)
    if len(nstr) == 1:
        nstr = "0" + nstr
    return nstr


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


def filter_obs_with_mask(mask, df):
    """
    Select observations lying within a boolean mask
    mask is boolean xarray data array
    df is pandas dataframe with lat, lon, etc.
    """

    reference_lat_grid = mask["lat"].values
    reference_lon_grid = mask["lon"].values

    # Query lats/lons
    query_lats = df["lat"].values
    query_lons = df["lon"].values

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


def count_obs_in_mask(mask, df):
    """
    Count the number of observations in a boolean mask
    mask is boolean xarray data array
    df is pandas dataframe with lat, lon, etc.
    """

    df_filtered = filter_obs_with_mask(mask, df)
    n_obs = len(df_filtered)

    return n_obs


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
    cbar_label=None,
    mask=None,
    only_ROI=False,
    state_vector_labels=None,
    last_ROI_element=None,
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
        cbar_label : colorbar label
        mask       : mask for region of interest, boolean dataarray
        only_ROI   : zero out data outside the region of interest, true or false
    """

    # Select map features
    oceans_50m = cartopy.feature.NaturalEarthFeature("physical", "ocean", "50m")
    lakes_50m = cartopy.feature.NaturalEarthFeature("physical", "lakes", "50m")
    states_provinces_50m = cartopy.feature.NaturalEarthFeature(
        "cultural", "admin_1_states_provinces_lines", "50m"
    )
    ax.add_feature(cartopy.feature.BORDERS, facecolor="none")
    ax.add_feature(oceans_50m, facecolor=[1, 1, 1], edgecolor="black")
    ax.add_feature(lakes_50m, facecolor=[1, 1, 1], edgecolor="black")
    ax.add_feature(states_provinces_50m, facecolor="none", edgecolor="black")

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
    if mask is not None:
        mask.plot.contour(levels=1, colors="k", linewidths=4, ax=ax)

    # Remove duplicated axis labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0)
    gl.right_labels = False
    gl.top_labels = False

    # Title
    if title:
        ax.set_title(title)


def filter_tropomi(tropomi_data, xlim, ylim, startdate, enddate):
    """
    Description:
        Filter out any data that does not meet the following
        criteria: We only consider data within lat/lon/time bounds,
        with QA > 0.5, and with safe surface albedo values
    Returns:
        numpy array with satellite indices for filtered tropomi data.
    """
    return np.where(
        (tropomi_data["longitude"] > xlim[0])
        & (tropomi_data["longitude"] < xlim[1])
        & (tropomi_data["latitude"] > ylim[0])
        & (tropomi_data["latitude"] < ylim[1])
        & (tropomi_data["time"] >= startdate)
        & (tropomi_data["time"] <= enddate)
        & (tropomi_data["qa_value"] >= 0.5)
        & (tropomi_data["swir_albedo"] > 0.05)
        & (tropomi_data["blended_albedo"] < 0.85)
    )

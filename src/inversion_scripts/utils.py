import os
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
    point_sources=None,
    cbar_label=None,
    mask=None,
    only_ROI=False,
    state_vector_labels=None,
    last_ROI_element=None,
    is_regional=True
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
    """

    # Select map features
    if is_regional:
        oceans_50m = cartopy.feature.NaturalEarthFeature("physical", "ocean", "50m")
        lakes_50m = cartopy.feature.NaturalEarthFeature("physical", "lakes", "50m")
        states_provinces_50m = cartopy.feature.NaturalEarthFeature(
            "cultural", "admin_1_states_provinces_lines", "50m"
        )
        ax.add_feature(cartopy.feature.BORDERS, facecolor="none")
        ax.add_feature(oceans_50m, facecolor=[1, 1, 1], edgecolor="black")
        ax.add_feature(lakes_50m, facecolor=[1, 1, 1], edgecolor="black")
        ax.add_feature(states_provinces_50m, facecolor="none", edgecolor="black")
    else:
        ax.coastlines(resolution='110m')

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
    
    # Marks any specified high-resolution coordinates on the preview observation density map
    if point_sources:
        for coord in point_sources:
            ax.plot(coord[1], coord[0], marker="x", markeredgecolor="black")
        point = Line2D([0], [0], label='point source', marker='x', markersize=10, markeredgecolor='black', markerfacecolor='k', linestyle='')
        ax.legend(handles=[point])


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


def filter_tropomi(tropomi_data, xlim, ylim, startdate, enddate):
    """
    Description:
        Filter out any data that does not meet the following
        criteria: We only consider data within lat/lon/time bounds,
        with QA > 0.5 and that don't cross the antimeridian.
        Also, we filter out water pixels and south of 60S.
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
        & (tropomi_data["longitude_bounds"].ptp(axis=2) < 100)
        & (tropomi_data["surface_classification"] != 1)
        & (tropomi_data["latitude"] > -60)
    )

def filter_blended(blended_data, xlim, ylim, startdate, enddate):
    """
    Description:
        Filter out any data that does not meet the following
        criteria: We only consider data within lat/lon/time bounds,
        that don't cross the antimeridian, and we filter out all
        coastal pixels (surface classification 3) and inland water
        pixels with a poor fit (surface classifcation 2, 
        SWIR chi-2 > 20000) (recommendation from Balasus et al. 2023).
        Also, we filter out water pixels and south of 60S.
    Returns:
        numpy array with satellite indices for filtered tropomi data.
    """
    return np.where(
        (blended_data["longitude"] > xlim[0])
        & (blended_data["longitude"] < xlim[1])
        & (blended_data["latitude"] > ylim[0])
        & (blended_data["latitude"] < ylim[1])
        & (blended_data["time"] >= startdate)
        & (blended_data["time"] <= enddate)
        & (blended_data["longitude_bounds"].ptp(axis=2) < 100)
        & ~((blended_data["surface_classification"] == 3) | ((blended_data["surface_classification"] == 2) & (blended_data["chi_square_SWIR"][:] > 20000)))
        & (blended_data["surface_classification"] != 1)
        & (blended_data["latitude"] > -60)
    )


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
    s_super = np.sqrt(
        sO**2 * (((1 - r_retrieval) / p) + r_retrieval) + s_transport**2
    )
    return s_super

def get_posterior_emissions(prior, scale):
    """
    Function to calculate the posterior emissions from the prior 
    and the scale factors. Properly accounting for no optimization 
    of the soil sink.
    Args:
        prior  : xarray dataset
            prior emissions
        scales : xarray dataset scale factors
    Returns:
        posterior : xarray dataset
            posterior emissions
    """
    # we do not optimize soil absorbtion in the inversion. This 
    # means that we need to keep the soil sink constant and properly 
    # account for it in the posterior emissions calculation.
    # To do this, we:
    
    # make a copy of the original soil sink
    prior_soil_sink = prior["EmisCH4_SoilAbsorb"].copy()
    
    # remove the soil sink from the prior total before applying scale factors
    prior["EmisCH4_Total"] = prior["EmisCH4_Total"] - prior_soil_sink
    
    # scale the prior emissions for all sectors using the scale factors
    posterior = prior * scale["ScaleFactor"]
    
    # But reset the soil sink to the original value
    posterior["EmisCH4_SoilAbsorb"] = prior_soil_sink
    
    # Add the original soil sink back to the total emissions
    posterior["EmisCH4_Total"] = posterior["EmisCH4_Total"] + prior_soil_sink
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
        # Extract the date part from the filename
        date_str = file.split('.')[1]
        file_date = datetime.strptime(date_str, "%Y%m%d%H%M")

        # Check if the file date is within the specified range
        if start_date <= file_date <= end_date:
            filtered_files.append(file)
    
    return filtered_files
    
def get_mean_emissions(start_date, end_date, prior_cache_path):
    """
    Calculate the mean emissions for the specified date range.
    """
    # find all prior files in the specified date range
    prior_files = [f for f in os.listdir(prior_cache_path) if "HEMCO_sa_diagnostics" in f]
    prior_files = filter_prior_files(prior_files, str(start_date), str(end_date))
    hemco_diags = [xr.load_dataset(os.path.join(prior_cache_path, f)) for f in prior_files]
    
    # concatenate all datasets and aggregate into the mean prior 
    # emissions for the specified date range
    prior_ds = xr.concat(hemco_diags, dim="time")
    return prior_ds.mean(dim=["time"])

def get_period_mean_emissions(prior_cache_path, period, periods_csv_path):
    """
    Calculate the mean emissions for the specified kalman period.
    """
    period_df = pd.read_csv(periods_csv_path)
    period_df = period_df[period_df['period_number'] == period]
    period_df.reset_index(drop=True, inplace=True)
    start_date = str(period_df.loc[0,"Starts"])
    end_date = str(period_df.loc[0,"Ends"])
    return get_mean_emissions(start_date, end_date, prior_cache_path)
    
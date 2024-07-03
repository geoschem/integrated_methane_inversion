import numpy as np
import xarray as xr
from shapely.geometry.polygon import Polygon
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from pyproj import Geod
import cartopy
import cartopy.crs as ccrs
import pickle
import csv
from matplotlib.lines import Line2D


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
        raise ValueError(f"{species} is not recognized. Please add a line to "
                         "mixing_ratio_conv_factor in src/inversion_scripts/utils.py")


def species_molar_mass(species):
    if species == "CH4":
        M = 0.01604  # Molar mass of methane [kg/mol]
    elif species == "CO2":
        M = 0.04401
    else:
        raise ValueError(f"{species} is not recognized. Please add a line to "
                         "species_molar_mass in src/inversion_scripts/utils.py")


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
            dat["CH4"] = tropomi_data["methane_mixing_ratio_bias_corrected"].values[0, :, :]
            dat["qa_value"] = tropomi_data["qa_value"].values[0, :, :]
            dat["longitude"] = tropomi_data["longitude"].values[0, :, :]
            dat["latitude"] = tropomi_data["latitude"].values[0, :, :]

            utc_str = tropomi_data["time_utc"].values[0,:]
            utc_str = np.array([d.replace("Z","") for d in utc_str]).astype("datetime64[ns]")
            dat["time"] = np.repeat(utc_str[:, np.newaxis], dat["CH4"].shape[1], axis=1)

        # Store column averaging kernel, SWIR and NIR surface albedo
        with xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS") as tropomi_data:
            dat["column_AK"] = tropomi_data["column_averaging_kernel"].values[0, :, :, ::-1]
            dat["swir_albedo"] = tropomi_data["surface_albedo_SWIR"].values[0, :, :]
            dat["nir_albedo"] = tropomi_data["surface_albedo_NIR"].values[0, :, :]
            dat["blended_albedo"] = 2.4 * dat["nir_albedo"] - 1.13 * dat["swir_albedo"]

        # Store methane prior profile, dry air subcolumns
        with xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/INPUT_DATA") as tropomi_data:
            dat["profile_apriori"] = tropomi_data["methane_profile_apriori"].values[0, :, :, ::-1]  # mol m-2
            dat["dry_air_subcolumns"] = tropomi_data["dry_air_subcolumns"].values[0, :, :, ::-1]  # mol m-2
            dat["surface_classification"] = (tropomi_data["surface_classification"].values[0, :, :].astype("uint8") & 0x03).astype(int)

            # Also get pressure interval and surface pressure for use below
            pressure_interval = (tropomi_data["pressure_interval"].values[0, :, :] / 100)  # Pa -> hPa
            surface_pressure = (tropomi_data["surface_pressure"].values[0, :, :] / 100)  # Pa -> hPa

        # Store latitude and longitude bounds for pixels
        with xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/GEOLOCATIONS") as tropomi_data:
            dat["longitude_bounds"] = tropomi_data["longitude_bounds"].values[0, :, :, :]
            dat["latitude_bounds"] = tropomi_data["latitude_bounds"].values[0, :, :, :]

        # Store vertical pressure profile
        n1 = dat["CH4"].shape[0]  # length of along-track dimension (scanline) of retrieval field
        n2 = dat["CH4"].shape[1]  # length of across-track dimension (ground_pixel) of retrieval field
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
    assert "BLND" in filename, f"BLND not in filename {filename}, but a blended function is being used"

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
            dat["profile_apriori"] = blended_data["methane_profile_apriori"].values[:, ::-1]
            dat["dry_air_subcolumns"] = blended_data["dry_air_subcolumns"].values[:, ::-1]
            dat["longitude_bounds"] = blended_data["longitude_bounds"].values[:]
            dat["latitude_bounds"] = blended_data["latitude_bounds"].values[:]
            dat["surface_classification"] = (blended_data["surface_classification"].values[:].astype("uint8") & 0x03).astype(int)
            dat["chi_square_SWIR"] = blended_data["chi_square_SWIR"].values[:]

            # Remove "Z" from time so that numpy doesn't throw a warning
            utc_str = blended_data["time_utc"].values[:]
            dat["time"] = np.array([d.replace("Z","") for d in utc_str]).astype("datetime64[ns]")

            # Need to calculate the pressure for the 13 TROPOMI levels (12 layer edges)
            pressure_interval = (blended_data["pressure_interval"].values[:] / 100)  # Pa -> hPa
            surface_pressure = (blended_data["surface_pressure"].values[:] / 100)    # Pa -> hPa
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
    filename,
    satellite_str,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
):
    # Read TROPOMI data
    assert satellite_str in ["BlendedTROPOMI", "TROPOMI", "Other"], "satellite_str  is not one of BlendedTROPOMI, TROPOMI, or Other"
    if satellite_str  == "BlendedTROPOMI":
        satellite = read_blended(filename)
    elif satellite_str  == "TROPOMI":
        satellite = read_tropomi(filename)
    else:
        satellite = ...
        print("Other data source is not currently supported --HON")

    # If empty, skip this file
    if satellite == None:
        print(f"Skipping {filename} due to file processing issue.")
        return satellite

    # Filter the data
    if satellite_str  == "BlendedTROPOMI":
        # Only going to consider blended data within lat/lon/time bounds and wihtout problematic coastal pixels
        sat_ind = filter_blended(satellite, xlim, ylim, gc_startdate, gc_enddate)
    elif satellite_str  == "TROPOMI":
        # Only going to consider TROPOMI data within lat/lon/time bounds and with QA > 0.5
        sat_ind = filter_tropomi(satellite, xlim, ylim, gc_startdate, gc_enddate)
    else:
        sat_ind = ...
        print("Other data source filtering is not currently supported --HON")

    return satellite, sat_ind

def get_posterior_emissions(prior, scale, species):
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
    prior_soil_sink = prior[f"Emis{species}_SoilAbsorb"].copy()
    
    # remove the soil sink from the prior total before applying scale factors
    prior[f"Emis{species}_Total"] = prior[f"Emis{species}_Total"] - prior_soil_sink
    
    # scale the prior emissions for all sectors using the scale factors
    posterior = prior * scale["ScaleFactor"]
    
    # But reset the soil sink to the original value
    posterior[f"Emis{species}_SoilAbsorb"] = prior_soil_sink
    
    # Add the original soil sink back to the total emissions
    posterior[f"Emis{species}_Total"] = posterior[f"Emis{species}_Total"] + prior_soil_sink
    return posterior

def get_strdate(current_time, date_threshold):
    # round observation time to nearest hour
    strdate = current_time.round("60min").strftime("%Y%m%d_%H")
    # Unless it equals the date threshold (hour 00 after the inversion period)
    if strdate == date_threshold:
        strdate = current_time.floor("60min").strftime("%Y%m%d_%H")
 
    return strdate

import glob
import os
import re
import datetime
import numpy as np
from joblib import Parallel, delayed
import xarray as xr
from netCDF4 import Dataset
import sys
import pandas as pd

import yaml

with open("config_boundary_conditions.yml", "r") as f:
    config = yaml.safe_load(f)

sys.path.insert(0, "../../")
from src.inversion_scripts.operators.operator_utilities import nearest_loc
from src.inversion_scripts.operators.TROPOMI_operator import apply_tropomi_operator
from src.inversion_scripts.utils import save_obj, load_obj


def get_TROPOMI_times(filename):
    """
    Function that parses the TROPOMI filenames to get the start and end times.
    Example input (str): S5P_RPRO_L2__CH4____20220725T152751_20220725T170921_24775_03_020400_20230201T100624.nc
    Example output (tuple): (np.datetime64('2022-07-25T15:27:51'), np.datetime64('2022-07-25T17:09:21'))
    """

    file_times = re.search(r"(\d{8}T\d{6})_(\d{8}T\d{6})", filename)
    assert (
        file_times is not None
    ), "check TROPOMI filename - wasn't able to find start and end times in the filename"
    start_TROPOMI_time = np.datetime64(
        datetime.datetime.strptime(file_times.group(1), "%Y%m%dT%H%M%S")
    )
    end_TROPOMI_time = np.datetime64(
        datetime.datetime.strptime(file_times.group(2), "%Y%m%dT%H%M%S")
    )

    return start_TROPOMI_time, end_TROPOMI_time


def apply_tropomi_operator_to_one_tropomi_file(filename):
    """
    Run apply_tropomi_operator from src/inversion_scripts/operators/TROPOMI_operator.py for a single TROPOMI file
    Example input (str): S5P_RPRO_L2__CH4____20220725T152751_20220725T170921_24775_03_020400_20230201T100624.nc
    """

    result = apply_tropomi_operator(
        filename=filename,
        BlendedTROPOMI=blendedTROPOMI,
        n_elements=False,  # Not relevant
        gc_startdate=start_time_of_interest,
        gc_enddate=end_time_of_interest,
        xlim=[-180, 180],
        ylim=[-90, 90],
        gc_cache=os.path.join(config["workDir"], "gc_run", "OutputDir"),
        build_jacobian=False,  # Not relevant
        period_i=False,  # Not relevant
        config=False,
    )  # Not relevant

    return result["obs_GC"], filename


def create_daily_means(satelliteDir, start_time_of_interest, end_time_of_interest):

    # List of all TROPOMI files that interesct our time period of interest
    TROPOMI_files = sorted(
        [
            file
            for file in glob.glob(os.path.join(satelliteDir, "*.nc"))
            if (
                start_time_of_interest
                <= get_TROPOMI_times(file)[0]
                <= end_time_of_interest
            )
            or (
                start_time_of_interest
                <= get_TROPOMI_times(file)[1]
                <= end_time_of_interest
            )
        ]
    )
    print(f"First TROPOMI file -> {TROPOMI_files[0]}")
    print(f"Last TROPOMI file  -> {TROPOMI_files[-1]}")

    # Using as many cores as you have, apply the TROPOMI operator to each file
    obsGC_and_filenames = Parallel(n_jobs=-1)(
        delayed(apply_tropomi_operator_to_one_tropomi_file)(filename)
        for filename in TROPOMI_files
    )

    # Read any of the GEOS-Chem files to get the lat/lon grid
    with xr.open_dataset(
        glob.glob(
            os.path.join(
                config["workDir"], "gc_run", "OutputDir", "GEOSChem.SpeciesConc*.nc4"
            )
        )[0]
    ) as data:
        LON = data["lon"].values
        LAT = data["lat"].values

    # List of all days in our time range of interest
    alldates = np.arange(
        start_time_of_interest, end_time_of_interest, dtype="datetime64[D]"
    )
    alldates = [day.astype(datetime.datetime).strftime("%Y%m%d") for day in alldates]

    # Initialize arrays for regridding
    daily_TROPOMI = np.zeros((len(LON), len(LAT), len(alldates)))
    daily_GC = np.zeros((len(LON), len(LAT), len(alldates)))
    daily_count = np.zeros((len(LON), len(LAT), len(alldates)))

    # Loop thorugh all of the files which now contain TROPOMI and the corresponding GC XCH4
    for obsGC, filename in obsGC_and_filenames:
        NN = obsGC.shape[0]
        if NN == 0:
            continue

        # For each TROPOMI observation, assign it to a GEOS-Chem grid cell
        for iNN in range(NN):

            # Which day are we on (this is not perfect right now because orbits can cross from one day to the next...
            # but it is the best we can do right now without changing apply_tropomi_operator)
            file_times = re.search(r"(\d{8}T\d{6})_(\d{8}T\d{6})", filename)
            assert (
                file_times is not None
            ), "check TROPOMI filename - wasn't able to find start and end times in the filename"
            date = datetime.datetime.strptime(
                file_times.group(1), "%Y%m%dT%H%M%S"
            ).strftime("%Y%m%d")
            time_ind = alldates.index(date)

            c_TROPOMI, c_GC, lon0, lat0 = obsGC[iNN, :4]
            ii = nearest_loc(lon0, LON, tolerance=5)
            jj = nearest_loc(lat0, LAT, tolerance=4)
            daily_TROPOMI[ii, jj, time_ind] += c_TROPOMI
            daily_GC[ii, jj, time_ind] += c_GC
            daily_count[ii, jj, time_ind] += 1

    # Normalize by how many observations got assigned to a grid cell to finish the regridding
    daily_count[daily_count == 0] = np.nan
    daily_TROPOMI = daily_TROPOMI / daily_count
    daily_GC = daily_GC / daily_count

    # Change dimensions
    regrid_TROPOMI = np.einsum(
        "ijl->lji", daily_TROPOMI
    )  # (lon, lat, time) -> (time, lat, lon)
    regrid_GC = np.einsum("ijl->lji", daily_GC)  # (lon, lat, time) -> (time, lat, lon)

    # Make a Dataset with variables of (TROPOMI_CH4, GC_CH4) and dims of (lon, lat, time)
    daily_means = xr.Dataset(
        {
            "TROPOMI_CH4": xr.DataArray(
                data=regrid_TROPOMI,
                dims=["time", "lat", "lon"],
                coords={"time": alldates, "lat": LAT, "lon": LON},
            ),
            "GC_CH4": xr.DataArray(
                data=regrid_GC,
                dims=["time", "lat", "lon"],
                coords={"time": alldates, "lat": LAT, "lon": LON},
            ),
        }
    )

    return daily_means


def calculate_bias(daily_means):

    bias = daily_means["GC_CH4"] - daily_means["TROPOMI_CH4"]

    # Smooth spatially
    bias = bias.rolling(
        lat=5,  # five lat grid boxes (10 degrees)
        lon=5,  # five lon grid boxes (12.5 degrees)
        center=True,  # five boxes includes the one we are cented on
        min_periods=25
        / 2,  # half (13) of the grid cells have a value to not output NaN
    ).mean(skipna=True)

    # Smooth temporally
    bias_15 = bias.rolling(
        time=15,  # average 15 days back in time (including the time we are centered on)
        min_periods=1,  # only one of the time values must have a value to not output NaN
    ).mean(skipna=True)

    bias_30 = bias.rolling(
        time=30,  # average 30 days back in time (including the time we are centered on)
        min_periods=1,  # only one of the time values must have a value to not output NaN
    ).mean(skipna=True)

    bias = bias_15.fillna(bias_30)  # fill in NaN values with the 30 day average

    # Create a dataarray with latitudinal average for each time step
    # We will fill the NaN values in bias with these averages
    nan_value_filler_2d = bias.copy()
    nan_value_filler_2d = (
        nan_value_filler_2d.where(
            nan_value_filler_2d.count("lon") >= 15
        )  # there needs to be 15 grid boxes
        .mean(dim=["lon"], skipna=True)  #  at this lat to define a mean
        .interpolate_na(dim="lat", method="linear")  # fill in "middle" NaN values
        .bfill(dim="lat")  # fill in NaN values towards -90 deg
        .ffill(dim="lat")  # fill in NaN values towards +90 deg
    )

    # Expand to 3 dimensions
    nan_value_filler_3d = bias.copy() * np.nan
    for i in range(len(daily_means["lon"].values)):
        nan_value_filler_3d[:, :, i] = nan_value_filler_2d

    # Use these values to fill NaNs
    bias = bias.fillna(nan_value_filler_3d)

    print(f"Average bias (GC-TROPOMI): {bias.mean().values:.2f} ppb\n")

    # If there are still NaNs (this will happen when TROPOMI data is missing), use 0.0 ppb as the bias but warn the user
    for t in range(len(bias["time"].values)):
        if np.any(np.isnan(bias[t, :, :].values)):
            print(f"WARNING -> using 0.0 ppb as bias for {bias['time'].values[t]}")
            bias[t, :, :] = bias[t, :, :].fillna(0)

    return bias


def write_bias_corrected_files(bias):

    # Get dates and convert the total column bias to mol/mol
    strdate = bias["time"].values
    bias_mol_mol = bias.values * 1e-9

    # Only write BCs for our date range
    files = sorted(
        glob.glob(
            os.path.join(
                config["workDir"],
                "gc_run",
                "OutputDir",
                "GEOSChem.BoundaryConditions*.nc4",
            )
        )
    )
    files = [
        f
        for f in files
        if (
            (
                np.datetime64(re.search(r"(\d{8})_(\d{4}z)", f).group(1))
                >= np.datetime64(config["startDate"])
            )
            & (
                np.datetime64(re.search(r"(\d{8})_(\d{4}z)", f).group(1))
                <= np.datetime64(config["endDate"])
            )
        )
    ]
    assert len(files) == len(
        strdate
    ), "ERROR -> bias dimension is not the same as number of boundary condition files"

    # For each file, remove the total column bias from each level of the GEOS-Chem boundary condition
    for filename in files:
        l = [
            index
            for index, date in enumerate(strdate)
            if date == re.search(r"(\d{8})_(\d{4}z)", filename).group(1)
        ]
        assert (
            len(l) == 1
        ), "ERROR -> there should only be bias per boundary condition file"
        index = l[0]
        bias_for_this_boundary_condition_file = bias_mol_mol[index, :, :]

        with xr.open_dataset(filename) as ds:
            original_data = ds["SpeciesBC_CH4"].values.copy()
            for t in range(original_data.shape[0]):
                for lev in range(original_data.shape[1]):
                    original_data[t, lev, :, :] -= bias_for_this_boundary_condition_file
            ds["SpeciesBC_CH4"].values = original_data
            if blendedTROPOMI:
                print(
                    f"Writing to {os.path.join(config['workDir'], 'blended-boundary-conditions', os.path.basename(filename))}"
                )
                ds.to_netcdf(
                    os.path.join(
                        config["workDir"],
                        "blended-boundary-conditions",
                        os.path.basename(filename),
                    )
                )
            else:
                print(
                    f"Writing to {os.path.join(config['workDir'], 'tropomi-boundary-conditions', os.path.basename(filename))}"
                )
                ds.to_netcdf(
                    os.path.join(
                        config["workDir"],
                        "tropomi-boundary-conditions",
                        os.path.basename(filename),
                    )
                )


if __name__ == "__main__":

    # Arguments from run_boundary_conditions.sh
    blendedTROPOMI = sys.argv[1] == "True"  # use blended data?
    satelliteDir = sys.argv[2]  # where is the satellite data?
    # Start of GC output (+1 day except 1 Apr 2018 because we ran 1 day extra at the start to account for data not being written at t=0)
    start_time_of_interest = np.datetime64(
        datetime.datetime.strptime(sys.argv[3], "%Y%m%d")
    )
    if start_time_of_interest != np.datetime64("2018-04-01T00:00:00"):
        start_time_of_interest += np.timedelta64(1, "D")
    # End of GC output
    end_time_of_interest = np.datetime64(
        datetime.datetime.strptime(sys.argv[4], "%Y%m%d")
    )
    print(f"\nwrite_boundary_conditions.py output for blendedTROPOMI={blendedTROPOMI}")
    print(f"Using files at {satelliteDir}")

    """
    This script works in three parts, utilizing the GEOS-Chem output from run_boundary_conditions.sh
    (1) Make a gridded (2.0 x 2.5 x daily) field of TROPOMI/GEOS-Chem co-locations
        - for every TROPOMI observation, apply the TROPOMI operator to the GEOS-Chem fields
        - average the TROPOMI XCH4 and GEOS-Chem XCH4 to a (2.0 x 2.5) grid for each day
    (2) Make a gridded (2.0 x 2.5 x daily) field of the bias between TROPOMI and GEOS-Chem
        - subtract the TROPOMI and GEOS-Chem grids from part 1 to get a starting point for the bias
        - smooth this field spatially (5 lon grid boxes, 5 lat grid boxes) then temporally (15 days backwards)
            - use a temporal smoothing of 30 days back in time if there are no data in the past 15 days
        - fill NaN values with the latitudinal average at that time
            - for a latitudinal average to be defined, there must be >= 15 grid cells at that latitude
            - when a latitudinal average cannot be found, the closest latitudinal average is used
    (3) Write the boundary conditions
        - using the bias from Part 2, subtract the (GC-TROPOMI) bias from the GC boundary conditions
    """

    daily_means = create_daily_means(
        satelliteDir, start_time_of_interest, end_time_of_interest
    )
    bias = calculate_bias(daily_means)
    write_bias_corrected_files(bias)

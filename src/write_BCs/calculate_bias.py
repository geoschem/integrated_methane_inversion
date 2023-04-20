import numpy as np
import xarray as xr
import os

import yaml
with open("config_write_BCs.yml", "r") as f:
    config = yaml.safe_load(f)

#######################################################
# This script takes the preprocessed daily data 
# from tropomi/GC simulation and calculates the bias 
# between the two using the latitudinal averages and 
# spatial/temporal smoothing. 
# This script is based off of Lu Shen's 
# methods/ R script.
#######################################################

# settings to adjust default smoothing windows
# we use a +/- 15 day time window and 3x3
# smoothing for lat/lon
smoothing_lat_window = 3
smoothing_lon_window = 3
smoothing_time_window = 30

# replace extreme outliers by finding values above the given percentile
# using perc=.01 finds the highest negative bias and perc=.99 finds the 
# highest positive bias
def replace_outliers(data, perc=0.01):
    # calculate percentile 
    threshold = data.quantile(perc, skipna=True)
    sign = abs(threshold)/threshold

    # find outliers and replace them with max among remaining values 
    # .where replace outliers with nan
    mask = data.where(abs(data) <= abs(threshold))
    max_value = abs(mask).max().values

    # fill na values with the next highest values below the threshold
    mask = mask.fillna(max_value*sign)

    # remove the quantile coordinate
    data = mask.drop_vars('quantile')
    
    return data

# smooth out 3d dataarray over time, lat, lon using desired windows
# requires at least 20 of the gridcells to be non-nan values to 
# calculate the mean (specified by min_periods)
def smooth_3D_da(
    da,
    t_window=smoothing_time_window,
    lat_window=smoothing_lat_window,
    lon_window=smoothing_lon_window,
):
    return da.rolling(
        time=t_window,
        lat=lat_window,
        lon=lon_window,
        min_periods=20,
        center=True,
    ).mean(skipna=True)

if __name__ == "__main__":

    # access the preprocessed CH4 data
    filepath = os.path.join(config["workdir"], "step2", "Daily_CH4.nc")
    with xr.open_dataset(filepath) as daily_CH4:
        TROPOMI_lon = daily_CH4["lon"]
        TROPOMI_lat = daily_CH4["lat"]
        GC_CH4 = daily_CH4["GC"]
        TROPOMI_CH4 = daily_CH4["CH4"]
        date = daily_CH4["date"]

    # smooth the background GC and TROPOMI data
    GC_bkgd = smooth_3D_da(GC_CH4)
    TROPOMI_bkgd = smooth_3D_da(TROPOMI_CH4)

    # calculate bias between GC background CH4 and
    # TROPOMI observational background CH4
    bias_4x5 = GC_bkgd - TROPOMI_bkgd

    # build a smoothed dataset to fill in nan values of raw bias
    bias_avg_base = bias_4x5

    # create a dataarray with latitudinal average for each time step
    # we use nearest neighbor interpolation to fill in data gaps
    # for the edges we use backfill and forwardfill which takes 
    # the nearest lat and infills all the way to the edge
    lat_base = (
        bias_avg_base.mean(dim=["lon"], skipna=True)
        .interpolate_na(dim="lat", method="nearest")
        .bfill(dim="lat")
        .ffill(dim="lat")
        .rolling(
            time=smoothing_time_window,
            lat=smoothing_lat_window,
            center=True,
            min_periods=1,
        )
        .mean(skipna=True)
    )

    # expand the latitudinal averages into a 3D dataarray for easier infilling
    lat_base_full = bias_4x5.copy()
    lat_base_full = xr.where(lat_base_full != np.nan, np.nan, np.nan)
    for i in range(len(TROPOMI_lon)):
        lat_base_full[:, :, i] = lat_base

    # infill nan values in bias with latitudinal averages
    # for each corresponding day
    bias_avg_base = bias_avg_base.fillna(lat_base_full)

    # infill nan values of raw bias data with the
    # average background data, replace any extreme outliers with more 
    # reasonable maximums, and smooth the final product over lat and lon
    bias_4x5_new = replace_outliers(bias_4x5.fillna(bias_avg_base), perc=.01)
    bias_4x5_new = bias_4x5_new.rolling(
            lat=smoothing_lat_window,
            lon=smoothing_lon_window,
            min_periods = 1,
            center=True,
        ).mean()

    print(f"Average bias in ppb (GC-tropomi): {bias_4x5_new.mean()}")

    # create dataset and export to netcdf file
    ds = bias_4x5_new.to_dataset(name="Bias")
    ds = ds.assign_coords({"time": date.values})
    ds.to_netcdf(os.path.join(config["workdir"], "step3", "Bias_4x5_dk_2_updated.nc"))
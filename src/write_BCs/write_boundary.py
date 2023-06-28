import numpy as np
import re
import xarray as xr
import glob
import os

import yaml
with open("config_write_BCs.yml", "r") as f:
    config = yaml.safe_load(f)

if __name__ == "__main__":

    outputDir = os.path.join(config["workdir"], "smoothed-boundary-conditions")

    # Read smooth biases calculated in calculate_bias.py
    with xr.open_dataset(os.path.join(config["workdir"], "step3", "Bias_4x5_dk_2_updated.nc")) as file1:
        all_Bias = file1["Bias"].values * 1e-9
        strdate = file1["time"].values

    # Only write BCs for our date range (otherwise, files would be written for the last few days, which are invalid due to a lack of TROPOMI data on the +15 day end)
    os.chdir(os.path.join(config["workdir"],"runGCC1402","OutputDir"))
    files = sorted(glob.glob("GEOSChem.BoundaryConditions*.nc4"))
    files = [f for f in files if ((np.datetime64(f[28:36]) >= np.datetime64(config["startdate"])) &
                                  (np.datetime64(f[28:36]) <= np.datetime64(config["enddate"])))]

    # For each file, remove the bias calculated in calculate_bias.py from each level of the GEOS-Chem boundary conditions
    for filename in files:
        print(filename)
        temp = int(re.split("\.|_", filename)[2])
        ind = [index for index, date in enumerate(strdate) if date == temp]
        if len(ind) == 0:
            print("skipping file")
            continue
        ind = ind[0]
        Bias = all_Bias[ind, :, :]

        with xr.open_dataset(filename) as file2: 
            orig_data = file2["SpeciesBC_CH4"].values.copy()
            for t in range(orig_data.shape[0]):
                for lev in range(47):
                    orig_data[t, lev, :, :] = orig_data[t, lev, :, :] - Bias
            file2["SpeciesBC_CH4"].values = orig_data
            file2.to_netcdf(f"{outputDir}/{filename}")
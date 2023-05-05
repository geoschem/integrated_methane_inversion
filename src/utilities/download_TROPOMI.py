#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Note: requries aws cli to be installed and configured
# Description: Download TROPOMI data from meeo S3 bucket for desired dates.
#              Function can be called from another script or run as a 
#              directly as a script.               
# Example Usage as a script:
#     $ python download_TROPOMI.py 20190101 20190214 TROPOMI_data

import os
import sys
import subprocess
import datetime
import numpy as np


def download_TROPOMI(startdate, enddate, Sat_datadir):
    """
    Download TROPOMI data from AWS for desired dates.

    Arguments
        startdate    [np.datetime64]  : Start date of download range
        enddate      [np.datetime64]  : End date of download range
        Sat_datadir  [str]            : TROPOMI data directory for storing data

    """
    # offline: 07/26/2022 to present
    # s3://meeo-s5p/OFFL/L2__CH4___/YYYY/MM/DD/*.nc
    # reprocessed: 04/30/2018 to 07/25/2022
    # s3://meeo-s5p/RPRO/L2__CH4___/YYYY/MM/DD/*.nc
    # --no-sign-request

    # list of all dates between startdate and enddate excluding enddate
    download_dates = np.arange(startdate, enddate, dtype="datetime64[D]")

    # we create a bash script to download the data using the aws cli
    # This is supposedly faster than downloading data with boto3 via python
    DATA_DOWNLOAD_SCRIPT = "./auto_generated_download_script.sh"
    cmd_prefix = "aws s3 sync --only-show-errors "
    remote_root = "s3://meeo-s5p/"

    with open(DATA_DOWNLOAD_SCRIPT, "w") as f:
        print("#!/bin/bash\n", file=f)

        for date in download_dates:
            # extract year, month, day strings from date
            year, month, day = str(date).split("-")

            # determine download directory based on date
            # data starts on 2018-04-30
            # use reprocessed data for dates before 2022-07-26
            # this ensures use of the v02.04.00/v02.05.00 product
            if date < np.datetime64("2018-04-30"):
                print(
                    f"Skipping TROPOMI data download for {date}, : no data from this date"
                )
                continue
            elif date < np.datetime64("2022-07-26"):
                subdir = "RPRO"
                version = "*_020400_*"
                # handle corner case where a single orbit file would not be downloaded
                if date == np.datetime64("2022-07-25"):
                    download_str = ("aws s3 cp --only-show-errors "
                    + f"{remote_root}OFFL/L2__CH4___/{year}/{month}/{day}/"
                    + "S5P_OFFL_L2__CH4____20220725T235519_20220726T013649_24780_03_020400_20220727T154716.nc"
                    + f" {Sat_datadir}")
                    f.write(download_str)
                    f.write("\n")
            # use offline data for dates after 2022-07-25
            # get v020400 before 2023-03-12 and v020500 after (they are the same)
            else:
                subdir = "OFFL"
                if date < np.datetime64("2023-03-12"):
                    version = "*_020400_*"
                else:
                    version = "*_020500_*"
                    # corner case for one file on 2023-03-12 that is v020400 while the rest are v020500
                    if date == np.datetime64("2023-03-12"):
                        download_str = ("aws s3 cp --only-show-errors "
                        + f"{remote_root}OFFL/L2__CH4___/{year}/{month}/{day}/"
                        + "S5P_OFFL_L2__CH4____20230312T013758_20230312T031928_28030_03_020400_20230313T172452.nc"
                        + f" {Sat_datadir}")
                        f.write(download_str)
                        f.write("\n")
                    
            # build download string
            download_str = cmd_prefix + remote_root
            download_str = (
                f'{download_str}{subdir}/L2__CH4___/{year}/{month}/{day}/ {Sat_datadir} --exclude "*" --include "{version}"'
            )

            # write download string to file
            f.write(download_str)
            f.write("\n")

    # Run the data download script
    # Remove the file afterwards
    os.chmod(DATA_DOWNLOAD_SCRIPT, 0o755)
    status = subprocess.call(DATA_DOWNLOAD_SCRIPT)
    os.remove(DATA_DOWNLOAD_SCRIPT)


if __name__ == "__main__":
    startday = sys.argv[1]
    endday = sys.argv[2]
    Sat_datadir = sys.argv[3]

    # Reformat start and end days for datetime in configuration
    start = f"{startday[0:4]}-{startday[4:6]}-{startday[6:8]} 00:00:00"
    end = f"{endday[0:4]}-{endday[4:6]}-{endday[6:8]} 23:59:59"

    # Convert to datetime64
    GC_startdate = np.datetime64(datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S"))
    GC_enddate = np.datetime64(datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S"))

    download_TROPOMI(GC_startdate, GC_enddate, Sat_datadir)

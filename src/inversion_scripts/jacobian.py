#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import numpy as np
import re
import os
import datetime
from utils import save_obj
from operators.TROPOMI_operator import apply_average_tropomi_operator


if __name__ == "__main__":
    import sys

    startday = sys.argv[1]
    endday = sys.argv[2]
    lonmin = float(sys.argv[3])
    lonmax = float(sys.argv[4])
    latmin = float(sys.argv[5])
    latmax = float(sys.argv[6])
    n_elements = int(sys.argv[7])
    tropomi_cache = sys.argv[8]
    isPost = sys.argv[9]

    # Reformat start and end days for datetime in configuration
    start = f"{startday[0:4]}-{startday[4:6]}-{startday[6:8]} 00:00:00"
    end = f"{endday[0:4]}-{endday[4:6]}-{endday[6:8]} 23:59:59"

    # Configuration
    workdir = "."
    sensi_cache = f"{workdir}/data_sensitivities"
    if isPost.lower() == "false":
        build_jacobian = True
        gc_cache = f"{workdir}/data_geoschem"
        outputdir = f"{workdir}/data_converted"
    else:  # if sampling posterior simulation
        build_jacobian = False
        gc_cache = f"{workdir}/data_geoschem_posterior"
        outputdir = f"{workdir}/data_converted_posterior"
    xlim = [lonmin, lonmax]
    ylim = [latmin, latmax]
    gc_startdate = np.datetime64(datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S"))
    gc_enddate = np.datetime64(
        datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")
        - datetime.timedelta(days=1)
    )
    print("Start:", start)
    print("End:", end)

    # Get TROPOMI data filenames for the desired date range
    allfiles = glob.glob(f"{tropomi_cache}/*.nc")
    sat_files = []
    for index in range(len(allfiles)):
        filename = allfiles[index]
        shortname = re.split("\/", filename)[-1]
        shortname = re.split("\.", shortname)[0]
        strdate = re.split("\.|_+|T", shortname)[4]
        strdate = datetime.datetime.strptime(strdate, "%Y%m%d")
        if (strdate >= gc_startdate) and (strdate <= gc_enddate):
            sat_files.append(filename)
    sat_files.sort()
    print("Found", len(sat_files), "TROPOMI data files.")

    # Map GEOS-Chem to TROPOMI observation space
    # Also return Jacobian matrix if build_jacobian=True
    for filename in sat_files:

        # Check if TROPOMI file has already been processed
        print("========================")
        shortname = re.split("\/", filename)[-1]
        print(shortname)
        date = re.split("\.", shortname)[0]

        # If not yet processed, run apply_average_tropomi_operator()
        if not os.path.isfile(f"{outputdir}/{date}_GCtoTROPOMI.pkl"):
            print("Applying TROPOMI operator...")
            output = apply_average_tropomi_operator(
                filename,
                n_elements,
                gc_startdate,
                gc_enddate,
                xlim,
                ylim,
                gc_cache,
                build_jacobian,
                sensi_cache,
            )
            if output == None:
                continue

        if output["obs_GC"].shape[0] > 0:
            print("Saving .pkl file")
            save_obj(output, f"{outputdir}/{date}_GCtoTROPOMI.pkl")

    print(f"Wrote files to {outputdir}")

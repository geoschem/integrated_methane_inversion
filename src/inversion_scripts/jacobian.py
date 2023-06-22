#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import glob
import numpy as np
import re
import os
import datetime
from utils import save_obj
from operators.TROPOMI_operator import (
    apply_average_tropomi_operator,
    apply_tropomi_operator,
)


def apply_operator(operator, params):
    """
    Run the chosen operator based on selected instrument

    Arguments
        operator [str]    : Data conversion operator to use
        params   [dict]   : parameters to run the given operator
    Returns
        output   [dict]   : Dictionary with:
                            - obs_GC : GEOS-Chem and TROPOMI methane data
                            - TROPOMI methane
                            - GEOS-Chem methane
                            - TROPOMI lat, lon
                            - TROPOMI lat index, lon index
                              If build_jacobian=True, also include:
                                - K      : Jacobian matrix
    """
    if operator == "TROPOMI_average":
        return apply_average_tropomi_operator(
            params["filename"],
            params["blended"],
            params["n_elements"],
            params["gc_startdate"],
            params["gc_enddate"],
            params["xlim"],
            params["ylim"],
            params["gc_cache"],
            params["build_jacobian"],
            params["sensi_cache"],
        )
    elif operator == "TROPOMI":
        return apply_tropomi_operator(
            params["filename"],
            params["blended"],
            params["n_elements"],
            params["gc_startdate"],
            params["gc_enddate"],
            params["xlim"],
            params["ylim"],
            params["gc_cache"],
            params["build_jacobian"],
            params["sensi_cache"],
        )
    else:
        raise ValueError("Error: invalid operator selected.")


if __name__ == "__main__":

    startday = sys.argv[1]
    endday = sys.argv[2]
    lonmin = float(sys.argv[3])
    lonmax = float(sys.argv[4])
    latmin = float(sys.argv[5])
    latmax = float(sys.argv[6])
    n_elements = int(sys.argv[7])
    tropomi_cache = sys.argv[8]
    blended = sys.argv[9] == "true"
    isPost = sys.argv[10]

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
        vizdir = f"{workdir}/data_visualization"
    else:  # if sampling posterior simulation
        build_jacobian = False
        gc_cache = f"{workdir}/data_geoschem_posterior"
        outputdir = f"{workdir}/data_converted_posterior"
        vizdir = f"{workdir}/data_visualization_posterior"
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

            output = apply_operator(
                "TROPOMI_average",
                {
                    "filename": filename,
                    "blended": blended,
                    "n_elements": n_elements,
                    "gc_startdate": gc_startdate,
                    "gc_enddate": gc_enddate,
                    "xlim": xlim,
                    "ylim": ylim,
                    "gc_cache": gc_cache,
                    "build_jacobian": build_jacobian,
                    "sensi_cache": sensi_cache,
                },
            )

            # we also save out the unaveraged tropomi operator for visualization purposes
            viz_output = apply_operator(
                "TROPOMI",
                {
                    "filename": filename,
                    "blended": blended,
                    "n_elements": n_elements,
                    "gc_startdate": gc_startdate,
                    "gc_enddate": gc_enddate,
                    "xlim": xlim,
                    "ylim": ylim,
                    "gc_cache": gc_cache,
                    "build_jacobian": build_jacobian,
                    "sensi_cache": sensi_cache,
                },
            )

            if output == None:
                continue

        if output["obs_GC"].shape[0] > 0:
            print("Saving .pkl file")
            save_obj(output, f"{outputdir}/{date}_GCtoTROPOMI.pkl")
            save_obj(viz_output, f"{vizdir}/{date}_GCtoTROPOMI.pkl")

    print(f"Wrote files to {outputdir}")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import glob
import numpy as np
import re
import os
import datetime
from src.inversion_scripts.utils import save_obj
from src.inversion_scripts.operators.satellite_operator import (
    apply_average_satellite_operator,
    apply_satellite_operator,
)
from joblib import Parallel, delayed


def apply_operator(operator, params):
    """
    Run the chosen operator based on selected instrument

    Arguments
        operator [str]    : Data conversion operator to use
        params   [dict]   : parameters to run the given operator
    Returns
        output   [dict]   : Dictionary with:
                            - obs_GC : GEOS-Chem and satellite column data
                            - satellite columns
                            - GEOS-Chem columns
                            - satellite lat, lon
                            - satellite lat index, lon index
                              If build_jacobian=True, also include:
                                - K      : Jacobian matrix
    """
    if operator == "satellite_average":
        return apply_average_satellite_operator(
            params["filename"],
            params["species"],
            params["satellite_product"],
            params["n_elements"],
            params["gc_startdate"],
            params["gc_enddate"],
            params["xlim"],
            params["ylim"],
            params["gc_cache"],
            params["build_jacobian"],
            params["sensi_cache"],
        )
    elif operator == "satellite":
        return apply_satellite_operator(
            params["filename"],
            params["species"],
            params["satellite_product"],
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
    species = sys.argv[8]
    satellite_cache = sys.argv[9]
    satellite_product = sys.argv[10]
    isPost = sys.argv[11]
    build_jacobian = sys.argv[12]

    # Reformat start and end days for datetime in configuration
    start = f"{startday[0:4]}-{startday[4:6]}-{startday[6:8]} 00:00:00"
    end = f"{endday[0:4]}-{endday[4:6]}-{endday[6:8]} 23:59:59"

    # Configuration
    workdir = "."
    sensi_cache = f"{workdir}/data_sensitivities"
    if build_jacobian.lower() == "true":
        build_jacobian = True
    else:
        build_jacobian = False
    if isPost.lower() == "false":  # if sampling prior simulation
        gc_cache = f"{workdir}/data_geoschem"
        outputdir = f"{workdir}/data_converted"
        vizdir = f"{workdir}/data_visualization"
    else:  # if sampling posterior simulation
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

    # Get satellite data filenames for the desired date range
    allfiles = glob.glob(f"{satellite_cache}/*.nc")
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
    print("Found", len(sat_files), "satellite data files.")

    # Map GEOS-Chem to satellite observation space
    # Also return Jacobian matrix if build_jacobian=True
    def process(filename):

        # Check if satellite file has already been processed
        print("========================")
        shortname = re.split("\/", filename)[-1]
        print(shortname)
        date = re.split("\.", shortname)[0]

        # If not yet processed, run apply_average_satellite_operator()
        if not os.path.isfile(f"{outputdir}/{date}_GCtoSatellite.pkl"):
            print("Applying satellite operator...")

            output = apply_operator(
                "satellite_average",
                {
                    "filename": filename,
                    "species" : species,
                    "satellite_product": satellite_product,
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

            # we also save out the unaveraged satellite operator for visualization purposes
            viz_output = apply_operator(
                "satellite",
                {
                    "filename": filename,
                    "species" : species,
                    "satellite_product": satellite_product,
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
                return 0
        else:
            return 0

        if output["obs_GC"].shape[0] > 0:
            print("Saving .pkl file")
            save_obj(output, f"{outputdir}/{date}_GCtoSatellite.pkl")
            save_obj(viz_output, f"{vizdir}/{date}_GCtoSatellite.pkl")
        return 0

    results = Parallel(n_jobs=-1)(delayed(process)(filename) for filename in sat_files)
    print(f"Wrote files to {outputdir}")

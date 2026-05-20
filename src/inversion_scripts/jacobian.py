#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import glob
import numpy as np
import re
import os
import datetime
import gc
from src.inversion_scripts.utils import save_obj
from src.utilities.config_utils import load_config
from src.inversion_scripts.operators.satellite_operator import (
    apply_operator,
    get_virtual_satellite,
    get_virtual_satellite_pert_and_base,
    superobservations,
)
from src.inversion_scripts.utils import (
    check_is_OH_element,
    check_is_BC_element,
)
from joblib import Parallel, delayed


def construct_jacobian(
    obs_mapped_to_gc: np.ndarray,
    n_elements,
    gc_cache,
    period_i,
    config,
):
    """
    Construct the Jacobian matrix from the perturbation runs

    Arguments
        obs_mapped_to_gc      : satellite observations mapped to GC gridcells
        n_elements           : number of state vector elements
        gc_cache             : path to GEOS-Chem output data
        period_i             : kalman filter period
        config               : inversion configuration dictionary
    
    Returns
         K      : Jacobian matrix
    """
    # Initialize Jacobian K
    n_gridcells = len(obs_mapped_to_gc)
    jacobian_K = np.empty([n_gridcells, n_elements], dtype=np.float32)
    jacobian_K.fill(np.nan)

    pertf = os.path.expandvars(
        f'{config["OutputPath"]}/{config["RunName"]}/'
        f"archive_perturbation_sfs/pert_sf_{period_i}.npz"
    )

    emis_perturbations_dict = np.load(pertf, mmap_mode='r')
    emis_perturbations = emis_perturbations_dict["effective_pert_sf"]

    # Calculate sensitivities and save in K matrix
    # determine which elements are for emis,
    # BCs, and OH
    oh_indices = []
    bc_indices = []
    emis_indices = []

    for e in range(n_elements):
        i_elem = e + 1
        # booleans for whether this element is a
        # BC element or OH element
        is_OH_element = check_is_OH_element(
            i_elem, n_elements, config["OptimizeOH"], config["isRegional"]
        )

        is_BC_element = check_is_BC_element(
            i_elem,
            n_elements,
            config["OptimizeOH"],
            config["OptimizeBCs"],
            is_OH_element,
            config["isRegional"],
        )

        if is_OH_element:
            oh_indices.append(e)
        elif is_BC_element:
            bc_indices.append(e)
        else:
            emis_indices.append(e)
    
    all_strdate = [gridcell["time"] for gridcell in obs_mapped_to_gc]
    all_strdate = list(set(all_strdate))

    for strdate in all_strdate:
        gridcell_dict = obs_mapped_to_gc[obs_mapped_to_gc["time"] == strdate]
        sel_idx = np.where(obs_mapped_to_gc["time"] == strdate)[0]
        virtual_satellite = get_virtual_satellite(
            strdate, gc_cache, gridcell_dict, n_elements, config
        )
        virtual_satellite_pert, virtual_satellite_base = get_virtual_satellite_pert_and_base(
            strdate, gc_cache, gridcell_dict, n_elements, config
        )

        pert_jacobian_xspecies = virtual_satellite_pert # (n_superobs, n_element)
        emis_base_xspecies = virtual_satellite_base # emis_base and BC_base is "RunName_0001" and "SpeciesConcVV_species"
        oh_base_xspecies = virtual_satellite # OH base is "RunName_0000"

        # get perturbations and calculate sensitivities
        perturbations = np.ones((len(gridcell_dict), n_elements), dtype=np.float32)

        # fill pert base array with values
        # array contains 1 entry for each state vector element
        # fill array with nans
        base_xspecies = np.full((len(gridcell_dict), n_elements), np.nan, dtype=np.float32)
        # fill emission elements with the base value
        base_xspecies[:,emis_indices] = np.repeat(emis_base_xspecies,
                                                np.asarray(emis_indices).size, axis=1)

        # emissions perturbations
        perturbations[:,emis_indices] = np.repeat(emis_perturbations[None,:],
                                                    len(gridcell_dict), axis=0)

        # OH perturbations
        if config["OptimizeOH"]:
            # fill OH elements with the OH base value
            base_xspecies[:,oh_indices] = np.repeat(oh_base_xspecies[:,None],
                                                np.asarray(oh_indices).size, axis=1)
            # update perturbations array to include OH perturbations
            perturbations[:,oh_indices] = float(config["PerturbValueOH"]) - 1.0

        # BC perturbations
        if config["OptimizeBCs"]:
            # fill BC elements with the base value, which is same as emis value
            base_xspecies[:,bc_indices] = np.repeat(emis_base_xspecies,
                                                np.asarray(bc_indices).size, axis=1)

            # compute BC perturbation for jacobian construction
            perturbations[:,bc_indices] = config["PerturbValueBCs"]

        # calculate sensitivities
        jacobian_K[sel_idx,:] = ((pert_jacobian_xspecies - base_xspecies) / perturbations).astype(np.float32)

    return jacobian_K

if __name__ == "__main__":

    workdir = sys.argv[1]
    config = load_config(sys.argv[2])
    startday = sys.argv[3]
    endday = sys.argv[4]
    lonmin = float(sys.argv[5])
    lonmax = float(sys.argv[6])
    latmin = float(sys.argv[7])
    latmax = float(sys.argv[8])
    n_elements = int(sys.argv[9])
    species = sys.argv[10]
    satellite_cache = sys.argv[11]
    satellite_product = sys.argv[12]
    use_water_obs = sys.argv[13]
    isPost = sys.argv[14]
    period_i = int(sys.argv[15])
    build_jacobian = sys.argv[16]
    viz_prior = sys.argv[17]

    # Reformat start and end days for datetime in configuration
    start = f"{startday[0:4]}-{startday[4:6]}-{startday[6:8]} 00:00:00"
    end = f"{endday[0:4]}-{endday[4:6]}-{endday[6:8]} 23:59:59"

    # Configuration
    if build_jacobian.lower() == "true":
        build_jacobian = True
    else:
        build_jacobian = False
    if isPost.lower() == "false":  # if sampling prior simulation
        gc_cache = f"{workdir}/data_geoschem"
        outputdir = f"{workdir}/data_converted"
        vizdir = f"{workdir}/data_visualization"

        # for lognormal, we also sample the prior simulation in a
        # separate call to jacobian.py solely for visualization purposes
        if viz_prior.lower() == "true":
            gc_cache = f"{gc_cache}_prior"
            outputdir = f"{outputdir}_prior"
            vizdir = f"{vizdir}_prior"

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
    print("Start:", gc_startdate)
    print("End:", gc_enddate)

    # Get satellite data filenames for the desired date range
    allfiles = glob.glob(f"{satellite_cache}/*.nc")
    sat_files = []
    for index in range(len(allfiles)):
        filename = allfiles[index]
        shortname = re.split(r"\/", filename)[-1]
        shortname = re.split(r"\.", shortname)[0]
        strdate = re.split(r"\.|_+|T", shortname)[4]
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
        shortname = re.split(r"\/", filename)[-1]
        print(shortname)
        date = re.split(r"\.", shortname)[0]

        # If not yet processed, run apply_average_satellite_operator()
        if not os.path.isfile(f"{outputdir}/{date}_GCtoSatellite.pkl"):
            print("Applying satellite operator...")

            # Compute super-observations for this satellite file
            obs_mapped_to_gc = superobservations(
                filename,
                species,
                satellite_product,
                n_elements,
                gc_startdate,
                gc_enddate,
                xlim,
                ylim,
                gc_cache,
                period_i,
                config,
                use_water_obs=use_water_obs,
            )
            if obs_mapped_to_gc is None:
                return 0

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
                    "period_i": period_i,
                    "use_water_obs": use_water_obs,
                },
                obs_mapped_to_gc,
                config,
            )
            if build_jacobian:
                jacobian = construct_jacobian(
                    obs_mapped_to_gc,
                    n_elements,
                    gc_cache,
                    period_i,
                    config,
                )
                output['K'] = jacobian

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
                    "period_i": period_i,
                    "use_water_obs": use_water_obs,
                },
                obs_mapped_to_gc,
                config,
            )

            if output is None:
                return 0
        else:
            return 0

        if output["obs_GC"].shape[0] > 0:
            print("Saving .pkl file")
            save_obj(output, f"{outputdir}/{date}_GCtoSatellite.pkl")
            save_obj(viz_output, f"{vizdir}/{date}_GCtoSatellite.pkl")

        #Clean up to reduce memory use
        del output, viz_output
        gc.collect()

        return 0

    results = Parallel(n_jobs=-1)(delayed(process)(filename) for filename in sat_files)
    print(f"Wrote files to {outputdir}")

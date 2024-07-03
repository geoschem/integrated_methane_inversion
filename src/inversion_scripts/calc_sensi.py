import os
import glob
import math
import datetime
import numpy as np
import xarray as xr
from joblib import Parallel, delayed


def zero_pad_num(n):
    nstr = str(n)
    if len(nstr) == 1:
        nstr = "000" + nstr
    if len(nstr) == 2:
        nstr = "00" + nstr
    if len(nstr) == 3:
        nstr = "0" + nstr
    return nstr


def check_is_OH_element(sv_elem, nelements, opt_OH):
    """
    Determine if the current state vector element is the OH element
    """
    return opt_OH and (sv_elem == nelements)


def check_is_BC_element(sv_elem, nelements, opt_OH, opt_BC, is_OH_element):
    """
    Determine if the current state vector element is a boundary condition element
    """
    return (
        not is_OH_element
        and opt_BC
        and (
            (opt_OH and (sv_elem > (nelements - 5)))
            or ((not opt_OH) and (sv_elem > (nelements - 4)))
        )
    )


def test_GC_output_for_BC_perturbations(e, nelements, sensitivities, opt_OH):
    """
    Ensures that CH4 boundary condition perturbation in GEOS-Chem is working as intended
    sensitivities = (pert-base)/perturbationBC which should equal 1e-9 inside the perturbed borders
    example: the north boundary is perturbed by 10 ppb
             pert-base=10e-9 mol/mol in the 3 grid cells that have been perturbed
             perturbationBC=10 ppb
             sensitivities = (pert-base)/perturbationBC = 1e-9
    """
    # if optimizing OH, adjust which elements
    # we are checking by 1
    if opt_OH:
        e_pad = 1
    else:
        e_pad = 0

    if e == (nelements - e_pad - 4):  # North boundary
        check = np.mean(sensitivities[:, -3:, 3:-3])
    elif e == (nelements - e_pad - 3):  # South boundary
        check = np.mean(sensitivities[:, 0:3, 3:-3])
    elif e == (nelements - e_pad - 2):  # East boundary
        check = np.mean(sensitivities[:, :, -3:])
    elif e == (nelements - e_pad - 1):  # West boundary
        check = np.mean(sensitivities[:, :, 0:3])
    else:
        msg = (
            "GC CH4 perturb not working... "
            "Check OH and BC optimization options. "
            "Ensure perturbations are >0 if optimizing "
            "BCs and/or OH."
        )
        raise RuntimeError(msg)
    assert (
        abs(check - 1e-9) < 1e-11
    ), f"GC CH4 perturb not working... perturbation is off by {abs(check - 1e-9)} mol/mol/ppb"


def calc_sensi(
    nelements,
    ntracers,
    perturbation,
    startday,
    endday,
    run_dirs_pth,
    run_name,
    sensi_save_pth,
    perturbationBC,
    perturbationOH,
):
    """
    Loops over output data from GEOS-Chem perturbation simulations to compute sensitivities
    for the Jacobian matrix.

    Arguments
        nelements      [int]   : Number of state vector elements
        ntracers       [int]   : Number of Jacobian tracers in simulations
        perturbation   [str]   : Path to perturbation array file
        startday       [str]   : First day of inversion period; formatted YYYYMMDD
        endday         [str]   : Last day of inversion period; formatted YYYYMMDD
        run_dirs_pth   [str]   : Path to directory containing GC Jacobian run directories
        run_name       [str]   : Simulation run name; e.g. 'CH4_Jacobian'
        sensi_save_pth [str]   : Path to save the sensitivity data
        perturbationBC [float] : Size of BC perturbation in ppb (eg. 10.0)
        perturbationOH [float] : Size of OH perturbation in ppb (eg. 1.5)

    Resulting 'sensi' files look like:

        <xarray.Dataset>
        Dimensions:  (grid: 1207, lat: 105, lev: 47, lon: 87)
        Coordinates:
        * lon      (lon) float64 -107.8 -107.5 -107.2 -106.9 ... -81.56 -81.25 -80.94
        * lat      (lat) float64 10.0 10.25 10.5 10.75 11.0 ... 35.25 35.5 35.75 36.0
        * lev      (lev) int32 1 2 3 4 5 6 7 8 9 10 ... 38 39 40 41 42 43 44 45 46 47
        * grid     (grid) int32 1 2 3 4 5 6 7 8 ... 1201 1202 1203 1204 1205 1206 1207
        Data variables:
            sensi    (grid, lev, lat, lon) float32 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0

    Pseudocode summary:

        for each day:
            load the base run SpeciesConc file
            nlon = count the number of longitudes
            nlat = count the number of latitudes
            nlev = count the number of vertical levels
            for each hour:
                base = extract the base run data for the hour
                sensi = np.empty((nelements, nlev, nlat, nlon))
                sensi.fill(np.nan)
                for each state vector element:
                    load the SpeciesConc .nc file for the element and day
                    pert = extract the data for the hour
                    sens = pert - base
                    sensi[element,:,:,:] = sens
                save sensi as netcdf with appropriate coordinate variables
    """
    # subtract by 1 because here we assume .5 is a +50% perturbation
    # perturbation = perturbation - 1

    # Make date range
    days = []
    dt = datetime.datetime.strptime(startday, "%Y%m%d")
    dt_max = datetime.datetime.strptime(endday, "%Y%m%d")
    while dt < dt_max:
        dt_str = str(dt)[0:10].replace("-", "")
        days.append(dt_str)
        delta = datetime.timedelta(days=1)
        dt += delta
    # count number of perturbation simulations in run_dirs_pth
    # we subtract 1 because of the prior simulation
    pattern = os.path.join(run_dirs_pth, "*_[0-9][0-9][0-9][0-9]")
    nruns = len([d for d in glob.glob(pattern) if os.path.isdir(d)]) - 1

    # Loop over model data to get sensitivities
    hours = range(24)
    elements = range(nelements)

    # whether we have OH and BC perturbations
    opt_OH = True if (perturbationOH > 0.0) else False
    opt_BC = True if (perturbationBC > 0.0) else False

    # Dictionary that stores mapping of state vector elements to
    # perturbation simulation numbers
    pert_simulations_dict = {}
    for e in elements:
        # State vector elements are numbered 1..nelements
        sv_elem = e + 1
        
        is_OH_element = check_is_OH_element(sv_elem, nelements, opt_OH)
        # Determine which run directory to look in
        if is_OH_element:
            run_number = nruns
        elif check_is_BC_element(sv_elem, nelements, opt_OH, opt_BC, is_OH_element):
            num_back = nelements % sv_elem
            run_number = nruns - num_back
        else:
            run_number = math.ceil(sv_elem / ntracers)

        run_num = str(run_number).zfill(4)

        # add the element to the dictionary for the relevant simulation number
        if run_num not in pert_simulations_dict:
            pert_simulations_dict[run_num] = [sv_elem]
        else:
            pert_simulations_dict[run_num].append(sv_elem)
    
    # For each day
    for d in days:
        # Load the base run SpeciesConc file
        prior = xr.load_dataset(
            f"{run_dirs_pth}/{run_name}_0000/OutputDir/GEOSChem.SpeciesConc.{d}_0000z.nc4"
        )
        bc_base = xr.load_dataset(
            f"{run_dirs_pth}/{run_name}_0001/OutputDir/GEOSChem.SpeciesConc.{d}_0000z.nc4"
        )
        # Count nlat, nlon, nlev
        nlon = len(prior["lon"])  # 52
        nlat = len(prior["lat"])  # 61
        nlev = len(prior["lev"])  # 47

        # For each hour
        def process(h):
            # Initialize sensitivities array
            sensi = np.empty((nelements, nlev, nlat, nlon))
            sensi.fill(np.nan)

            # Loop through each perturbation simulation
            for run_num in pert_simulations_dict.keys():
                # Load the SpeciesConc file for the current pert simulation and day
                pert_data = xr.open_dataset(
                    f"{run_dirs_pth}/{run_name}_{run_num}/OutputDir/GEOSChem.SpeciesConc.{d}_0000z.nc4"
                )

                # For each state vector element in the loaded simulation
                for sv_elem in pert_simulations_dict[run_num]:
                    e_idx = sv_elem - 1
                    elem = str(sv_elem).zfill(4)

                    # booleans for whether this element is a BC element or OH element
                    is_OH_element = check_is_OH_element(sv_elem, nelements, opt_OH)

                    is_BC_element = check_is_BC_element(
                        sv_elem, nelements, opt_OH, opt_BC, is_OH_element
                    )

                    # Get the data for the current hour
                    key = (
                        "SpeciesConcVV_CH4"
                        if is_OH_element or is_BC_element
                        else f"SpeciesConcVV_CH4_{elem}"
                    )
                    pert = pert_data[key][h, :, :, :]

                    # Compute and store the sensitivities
                    if is_OH_element:
                        # Get the base run data for the hour
                        oh_base = prior["SpeciesConcVV_CH4"][h, :, :, :]
                        # calculate OH sensitivities
                        sensitivities = (pert.values - oh_base.values) / perturbationOH
                    elif is_BC_element:
                        pert_base = bc_base[f"SpeciesConcVV_CH4"][h, :, :, :]
                        # calculate BC sensitivities
                        sensitivities = (
                            pert.values - pert_base.values
                        ) / perturbationBC

                        # because we take the first hour on the first day from spinup
                        if h != 0:
                            test_GC_output_for_BC_perturbations(
                                e_idx, nelements, sensitivities, opt_OH
                            )
                    else:
                        pert_base = pert_data["SpeciesConcVV_CH4"][h, :, :, :]
                        # Calculate emission perturbations
                        sensitivities = (pert.values - pert_base.values) / perturbation[
                            e_idx
                        ]
                    sensi[e_idx, :, :, :] = sensitivities

                # close pert data file to reduce memory load
                pert_data.close()
            # Save sensi as netcdf with appropriate coordinate variables
            sensi = xr.DataArray(
                sensi,
                coords=(
                    np.arange(1, nelements + 1),
                    np.arange(1, nlev + 1),
                    prior.lat,
                    prior.lon,
                ),
                dims=["element", "lev", "lat", "lon"],
                name="Sensitivities",
            )
            sensi = sensi.to_dataset()
            sensi.to_netcdf(
                f"{sensi_save_pth}/sensi_{d}_{str(h).zfill(2)}.nc",
                encoding={
                    v: {"zlib": True, "complevel": 9} for v in sensi.data_vars
                },
            )

        results = Parallel(n_jobs=-1)(delayed(process)(hour) for hour in hours)
    print(f"Saved GEOS-Chem sensitivity files to {sensi_save_pth}")


if __name__ == "__main__":
    import sys

    nelements = int(sys.argv[1])
    ntracers = int(sys.argv[2])
    perturbation = sys.argv[3]
    startday = sys.argv[4]
    endday = sys.argv[5]
    run_dirs_pth = sys.argv[6]
    run_name = sys.argv[7]
    sensi_save_pth = sys.argv[8]
    perturbationBC = float(sys.argv[9])
    perturbationOH = float(sys.argv[10])

    perturbation = np.load(perturbation)
    calc_sensi(
        nelements,
        ntracers,
        perturbation,
        startday,
        endday,
        run_dirs_pth,
        run_name,
        sensi_save_pth,
        perturbationBC,
        perturbationOH,
    )

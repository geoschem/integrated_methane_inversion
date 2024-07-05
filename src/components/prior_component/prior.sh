#!/bin/bash

# Functions available in this file include:
#   - run_prior
#   - exclude_soil_sink

# Description: Run a HEMCO standalone simulation to generate prior emissions
#
# Usage:
#   run_prior
run_prior() {
    prior_start=$(date +%s)
    if [[ -d ${RunDirs}/prior_run ]]; then
        printf "\nERROR: ${PriorDir} already exists. Please remove or set 'DoPriorEmis: false' in config.yml.\n"
        exit 9999
    fi

    printf "\n=== GENERATING PRIOR EMISSIONS ===\n"

    cd ${GCClassicPath}/src/HEMCO/run

    # Commands to feed to createRunDir.sh
    # HEMCO standalone run directories are created for the global domain by default
    if [[ "$Met" == "MERRA2" || "$Met" == "MERRA-2" || "$Met" == "merra2" ]]; then
        metNum="1"
    elif [[ "$Met" == "GEOSFP" || "$Met" == "GEOS-FP" || "$Met" == "geosfp" ]]; then
        metNum="2"
    else
        printf "\nERROR: Meteorology field ${Met} is not supported by the IMI. "
        printf "\n Options are GEOSFP or MERRA2.\n"
        exit 1
    fi
    if [ "$Res" = "4.0x5.0" ]; then
        resnum="1"
    elif [ "$Res" == "2.0x2.5" ]; then
        resnum="2"
    elif [ "$Res" == "0.5x0.625" ]; then
        resnum="3"
    elif [ "$Res" == "0.25x0.3125" ]; then
        resnum="4"
    else
        printf "\nERROR: Grid resolution ${Res} is not supported by the IMI. "
        printf "\n Options are 0.25x0.3125, 0.5x0.625, 2.0x2.5, or 4.0x5.0.\n"
        exit 1
    fi
    HEMCOconfig=${RunTemplate}/HEMCO_Config.rc

    priorDir="prior_run"
    cmd="${metNum}\n${resnum}\n${HEMCOconfig}\n${RunDirs}\n${priorDir}\nn\n"

    # Create HEMCO standalone directory
    printf ${cmd} | ./createRunDir.sh >>createHemcoDir.log 2>&1
    rm -f createHemcoDir.log
    printf "\nCreated ${RunDirs}/prior_run\n"

    cd ${RunDirs}/prior_run

    # Modify HEMCO files based on settings in config.yml
    sed -i -e "/DiagnFreq:           00000100 000000/d" \
        -e "/Negative values:     0/d" HEMCO_sa_Config.rc
    sed -i -e "s/METEOROLOGY            :       true/METEOROLOGY            :       false/g" \
        -e "s|DiagnFreq:                   Monthly|DiagnFreq:                   Daily|g" HEMCO_Config.rc
    sed -i -e "/#SBATCH -c 8/d" runHEMCO.sh
    sed -i -e "/#SBATCH -t 0-12:00/d" runHEMCO.sh
    sed -i -e "/#SBATCH -p huce_intel/d" runHEMCO.sh
    sed -i -e "/#SBATCH --mem=15000/d" runHEMCO.sh
    sed -i '/.*hemco_standalone.*/a\
retVal=$?\
if [ $retVal -ne 0 ]; then\
    rm -f .error_status_file.txt\
    echo "Error Status: $retVal" > .error_status_file.txt\
    echo "prior run exited with error code: $retVal"\
    exit $retVal\
fi' runHEMCO.sh

    sed -i -e "/retVal=\$?/a retVal=\$?" runHEMCO.sh
    sa_XMIN=$(python ${InversionPath}/src/components/prior_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc XMIN)
    sa_XMAX=$(python ${InversionPath}/src/components/prior_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc XMAX)
    sa_YMIN=$(python ${InversionPath}/src/components/prior_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc YMIN)
    sa_YMAX=$(python ${InversionPath}/src/components/prior_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc YMAX)
    sa_NX=$(python ${InversionPath}/src/components/prior_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc NX)
    sa_NY=$(python ${InversionPath}/src/components/prior_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc NY)
    sa_YEDGE=$(python ${InversionPath}/src/components/prior_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc YEDGE)
    sa_YMID=$(python ${InversionPath}/src/components/prior_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc YMID)
    sed -i -e "s/XMIN.*/XMIN: ${sa_XMIN}/" \
        -e "s/XMAX.*/XMAX: ${sa_XMAX}/" \
        -e "s/YMIN.*/YMIN: ${sa_YMIN}/" \
        -e "s/YMAX.*/YMAX: ${sa_YMAX}/" \
        -e "s/NX.*/NX: ${sa_NX}/" \
        -e "s/NY.*/NY: ${sa_NY}/" \
        -e "s/YEDGE.*/YEDGE: ${sa_YEDGE}/" \
        -e "s/YMID.*/YMID: ${sa_YMID}/" HEMCO_sa_Grid.${gridFile}.rc
    mv runHEMCO.sh ${RunName}_Prior.run

    # Compile HEMCO and store executable in template run directory
    printf "\nCompiling HEMCO...\n"
    cd build
    cmake ${InversionPath}/GCClassic/src/HEMCO >>build_hemco.log 2>&1
    cmake . -DRUNDIR=.. >>build_hemco.log 2>&1
    make -j install >>build_hemco.log 2>&1
    cd ..
    if [[ -f hemco_standalone ]]; then
        mkdir HEMCO_build_info
        mv build/CMakeCache.txt HEMCO_build_info
        rm -rf build
    else
        printf "\nHEMCO build failed! \n\nSee ${RunTemplate}/build/build_hemco.log for details\n"
        exit 999
    fi
    printf "\nDone compiling HEMCO \n\nSee ${RunDirs}/prior_run/HEMCO_build_info for details\n\n"

    printf "\nSubmitting prior emissions hemco simulation\n\n"

    run_hemco_sa $StartDate $EndDate

    printf "\nDone prior emissions hemco simulation\n\n"

    printf "\n=== DONE GENERATING PRIOR EMISSIONS ===\n"
    prior_end=$(date +%s)
}

# Description: Run HEMCO standalone simulation
#     to generate prior emissions for given dates
# Usage: run_hemco_sa <hemco_start> <hemco_end>
run_hemco_sa() {
    hemco_start=$1
    hemco_end=$2
    set -e

    pushd ${RunDirs}/prior_run
    # replace start and end times in HEMCO_sa_Time.rc
    sed -i -e "s|START.*|START: ${hemco_start:0:4}-${hemco_start:4:2}-${hemco_start:6:2} 00:00:00|g" \
        -e "s|END.*|END: ${hemco_end:0:4}-${hemco_end:4:2}-${hemco_end:6:2} 00:00:00|g" HEMCO_sa_Time.rc

    # Submit job to job scheduler
    sbatch --mem $RequestedMemory \
        -c $RequestedCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -W ${RunName}_Prior.run
    wait

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    # Remove soil absorption uptake from total emissions
    pushd OutputDir
    for file in HEMCO_sa_diagnostics*.nc; do
        exclude_soil_sink $file $file
    done
    popd
    popd
    set +e
}

# Description: Create new netCDF file with EmisCH4_Total_ExclSoilAbs
#   which removes the soil sink from the total emissions
# Usage:
#   exclude_soil_sink <src-file> <target-file>
exclude_soil_sink() {
    python -c "import sys; import xarray; import numpy as np; \
    emis = xarray.load_dataset(sys.argv[1]); \
    emis['EmisCH4_Total_ExclSoilAbs'] = emis['EmisCH4_Total'] - emis['EmisCH4_SoilAbsorb']; \
    emis['EmisCH4_Total_ExclSoilAbs'].attrs = emis['EmisCH4_Total'].attrs; \
    emis.to_netcdf(sys.argv[2])" $1 $2
}

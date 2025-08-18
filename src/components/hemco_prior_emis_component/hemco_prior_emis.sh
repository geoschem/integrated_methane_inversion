#!/bin/bash

# Functions available in this file include:
#   - run_hemco_prior_emis
#   - exclude_soil_sink
#   - for GCHP: set_prior_gchp; run_prior_gchp

# Description: Run a HEMCO standalone simulation to generate prior emissions
#
# Usage:
#   run_hemco_prior_emis
run_hemco_prior_emis() {
    hemco_prior_emis_start=$(date +%s)
    
    HEMCOdir="hemco_prior_emis"
    if [[ -d ${RunDirs}/${HEMCOdir} ]]; then
        printf "\nERROR: ${RunDirs}/${HEMCOdir} already exists. Please remove or set 'DoHemcoPriorEmis: false' in config.yml.\n"
        exit 9999
    fi

    ### Perform dry run if requested
    if "$HemcoPriorEmisDryRun"; then
        pushd ${RunDirs}/template_run
        printf "\nExecuting dry-run for HEMCO prior emissions run...\n"
        ../GEOSChem_build/gcclassic --dryrun &> log.dryrun
        # prevent restart file from getting downloaded
        sed -i '/GEOSChem.Restart/d' log.dryrun
        # prevent download of GEOS met fields
        sed -i "/GEOS_${Res}/d" log.dryrun
        python download_gc_data.py log.dryrun aws
        popd
    fi

    printf "\n=== GENERATING PRIOR EMISSIONS WITH HEMCO STANDALONE ===\n"

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
    elif [ "$Res" == "0.125x0.15625" ]; then
       resnum="5"
    else
       printf "\nERROR: Grid resolution ${Res} is not supported by the IMI."
       printf "\n Options are 0.125x0.15625, 0.25x0.3125, 0.5x0.625, 2.0x2.5, or 4.0x5.0.\n"
       exit 1
    fi
    HEMCOconfig=${RunTemplate}/HEMCO_Config.rc

    # Create HEMCO standalone directory
    cmd="${metNum}\n${resnum}\n${HEMCOconfig}\n${RunDirs}\n${HEMCOdir}\nn\n"
    printf ${cmd} | ./createRunDir.sh >>createHemcoDir.log 2>&1
    rm -f createHemcoDir.log
    printf "\nCreated ${RunDirs}/${HEMCOdir}\n"

    cd ${RunDirs}/${HEMCOdir}

    # Modify HEMCO files based on settings in config.yml
    sed -i -e "/DiagnFreq:           00000100 000000/d" \
        -e "/Negative values:     0/d" HEMCO_sa_Config.rc
    sed -i -e "s/METEOROLOGY            :       true/METEOROLOGY            :       false/g" \
        -e "s|DiagnFreq:                   End|DiagnFreq:                   Daily|g" HEMCO_Config.rc
    sed -i -e "/#SBATCH -c 8/d" runHEMCO.sh
    sed -i -e "/#SBATCH -t 0-12:00/d" runHEMCO.sh
    sed -i -e "/#SBATCH -p huce_intel/d" runHEMCO.sh
    sed -i -e "/#SBATCH --mem=15000/d" runHEMCO.sh
    sed -i '/.*hemco_standalone.*/a\
retVal=$?\
if [ $retVal -ne 0 ]; then\
    rm -f .error_status_file.txt\
    echo "Error Status: $retVal" > .error_status_file.txt\
    echo "HEMCO prior emis run exited with error code: $retVal"\
    exit $retVal\
fi' runHEMCO.sh

    sa_XMIN=$(python ${InversionPath}/src/components/hemco_prior_emis_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc XMIN)
    sa_XMAX=$(python ${InversionPath}/src/components/hemco_prior_emis_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc XMAX)
    sa_YMIN=$(python ${InversionPath}/src/components/hemco_prior_emis_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc YMIN)
    sa_YMAX=$(python ${InversionPath}/src/components/hemco_prior_emis_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc YMAX)
    sa_NX=$(python ${InversionPath}/src/components/hemco_prior_emis_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc NX)
    sa_NY=$(python ${InversionPath}/src/components/hemco_prior_emis_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc NY)
    sa_YEDGE=$(python ${InversionPath}/src/components/hemco_prior_emis_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc YEDGE)
    sa_YMID=$(python ${InversionPath}/src/components/hemco_prior_emis_component/get_hemco_grid_vars.py ${RunDirs}/StateVector.nc YMID)
    sed -i -e "s/XMIN.*/XMIN: ${sa_XMIN}/" \
        -e "s/XMAX.*/XMAX: ${sa_XMAX}/" \
        -e "s/YMIN.*/YMIN: ${sa_YMIN}/" \
        -e "s/YMAX.*/YMAX: ${sa_YMAX}/" \
        -e "s/NX.*/NX: ${sa_NX}/" \
        -e "s/NY.*/NY: ${sa_NY}/" \
        -e "s/YEDGE.*/YEDGE: ${sa_YEDGE}/" \
        -e "s/YMID.*/YMID: ${sa_YMID}/" HEMCO_sa_Grid.${gridFile}.rc
    mv runHEMCO.sh ${RunName}_HEMCO_Prior_Emis.run

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
    printf "\nDone compiling HEMCO \n\nSee ${RunDirs}/${HEMCOdir}/HEMCO_build_info for details\n\n"

    printf "\nSubmitting HEMCO prior emissions simulation\n\n"

    run_hemco_sa $StartDate $EndDate

    printf "\nDone HEMCO prior emissions simulation\n\n"

    printf "\n=== DONE GENERATING PRIOR EMISSIONS WITH HEMCO STANDALONE ===\n"
    hemco_prior_emis_end=$(date +%s)
}

# Description: Run HEMCO standalone simulation
#     to generate prior emissions for given dates
# Usage: run_hemco_sa <hemco_start> <hemco_end>
run_hemco_sa() {
    hemco_start=$1
    hemco_end=$2
    set -e

    pushd ${RunDirs}/${HEMCOdir}
    # replace start and end times in HEMCO_sa_Time.rc
    sed -i -e "s|START.*|START: ${hemco_start:0:4}-${hemco_start:4:2}-${hemco_start:6:2} 00:00:00|g" \
        -e "s|END.*|END: ${hemco_end:0:4}-${hemco_end:4:2}-${hemco_end:6:2} 00:00:00|g" HEMCO_sa_Time.rc

    rm -f .error_status_file.txt
    # Submit job to job scheduler
    sbatch --mem $RequestedMemory \
        -c $RequestedCPUs \
        -t $RequestedTime \
        -o ${RunName}_HEMCO_Prior_Emis.log \
        -p $SchedulerPartition \
        -W ${RunName}_HEMCO_Prior_Emis.run
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
    python -c '
import sys, warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
import xarray as xr
import numpy as np

emis = xr.load_dataset(sys.argv[1], decode_times=False)
emis["EmisCH4_Total_ExclSoilAbs"] = emis["EmisCH4_Total"] - emis["EmisCH4_SoilAbsorb"]
emis["EmisCH4_Total_ExclSoilAbs"].attrs = emis["EmisCH4_Total"].attrs.copy()
emis["EmisCH4_Total_ExclSoilAbs"].encoding = emis["EmisCH4_Total"].encoding.copy()

if "time" in emis.coords:
    original_units = emis["time"].attrs.get("units", "")
    if "since " in original_units:
        idx = original_units.index("since ") + len("since ")
        date_part = original_units[idx:idx+10]
        unit_part = original_units.split(" since")[0]
        new_units = f"{unit_part} since {date_part} 00:00:00"
        emis["time"].attrs.clear()
        emis["time"].attrs["standard_name"] = "time"
        emis["time"].attrs["long_name"] = "Time"
        emis["time"].attrs["units"] = new_units
        emis["time"].attrs["calendar"] = "standard"
        emis["time"].attrs["axis"] = "T"

keep_vars = [
    var for var in emis.data_vars
    if any(kw in var.lower() for kw in ["time", "lon", "lat", "area"]) or var.startswith("Emis")
]
emis = emis[keep_vars]

# Define the attributes to rename
attrs_to_rename = ["stretch_factor", "target_lat", "target_lon"]

for attr in attrs_to_rename:
    if attr in emis.attrs:
        # Get the value
        value = emis.attrs.pop(attr)
        # Set new attribute with uppercase name
        emis.attrs[attr.upper()] = value
emis.to_netcdf(sys.argv[2])
' "$1" "$2"
}

# Description: Setup Spinup Directory
# Usage:
#   setup_prior_gchp
setup_prior_gchp() {
    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml\n"
        exit 9999
    fi

    printf "\n=== CREATING GCHP Prior RUN DIRECTORY ===\n"

    cd ${RunDirs}

    # Make the directory
    runDir="hemco_prior_emis"
    if [[ -d ${RunDirs}/${runDir} ]]; then
        printf "\n${RunDirs}/${runDir} already exists. Skipping creation.\n"
        return
    fi
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/* ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable
    ln -nsf ../GEOSChem_build/gchp .
    sed -i -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
        -e "s/^TOTAL_CORES=.*/TOTAL_CORES=${TOTAL_CORES}/" \
        -e "s/^NUM_NODES=.*/NUM_NODES=${NUM_NODES}/" \
        -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=${NUM_CORES_PER_NODE}/" \
        setCommonRunSettings.sh

    # regrid restart file to GCHP resolution
    TROPOMIBC=${RestartFilePrefix}${StartDate}_0000z.nc4
    TemplatePrefix="${RunDirs}/${runDir}/Restarts/GEOSChem.Restart.20190101_0000z"
    FilePrefix="GEOSChem.Restart.${StartDate}_0000z"
    cd ../CS_grids
    TROPOMIBC72="temp_tropomi-bc.nc4"
    python ${InversionPath}/src/utilities/regrid_vertgrid_47-to-72.py $TROPOMIBC $TROPOMIBC72
    regrid_tropomi-BC-restart_gcc2gchp ${TROPOMIBC72} ${TemplatePrefix} ${FilePrefix} ${CS_RES} ${STRETCH_GRID} ${STRETCH_FACTOR} ${TARGET_LAT} ${TARGET_LON}
    RestartFile="${RunDirs}/CS_grids/${FilePrefix}.c${CS_RES}.nc4"
    cd ../${runDir}
    ln -nsf $RestartFile Restarts/${FilePrefix}.c${CS_RES}.nc4

    # a temporary fix for GCHP: get day+1 emissions for running GCHP
    RunDuration=$(get_run_duration "$StartDate" "$EndDate")
    NextEndDate=$(date -d "$EndDate +1 day" +%Y%m%d)
    NextRunDuration=$(get_run_duration "$StartDate" "$NextEndDate")
    sed -i -e "s/Run_Duration=\"[0-9]\{8\} 000000\"/Run_Duration=\"${NextRunDuration} 000000\"/" \
        -e "s/Do_Chemistry=.*/Do_Chemistry=false/" \
        -e "s/Do_Advection=.*/Do_Advection=false/" \
        -e "s/Do_Cloud_Conv=.*/Do_Cloud_Conv=false/" \
        -e "s/Do_PBL_Mixing=.*/Do_PBL_Mixing=false/" \
        -e "s/Do_Non_Local_Mixing=.*/Do_Non_Local_Mixing=false/" \
        -e "s/Do_DryDep=.*/Do_DryDep=false/" \
        -e "s/Do_WetDep=.*/Do_WetDep=false/" \
        -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
        setCommonRunSettings.sh

    # get daily emissions output only
    sed -i -e 's/#'\''Emissions/'\''Emissions/g' \
        -e "s/'SpeciesConc',/#'SpeciesConc',/" HISTORY.rc
    sed -i -e "s/Emissions.frequency:[[:space:]]*010000/Emissions.frequency:        240000/" \
        -e "s/Emissions.duration:[[:space:]]*010000/Emissions.duration:         240000/" \
        HISTORY.rc
    # Create run script from template
    sed -e "s:namename:${RunName}_HEMCO_Prior_Emis:g" \
        -e "s:##:#:g" gchp_ch4_run.template >${RunName}_HEMCO_Prior_Emis.run
    chmod 755 ${RunName}_HEMCO_Prior_Emis.run
    rm -f gchp_ch4_run.template

    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING GCHP Prior RUN DIRECTORY ===\n"
}

# Description: Run Prior Directory
# Usage:
#   run_prior_gchp
run_prior_gchp() {
    hemco_start=$1
    hemco_end=$2
    
    prior_start=$(date +%s)
    printf "\n=== SUBMITTING GCHP Prior SIMULATION ===\n"

    cd ${RunDirs}/hemco_prior_emis

    echo "$hemco_start 000000" > cap_restart
    # a temporary fix for GCHP: get day+1 emissions for running GCHP
    RunDuration=$(get_run_duration "$StartDate" "$EndDate")
    NextEndDate=$(date -d "$EndDate +1 day" +%Y%m%d)
    NextRunDuration=$(get_run_duration "$StartDate" "$NextEndDate")
    sed -i -e "s/Run_Duration=\"[0-9]\{8\} 000000\"/Run_Duration=\"${NextRunDuration} 000000\"/" \
        -e "s/Do_Chemistry=.*/Do_Chemistry=false/" \
        -e "s/Do_Advection=.*/Do_Advection=false/" \
        -e "s/Do_Cloud_Conv=.*/Do_Cloud_Conv=false/" \
        -e "s/Do_PBL_Mixing=.*/Do_PBL_Mixing=false/" \
        -e "s/Do_Non_Local_Mixing=.*/Do_Non_Local_Mixing=false/" \
        -e "s/Do_DryDep=.*/Do_DryDep=false/" \
        -e "s/Do_WetDep=.*/Do_WetDep=false/" \
        -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
        setCommonRunSettings.sh

    # Submit job to job scheduler
    sbatch --mem $RequestedMemory \
        -N $NUM_NODES \
        -n $TOTAL_CORES \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -W ${RunName}_HEMCO_Prior_Emis.run
    wait

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    # Remove soil absorption uptake from total emissions
    cd OutputDir
    for file in GEOSChem.Emissions*.nc4; do
        exclude_soil_sink $file $file
    done
    set +e

    # Navigate back to top-level directory
    cd ${RunDirs}

    printf "\n=== DONE GCHP Prior SIMULATION ===\n"
    spinup_end=$(date +%s)
}
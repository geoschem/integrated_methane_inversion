#!/bin/bash

# Functions available in this file include:
#   - setup_spinup
#   - run_spinup

# Description: Setup Spinup Directory
# Usage:
#   setup_spinup
setup_spinup() {
    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml\n"
        exit 9999
    fi

    printf "\n=== CREATING SPINUP RUN DIRECTORY ===\n"

    cd ${RunDirs}

    # Define the run directory name
    SpinupName="${RunName}_Spinup"

    # Make the directory
    runDir="spinup_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/* ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable
    if "$UseGCHP"; then
        sed -i -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
            -e "s/^TOTAL_CORES=.*/TOTAL_CORES=${TOTAL_CORES}/" \
            -e "s/^NUM_NODES=.*/NUM_NODES=${NUM_NODES}/" \
            -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=${NUM_CORES_PER_NODE}/" \
            setCommonRunSettings.sh
        ln -nsf ../GEOSChem_build/gchp .
    else
        ln -nsf ../GEOSChem_build/gcclassic .
    fi

    # Link to restart file
    if "$UseGCHP"; then
        # regrid restart file to GCHP resolution
        TROPOMIBC="${RestartFilePrefix}${SpinupStart}_0000z.nc4"
        TemplatePrefix="${RunDirs}/${runDir}/Restarts/GEOSChem.Restart.20190101_0000z"
        FilePrefix="GEOSChem.Restart.${SpinupStart}_0000z"
        cd ../CS_grids
        TROPOMIBC72="temp_tropomi-bc.nc4"
        python ${InversionPath}/src/utilities/regrid_vertgrid_47-to-72.py $TROPOMIBC $TROPOMIBC72
        regrid_tropomi-BC-restart_gcc2gchp ${TROPOMIBC72} ${TemplatePrefix} ${FilePrefix} ${CS_RES} ${STRETCH_GRID} ${STRETCH_FACTOR} ${TARGET_LAT} ${TARGET_LON}
        RestartFile="${RunDirs}/CS_grids/${FilePrefix}.c${CS_RES}.nc4"
        cd ../${runDir}
        ln -nsf $RestartFile Restarts/GEOSChem.Restart.${SpinupStart}_0000z.c${CS_RES}.nc4
    else
        RestartFile=${RestartFilePrefix}${SpinupStart}_0000z.nc4
        ln -nsf $RestartFile Restarts/GEOSChem.Restart.${SpinupStart}_0000z.nc4
        if "$UseBCsForRestart"; then
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n"
        fi
    fi

    # Update settings in geoschem_config.yml
    if "$UseGCHP"; then
        # Convert months into years and remaining months
        years=$(( ${SpinupMonths} / 12 ))
        months=$(( ${SpinupMonths} % 12 ))
        days=0

        # Format to YYYYMMDD
        SpinupDuration=$(printf "%04d%02d%02d" $years $months $days)
        sed -i -e "s/Run_Duration=\"[0-9]\{8\} 000000\"/Run_Duration=\"${SpinupDuration} 000000\"/" \
            -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
            -e "s/^TOTAL_CORES=.*/TOTAL_CORES=${TOTAL_CORES}/" \
            -e "s/^NUM_NODES=.*/NUM_NODES=${NUM_NODES}/" \
            -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=${NUM_CORES_PER_NODE}/" \
            setCommonRunSettings.sh
        echo "$SpinupStart 000000" > cap_restart
    else
        sed -i -e "s|${StartDate}|${SpinupStart}|g" \
            -e "s|${EndDate}|${SpinupEnd}|g" geoschem_config.yml
    fi

    # Turn on LevelEdgeDiags output
    if "$HourlyCH4"; then
        if "$UseGCHP"; then
            sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                -e 's/LevelEdgeDiags.frequency:.*/LevelEdgeDiags.frequency:      010000/g' \
                -e 's/LevelEdgeDiags.duration:.*/LevelEdgeDiags.duration:    240000/g' HISTORY.rc
        else
            sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
                -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' HISTORY.rc
        fi
    fi

    # Create run script from template
    if "$UseGCHP"; then
        sed -e "s:namename:${SpinupName}:g" \
            -e "s:##:#:g" gchp_ch4_run.template >${SpinupName}.run
        chmod 755 ${SpinupName}.run
        rm -f gchp_ch4_run.template
    else
        sed -e "s:namename:${SpinupName}:g" \
            -e "s:##:#:g" ch4_run.template >${SpinupName}.run
        chmod 755 ${SpinupName}.run
        rm -f ch4_run.template
    fi

    ### Perform dry run if requested
    if [ "$UseGCHP" != "true" ]; then
        if "$SpinupDryrun"; then
            printf "\nExecuting dry-run for spinup run...\n"
            ./gcclassic --dryrun &>log.dryrun
            # prevent restart file from getting downloaded since
            # we don't want to overwrite the one we link to above
            sed -i '/GEOSChem.Restart/d' log.dryrun
            python download_gc_data.py log.dryrun aws
        fi
    fi
    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING SPINUP RUN DIRECTORY ===\n"
}

# Description: Run Spinup Directory
# Usage:
#   run_spinup
run_spinup() {
    spinup_start=$(date +%s)
    printf "\n=== SUBMITTING SPINUP SIMULATION ===\n"

    cd ${RunDirs}/spinup_run

    # Submit job to job scheduler
    if "$UseGCHP"; then
        sbatch --mem $RequestedMemory \
            -N $NUM_NODES \
            -n $TOTAL_CORES \
            -t $RequestedTime \
            -p $SchedulerPartition \
            -W ${RunName}_Spinup.run
    else
        sbatch --mem $RequestedMemory \
            -c $RequestedCPUs \
            -t $RequestedTime \
            -p $SchedulerPartition \
            -W ${RunName}_Spinup.run
    fi
    wait

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    printf "\n=== DONE SPINUP SIMULATION ===\n"
    spinup_end=$(date +%s)
}

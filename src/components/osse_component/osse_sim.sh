#!/bin/bash

# Functions available in this file include:
#   - setup_osse
#   - run_osse

# Description: Setup OSSE Directory
# Usage:
#   setup_osse
setup_osse() {
    set -e
    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml\n"
        exit 9999
    fi

    printf "\n=== CREATING OSSE OBS RUN DIRECTORY ===\n"

    cd ${RunDirs}

    # Define the run directory name
    OSSEName="${RunName}_OSSE_Observations"

    # Make the directory
    runDir="osse_observations_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/* ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable
    ln -sfn ../GEOSChem_build/gcclassic .

    # Update settings in geoschem_config.yml
    # sed -i -e "s|${StartDate}|${SpinupStart}|g" \
    #     -e "s|${EndDate}|${SpinupEnd}|g" geoschem_config.yml

    # Turn on LevelEdgeDiags output
    if "$HourlyCH4"; then
        sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
            -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
            -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
            -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
    fi

    # Create run script from template
    sed -e "s:namename:${OSSEName}:g" \
        -e "s:##:#:g" ch4_run.template >${OSSEName}.run
    chmod 755 ${OSSEName}.run
    rm -f ch4_run.template

    ### Perform dry run if requested
    if "$ProductionDryRun"; then
        printf "\nExecuting dry-run for OSSE run...\n"
        ./gcclassic --dryrun &>log.dryrun
        # prevent restart file from getting downloaded since
        # we don't want to overwrite the one we link to above
        sed -i '/GEOSChem.Restart/d' log.dryrun
        python download_gc_data.py log.dryrun aws
    fi

    # Create random scaling factors
    if "$CreateAutomaticScaleFactorFileOSSE"; then
        printf "\nCreating random scaling factors for OSSE run...\n"
        python ${InversionPath}/src/components/osse_component/make_random_sf.py \
        ${RunDirs}/StateVector.nc \
        ${InversionPath}/${ConfigFile}
    else
        # Copy custom scale factor to $OSSEDirs directory for later use
        printf "\nCopying scale factor file\n"
        cp -v $ScaleFactorFileOSSE ${RunDirs}/osse_observations_run/ScaleFactors.nc
        python ${InversionPath}/src/components/osse_component/make_random_sf.py \
        ${RunDirs}/StateVector.nc \
        ${InversionPath}/${ConfigFile}
    fi

    # Link to restart file
    RestartFileFromSpinup=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        ln -sfn $RestartFileFromSpinup Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
    else
        RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
        ln -sfn $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
        if "$UseBCsForRestart"; then
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n"
        fi
    fi

    # Apply random scale factors to emissions and use precomputed emissions
    # But exclude soil absorption from the application of the scale factors
    # Read in soil absorption separately with MeMo_SOIL_ABSORPTION
    sed -i -e "s|\.\./\.\.|\.\.|g" \
        -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" \
        -e "s|--> Emis_PosteriorSF       :       false|--> Emis_PosteriorSF       :       true|g" \
        -e "s|--> UseTotalPriorEmis      :       false|--> UseTotalPriorEmis      :       true|g" \
        -e "/(((MeMo_SOIL_ABSORPTION/i ))).not.UseTotalPriorEmis" \
        -e "/)))MeMo_SOIL_ABSORPTION/a (((.not.UseTotalPriorEmis" \
        -e "s|gridded_posterior.nc|ScaleFactors.nc|g" \
        -e "s|GFED                   : on|GFED                   : off|g" HEMCO_Config.rc


    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING OSSE RUN DIRECTORY ===\n"
    set +e
}

# Description: Run OSSE Directory
# Usage:
#   run_osse
run_osse() {
    set -e
    osse_start=$(date +%s)
    printf "\n=== SUBMITTING OSSE SIMULATION ===\n"

    cd ${RunDirs}/osse_observations_run

    # Submit job to job scheduler
    sbatch --mem $RequestedMemory \
        -c $RequestedCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -W ${OSSEName}.run
    wait

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    printf "\n=== DONE OSSE SIMULATION ===\n"
    osse_end=$(date +%s)
    set +e
}

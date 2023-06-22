#!/bin/bash

# Functions available in this file include:
#   - setup_posterior 
#   - run_posterior 

# Description: Setup posterior GCClassic run directory
# Usage:
#   setup_posterior
setup_posterior() {
    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml" 
        exit 9999
    fi

    printf "\n=== CREATING POSTERIOR RUN DIRECTORY ===\n"
    
    cd ${RunDirs}
    
    # Define the run directory name
    PosteriorName="${RunName}_Posterior"

    # Make the directory
    runDir="posterior_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/*  ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFileFromSpinup=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        ln -s $RestartFileFromSpinup Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
    else
        RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
        ln -s $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
        if "$UseBCsForRestart"; then
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n" 
        fi
    fi
    
    # Update settings in geoschem_config.yml
    sed -i "/analytical_inversion/{N;s/activate: true/activate: false/}" geoschem_config.yml
    sed -i "s/use_emission_scale_factor: false/use_emission_scale_factor: true/g" geoschem_config.yml
    
    # Update settings in HEMCO_Config.rc
    sed -i -e "s|--> Emis_ScaleFactor       :       false|--> Emis_ScaleFactor       :       true|g" \
           -e "s|gridded_posterior.nc|${RunDirs}/inversion/gridded_posterior.nc|g" HEMCO_Config.rc

    # Turn on LevelEdgeDiags output
    if "$HourlyCH4"; then
        sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
               -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
               -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
               -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
    fi

    ### Turn on observation operators if requested, for posterior run
    activate_observations

    # Create run script from template
    sed -e "s:namename:${PosteriorName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PosteriorName}.run
    chmod 755 ${PosteriorName}.run
    rm -f ch4_run.template

    ### Perform dry run if requested
    if "$PosteriorDryRun"; then
        printf "\nExecuting dry-run for posterior run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
    fi
    
    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING POSTERIOR RUN DIRECTORY ===\n"
}


# Description: Run posterior simulation and process output
# Usage:
#   setup_posterior
run_posterior() {
    posterior_start=$(date +%s)
    cd ${RunDirs}/posterior_run
    
    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Submit job to job scheduler
    printf "\n=== SUBMITTING POSTERIOR SIMULATION ===\n"
    sbatch --mem $SimulationMemory \
           -c $SimulationCPUs \
           -t $RequestedTime \
           -p $SchedulerPartition \
           -W ${RunName}_Posterior.run; wait;
    
    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO
    
    printf "\n=== DONE POSTERIOR SIMULATION ===\n"

    cd ${RunDirs}/inversion

    # Fill missing data (first hour of simulation) in posterior output
    PosteriorRunDir="${RunDirs}/posterior_run"
    PrevDir="${RunDirs}/spinup_run"
    printf "\n=== Calling postproc_diags.py for posterior ===\n"
    python postproc_diags.py $RunName $PosteriorRunDir $PrevDir $StartDate; wait
    printf "\n=== DONE -- postproc_diags.py ===\n"

    # Build directory for hourly posterior GEOS-Chem output data
    mkdir -p data_converted_posterior
    mkdir -p data_visualization_posterior
    mkdir -p data_geoschem_posterior
    GCsourcepth="${PosteriorRunDir}/OutputDir"
    GCDir="./data_geoschem_posterior"
    printf "\n=== Calling setup_gc_cache.py for posterior ===\n"
    python setup_gc_cache.py $StartDate $EndDate $GCsourcepth $GCDir; wait
    printf "\n=== DONE -- setup_gc_cache.py ===\n"

    # Sample GEOS-Chem atmosphere with TROPOMI
    LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
    LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
    LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
    LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)
    nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc)
    FetchTROPOMI="False"
    isPost="True"

    printf "\n=== Calling jacobian.py to sample posterior simulation (without jacobian sensitivity analysis) ===\n"
    python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $isPost; wait
    printf "\n=== DONE sampling the posterior simulation ===\n\n"
    posterior_end=$(date +%s)
}

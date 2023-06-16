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
    cp -r ${RunTemplate}/*  ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFile=${RestartFilePrefix}${SpinupStart}_0000z.nc4
    ln -s $RestartFile Restarts/GEOSChem.Restart.${SpinupStart}_0000z.nc4
    if "$UseBCsForRestart"; then
        sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
	printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n" 
    fi
    
    # Update settings in geoschem_config.yml
    sed -i -e "s|${StartDate}|${SpinupStart}|g" \
           -e "s|${EndDate}|${SpinupEnd}|g" geoschem_config.yml
    sed -i "/analytical_inversion/{N;s/activate: true/activate: false/}" geoschem_config.yml

    # Turn on LevelEdgeDiags output
    if "$HourlyCH4"; then
        sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
               -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
               -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
               -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
    fi

    # Create run script from template
    sed -e "s:namename:${SpinupName}:g" \
        -e "s:##:#:g" ch4_run.template > ${SpinupName}.run
    chmod 755 ${SpinupName}.run
    rm -f ch4_run.template

    ### Perform dry run if requested
    if "$SpinupDryrun"; then
        printf "\nExecuting dry-run for spinup run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
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

    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Submit job to job scheduler
    sbatch --mem $SimulationMemory \
    -c $SimulationCPUs \
    -t $RequestedTime \
    -p $SchedulerPartition \
    -W ${RunName}_Spinup.run; wait;

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    printf "\n=== DONE SPINUP SIMULATION ===\n"
    spinup_end=$(date +%s)
}

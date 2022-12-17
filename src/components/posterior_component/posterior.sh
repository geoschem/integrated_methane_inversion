#!/bin/bash

# Functions available in this file include:
#   - setup_posterior 

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
    RestartFileFromSpinup=../../spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
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

    # Create run script from template
    sed -e "s:namename:${PosteriorName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PosteriorName}.run
    chmod 755 ${PosteriorName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${PosteriorName}.run
    fi

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
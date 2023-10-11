#!/bin/bash

# Functions available in this file include:
#   - run_preview 

# Description: Run the IMI Preview
#   The IMI Preview estimates the quality and cost given the set config file
# Usage:
#   run_preview
run_preview() {
    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml\n" 
        exit 9999
    fi

    printf "\n=== CREATING IMI PREVIEW RUN DIRECTORY ===\n"

    cd ${RunDirs}
    
    # Define the preview run name
    PreviewName="${RunName}_Preview"

    # Make the directory
    runDir="preview_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/*  ${runDir}
    cd $runDir

    # Remove old error status file if present
    rm -f .error_status_file.txt
    
    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFilePreview=${RestartFilePreviewPrefix}${StartDate}_0000z.nc4
    ln -s $RestartFilePreview Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
    if "$UseBCsForRestart"; then
        sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
    fi

    # End date for the preview simulation
    PreviewEnd=$(date --date="${StartDate} +1 day" +%Y%m%d)

    # Update settings in geoschem_config.yml
    sed -i -e "s|${EndDate}|${PreviewEnd}|g" geoschem_config.yml
    sed -i "/analytical_inversion/{N;s/activate: true/activate: false/}" geoschem_config.yml

    # Update settings in HEMCO_Config.rc
    sed -i -e "s|DiagnFreq:                   Monthly|DiagnFreq:                   End|g" HEMCO_Config.rc

    # Update for Kalman filter option
    if "$KalmanMode"; then
        sed -i -e "s|use_emission_scale_factor: true|use_emission_scale_factor: false|g" geoschem_config.yml
        sed -i -e "s|--> Emis_ScaleFactor       :       true|--> Emis_ScaleFactor       :       false|g" HEMCO_Config.rc
    fi

    # Create run script from template
    sed -e "s:namename:${PreviewName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PreviewName}.run
    chmod 755 ${PreviewName}.run
    rm -f ch4_run.template

    ### Perform dry run if requested
    if "$PreviewDryRun"; then
        printf "\nExecuting dry-run for preview run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
    fi

    printf "\n=== DONE CREATING PREVIEW RUN DIRECTORY ===\n"

    ##===============##
    ##  Run preview  ##
    ##===============##

    printf "\n=== RUNNING IMI PREVIEW ===\n"

    # Submit preview GEOS-Chem job to job scheduler
    printf "\nRunning preview GEOS-Chem simulation... "
    if "$UseSlurm"; then
        sbatch --mem $SimulationMemory \
               -c $SimulationCPUs \
               -t $RequestedTime \
               -p $SchedulerPartition \
               -W ${RunName}_Preview.run; wait;
    else
        ./${RunName}_Preview.run
    fi

    # Specify inputs for preview script
    config_path=${InversionPath}/${ConfigFile}
    state_vector_path=${RunDirs}/StateVector.nc
    preview_dir=${RunDirs}/${runDir}
    tropomi_cache=${RunDirs}/data_TROPOMI
    preview_file=${InversionPath}/src/inversion_scripts/imi_preview.py

    # Run preview script
    # If running end to end script with sbatch then use
    # sbatch to take advantage of multiple cores
    printf "\nCreating preview plots and statistics... "
    if "$UseSlurm"; then
        chmod +x $preview_file
        sbatch --mem $SimulationMemory \
        -c $SimulationCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -W $preview_file $InversionPath $config_path $state_vector_path $preview_dir $tropomi_cache; wait;
    else
        python $preview_file $InversionPath $config_path $state_vector_path $preview_dir $tropomi_cache
    fi
    printf "\n=== DONE RUNNING IMI PREVIEW ===\n"

    # check if sbatch commands exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    # Navigate back to top-level directory
    cd ..
}


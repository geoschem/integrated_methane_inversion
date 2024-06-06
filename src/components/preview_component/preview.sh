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

    # First run the Preview if necessary to get prior emissions
    # needed for prepare_sf.py
    if [[ ! -d ${RunDirs}/prior_run/OutputDir ]]; then
        printf "\Prior Dir not detected. Running HEMCO for prior emissions as a prerequisite for IMI Preview.\n"
        run_prior
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
    
    # Link to GEOS-Chem executable
    ln -s ../GEOSChem_build/gcclassic .

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


    # Create run script from template
    sed -e "s:namename:${PreviewName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PreviewName}.run
    chmod 755 ${PreviewName}.run
    rm -f ch4_run.template

    ### Perform dry run if requested
    if "$PreviewDryRun"; then
        printf "\nExecuting dry-run for preview run...\n"
        ./gcclassic --dryrun &> log.dryrun
        # prevent restart file from getting downloaded since
        # we don't want to overwrite the one we link to above
        sed -i '/GEOSChem.Restart/d' log.dryrun
        ./download_data.py log.dryrun aws
    fi

    printf "\n=== DONE CREATING PREVIEW RUN DIRECTORY ===\n"

    ##===============##
    ##  Run preview  ##
    ##===============##

    printf "\n=== RUNNING IMI PREVIEW ===\n"

    # Specify inputs for preview script
    config_path=${InversionPath}/${ConfigFile}
    state_vector_path=${RunDirs}/StateVector.nc
    preview_dir=${RunDirs}/${runDir}
    tropomi_cache=${RunDirs}/satellite_data
    preview_file=${InversionPath}/src/inversion_scripts/imi_preview.py

    # Run preview script
    # If running end to end script with sbatch then use
    # sbatch to take advantage of multiple cores
    printf "\nCreating preview plots and statistics... "
    if "$UseSlurm"; then
        chmod +x $preview_file
        sbatch --mem $RequestedMemory \
        -c $RequestedCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -o imi_output.tmp \
        -W $preview_file $InversionPath $ConfigPath $state_vector_path $preview_dir $tropomi_cache; wait;
        cat imi_output.tmp >> ${InversionPath}/imi_output.log
        rm imi_output.tmp
    else
        python $preview_file $InversionPath $ConfigPath $state_vector_path $preview_dir $tropomi_cache
    fi
    printf "\n=== DONE RUNNING IMI PREVIEW ===\n"

    # check if sbatch commands exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    # Navigate back to top-level directory
    cd ..
}


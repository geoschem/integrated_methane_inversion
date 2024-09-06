#!/bin/bash

# Functions available in this file include:
#   - run_preview

# Description: Run the IMI Preview
#   The IMI Preview estimates the quality and cost given the set config file
# Usage:
#   run_preview
run_preview() {

    # First run the HEMCO standalone if necessary to get prior emissions
    # needed for prepare_sf.py
    if [[ ! -d ${RunDirs}/hemco_prior_emis/OutputDir ]]; then
        printf "\hemco_prior_emis directory not detected. Running HEMCO for prior emissions as a prerequisite for IMI Preview.\n"
        run_hemco_prior_emis
    fi

    # Make the preview directory
    runDir="preview"
    mkdir -p -v ${RunDirs}/${runDir}

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
        rm -f .preview_error_status.txt
        chmod +x $preview_file
        sbatch --mem $RequestedMemory \
            -c $RequestedCPUs \
            -t $RequestedTime \
            -p $SchedulerPartition \
            -o imi_output.tmp \
            -W $preview_file $InversionPath $ConfigPath $state_vector_path $preview_dir $tropomi_cache
        wait
        cat imi_output.tmp >>${InversionPath}/imi_output.log
        rm imi_output.tmp
        # check for any errors
        [ ! -f ".preview_error_status.txt" ] || imi_failed $LINENO
    else
        python $preview_file $InversionPath $ConfigPath $state_vector_path $preview_dir $tropomi_cache
    fi
    printf "\n=== DONE RUNNING IMI PREVIEW ===\n"

    # check if sbatch commands exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    # Navigate back to top-level directory
    cd ..
}

#!/bin/bash

# Description:
#   This script will cleanup the directories and files created by 
#   the imi inversion workflow. Note that some files, like the 
#   calculated sensitivities may be necessary for adding additional 
#   observation operators, so be careful when deleting files.
# Usage:
#   cleanup_script.sh [config.yml]

printf "\n###########################################"
printf "\n#            IMI Cleanup Script           #"
printf "\n###########################################\n\n"

# check if config file is provided
if [[ $# == 1 ]] ; then
    ConfigFile=$1
else
    printf "\nError no config file provided. Usage: cleanup_script.sh [config.yml]\n"
    exit 1
fi

# parse config file
source parse_yaml.sh
eval $(parse_yaml $1)

# Path for IMI rundir based on config file
RunDir="${OutputPath}/${RunName}"

# check to see if user wants to do a dry run
printf "\nDo you want to do a dry-run (y/n)?" 
printf "\nThis will print out all files to be deleted but not any deletions.\n"
read response
if [[ "x${response}" == "xy" ]]; then
    printf "\nEntering dry-run mode. No deletions will take place.\n"
    rm_command="echo rm"
else
    printf "\nEntering deletion mode. Caution: All deletions permanent.\n"
    rm_command=rm
fi

# Ask user if they want to remove data_sensitivities directory
printf "\nDo you want to remove the data_sensitivity directories and files (y/n)?\n"
printf "\nNote: removing these files means you will need to recalculate sensitivities for jacobian perturbations."
printf " This is recommended if you do not plan on adding additional observation operators.\n"
read response

# Remove data_sensitivities directories
if [[ "x${response}" == "xy" ]]; then
    # Check if directory exists before removing them
    kf_inversion_dirs="${RunDir}/kf_inversions/"
    inversion_sensi_dir="${RunDir}/inversion/data_sensitivities"
    # look for KF sensitivities
    if [ -d "${kf_inversion_dirs}" ]; then
        printf "\nRemoving data_sensitivity directories: ${kf_inversion_dirs}/period*/data_sensitivities\n"
        find ${kf_inversion_dirs} -type d -path '*period*/data_sensitivities' -exec ${rm_command} -rf {}/* \;
    # look for inversion sensitivities
    elif [ -d "${inversion_sensi_dir}" ]; then
        printf "\nRemoving data_sensitivity directories: ${inversion_sensi_dir}\n"
        ${rm_command} -r ${inversion_sensi_dir}
    else
        printf "\nNo data_sensitivities directories found in ${RunDir}\n"
    fi
else 
    printf "\nSkipping removal of data_sensitivities directories.\n"
fi

# Ask user if they want to remove the jacobian output directories
printf "\nDo you want to remove the jacobian simulation output directories and files (y/n)?\n"
printf "\nNote: removing these files means you will be unable to recalculate sensitivities " 
printf "without rerunning all simulations."
printf "This is recommended if you have already calculated the sensitivities.\n"
read response

# Remove jacobian run output directories
if [[ "x${response}" == "xy" ]]; then
    # Check if directory exists before removing them
    jacobian_dir="${RunDir}/jacobian_runs"
    if [ -d "${jacobian_dir}/${RunName}_0001/OutputDir" ]; then
        printf "\nRemoving jacobian run outputs (except for the prior and background sim): ${jacobian_dir}/{RunName}_****/OutputDir\n"
        find ${jacobian_dir}/ -type d -name "${RunName}_*" ! -name "${RunName}_background" ! -name "${RunName}_0000" -exec ${rm_command} -rf {}/OutputDir/* \;
    else
        printf "\nNo jacobian_run OutputDirs found in ${jacobian_dir}\n"
    fi
else 
    printf "\nSkipping removal of jacobian simulation output directories\n"
fi

# Finish cleanup
printf "\nCleanup complete\n"
exit 0

#!/bin/bash

# Common shell function for the IMI
# Functions available in this file include:
#   - imi_failed 
#   - run_inversion 


# Description: Print error message for if the IMI fails
#   Copy output file to output directory if it exists
# Usage:
#   imi_failed
imi_failed() {
    file=`basename "$0"`
    printf "\nFATAL ERROR on line number ${1} of ${file}: IMI exiting."
    if [ -d "${OutputPath}/${RunName}" ]; then
        cp "${InversionPath}/imi_output.log" "${OutputPath}/${RunName}/imi_output.log"
    fi
    exit 1
}

# Description: replace sbatch, time, memory, and cpu resources
#   in a run script
# Usage:
#   replace_sbatch_resources 3 10000 runscript.run
replace_sbatch_resources() {
    requested_cpus=$1
    requested_mem=$2
    file=$3
    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" $file
    fi
    sed -i -e "s:{SIMULATION_CPUS}:${JacobianCPUs}:g" \
           -e "s:{SIMULATION_MEMORY}:${JacobianMemory}:g" ${file}
}
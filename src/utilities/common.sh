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

# Description: remove sbatch memory and cpu resources
# remove time too on AWS
# Usage:
#   remove_sbatch_headers runscript.run
remove_sbatch_headers() {
    file=$1
    sed -i -e "/#SBATCH -n/d" $file
    sed -i -e "/#SBATCH --mem/d" $file
    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" $file
    fi
}
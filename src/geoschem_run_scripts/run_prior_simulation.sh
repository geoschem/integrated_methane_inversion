#!/bin/bash
#SBATCH -J {RunName}
#SBATCH -N 1

### Run directory
RUNDIR=$(pwd -P)

### Get current task ID
xstr="0000"

# This checks for the presence of the error status file. If present, this indicates 
# a previous prior sim exited with an error, so this prior will not run
FILE=.error_status_file.txt
if test -f "$FILE"; then
    echo "$FILE exists. Exiting."
    echo "prior simulation: ${xstr} exited without running."
    exit 1
fi

### Run GEOS-Chem in the directory corresponding to the cluster Id
cd  ${RUNDIR}/{RunName}_${xstr}
./{RunName}_${xstr}.run

# save the exit code of the prior simulation cmd
retVal=$?

# Check whether the prior sim finished successfully. If not, write to a hidden file. 
# The presence of the .error_status_file.txt indicates whether an error ocurred. 
# This is needed because scripts that set off sbatch jobs have no knowledge of 
# whether the job finished successfully.
if [ $retVal -ne 0 ]; then
    rm -f .error_status_file.txt
    echo "Error Status: $retVal" > ../.error_status_file.txt
    echo "prior simulation: ${xstr} exited with error code: $retVal"
    echo "Check the log file in the ${RUNDIR}/{RunName}_${xstr} directory for more details."
    exit $retVal
fi

echo "finished prior simulation: ${xstr}"

exit 0

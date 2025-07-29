#!/bin/bash
#SBATCH -J {RunName}

### Run directory
RUNDIR=$(pwd -P)

### Get current task ID
xstr="0000"

# rerunning prior simulation 
rm -f .error_status_file.txt

### Run GEOS-Chem in the directory corresponding to the cluster Id
cd  ${RUNDIR}/{RunName}_${xstr}
if {UseGCHP}; then
    ./cleanRunDir.sh
    echo "{StartDate} 000000" > cap_restart
    sed -i -e "s/Run_Duration=\"[0-9]\{8\} 000000\"/Run_Duration=\"{RunDuration} 000000\"/" \
        setCommonRunSettings.sh
fi
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

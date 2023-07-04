#!/bin/bash
#SBATCH -J {RunName}
#SBATCH -N 1

### Run directory
RUNDIR=$(pwd -P)

### Get current task ID
x=${SLURM_ARRAY_TASK_ID}

### Add zeros to the cluster Id
if [ $x -lt 10 ]; then
    xstr="000${x}"
elif [ $x -lt 100 ]; then
    xstr="00${x}"
elif [ $x -lt 1000 ]; then
    xstr="0${x}"
else
    xstr="${x}"
fi

output_log_file={InversionPath}/imi_output.log

# This checks for the presence of the error status file. If present, this indicates 
# a prior jacobian exited with an error, so this jacobian will not run
FILE=.error_status_file.txt
if test -f "$FILE"; then
    echo "$FILE exists. Exiting."
    echo "jacobian simulation: ${xstr} exited without running." >> $output_log_file
    exit 1
fi

### Run GEOS-Chem in the directory corresponding to the cluster Id
cd  ${RUNDIR}/{RunName}_${xstr}
./{RunName}_${xstr}.run

# save the exit code of the jacobian simulation cmd
retVal=$?

# Check whether the jacobian finished successfully. If not, write to a hidden file. 
# The presence of the .error_status_file.txt indicates whether an error ocurred. 
# This is needed because scripts that set off sbatch jobs have no knowledge of 
# whether the job finished successfully.
if [ $retVal -ne 0 ]; then
    rm -f .error_status_file.txt
    echo "Error Status: $retVal" > ../.error_status_file.txt
    echo "jacobian simulation: ${xstr} exited with error code: $retVal" >> $output_log_file
    echo "Check the log file in the ${RUNDIR}/{RunName}_${xstr} directory for more details." >> $output_log_file
    exit $retVal
fi

echo "finished jacobian simulation: ${xstr}" >> $output_log_file

exit 0

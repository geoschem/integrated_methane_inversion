#!/bin/bash
#SBATCH -J {RunName}

is_valid_nc() {
    local file="$1"
    ncks -m "$file" > /dev/null 2>&1
}

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

# This checks for the presence of the error status file. If present, this indicates
# a prior jacobian exited with an error, so this jacobian will not run
FILE=.error_status_file.txt
if test -f "$FILE"; then
    echo "$FILE exists. Exiting."
    echo "jacobian simulation: ${xstr} exited without running."
    exit 1
fi

if {ReDoJacobian}; then

    cd ${RUNDIR}/{RunName}_${xstr}

    # check for last conc file
    # it has only 1 timestep, so we
    # only need check its existence
    # rerun if it is not there
    yyyymmdd={EndDate}
    LastConcFile=$(date -d "${yyyymmdd} -1 day" +GEOSChem.SpeciesConc.%Y%m%d_0000z.nc4)
    next_day=$(date -d "${yyyymmdd} +1 day" +%Y%m%d)

    if {UseGCHP}; then
        date_val=$(cat cap_restart | sed 's/ .*//')
        if [[ "$date_val" == "$next_day" ]]; then
            echo "Not re-running jacobian simulation: ${xstr}"
            exit 0
        else
            ### Run GEOS-Chem in the directory corresponding to the cluster Id
            echo "Re-running jacobian simulation: ${xstr}"
            ./cleanRunDir.sh
            echo "{StartDate} 000000" > cap_restart
            ./{RunName}_${xstr}.run
            # save the exit code of the jacobian simulation cmd
            retVal=$?
        fi
    else
        if is_valid_nc "OutputDir/$LastConcFile"; then
            echo "Not re-running jacobian simulation: ${xstr}"
            exit 0
        else
            ### Run GEOS-Chem in the directory corresponding to the cluster Id
            echo "Re-running jacobian simulation: ${xstr}"
            ./{RunName}_${xstr}.run
            # save the exit code of the jacobian simulation cmd
            retVal=$?
        fi
    fi
else
    ### Run GEOS-Chem in the directory corresponding to the cluster Id
    cd  ${RUNDIR}/{RunName}_${xstr}
    ./{RunName}_${xstr}.run

    # save the exit code of the jacobian simulation cmd
    retVal=$?
fi

# Check whether the jacobian finished successfully. If not, write to a hidden file.
# The presence of the .error_status_file.txt indicates whether an error ocurred.
# This is needed because scripts that set off sbatch jobs have no knowledge of
# whether the job finished successfully.
if [ $retVal -ne 0 ]; then
    rm -f .error_status_file.txt
    echo "Error Status: $retVal" > ../.error_status_file.txt
    echo "jacobian simulation: ${xstr} exited with error code: $retVal"
    echo "Check the log file in the ${RUNDIR}/{RunName}_${xstr} directory for more details."
    exit $retVal
fi

echo "finished jacobian simulation: ${xstr}"

exit 0
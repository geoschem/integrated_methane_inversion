#!/bin/bash
#SBATCH -J {RunName}

is_valid_nc() {
    local file="$1"
    local yyyymmdd="$2"

    # Validate file structure
    if ! ncks -m "$file" > /dev/null 2>&1; then
        return 1
    fi

    # Extract last time string
    local last_time
    last_time=$(ncdump -t -v time "$file" | grep -oE '"[^"]+"' | tail -n 1 | tr -d '"')

    # Format YYYYMMDD to YYYY-MM-DD
    local expected_time="${yyyymmdd:0:4}-${yyyymmdd:4:2}-${yyyymmdd:6:2} 23"

    if [[ "$last_time" != "$expected_time" ]]; then
        return 1
    fi

    return 0
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
    # it has 24 timestep
    # check if it is valid and has 24 entries of time
    yyyymmdd={EndDate}
    last_date=$(date -d "${yyyymmdd} -1 day" +%Y%m%d)
    LastConcFile="GEOSChem.SpeciesConc.${last_date}_0000z.nc4"

    if is_valid_nc "OutputDir/$LastConcFile" "$last_date"; then
        echo "Not re-running jacobian simulation: ${xstr}"
        exit 0
    else
        ### Run GEOS-Chem in the directory corresponding to the cluster Id
        echo "Re-running jacobian simulation: ${xstr}"
        if {UseGCHP}; then
            ./cleanRunDir.sh
            echo "{StartDate} 000000" > cap_restart
        fi
        ./{RunName}_${xstr}.run
        # save the exit code of the jacobian simulation cmd
        retVal=$?
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
#!/bin/bash
#SBATCH -J {RunName}

is_valid_nc() {
    local file="$1"

    # Validate file structure
    if ! ncks -m "$file" > /dev/null 2>&1; then
        return 1
    fi

    # Extract the length of time dimension
    local time_len
    time_len=$(ncdump -h "$file" \
        | grep -- "time = UNLIMITED" \
        | sed -E 's/.*\(([0-9]+) currently\).*/\1/')

    # Check extraction worked and length matches
    if [[ -z "$time_len" || "$time_len" != "24" ]]; then
        return 1
    fi

    return 0
}

# Used both before submission and by workers. EndDate is exclusive and is
# rendered for the current Kalman period when this script is copied.
if {ReDoJacobian}; then
    if ! last_date=$(date -d "{EndDate} -1 day" +%Y%m%d); then
        echo "Cannot determine the last Jacobian output date." >&2
        exit 1
    fi
    LastConcFile="GEOSChem.SpeciesConc.${last_date}_0000z.nc4"
fi

jacobian_is_complete() {
    {ReDoJacobian} && is_valid_nc "$1/OutputDir/$LastConcFile"
}

# Print only the comma-separated array indices to stdout. No model runs or
# error-marker changes occur in this mode. Keep original IDs (do not renumber).
if [[ ${1:-} == --list-pending ]]; then
    pending_ids=""
    for ((candidate=$2; candidate<=$3; candidate++)); do
        # this line does not print to stdout, rather it defines the candidate_name variable
        printf -v candidate_name '{RunName}_%04d' "$candidate"
        if jacobian_is_complete "$candidate_name"; then
            echo "Not submitting completed Jacobian simulation: $candidate_name" >&2
        else
            pending_ids="${pending_ids:+${pending_ids},}${candidate}"
        fi
    done
    printf '%s\n' "$pending_ids"
    exit 0
fi

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
    if jacobian_is_complete .; then
        echo "Not re-running jacobian simulation: ${xstr}"
        exit 0
    else
        ### Run GEOS-Chem in the directory corresponding to the cluster Id
        echo "Re-running jacobian simulation: ${xstr}"
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

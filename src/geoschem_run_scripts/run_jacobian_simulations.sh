#!/bin/bash
#SBATCH -J {RunName}

is_valid_nc() {
    local file="$1"
    local yyyymmdd="$2"

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

fail() {   # usage: fail <code> <message>
    echo "$2"
    echo "Error Status: $1" > "${RUNDIR}/.error_status_file.txt"
    exit "$1"
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

end_yyyymmdd={EndDate}
start_yyyymmdd={StartDate}

cd "${RUNDIR}/{RunName}_${xstr}" || fail 1 "Cannot cd to run dir for jacobian simulation: ${xstr}"

if {ContinueIncompleteJacobians}; then

    # Check if last conc file is complete (has 24 timesteps)
    # if yes, exit as sim ran fully
    last_date=$(date -d "${end_yyyymmdd} -1 day" +%Y%m%d)
    LastConcFile="GEOSChem.SpeciesConc.${last_date}_0000z.nc4"
    if is_valid_nc "OutputDir/$LastConcFile"; then
        echo "Not re-running jacobian simulation: ${xstr}"
        exit 0
    fi

    # Find the latest usable restart in [start, end)
    RestartDate=""
    RestartTime=""
    RestartPath=""
    
    # looping over restarts, newest first
    for f in $(ls -1 Restarts/GEOSChem.Restart.????????_????z.nc4 2>/dev/null | sort -r); do
        rst=$(basename "$f")
        d=$(echo "$rst" | cut -d. -f3 | cut -d_ -f1)
        t=$(echo "$rst" | cut -d. -f3 | cut -d_ -f2 | cut -dz -f1)

	# skip if restart is before start date or after end date
        if [ "$d" -lt "$start_yyyymmdd" ] || [ "$d" -ge "$end_yyyymmdd" ]; then
            continue
        fi

	# skips restarts that are incomplete (from killed runs)
        if ! ncks -m "$f" > /dev/null 2>&1; then
            echo "skipping unreadable restart file: ${rst}"
            continue
        fi

	# loop exits on first success
        RestartDate=$d
        RestartTime="${t}00"
        RestartPath=$f
        break
    done

    ## Handling when restart is found and is NOT start-date (true redo case) 
    if [ -n "$RestartDate" ] && [ "$RestartDate" != "$start_yyyymmdd" ]; then
	
	# look for backup geoschem and hemco config files, create one if it doesn't exist
	# preserves original settings when restarting several times
        [ -f geoschem_config.yml.orig ] || cp -p geoschem_config.yml geoschem_config.yml.orig
        [ -f HEMCO_Config.rc.orig ]     || cp -p HEMCO_Config.rc HEMCO_Config.rc.orig

	# replacing start date in geoschem_config
        sed -i -E "s/^([[:space:]]*start_date:[[:space:]]*\[)[0-9]{8}(,[[:space:]]*)[0-9]+(\])/\1${RestartDate}\2${RestartTime}\3/" geoschem_config.yml
        
	# checking restart replace worked (nothing strange in formatting)
	check=$(grep -E '^[[:space:]]*start_date:' geoschem_config.yml | grep -oE '[0-9]{8}' | head -1)
        [ "$check" = "$RestartDate" ] || fail 1 "Failed to set start_date for jacobian simulation: ${xstr}"

        # Match the HEMCO restart variable prefix to what the file actually contains
        if ncdump -h "$RestartPath" | grep -q 'SpeciesRst_'; then
            prefix=SpeciesRst_
            sed -i '/GEOSChem\.Restart/ s/SpeciesBC_/SpeciesRst_/g' HEMCO_Config.rc
        else
            prefix=SpeciesBC_
            sed -i '/GEOSChem\.Restart/ s/SpeciesRst_/SpeciesBC_/g' HEMCO_Config.rc
        fi

	# checking replace worked
        grep 'GEOSChem\.Restart' HEMCO_Config.rc | grep -q "$prefix" \
            || fail 1 "Failed to set restart prefix in HEMCO_Config.rc for jacobian simulation: ${xstr}"

        echo "Resuming jacobian simulation: ${xstr} from ${RestartDate} (was ${start_yyyymmdd}), restart prefix ${prefix}"
    else
	## Handling when no restart is found or restart = startdate
    	# Resets to the original geoschem and hemco config state (a clean run)
        [ -f geoschem_config.yml.orig ] && cp -p geoschem_config.yml.orig geoschem_config.yml
        [ -f HEMCO_Config.rc.orig ]     && cp -p HEMCO_Config.rc.orig HEMCO_Config.rc
        echo "No mid-run restart for jacobian simulation: ${xstr}; running from ${start_yyyymmdd}"
    fi
fi

### Run GEOS-Chem in the directory corresponding to the cluster Id
echo "Running jacobian simulation: ${xstr}"
./{RunName}_${xstr}.run
retVal=$?

if [ "$retVal" -ne 0 ]; then
    echo "jacobian simulation: ${xstr} exited with error code: $retVal"
    echo "Check the log file in the ${RUNDIR}/{RunName}_${xstr} directory for more details."
    fail "$retVal" "Error Status: $retVal"
fi

echo "finished jacobian simulation: ${xstr}"
exit 0

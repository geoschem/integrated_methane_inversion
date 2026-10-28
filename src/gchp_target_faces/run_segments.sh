#!/bin/bash

# Run the Jacobian runs and the inversion stage segment by segment, for the
# target faces set in gchp_imi_common.sh, until a target end date. No
# questions are asked, so it can run unattended. Run it on the head node
# inside tmux:
#   tmux new -s segments
#   ./run_segments.sh <target_end YYYYMMDD> <segment_days>
# For example, with StartDate 20250101,
#   ./run_segments.sh 20250106 2
# runs segments ending 20250103, 20250105 and 20250106.
#
# For each segment:
#   1. Submit every Jacobian run that has not reached the segment end
#   2. Wait until every run has reached it, checking every poll_minutes.
#      Runs that fail, or stop before the segment end, are resubmitted once.
#   3. Run the inversion stage for the segment (run_gchp_imi_inversion.sh),
#      including pruning if prune_output is true
#
# It stops before starting the next segment if runs still have not reached
# the segment end after one resubmission, or if the inversion stage of any
# face is not DONE. Fix the cause and run it again with the same arguments:
# it continues after the last segment that every face has completed.
#
# target_end must not be after EndDate in the base IMI config. EndDate stays
# the end of the full inversion period.
#
# Progress is written to the screen and to <test_dir>/logs/run_segments.log.

set -e
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

#--------------------------------------------------------------------------
# Check arguments
#--------------------------------------------------------------------------
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <target_end YYYYMMDD> <segment_days>"
    exit 1
fi
target_end=$1
segment_days=$2

start_date=$(grep "^StartDate:" "${template_config}" | awk '{print $2}' | tr -d '"')
end_date=$(grep "^EndDate:" "${template_config}" | awk '{print $2}' | tr -d '"')

if [[ ! "${target_end}" =~ ^[0-9]{8}$ ]] || ! date -d "${target_end}" >/dev/null 2>&1; then
    echo "ERROR: Target end '${target_end}' is not a valid date YYYYMMDD"
    exit 1
elif [[ "${target_end}" > "${end_date}" ]]; then
    echo "ERROR: Target end ${target_end} is after EndDate ${end_date} in ${template_config}"
    exit 1
elif [[ ! "${segment_days}" =~ ^[0-9]+$ ]] || (( segment_days < 1 )); then
    echo "ERROR: Segment days '${segment_days}' must be a whole number of at least 1"
    exit 1
fi

mkdir -p "${test_dir}/logs"
log_file="${test_dir}/logs/run_segments.log"

#--------------------------------------------------------------------------
# Helper functions
#--------------------------------------------------------------------------

# Function to print a message with the time, to the screen and the log
# Argument: message
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') $1" | tee -a "${log_file}"
}

# Function to print the number of days from one date to a later one
# Arguments: first date, second date (YYYYMMDD)
days_between() {
    echo $(( ( $(date -d "$2" +%s) - $(date -d "$1" +%s) ) / 86400 ))
}

# Function to list the Jacobian run directories of all faces, one per line
all_run_dirs() {
    local i_face face_name run_dir
    for ((i_face=start_face; i_face<=num_faces; i_face++)); do
        face_name=$(face_name_for ${i_face})
        for run_dir in "${test_dir}/${face_name}/jacobian_runs/${face_name}_"*; do
            if [[ -d "${run_dir}" ]]; then
                echo "${run_dir}"
            fi
        done
    done
}

# Function to print the last segment end for which the inversion stage of
# every face is DONE, or StartDate if a face has none
last_done_segment() {
    local i_face face_name done_date oldest=""
    for ((i_face=start_face; i_face<=num_faces; i_face++)); do
        face_name=$(face_name_for ${i_face})
        done_date=$(grep "^DONE" "${test_dir}/${face_name}/inversion_status" 2>/dev/null \
                    | tail -n 1 | awk '{print $2}')
        if [[ -z "${done_date}" ]]; then
            echo "${start_date}"
            return
        fi
        if [[ -z "${oldest}" || "${done_date}" < "${oldest}" ]]; then
            oldest=${done_date}
        fi
    done
    echo "${oldest}"
}

# Function to submit every run that is READY for the segment end, and also
# every FAILED run if the second argument is true
# Arguments: segment end, true or false
submit_segment() {
    local seg_end=$1
    local retry_failed=$2
    local active run_dir status start days msg
    local num_submitted=0

    if ! active=$(active_run_dirs); then
        log "ERROR: squeue failed, so active runs are unknown. Nothing submitted."
        return 1
    fi

    while read -r run_dir; do
        [[ -n "${run_dir}" ]] || continue
        status=$(run_status "${run_dir}" "${seg_end}" "${active}")
        if [[ "${status}" == READY || ( "${status}" == FAILED && "${retry_failed}" == true ) ]]; then
            # The run's maximum days is exactly what it needs to reach the
            # segment end, since this script works out the dates itself
            start=$(cut -c 1-8 "${run_dir}/cap_restart")
            days=$(days_between "${start}" "${seg_end}")
            if ! msg=$(submit_run "${run_dir}" "${seg_end}" "${days}"); then
                log "ERROR: Submitting $(basename "${run_dir}") failed"
                return 1
            fi
            log "  $(basename "${run_dir}") (${status}): ${msg}"
            num_submitted=$((num_submitted + 1))
        fi
    done <<< "$(all_run_dirs)"

    log "Submitted ${num_submitted} run(s) to ${seg_end}"
}

# Function to wait until every run has reached the segment end. Runs that
# failed, or stopped before the segment end, are resubmitted once.
# Returns 1 if runs still have not reached it after that.
# Argument: segment end
wait_for_runs() {
    local seg_end=$1
    local resubmitted=false
    local seen_before=false
    local active run_dir status num_done num_active num_redo num_other

    while true; do
        if ! active=$(active_run_dirs); then
            log "squeue failed; checking again in ${poll_minutes} minutes"
            sleep $((poll_minutes * 60))
            continue
        fi

        num_done=0; num_active=0; num_redo=0; num_other=0
        while read -r run_dir; do
            [[ -n "${run_dir}" ]] || continue
            status=$(run_status "${run_dir}" "${seg_end}" "${active}")
            case "${status}" in
                DONE)         num_done=$((num_done + 1)) ;;
                ACTIVE)       num_active=$((num_active + 1)) ;;
                READY|FAILED) num_redo=$((num_redo + 1)) ;;
                *)            num_other=$((num_other + 1)) ;;
            esac
        done <<< "$(all_run_dirs)"

        log "Runs for ${seg_end}: ${num_done} DONE, ${num_active} ACTIVE, ${num_redo} failed or not finished, ${num_other} other"

        if (( num_other > 0 )); then
            log "ERROR: ${num_other} run(s) have no cap_restart. See ./check_gchp_imi_jacobian.sh ${seg_end} ${segment_days}"
            return 1
        elif (( num_active > 0 )); then
            seen_before=false
            sleep $((poll_minutes * 60))
        elif (( num_redo == 0 )); then
            return 0
        elif [[ "${seen_before}" == false ]]; then
            # A requeued job can briefly be missing from the queue. Check
            # once more before resubmitting, so a run never gets two jobs.
            seen_before=true
            sleep $((poll_minutes * 60))
        elif [[ "${resubmitted}" == true ]]; then
            log "ERROR: ${num_redo} run(s) did not reach ${seg_end} after one resubmission. See ./check_gchp_imi_jacobian.sh ${seg_end} ${segment_days}"
            return 1
        else
            log "Resubmitting ${num_redo} run(s) that failed or did not finish"
            submit_segment "${seg_end}" true || return 1
            resubmitted=true
            seen_before=false
        fi
    done
}

#--------------------------------------------------------------------------
# Run segments until every face is DONE through the target end
#--------------------------------------------------------------------------
log "Running ${segment_days}-day segments to ${target_end} for faces ${start_face} to ${num_faces}"

while true; do
    last_done=$(last_done_segment)
    if [[ ! "${last_done}" < "${target_end}" ]]; then
        log "Every face is DONE through ${target_end}. Finished."
        break
    fi

    seg_end=$(date -d "${last_done} +${segment_days} days" +%Y%m%d)
    if [[ "${seg_end}" > "${target_end}" ]]; then
        seg_end=${target_end}
    fi
    log "=== Segment ending ${seg_end} ==="

    # 1. Submit the Jacobian runs
    if ! submit_segment "${seg_end}" false; then
        log "Stopped. Fix the problem and run again with the same arguments."
        exit 1
    fi

    # 2. Wait for all runs to reach the segment end
    if ! wait_for_runs "${seg_end}"; then
        log "Stopped. Fix the problem and run again with the same arguments."
        exit 1
    fi

    # 3. Run the inversion stage, then check every face is DONE
    log "All runs reached ${seg_end}. Running the inversion stage."
    ./run_gchp_imi_inversion.sh "${seg_end}" >> "${log_file}" 2>&1 || true
    if [[ "$(last_done_segment)" < "${seg_end}" ]]; then
        log "ERROR: The inversion stage is not DONE for every face. See the end of ${log_file} and each face's inversion_status."
        log "Stopped. Fix the problem and run again with the same arguments."
        exit 1
    fi
    log "Inversion stage DONE for every face through ${seg_end}"
done

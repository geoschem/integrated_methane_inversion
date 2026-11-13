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
# Within a segment, faces are started in order of face number, but only
# while the space they will need fits under fsx_max_use_percent of FSx
# (see "Storage" in gchp_imi_common.sh). Every poll_minutes it checks each
# face:
#   - Not started: submit its Jacobian runs if its space fits. Faces are
#     never started out of order.
#   - Runs going: wait. Runs that failed, or stopped before the segment end,
#     are resubmitted once.
#   - All runs at the segment end: run its inversion stage (inversion_face.sh,
#     up to max_parallel_inversions faces at a time), which prunes its
#     Jacobian run output if prune_output is true. That frees its space for
#     the next faces.
#   - Inversion stage DONE: if release_done_faces is true, release its large
#     files that are safely in S3 (restarts, the hourly output kept for the
#     next inversion stage, overpass and .pkl files), so the next faces can
#     use the space. Files not exported yet are tried again at the next
#     check. They are loaded back from S3 when the next segment reads them.
# The next segment starts only when every face has finished this one. Then
# met fields and exported products that are no longer needed are released
# from FSx (they stay in S3), if release_met and release_products are true.
#
# If something fails (runs fail again after one resubmission, an inversion
# stage is not DONE, or a face cannot fit even on its own), no more faces
# are started; it waits for running inversion stages to end, then stops.
# Fix the cause and run it again with the same arguments: it continues from
# where things stand, after the last segment that every face completed.
#
# target_end must not be after EndDate in the base IMI config. EndDate stays
# the end of the full inversion period.
#
# Progress is written to the screen and to <test_dir>/logs/run_segments.log.

set -e
shopt -s nullglob
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

#--------------------------------------------------------------------------
# Check arguments and settings
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

# Function to print the dates from one date to another, inclusive, one per line
# Arguments: first date, last date (YYYYMMDD)
dates_from_to() {
    local d=$1
    while [[ ! "${d}" > "$2" ]]; do
        echo "${d}"
        d=$(date -d "${d} +1 days" +%Y%m%d)
    done
}

# Function to print a size in bytes as GB, for the log
# Argument: bytes
gb() {
    awk -v b="$1" 'BEGIN { printf "%.0f GB", b / 1e9 }'
}

# Function to test whether a face's inversion stage is DONE for a segment
# end (or a later one)
# Arguments: face name, segment end
face_done() {
    local done_date
    done_date=$(grep "^DONE" "${test_dir}/$1/inversion_status" 2>/dev/null \
                | tail -n 1 | awk '{print $2}')
    [[ -n "${done_date}" && ! "${done_date}" < "$2" ]]
}

# Function to print the last segment end for which the inversion stage of
# every face is DONE, or StartDate if a face has none
last_done_segment() {
    local i_face face_name done_date oldest=""
    end_face=$(( start_face + num_faces - 1 ))
    for ((i_face=start_face; i_face<=end_face; i_face++)); do
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

# Function to print the space, in bytes, a face still needs on FSx in this
# segment: its hourly output at the peak (segment days plus the day kept
# from the previous segment) and the restarts written during the segment,
# minus the hourly output already taking space on FSx. Files released, or
# not loaded from S3 yet, take no space and are not subtracted: they take
# space again when read. Never less than zero.
# Argument: face number
face_space_needed() {
    local face_name num_sv written
    face_name=$(face_name_for "$1")
    num_sv=$(num_sv_for "$1")
    written=$(du -sc --block-size=1 "${test_dir}/${face_name}/jacobian_runs/"*/OutputDir 2>/dev/null \
              | tail -n 1 | cut -f1)
    awk -v n="${num_sv}" -v d="${segment_days}" -v o="${gb_per_element_day}" \
        -v r="${gb_per_element_restart}" -v f="${storage_safety_factor}" \
        -v w="${written:-0}" \
        'BEGIN { need = n * (o * (d + 1) + r * d) * f * 1e9 - w
                 printf "%.0f", (need > 0 ? need : 0) }'
}

# Function to print how many Jacobian and inversion-stage jobs are running,
# starting (node booting) and pending in Slurm, for the status line.
# Jacobian jobs are named <face>_NNNN; inversion-stage jobs are named
# overpass_<face> and inversion_<face>.
job_counts() {
    squeue -h -u "${USER}" -o "%j %T" 2>/dev/null | awk -v p="${target_prefix}" '
        { kind = "" }
        $1 ~ "^" p "_T[0-9]+_[0-9]+$"                { kind = "jac" }
        $1 ~ "^(overpass|inversion)_" p "_T[0-9]+$"  { kind = "inv" }
        kind != "" && ($2 == "RUNNING" || $2 == "COMPLETING") { run[kind]++ }
        kind != "" && $2 == "CONFIGURING"                      { cf[kind]++ }
        kind != "" && $2 == "PENDING"                          { pd[kind]++ }
        END { printf "Jacobian jobs %d running, %d starting, %d pending; inversion jobs %d running, %d starting, %d pending",
              run["jac"], cf["jac"], pd["jac"], run["inv"], cf["inv"], pd["inv"] }'
}

# Function to submit the runs of one face that have not reached the segment
# end. Each run gets exactly the days it needs to reach it.
# Arguments: face name, segment end, true to also resubmit FAILED runs,
#            output of active_run_dirs
submit_face() {
    local face_name=$1
    local seg_end=$2
    local retry_failed=$3
    local active=$4
    local run_dir status start days msg

    for run_dir in "${test_dir}/${face_name}/jacobian_runs/${face_name}_"*; do
        [[ -d "${run_dir}" ]] || continue
        status=$(run_status "${run_dir}" "${seg_end}" "${active}")
        if [[ "${status}" == READY || ( "${status}" == FAILED && "${retry_failed}" == true ) ]]; then
            start=$(cut -c 1-8 "${run_dir}/cap_restart")
            days=$(days_between "${start}" "${seg_end}")
            if ! msg=$(submit_run "${run_dir}" "${seg_end}" "${days}"); then
                log "ERROR: Submitting $(basename "${run_dir}") failed"
                return 1
            fi
            log "  $(basename "${run_dir}") (${status}): ${msg}"
        fi
    done
}

# Function to release files that are safely in S3: archived, not changed
# since, and not already released. This frees their space on FSx; they are
# loaded back from S3 if anything reads them. Files not exported yet are
# left for later. Sets num_released, num_release_failed and
# num_release_waiting for the caller.
# Arguments: label for the log (empty: log only if something was released
#            or failed), then the files
release_archived() {
    local label=$1
    shift
    local states releasable num_releasable num_waiting failures num_failed

    num_released=0; num_release_failed=0; num_release_waiting=0
    if (( $# == 0 )); then
        return 0
    fi
    states=$(printf '%s\0' "$@" | xargs -0 lfs hsm_state 2>/dev/null || true)
    releasable=$(grep " archived" <<< "${states}" | grep -v -e " released" -e "dirty" \
                 | sed 's/: (0x[0-9a-fA-F]*).*//' || true)
    num_releasable=$(grep -c . <<< "${releasable}" || true)
    num_waiting=$(grep -v " archived" <<< "${states}" | grep -c . || true)

    # Release each file separately, so a file that cannot be released does
    # not stop the others. Use one sudo call per batch of files to be fast.
    failures=""
    if (( num_releasable > 0 )); then
        failures=$(tr '\n' '\0' <<< "${releasable}" \
            | xargs -0 sudo -n bash -c 'for f; do err=$(lfs hsm_release "$f" 2>&1) || echo "$f: could not release: $(echo ${err})"; done' _ \
            || true)
    fi
    num_failed=$(grep -c ": could not release" <<< "${failures}" || true)
    num_released=$(( num_releasable - num_failed ))
    num_release_failed=${num_failed}
    num_release_waiting=${num_waiting}
    if (( num_failed > 0 )); then
        log "  WARNING: ${num_failed} ${label:-files} could not be released. The first few:"
        head -n 5 <<< "${failures}" | while read -r line; do log "    ${line}"; done
    fi
    if [[ -n "${label}" ]]; then
        log "  ${label}: released ${num_released}, ${num_failed} failed, ${num_waiting} not exported to S3 yet"
    fi
}

# Function to release the large files of a face whose inversion stage is
# DONE for the segment: restarts, hourly output, overpass and .pkl files.
# Only regular files over 1 MB: small files free no space, and every job
# reading them would wait for them to load. Links are never followed.
# Returns 0 once nothing is left that still needs exporting.
# Argument: face name
release_face() {
    local face_dir="${test_dir}/$1"
    local -a dirs files
    local file
    dirs=( "${face_dir}"/jacobian_runs/*/Restarts "${face_dir}"/jacobian_runs/*/OutputDir
           "${face_dir}"/jacobian_runs/*/OverpassDiagnostics "${face_dir}"/inversion/data_converted )
    if (( ${#dirs[@]} == 0 )); then
        return 0
    fi
    files=()
    while IFS= read -r file; do
        files+=( "${file}" )
    done < <(find "${dirs[@]}" -maxdepth 1 -type f -size +1M 2>/dev/null)
    release_archived "" "${files[@]}"
    if (( num_released > 0 )); then
        log "$1: released ${num_released} files from FSx, ${num_release_waiting} not exported to S3 yet"
    fi
    (( num_release_waiting == 0 ))
}

# Function to release met fields and exported products no longer needed,
# for dates from 3 segments back up to 2 days before the segment end
# (later dates may still be read). Earlier dates were handled before.
# Argument: segment end
release_after_segment() {
    local seg_end=$1
    local first last d files

    last=$(date -d "${seg_end} -2 days" +%Y%m%d)
    first=$(date -d "${seg_end} -$((3 * segment_days + 2)) days" +%Y%m%d)
    if [[ "${first}" < "${start_date}" ]]; then
        first=$(date -d "${start_date} -1 days" +%Y%m%d)
    fi
    if [[ "${last}" < "${first}" ]]; then
        return 0
    fi
    log "Releasing files dated ${first} to ${last} that are no longer needed on FSx"

    if [[ "${release_met}" == true ]]; then
        files=()
        for d in $(dates_from_to "${first}" "${last}"); do
            files+=( "${met_dir}/${d:0:4}/${d:4:2}/GEOSFP.${d}."* )
        done
        release_archived "met fields" "${files[@]}"
    fi

    if [[ "${release_products}" == true ]]; then
        files=()
        for d in $(dates_from_to "${first}" "${last}"); do
            files+=( "${test_dir}"/*/jacobian_runs/*/OverpassDiagnostics/*.overpass."${d}"_*.nc4 )
            files+=( "${test_dir}"/*/inversion/data_converted/*__"${d}"T*.pkl )
        done
        release_archived "overpass and .pkl files" "${files[@]}"
    fi
}

#--------------------------------------------------------------------------
# Function to run one segment for all faces. Returns 1 if something failed.
# Arguments: segment start (last segment end every face completed),
#            segment end
#--------------------------------------------------------------------------
run_segment() {
    local seg_start=$1
    local seg_end=$2
    local -a inversion_pid=()   # by face number: PID of its inversion stage
    local -a started=()         # by face number: true once its runs were submitted
    local -a resubmitted=()     # by face number: true once runs were resubmitted
    local -a seen_stopped=()    # by face number: true if runs were stopped at the last check
    local -a face_released=()   # by face number: true once its large files were released
    local failed=false
    local i_face face_name active fs_size fs_used limit reserved need
    local blocked num_done num_inversion num_running num_waiting num_inv_now
    local run_dir status run_done run_active run_stopped

    while true; do
        if ! active=$(active_run_dirs); then
            log "squeue failed; checking again in ${poll_minutes} minutes"
            sleep $((poll_minutes * 60))
            continue
        fi
        read -r fs_size fs_used < <(df -B1 --output=size,used "${test_dir}" | tail -n 1)
        limit=$(( fs_size / 100 * fsx_max_use_percent ))

        # Inversion stages still running
        num_inv_now=0
	end_face=$(( start_face + num_faces - 1 ))
        for ((i_face=start_face; i_face<=end_face; i_face++)); do
            if [[ -n "${inversion_pid[i_face]}" ]] && kill -0 "${inversion_pid[i_face]}" 2>/dev/null; then
                num_inv_now=$((num_inv_now + 1))
            fi
        done

        reserved=0
        blocked=false
        num_done=0; num_inversion=0; num_running=0; num_waiting=0

        for ((i_face=start_face; i_face<=end_face; i_face++)); do
            face_name=$(face_name_for ${i_face})

            # Finished this segment
            if face_done "${face_name}" "${seg_end}"; then
                if [[ -n "${inversion_pid[i_face]}" ]]; then
                    wait "${inversion_pid[i_face]}" || true
                    inversion_pid[i_face]=""
                    log "${face_name}: inversion stage DONE"
                fi
                # Free its space for the faces still to run this segment.
                # Tried again at each check until every file is in S3.
                if [[ "${release_done_faces}" == true && "${face_released[i_face]}" != true ]]; then
                    if release_face "${face_name}"; then
                        face_released[i_face]=true
                    fi
                fi
                num_done=$((num_done + 1))
                continue
            fi

            # Inversion stage running, or just ended
            if [[ -n "${inversion_pid[i_face]}" ]]; then
                if kill -0 "${inversion_pid[i_face]}" 2>/dev/null; then
                    num_inversion=$((num_inversion + 1))
                    reserved=$(( reserved + $(face_space_needed ${i_face}) ))
                    continue
                fi
                wait "${inversion_pid[i_face]}" || true
                inversion_pid[i_face]=""
                if face_done "${face_name}" "${seg_end}"; then
                    log "${face_name}: inversion stage DONE"
                    num_done=$((num_done + 1))
                else
                    log "ERROR: ${face_name}: inversion stage not DONE. See ${test_dir}/logs/inversion_${face_name}.${seg_end}.log"
                    failed=true
                fi
                continue
            fi

            # State of the face's Jacobian runs
            run_done=0; run_active=0; run_stopped=0
            for run_dir in "${test_dir}/${face_name}/jacobian_runs/${face_name}_"*; do
                [[ -d "${run_dir}" ]] || continue
                status=$(run_status "${run_dir}" "${seg_end}" "${active}")
                case "${status}" in
                    DONE)   run_done=$((run_done + 1)) ;;
                    ACTIVE) run_active=$((run_active + 1)) ;;
                    *)      run_stopped=$((run_stopped + 1)) ;;
                esac
            done

            # All runs at the segment end: start its inversion stage
            if (( run_active == 0 && run_stopped == 0 && run_done > 0 )); then
                num_inversion=$((num_inversion + 1))
                reserved=$(( reserved + $(face_space_needed ${i_face}) ))
                if pgrep -f "inversion_face.sh ${face_name} ${seg_end}" >/dev/null; then
                    # Still running from an earlier run of this script that
                    # was stopped. Never start a second one for the same face.
                    continue
                fi
                if [[ "${failed}" == false ]] && (( num_inv_now < max_parallel_inversions )); then
                    log "${face_name}: all runs reached ${seg_end}; starting its inversion stage"
                    bash inversion_face.sh "${face_name}" "${seg_end}" \
                        > "${test_dir}/logs/inversion_${face_name}.${seg_end}.log" 2>&1 &
                    inversion_pid[i_face]=$!
                    num_inv_now=$((num_inv_now + 1))
                fi
                continue
            fi

            # Started: runs going, or some stopped before the segment end.
            # A face counts as started if this script submitted it, or any
            # of its runs is in the queue. After a restart of this script,
            # a face with stopped runs is simply started again.
            if (( run_active > 0 )) || [[ "${started[i_face]}" == true ]]; then
                num_running=$((num_running + 1))
                reserved=$(( reserved + $(face_space_needed ${i_face}) ))
                if (( run_active > 0 || run_stopped == 0 )); then
                    seen_stopped[i_face]=false
                elif [[ "${seen_stopped[i_face]}" != true ]]; then
                    # A requeued job can briefly be missing from the queue.
                    # Check once more before resubmitting, so a run never
                    # gets two jobs.
                    seen_stopped[i_face]=true
                elif [[ "${resubmitted[i_face]}" == true ]]; then
                    log "ERROR: ${face_name}: ${run_stopped} run(s) did not reach ${seg_end} after one resubmission. See ./check_gchp_imi_jacobian.sh ${seg_end} ${segment_days}"
                    failed=true
                elif [[ "${failed}" == false ]]; then
                    log "${face_name}: resubmitting ${run_stopped} run(s) that failed or did not finish"
                    submit_face "${face_name}" "${seg_end}" true "${active}" || failed=true
                    resubmitted[i_face]=true
                    seen_stopped[i_face]=false
                fi
                continue
            fi

            # Not started yet. Start it only if its space fits, in face order.
            num_waiting=$((num_waiting + 1))
            if [[ "${blocked}" == true || "${failed}" == true ]]; then
                continue
            fi
            need=$(face_space_needed ${i_face})
            if (( fs_used + reserved + need <= limit )); then
                log "${face_name}: starting its Jacobian runs (needs about $(gb ${need}))"
                started[i_face]=true
                # FAILED runs are included: they failed before this script
                # started, e.g. in an earlier attempt at this segment
                if submit_face "${face_name}" "${seg_end}" true "${active}"; then
                    reserved=$(( reserved + need ))
                    num_waiting=$((num_waiting - 1))
                    num_running=$((num_running + 1))
                else
                    failed=true
                fi
            else
                blocked=true
                if (( reserved == 0 )); then
                    log "ERROR: ${face_name} needs about $(gb ${need}), but only $(gb $(( limit - fs_used ))) is free under ${fsx_max_use_percent}% of FSx even with no other face running. Use shorter segments or a larger FSx."
                    failed=true
                fi
            fi
        done

        if (( num_done == num_faces )); then
            return 0
        fi
        if [[ "${failed}" == true ]] && (( num_inv_now == 0 )); then
            return 1
        fi
        log "Segment ${seg_end}: ${num_done} done, ${num_inversion} in inversion stage, ${num_running} running, ${num_waiting} waiting; FSx $(gb ${fs_used}) used, $(gb ${reserved}) reserved, limit $(gb ${limit}); $(job_counts)"
        sleep $((poll_minutes * 60))
    done
}

#--------------------------------------------------------------------------
# Run segments until every face is DONE through the target end
#--------------------------------------------------------------------------
end_face=$(( start_face + num_faces - 1 ))
log "Running ${segment_days}-day segments to ${target_end} for faces ${start_face} to ${end_face}"
if [[ "${prune_output}" != true ]]; then
    log "WARNING: prune_output is not true, so no space is freed after each face. Faces will stop fitting on FSx."
fi
if [[ "${release_met}" == true || "${release_products}" == true || "${release_done_faces}" == true ]] \
        && ! sudo -n true 2>/dev/null; then
    log "WARNING: sudo needs a password, so files cannot be released from FSx. Releasing is turned off."
    release_met=false
    release_products=false
    release_done_faces=false
fi

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

    if ! run_segment "${last_done}" "${seg_end}"; then
        log "Stopped. Fix the problem and run again with the same arguments."
        exit 1
    fi
    log "Every face is DONE through ${seg_end}"

    release_after_segment "${seg_end}"
done

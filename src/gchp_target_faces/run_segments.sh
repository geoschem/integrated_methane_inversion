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
# Segments end every segment_days days from StartDate (and at target_end).
# Each face works on the next segment end after the date its inversion
# stage is DONE through (StartDate if none). Faces move on independently:
# a face that finishes a segment starts its next one without waiting for
# the others, but never more than max_segments_ahead segments ahead of the
# slowest face (set in gchp_imi_common.sh; 0 means all faces move together,
# segment by segment). Faces that start behind, e.g. converted runs at
# different dates, catch up first.
#
# Every poll_minutes it checks each face, the faces furthest behind first:
#   - Not started: submit its Jacobian runs to its segment end if its space
#     fits under fsx_max_use_percent of FSx (see "Storage" in
#     gchp_imi_common.sh) and it is not too far ahead. Once a face does
#     not fit, no face after it is started, so faces behind always get
#     space first.
#   - Runs going: wait. Runs that failed, or stopped before the segment end,
#     are resubmitted once.
#   - All runs at the segment end: run its inversion stage (inversion_face.sh,
#     up to max_parallel_inversions faces at a time), which prunes its
#     Jacobian run output if prune_output is true. That frees its space for
#     other faces. When it is DONE, the face moves on to its next segment.
#   - Waiting (too far ahead, or no space yet): if release_done_faces is
#     true, release its large files that are safely in S3 (restarts, the
#     hourly output kept for the next inversion stage, overpass and .pkl
#     files), so other faces can use the space. Files not exported yet are
#     tried again at the next check. They are loaded back from S3 when read.
# Each time the slowest face moves on, met fields and exported products
# that no face needs any more are released from FSx (they stay in S3), if
# release_met and release_products are true.
#
# If something fails (runs fail again after one resubmission, an inversion
# stage is not DONE, or a face cannot fit even on its own), no more faces
# are started; it waits for running inversion stages to end, then stops.
# Fix the cause and run it again with the same arguments: each face
# continues from the date its inversion stage is DONE through.
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

# How many segments a face may get ahead of the slowest face. Not set: 0,
# all faces move together.
max_segments_ahead=${max_segments_ahead:-0}
if [[ ! "${max_segments_ahead}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: max_segments_ahead '${max_segments_ahead}' in gchp_imi_common.sh must be a whole number, 0 or more"
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

# Function to print the date a face's inversion stage is DONE through, or
# StartDate if it has no DONE line
# Argument: face name
face_done_date() {
    local done_date
    done_date=$(grep "^DONE" "${test_dir}/$1/inversion_status" 2>/dev/null \
                | tail -n 1 | awk '{print $2}')
    echo "${done_date:-${start_date}}"
}

# Function to print the next segment end after a date. Segments end every
# segment_days days from StartDate, and at target_end.
# Argument: date (YYYYMMDD)
next_segment_end() {
    local days end
    days=$(days_between "${start_date}" "$1")
    end=$(date -d "${start_date} +$(( (days / segment_days + 1) * segment_days )) days" +%Y%m%d)
    if [[ "${end}" > "${target_end}" ]]; then
        end=${target_end}
    fi
    echo "${end}"
}

# Function to print the number of a segment end: 1 for the first segment
# end after StartDate, 2 for the next, and so on (target_end counts as the
# segment it ends)
# Argument: segment end (YYYYMMDD)
segment_number() {
    local days
    days=$(days_between "${start_date}" "$1")
    echo $(( (days + segment_days - 1) / segment_days ))
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
# Run until every face is DONE through the target end. Each face has its
# own segment end; the arrays below hold its state for that segment end
# and are reset when the face moves on.
#--------------------------------------------------------------------------
end_face=$(( start_face + num_faces - 1 ))
log "Running ${segment_days}-day segments to ${target_end} for faces ${start_face} to ${end_face}, at most ${max_segments_ahead} segment(s) ahead of the slowest face"
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

declare -a seg_end_of=()      # by face number: segment end it is working on
declare -a done_of=()         # by face number: date its inversion stage is DONE through
declare -a inversion_pid=()   # by face number: PID of its inversion stage
declare -a started=()         # by face number: true once its runs were submitted
declare -a resubmitted=()     # by face number: true once runs were resubmitted
declare -a seen_stopped=()    # by face number: true if runs were stopped at the last check
declare -a face_released=()   # by face number: true once its large files were released
failed=false
released_through=""           # slowest face's DONE date at the last release of met fields
last_positions=""             # how many faces are at each segment end, at the last check

while true; do
    if ! active=$(active_run_dirs); then
        log "squeue failed; checking again in ${poll_minutes} minutes"
        sleep $((poll_minutes * 60))
        continue
    fi
    read -r fs_size fs_used < <(df -B1 --output=size,used "${test_dir}" | tail -n 1)
    limit=$(( fs_size / 100 * fsx_max_use_percent ))

    # Where each face is, and the slowest face
    slowest=""
    for ((i_face=start_face; i_face<=end_face; i_face++)); do
        done_of[i_face]=$(face_done_date "$(face_name_for ${i_face})")
        if [[ -z "${slowest}" || "${done_of[i_face]}" < "${slowest}" ]]; then
            slowest=${done_of[i_face]}
        fi
    done

    # Release met fields and products no face needs any more, each time the
    # slowest face has moved on
    if [[ "${slowest}" != "${released_through}" ]]; then
        if [[ -n "${released_through}" ]]; then
            log "Every face is DONE through ${slowest}"
        fi
        release_after_segment "${slowest}"
        released_through=${slowest}
    fi
    if [[ ! "${slowest}" < "${target_end}" ]]; then
        log "Every face is DONE through ${target_end}. Finished."
        break
    fi

    # Faces may start runs for segment ends up to this number
    slowest_end=$(next_segment_end "${slowest}")
    max_number=$(( $(segment_number "${slowest_end}") + max_segments_ahead ))

    # Inversion stages still running
    num_inv_now=0
    for ((i_face=start_face; i_face<=end_face; i_face++)); do
        if [[ -n "${inversion_pid[i_face]}" ]] && kill -0 "${inversion_pid[i_face]}" 2>/dev/null; then
            num_inv_now=$((num_inv_now + 1))
        fi
    done

    # Faces furthest behind first, then by face number
    order=$(for ((i_face=start_face; i_face<=end_face; i_face++)); do
                echo "$(next_segment_end "${done_of[i_face]}") ${i_face}"
            done | sort -k1,1 -k2,2n | awk '{print $2}')

    reserved=0
    blocked=false
    num_done=0; num_inversion=0; num_running=0; num_waiting=0
    positions=""

    for i_face in ${order}; do
        face_name=$(face_name_for ${i_face})

        # DONE through the target end
        if [[ ! "${done_of[i_face]}" < "${target_end}" ]]; then
            if [[ -n "${inversion_pid[i_face]}" ]]; then
                wait "${inversion_pid[i_face]}" || true
                inversion_pid[i_face]=""
                log "${face_name}: inversion stage DONE through ${done_of[i_face]}"
            fi
            if [[ "${release_done_faces}" == true && "${face_released[i_face]}" != true ]]; then
                if release_face "${face_name}"; then
                    face_released[i_face]=true
                fi
            fi
            num_done=$((num_done + 1))
            positions="${positions}finished"$'\n'
            continue
        fi

        # Its segment end. A new one at the first check, or when its
        # inversion stage has just become DONE.
        seg_end=$(next_segment_end "${done_of[i_face]}")
        positions="${positions}${seg_end}"$'\n'
        if [[ "${seg_end_of[i_face]}" != "${seg_end}" ]]; then
            if [[ -n "${inversion_pid[i_face]}" ]]; then
                wait "${inversion_pid[i_face]}" || true
                inversion_pid[i_face]=""
                log "${face_name}: inversion stage DONE through ${done_of[i_face]}"
            fi
            seg_end_of[i_face]=${seg_end}
            started[i_face]=""; resubmitted[i_face]=""; seen_stopped[i_face]=""
            face_released[i_face]=""
        fi

        # Inversion stage running, or just ended. It may have ended since
        # done_of was read at the start of this check, so read it again.
        if [[ -n "${inversion_pid[i_face]}" ]]; then
            if kill -0 "${inversion_pid[i_face]}" 2>/dev/null; then
                num_inversion=$((num_inversion + 1))
                reserved=$(( reserved + $(face_space_needed ${i_face}) ))
                continue
            fi
            wait "${inversion_pid[i_face]}" || true
            inversion_pid[i_face]=""
            if [[ "$(face_done_date "${face_name}")" < "${seg_end}" ]]; then
                log "ERROR: ${face_name}: inversion stage not DONE. See ${test_dir}/logs/inversion_${face_name}.${seg_end}.log"
                failed=true
            else
                # It moves on to its next segment at the next check
                log "${face_name}: inversion stage DONE through ${seg_end}"
                num_done=$((num_done + 1))
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

        # All runs at the segment end: start its inversion stage. Also when
        # the face is ahead, since it needs no new runs.
        if (( run_active == 0 && run_stopped == 0 && run_done > 0 )); then
            num_inversion=$((num_inversion + 1))
            reserved=$(( reserved + $(face_space_needed ${i_face}) ))
            if pgrep -f "inversion_face.sh ${face_name} ${seg_end}" >/dev/null; then
                # Still running from an earlier run of this script that
                # was stopped. Never start a second one for the same face.
                continue
            fi
            if [[ ! "$(face_done_date "${face_name}")" < "${seg_end}" ]]; then
                # That earlier one has just finished: an inversion stage
                # writes DONE before it exits. The face moves on at the
                # next check.
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

        # Not started. Too far ahead of the slowest face (counted as done),
        # or waiting for space: release its files meanwhile.
        too_far=false
        if (( $(segment_number "${seg_end}") > max_number )); then
            too_far=true
        fi
        if [[ "${too_far}" == true || "${blocked}" == true || "${failed}" == true ]]; then
            if [[ "${too_far}" == true ]]; then
                num_done=$((num_done + 1))
            else
                num_waiting=$((num_waiting + 1))
            fi
            if [[ "${release_done_faces}" == true && "${face_released[i_face]}" != true ]]; then
                if release_face "${face_name}"; then
                    face_released[i_face]=true
                fi
            fi
            continue
        fi

        # Start it if its space fits. Otherwise no face after it is started.
        need=$(face_space_needed ${i_face})
        if (( fs_used + reserved + need <= limit )); then
            log "${face_name}: starting its Jacobian runs to ${seg_end} (needs about $(gb ${need}))"
            started[i_face]=true
            # FAILED runs are included: they failed before this script
            # started, e.g. in an earlier attempt at this segment
            if submit_face "${face_name}" "${seg_end}" true "${active}"; then
                reserved=$(( reserved + need ))
                num_running=$((num_running + 1))
            else
                failed=true
            fi
        else
            blocked=true
            num_waiting=$((num_waiting + 1))
            if (( reserved == 0 )); then
                log "ERROR: ${face_name} needs about $(gb ${need}), but only $(gb $(( limit - fs_used ))) is free under ${fsx_max_use_percent}% of FSx even with no other face running. Use shorter segments or a larger FSx."
                failed=true
            fi
            if [[ "${release_done_faces}" == true && "${face_released[i_face]}" != true ]]; then
                if release_face "${face_name}"; then
                    face_released[i_face]=true
                fi
            fi
        fi
    done

    if [[ "${failed}" == true ]] && (( num_inv_now == 0 )); then
        log "Stopped. Fix the problem and run again with the same arguments."
        exit 1
    fi

    # How many faces are at each segment end, when that changes
    positions=$(sort <<< "${positions}" | uniq -c | awk 'NF == 2 { printf "%s%s %d", sep, $2, $1; sep = ", " }')
    if [[ "${positions}" != "${last_positions}" ]]; then
        log "Faces by segment end: ${positions}"
        last_positions=${positions}
    fi
    log "Segment ${slowest_end}: ${num_done} done, ${num_inversion} in inversion stage, ${num_running} running, ${num_waiting} waiting; FSx $(gb ${fs_used}) used, $(gb ${reserved}) reserved, limit $(gb ${limit}); $(job_counts)"
    sleep $((poll_minutes * 60))
done

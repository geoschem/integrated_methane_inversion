#!/bin/bash

# One-time conversion of target face directories set up and run on another
# cluster with the IMI's standard GCHP run script (gchp_ch4_run.template),
# e.g. Dandan's runs, so they can be continued here with these scripts and
# the spot run script (gchp_ch4_run_spot.template).
#
# Usage:
#   ./convert_rundirs.sh            # dry run: print what would change
#   ./convert_rundirs.sh --apply    # make the changes
#
# Works on the faces set in gchp_imi_common.sh (test_dir, config_dir,
# config_name_style, start_face, num_faces). Copy the face directories into
# test_dir first; never convert the only copy. The paths of the other
# cluster are set below under "Paths on the other cluster".
#
# For each face, this:
#   1. Checks that it is safe to convert, and changes nothing if any check
#      fails (see "Before running" below)
#   2. Turns files that hold a link target back into links. Links copied
#      through S3 can come back as small files containing the target path.
#   3. Points links and paths in HEMCO_Config.rc, rundir_vars.txt,
#      submit_jacobian_simulations_array.sh and the face config to this
#      cluster's paths
#   4. Makes the face config in config_dir from the base config
#      (template_config), so it has every key the current IMI and these
#      scripts need and this cluster's settings. Only the face-specific
#      values (face_keys below) are taken from the face's own config, i.e.
#      the copy the IMI saved in the face directory. A face config in
#      config_dir that already has every key of the base config is kept.
#      OutputPath is set to test_dir and UseSpotRunScript to true.
#   5. In each Jacobian run directory:
#      a. Renames dated checkpoints
#           Restarts/gcchem_internal_checkpoint.YYYYMMDD_HHMMz.nc4
#         to standard restart names
#           Restarts/GEOSChem.Restart.YYYYMMDD_HHMMz.c<N>.nc4
#      b. Sets cap_restart to the date of the newest restart file
#      c. Keeps the old run script as <run name>.run.old and replaces it
#         with the spot run script
#      d. Links gchp.env and removeLocalRestarts.sh
#      e. Appends a CONVERTED line to job_status
#   6. Records how far its inversion stage got, from its inversion and
#      overpass markers, as a DONE line in inversion_status, so
#      run_segments.sh continues from there
#   7. Reports links that still point to nothing, and dates of hourly
#      output missing that the next inversion stage will need
#
# Before running with --apply:
#   - Cancel all jobs for these run directories, and stop anything else
#     that submits jobs or renames checkpoints in them (e.g. cron jobs).
#     This script stops if a job is in the Slurm queue or crontab still
#     runs rename_final_checkpoints or submit_lagging_jacobian_runs.
#   - Rename undated final checkpoints (Restarts/gcchem_internal_checkpoint),
#     e.g. with rename_final_checkpoints.sh --apply, which dates them from
#     the run logs. This script stops if any are left.
#   - Check in the dry run that StartDate and EndDate in the face configs
#     are the full inversion period. They must match the base config, since
#     the face configs are made from it.
#
# Faces whose runs are all converted (CONVERTED in job_status) are skipped,
# so this can be rerun.

set -e
shopt -s nullglob
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

#--------------------------------------------------------------------------
# Paths on the other cluster, and what they are on this cluster
#--------------------------------------------------------------------------
old_output_dir=/fsx_output/output     # its OutputPath  -> test_dir
old_imi_dir=/fsx_output/imi-gchp      # its IMI clone   -> imi_dir
old_input_dir=/fsx_input              # its input data  -> new_input_dir
new_input_dir="${HOME}/fsx/dzhang-imi-gchp"

# Text files that contain paths and are read by GCHP, HEMCO or the IMI
path_files="HEMCO_Config.rc rundir_vars.txt submit_jacobian_simulations_array.sh"

# Config keys taken from the face's own config. Everything else comes from
# the base config.
face_keys="RunName TARGET_LON TARGET_LAT StateVectorFile NumJacobianTracers CreateAutomaticRectilinearStateVectorFile"

if [[ $# -eq 1 && "$1" == "--apply" ]]; then
    apply=true
elif [[ $# -eq 0 ]]; then
    apply=false
else
    echo "Usage: $0 [--apply]"
    exit 1
fi

spot_template="${imi_dir}/src/geoschem_run_scripts/gchp_ch4_run_spot.template"

# Function to test whether changes are being made (--apply) or not
applying() {
    [[ "${apply}" == true ]]
}

# Function to print a path with the other cluster's locations replaced by
# this cluster's
# Argument: path
map_path() {
    local p=$1
    p=${p/#${old_output_dir}/${test_dir}}
    p=${p/#${old_imi_dir}/${imi_dir}}
    p=${p/#${old_input_dir}/${new_input_dir}}
    echo "${p}"
}

# Function to replace the other cluster's paths with this cluster's in a
# text file
# Argument: file
map_paths_in_file() {
    sed -i -e "s|${old_output_dir}|${test_dir}|g" \
           -e "s|${old_imi_dir}|${imi_dir}|g" \
           -e "s|${old_input_dir}|${new_input_dir}|g" "$1"
}

# Function to test whether a file holds a link target: a small one-line
# text file whose content starts with /
# Argument: file
is_link_file() {
    [[ -f "$1" && ! -L "$1" ]] || return 1
    (( $(wc -c < "$1") < 4096 )) || return 1
    (( $(wc -l < "$1") <= 1 )) || return 1
    grep -Iq . "$1" || return 1
    [[ "$(head -c 1 "$1")" == / ]]
}

# Function to list the files of a face that hold a link target, skipping
# model output and inversion products
# Argument: face directory
link_files_of() {
    local file
    while IFS= read -r file; do
        if is_link_file "${file}"; then
            echo "${file}"
        fi
    done < <(find "$1" -type f -size -4k -not -path '*/OutputDir/*' \
             -not -path '*/OverpassDiagnostics/*' -not -path '*/inversion/*')
}

# Function to print the segment end the face's inversion stage completed,
# from its inversion and overpass markers (the earlier of the two), or
# nothing if it has none
# Argument: face directory
marker_date() {
    local inv_files=( "$1"/inversion_complete.*_S* )
    local ovp_files=( "$1"/overpass_complete.*_S* )
    local inv ovp
    inv=${inv_files[0]##*_S}
    ovp=${ovp_files[0]##*_S}
    if [[ -n "${inv}" && -n "${ovp}" && "${ovp}" < "${inv}" ]]; then
        echo "${ovp}"
    else
        echo "${inv}"
    fi
}

# Function to test whether a config has every key the base config has
# Argument: config
has_all_keys() {
    [[ -z "$(comm -23 <(grep -oE '^[A-Za-z0-9_]+:' "${template_config}" | sort -u) \
                      <(grep -oE '^[A-Za-z0-9_]+:' "$1" | sort -u))" ]]
}

# Function to make a face config from the base config, with the
# face-specific values (face_keys) from the face's own config
# Arguments: face's own config, new config (may be the same file)
make_face_config_from() {
    local saved key value
    saved=$(mktemp)
    cp "$1" "${saved}"
    cp "${template_config}" "$2"
    for key in ${face_keys}; do
        value=$(grep "^${key}:" "${saved}" | head -n 1 | cut -d: -f2- | sed 's/^ *//')
        set_config "$2" "${key}" "${value}"
    done
    rm -f "${saved}"
    set_config "$2" OutputPath "\"${test_dir}\""
    map_paths_in_file "$2"
}

# Function to print the newest restart date and time (YYYYMMDD_HHMM) in a
# Restarts directory, under either restart name, or nothing if there is none
# Argument: Restarts directory
newest_restart() {
    local files=( "$1"/gcchem_internal_checkpoint.*z.nc4 "$1"/GEOSChem.Restart.*z.c*.nc4 )
    if (( ${#files[@]} > 0 )); then
        printf '%s\n' "${files[@]##*/}" | grep -oE '[0-9]{8}_[0-9]{4}' | sort | tail -n 1
    fi
}

# Function to print the dates from one date to the day before another,
# one per line
# Arguments: first date, end date (YYYYMMDD, not included)
dates_before() {
    local d=$1
    while [[ "${d}" < "$2" ]]; do
        echo "${d}"
        d=$(date -d "${d} +1 days" +%Y%m%d)
    done
}

if applying; then
    echo "Converting face directories in ${test_dir}"
else
    echo "DRY RUN: showing what would change in ${test_dir}. Nothing is changed."
fi

#--------------------------------------------------------------------------
# Pass 1: check that it is safe to convert, before changing anything
#--------------------------------------------------------------------------
problems=0

# Function to report a problem that prevents the conversion
# Argument: message
problem() {
    echo "  PROBLEM: $1"
    problems=$((problems + 1))
}

if crontab -l 2>/dev/null | grep -v "^ *#" | \
        grep -qE "rename_final_checkpoints|submit_lagging_jacobian_runs"; then
    problem "crontab still runs rename_final_checkpoints or submit_lagging_jacobian_runs. Remove those entries first (crontab -e)."
fi
if [[ ! -f "${spot_template}" ]]; then
    problem "Spot run script ${spot_template} does not exist"
fi
if [[ ! -d "${new_input_dir}" ]]; then
    problem "Input data directory ${new_input_dir} does not exist"
fi

# Config paths are relative to imi_dir
cd "${imi_dir}"
active=$(active_run_dirs)
periods=""
end_face=$(( start_face + num_faces - 1 ))
for ((i_face=start_face; i_face<=end_face; i_face++)); do
    face_name=$(face_name_for ${i_face})
    face_dir="${test_dir}/${face_name}"
    config=$(config_for "${face_name}")
    face_copy="${face_dir}/config_${face_name}.yml"

    if [[ ! -d "${face_dir}/jacobian_runs" ]]; then
        problem "${face_dir}/jacobian_runs does not exist"
        continue
    fi
    if [[ -f "${config}" ]]; then
        source_config=${config}
    elif [[ -f "${face_copy}" ]]; then
        source_config=${face_copy}
    else
        problem "No config: neither ${imi_dir}/${config} nor ${face_copy} exists"
        continue
    fi

    output_path=$(grep "^OutputPath:" "${source_config}" | awk '{print $2}' | tr -d '"')
    if [[ "$(map_path "${output_path%/}")" != "${test_dir}" ]]; then
        problem "OutputPath in ${source_config} is ${output_path}, which is neither ${test_dir} nor ${old_output_dir}"
    fi
    periods="${periods}$(grep -E "^(StartDate|EndDate):" "${source_config}" | tr '\n' ' ')"$'\n'
    for key in StartDate EndDate; do
        if [[ "$(grep "^${key}:" "${source_config}")" != "$(grep "^${key}:" "${template_config}")" ]]; then
            problem "${key} in ${source_config} differs from the base config ${template_config}. The face config is made from the base config, so they must match."
        fi
    done
    for key in ${face_keys}; do
        if ! grep -q "^${key}:" "${source_config}"; then
            problem "${source_config} has no ${key}"
        fi
    done

    for run_dir in "${face_dir}/jacobian_runs/${face_name}_"*; do
        [[ -d "${run_dir}" ]] || continue
        if grep -qxF "$(cd -P "${run_dir}" && pwd)" <<< "${active}"; then
            problem "$(basename "${run_dir}") has a job in the Slurm queue. Cancel it first."
        fi
        if [[ -e "${run_dir}/Restarts/gcchem_internal_checkpoint" ]]; then
            problem "$(basename "${run_dir}") has an undated checkpoint Restarts/gcchem_internal_checkpoint. Date it first (rename_final_checkpoints.sh --apply)."
        fi
    done
done

echo "StartDate and EndDate in the face configs (must be the full inversion period):"
sort -u <<< "${periods}" | sed '/^$/d; s/^/  /'

if (( problems > 0 )); then
    echo "Found ${problems} problem(s). Nothing was changed."
    exit 1
fi

#--------------------------------------------------------------------------
# Function to convert one Jacobian run directory. Without --apply it only
# prints what it would change.
# Argument: run directory
#--------------------------------------------------------------------------
convert_run_dir() {
    local run_dir=$1
    local run_name file stamp new newest old_cap new_cap N

    run_name=$(basename "${run_dir}")
    cd "${run_dir}"

    if grep -q "^CONVERTED" job_status 2>/dev/null; then
        echo "  ${run_name}: already converted, skipping"
        return
    fi
    echo "  ${run_name}:"

    # Grid resolution, used in restart file names (as in the run script)
    N=$(grep "CS_RES=" setCommonRunSettings.sh | cut -c 8- | xargs)

    # a. Rename dated checkpoints to standard restart names
    for file in Restarts/gcchem_internal_checkpoint.*z.nc4; do
        [[ -e "${file}" ]] || continue
        stamp=${file#Restarts/gcchem_internal_checkpoint.}
        stamp=${stamp%z.nc4}
        new="Restarts/GEOSChem.Restart.${stamp}z.c${N}.nc4"
        if [[ -e "${new}" ]]; then
            echo "    WARNING: ${new} already exists; leaving ${file} as it is"
            continue
        fi
        echo "    rename ${file} -> ${new}"
        if applying; then
            mv "${file}" "${new}"
        fi
    done

    # b. Set cap_restart to the newest restart date, under either name
    newest=$(newest_restart Restarts)
    if [[ -z "${newest}" ]]; then
        echo "    WARNING: no restart files found; cap_restart not changed"
    else
        old_cap=$(cat cap_restart 2>/dev/null || echo "none")
        new_cap="${newest:0:8} ${newest:9:4}00"
        echo "    cap_restart: ${old_cap} -> ${new_cap}"
        if applying; then
            echo "${new_cap}" > cap_restart
        fi
    fi

    # c. Replace the run script, keeping the old one
    echo "    replace ${run_name}.run with the spot run script (old one kept as ${run_name}.run.old)"
    if applying; then
        if [[ -f "${run_name}.run" && ! -e "${run_name}.run.old" ]]; then
            mv "${run_name}.run" "${run_name}.run.old"
        fi
        cp "${spot_template}" "${run_name}.run"
        chmod 755 "${run_name}.run"
    fi

    # d. Link the environment file and the restart cleanup script
    echo "    link gchp.env and removeLocalRestarts.sh"
    if applying; then
        ln -nsf "${env_file}" gchp.env
        ln -nsf "${cleanup_restarts_script}" removeLocalRestarts.sh
    fi

    # e. Record the conversion
    if applying; then
        echo "CONVERTED $(date): cap_restart set to ${new_cap:-unchanged}" >> job_status
    fi
}

#--------------------------------------------------------------------------
# Pass 2: convert each face
#--------------------------------------------------------------------------
num_broken=0
num_missing_days=0
for ((i_face=start_face; i_face<=end_face; i_face++)); do
    face_name=$(face_name_for ${i_face})
    face_dir="${test_dir}/${face_name}"
    status_file="${face_dir}/inversion_status"
    cd "${imi_dir}"
    config=$(config_for "${face_name}")
    echo "${face_name}:"

    # 2. Turn files that hold a link target back into links
    num_links=0
    while IFS= read -r file; do
        [[ -n "${file}" ]] || continue
        target=$(map_path "$(head -n 1 "${file}")")
        if (( num_links < 5 )); then
            echo "  link ${file#${face_dir}/} -> ${target}"
        fi
        num_links=$((num_links + 1))
        if applying; then
            ln -nsf "${target}" "${file}"
        fi
    done <<< "$(link_files_of "${face_dir}")"
    if (( num_links > 5 )); then
        echo "  ... and $((num_links - 5)) more files turned back into links"
    fi

    # Links that already are links but point to the other cluster
    while IFS= read -r file; do
        [[ -n "${file}" ]] || continue
        target=$(readlink "${file}")
        if [[ "$(map_path "${target}")" != "${target}" ]]; then
            echo "  relink ${file#${face_dir}/} -> $(map_path "${target}")"
            if applying; then
                ln -nsf "$(map_path "${target}")" "${file}"
            fi
        fi
    done <<< "$(find "${face_dir}" -type l -not -path '*/OutputDir/*')"

    # 3. Paths in text files
    for name in ${path_files} "config_${face_name}.yml"; do
        while IFS= read -r file; do
            [[ -n "${file}" ]] || continue
            if grep -qE "${old_output_dir}|${old_imi_dir}|${old_input_dir}" "${file}"; then
                echo "  paths in ${file#${face_dir}/}"
                if applying; then
                    map_paths_in_file "${file}"
                fi
            fi
        done <<< "$(find "${face_dir}" -name "${name}" -type f -not -path '*/OutputDir/*')"
    done

    # 4. Face config, made from the base config unless already complete
    if [[ -f "${config}" ]] && has_all_keys "${config}"; then
        echo "  ${config} has every key of the base config; kept"
    else
        if [[ -f "${config}" ]]; then
            face_config=${config}
        else
            face_config="${face_dir}/config_${face_name}.yml"
        fi
        echo "  make ${config} from the base config, with these values from ${face_config}:"
        for key in ${face_keys}; do
            echo "    $(grep "^${key}:" "${face_config}" | head -n 1)"
        done
        if applying; then
            mkdir -p "$(dirname "${config}")"
            make_face_config_from "${face_config}" "${config}"
        fi
    fi
    echo "  set UseSpotRunScript: true in ${config}"
    if applying; then
        set_config "${config}" UseSpotRunScript true
    fi

    # 5. Jacobian run directories
    for run_dir in "${face_dir}/jacobian_runs/${face_name}_"*; do
        [[ -d "${run_dir}" ]] || continue
        convert_run_dir "${run_dir}"
    done

    # 6. How far its inversion stage got
    done_date=$(marker_date "${face_dir}")
    if [[ -z "${done_date}" ]]; then
        echo "  no inversion marker: run_segments.sh will start this face from StartDate"
    elif grep -q "^DONE" "${status_file}" 2>/dev/null; then
        echo "  inversion_status already has a DONE line; not changed"
    else
        echo "  inversion stage done through ${done_date} (from its markers): record in inversion_status"
        if applying; then
            echo "DONE ${done_date} $(date): converted, from inversion and overpass markers" >> "${status_file}"
        fi
    fi

    # 7. Checks. Links still pointing to nothing (after --apply)
    if applying; then
        while IFS= read -r file; do
            [[ -n "${file}" ]] || continue
            echo "  WARNING: link points to nothing: ${file#${face_dir}/} -> $(readlink "${file}")"
            num_broken=$((num_broken + 1))
        done <<< "$(find "${face_dir}" -type l -not -path '*/OutputDir/*' ! -exec test -e {} \; -print)"
    fi

    # Hourly output the next inversion stage will read: from the day before
    # the marker date up to each run's newest restart
    if [[ -n "${done_date}" ]]; then
        first_needed=$(date -d "${done_date} -1 days" +%Y%m%d)
        for run_dir in "${face_dir}/jacobian_runs/${face_name}_"*; do
            [[ -d "${run_dir}" ]] || continue
            newest=$(newest_restart "${run_dir}/Restarts" | cut -c 1-8)
            [[ -n "${newest}" ]] || continue
            missing=""
            for d in $(dates_before "${first_needed}" "${newest}"); do
                files=( "${run_dir}"/OutputDir/GEOSChem.*SpeciesConc.${d}_0000z.nc4 )
                if (( ${#files[@]} == 0 )); then
                    missing="${missing} ${d}"
                fi
            done
            if [[ -n "${missing}" ]]; then
                echo "  WARNING: $(basename "${run_dir}") has no hourly output for:${missing}"
                num_missing_days=$((num_missing_days + 1))
            fi
        done
    fi
done

echo ""
if (( num_missing_days > 0 )); then
    echo "WARNING: ${num_missing_days} run(s) are missing hourly output the next inversion stage needs. Copy it from the original output, or the inversion stage will fail."
fi
if applying; then
    if (( num_broken > 0 )); then
        echo "WARNING: ${num_broken} link(s) point to nothing. Check the paths at the top of this script."
    fi
    echo "Done. Check the runs with: ./check_gchp_imi_jacobian.sh <segment end> <max days>"
else
    echo "Dry run done. Nothing was changed. Run with --apply to make these changes."
fi

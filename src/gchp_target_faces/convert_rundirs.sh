#!/bin/bash

# One-time conversion of Jacobian run directories that were set up and run
# with the IMI's standard GCHP run script (gchp_ch4_run.template), e.g. by
# Dandan's scripts, so they can be continued with these scripts and the spot
# run script (gchp_ch4_run_spot.template).
#
# Usage:
#   ./convert_rundirs.sh            # dry run: print what would change
#   ./convert_rundirs.sh --apply    # make the changes
#
# For the faces set in gchp_imi_common.sh (test_dir, config_dir,
# config_name_style, start_face, num_faces), this:
#   1. Checks that it is safe to convert (see below), and changes nothing
#      if any check fails
#   2. Sets UseSpotRunScript to true in each face config
#   3. In each Jacobian run directory:
#      a. Renames dated checkpoints
#           Restarts/gcchem_internal_checkpoint.YYYYMMDD_HHMMz.nc4
#         to standard restart names
#           Restarts/GEOSChem.Restart.YYYYMMDD_HHMMz.c<N>.nc4
#      b. Sets cap_restart to the date of the newest restart file
#      c. Keeps the old run script as <run name>.run.old and replaces it
#         with the spot run script
#      d. Links gchp.env and removeLocalRestarts.sh
#      e. Appends a CONVERTED line to job_status
#
# Before running with --apply:
#   - Remove cron jobs that rename checkpoints or submit runs in these
#     directories. This script stops if crontab still runs
#     rename_final_checkpoints or submit_lagging_jacobian_runs.
#   - Cancel all jobs for these run directories. This script stops if any
#     are in the Slurm queue.
#   - Rename undated final checkpoints (Restarts/gcchem_internal_checkpoint),
#     e.g. with rename_final_checkpoints.sh --apply, which dates them from
#     the run logs. This script stops if any are left.
#   - Check in the dry run that StartDate and EndDate in the face configs
#     are the full inversion period.
#
# Run directories already converted (CONVERTED in job_status) are skipped,
# so this can be rerun.

set -e
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

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

if applying; then
    echo "Converting run directories in ${test_dir}"
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

# Config paths are relative to imi_dir
cd "${imi_dir}"
active=$(active_run_dirs)
periods=""
for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    face_dir="${test_dir}/${face_name}"
    config=$(config_for "${face_name}")

    if [[ ! -d "${face_dir}/jacobian_runs" ]]; then
        problem "${face_dir}/jacobian_runs does not exist"
        continue
    fi
    if [[ ! -f "${config}" ]]; then
        problem "Config file ${imi_dir}/${config} does not exist"
        continue
    fi

    output_path=$(grep "^OutputPath:" "${config}" | awk '{print $2}' | tr -d '"')
    if [[ "${output_path%/}" != "${test_dir}" ]]; then
        problem "OutputPath in ${config} is ${output_path}, not ${test_dir}"
    fi
    periods="${periods}$(grep -E "^(StartDate|EndDate):" "${config}" | tr '\n' ' ')"$'\n'

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
    newest=$(ls -1 Restarts/gcchem_internal_checkpoint.*z.nc4 \
                   Restarts/GEOSChem.Restart.*z.c${N}.nc4 2>/dev/null \
             | grep -oE '[0-9]{8}_[0-9]{4}' | sort | tail -n 1)
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
for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    cd "${imi_dir}"
    config=$(config_for "${face_name}")
    echo "${face_name}:"

    # The key may be missing in configs made before it existed
    echo "  set UseSpotRunScript: true in ${config}"
    if applying; then
        if grep -q "^UseSpotRunScript:" "${config}"; then
            set_config "${config}" UseSpotRunScript true
        else
            echo "UseSpotRunScript: true" >> "${config}"
        fi
    fi

    for run_dir in "${test_dir}/${face_name}/jacobian_runs/${face_name}_"*; do
        [[ -d "${run_dir}" ]] || continue
        convert_run_dir "${run_dir}"
    done
done

echo ""
if applying; then
    echo "Done. Check the runs with: ./check_gchp_imi_jacobian.sh <segment end> <max days>"
else
    echo "Dry run done. Nothing was changed. Run with --apply to make these changes."
fi

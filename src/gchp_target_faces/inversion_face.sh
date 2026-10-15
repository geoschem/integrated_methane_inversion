#!/bin/bash

# Run the IMI inversion stage for one target face. With DisableRun0000 and
# ArchiveJacobiansOnly true, this builds and archives the Jacobian from the
# Jacobian run output. No posterior is calculated.
#
# Run by run_gchp_imi_inversion.sh on the head node. Usage:
#   bash inversion_face.sh <face_name>
#
# run_imi.sh submits the processing job itself (sbatch -W), using
# InvSchedulerPartition, InversionCPUs, InversionMemory and InversionTime
# from the IMI config, and waits for it to finish.
#
# Progress is written to file <face dir>/inversion_status. The last line
# is DONE on success, or FAILED on error.

set -e
source "$(dirname "$0")/gchp_imi_common.sh"

face_name=$1
face_dir="${test_dir}/${face_name}"
config="${config_dir}/config_${face_name}.yml"   # relative to imi_dir
status_file="${face_dir}/inversion_status"

echo "STARTED $(date) on $(hostname)" > "${status_file}"
trap 'echo "FAILED $(date) at line ${LINENO} of inversion_face.sh" >> "${status_file}"' ERR

# run_imi.sh must be run from the IMI directory, so go there now
cd "${imi_dir}"

#--------------------------------------------------------------------------
# Step 1: Check that every Jacobian run reached the end date. Otherwise
# the IMI would quietly use the earliest end date of all runs.
#--------------------------------------------------------------------------
end_date=$(grep "^EndDate:" "${config}" | awk '{print $2}' | tr -d '"')
active=$(active_run_dirs)
for run_dir in "${face_dir}/jacobian_runs/${face_name}_"*; do
    [[ -d "${run_dir}" ]] || continue
    status=$(run_status "${run_dir}" "${end_date}" "${active}")
    if [[ "${status}" != DONE ]]; then
        echo "ERROR: $(basename "${run_dir}") is ${status}, not DONE for end date ${end_date}"
        echo "FAILED $(date): Jacobian runs not finished" >> "${status_file}"
        exit 1
    fi
done

#--------------------------------------------------------------------------
# Step 2: Set up the inversion directory and run the inversion
#--------------------------------------------------------------------------
echo "inversion_face.sh: Running inversion stage for ${face_name}"
set_config "${config}" SetupTemplateRundir       false
set_config "${config}" SetupHemcoPriorEmisRunDir false
set_config "${config}" SetupJacobianRuns         false
set_config "${config}" DoSatelliteDiagnostics    true
set_config "${config}" SetupInversion            true
set_config "${config}" DoInversion               true
set_config "${config}" ArchiveJacobiansOnly      true
set_config "${config}" PruneOutputDir            "${prune_output}"
./run_imi.sh "${config}"

echo "DONE $(date)" >> "${status_file}"

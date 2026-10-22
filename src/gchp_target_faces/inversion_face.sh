#!/bin/bash

# Run the IMI processing stages for one target face, for one segment:
#   1. Satellite diagnostics (overpass files)
#   2. Inversion stage. With DisableRun0000 and ArchiveJacobiansOnly true,
#      this builds the Jacobian (.pkl files) from the Jacobian run output.
#      No posterior is calculated.
#   3. If prune_output is true in gchp_imi_common.sh, delete Jacobian run
#      output that both stages are done with
#   4. Checks that the outputs are complete; DONE is written only if so
#
# Run by run_gchp_imi_inversion.sh on the head node. Usage:
#   bash inversion_face.sh <face_name> <segment_end YYYYMMDD>
#
# All Jacobian runs of the face must have reached segment_end.
#
# IMPORTANT: StartDate and EndDate in the face config must stay the full
# inversion period for every segment. Do not set EndDate to the segment
# end. The IMI works out how far to process from the restart files of the
# Jacobian runs, and never remakes a .pkl file it has already made. A .pkl
# made with a shorter EndDate would miss observations near the segment end.
#
# The last segment should end one day after EndDate. Otherwise the overpass
# diagnostics for the last day are missing values (NaN) for the western
# hemisphere, where 13:30 local time is already the next UTC day.
#
# run_imi.sh submits each processing job itself (sbatch -W), using
# InvSchedulerPartition, InversionCPUs, InversionMemory and InversionTime
# from the IMI config, and waits for it to finish.
#
# Progress is appended to file <face dir>/inversion_status, one line per
# event. The last line is "DONE <segment_end> ..." on success, or
# "FAILED <segment_end> ..." on error.

set -e
source "$(dirname "$0")/gchp_imi_common.sh"

if [[ $# -ne 2 ]]; then
    echo "Usage: bash inversion_face.sh <face_name> <segment_end YYYYMMDD>"
    exit 1
fi
face_name=$1
segment_end=$2
face_dir="${test_dir}/${face_name}"
config=$(config_for "${face_name}")   # relative to imi_dir
status_file="${face_dir}/inversion_status"

if [[ ! -d "${face_dir}" ]]; then
    echo "ERROR: Face directory ${face_dir} does not exist"
    exit 1
fi

echo "STARTED ${segment_end} $(date) on $(hostname)" >> "${status_file}"
trap 'echo "FAILED ${segment_end} $(date) at line ${LINENO} of inversion_face.sh" >> "${status_file}"' ERR

# Function to print an error, record it in inversion_status, and stop
# Argument: error message
fail() {
    echo "ERROR: $1"
    echo "FAILED ${segment_end} $(date): $1" >> "${status_file}"
    exit 1
}

# run_imi.sh must be run from the IMI directory, so go there now
cd "${imi_dir}"

#--------------------------------------------------------------------------
# Step 1: Check the segment end, and that every Jacobian run reached it.
# Otherwise the IMI would quietly use the earliest date of all runs.
# The segment end may be up to one day after EndDate in the config, so the
# overpass diagnostics for the last day have the next day's model output.
#--------------------------------------------------------------------------
if [[ ! -f "${config}" ]]; then
    fail "Config file ${imi_dir}/${config} does not exist"
fi
end_date=$(grep "^EndDate:" "${config}" | awk '{print $2}' | tr -d '"')
if [[ ! "${segment_end}" =~ ^[0-9]{8}$ ]] || ! date -d "${segment_end}" >/dev/null 2>&1; then
    fail "Segment end '${segment_end}' is not a valid date YYYYMMDD"
elif [[ "${segment_end}" > "$(date -d "${end_date} +1 day" +%Y%m%d)" ]]; then
    fail "Segment end ${segment_end} is more than one day after EndDate ${end_date} in ${config}"
fi

# Do not allow segment end to be final EndDate since it would cause the last day
# to be marked with NaN in it, and later segment runs would not fix it.
if [[ "${segment_end}" == "${end_date}" ]]; then
    fail "Segment end equals EndDate ${end_date}. End the last segment one day later, $(date -d "${end_date} +1 day" +%Y%m%d), so the last day's overpass diagnostics are complete."
fi

active=$(active_run_dirs)
for run_dir in "${face_dir}/jacobian_runs/${face_name}_"*; do
    [[ -d "${run_dir}" ]] || continue
    status=$(run_status "${run_dir}" "${segment_end}" "${active}")
    if [[ "${status}" != DONE ]]; then
        fail "$(basename "${run_dir}") is ${status}, not DONE for segment end ${segment_end}"
    fi
done

#--------------------------------------------------------------------------
# Step 2: Run satellite diagnostics and the inversion stage, then prune
# if requested
#--------------------------------------------------------------------------
echo "inversion_face.sh: Processing ${face_name} through segment end ${segment_end}"
set_config "${config}" SetupTemplateRundir       false
set_config "${config}" SetupHemcoPriorEmisRunDir false
set_config "${config}" SetupJacobianRuns         false
set_config "${config}" DoSatelliteDiagnostics    true
set_config "${config}" SetupInversion            true
set_config "${config}" DoInversion               true
set_config "${config}" InvSchedulerPartition     "\"${inversion_partition}\""
set_config "${config}" ArchiveJacobiansOnly      true
set_config "${config}" PruneOutputDir            "${prune_output}"
./run_imi.sh "${config}"

#--------------------------------------------------------------------------
# Step 3: Check that the outputs are complete before recording DONE in the
# status file. This runs prune_outputdir.py without --execute and so deletes
# nothing. It checks that both stages' markers agree, every run has every
# overpass file, and every .pkl in the manifest exists, and exits with
# code 2 if not.
#--------------------------------------------------------------------------
echo "inversion_face.sh: Checking outputs"
conda_env=$(grep "^CondaEnv:" "${config}" | awk '{print $2}' | tr -d '"')
if ! conda run -n "${conda_env}" python src/utilities/prune_outputdir.py "${config}"; then
    fail "Output check failed for segment end ${segment_end}. See the log above."
fi

echo "DONE ${segment_end} $(date)" >> "${status_file}"

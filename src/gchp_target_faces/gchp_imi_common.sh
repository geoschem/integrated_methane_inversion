#!/bin/bash
# Settings and helper functions shared by the GCHP IMI setup scripts.

# *** Edit settings here only ***

#--------------------------------------------------------------------------
# Test settings
#--------------------------------------------------------------------------
test_name=dzhang_s3_output_cleanup
target_prefix=Global_1yr_2025_C36S10 # Do not change
start_face=1
num_faces=3       # For all faces use 220
max_num_tracers=200 # Do not change!

#--------------------------------------------------------------------------
# Slurm resources
#--------------------------------------------------------------------------

# For target face dir setup (except first face which will be done locally)
# IMPORTANT: use an on-demand partition, not spot, to prevent termination.
# The IMI does not allow setting up again if the directory already exists.
setup_partition=ondem4cpu # Put partition name here
face_cores=4
face_time=0-0:30

# For Jacobian runs. These override the #SBATCH lines in the run script
# template. The number of cores is set by TOTAL_CORES in the IMI config.
run_partition=spot-24or48
run_time=0-2:00

# How often run_segments.sh checks on the runs, in minutes
poll_minutes=10

#--------------------------------------------------------------------------
# Inversion stage. Each face runs run_imi.sh on the head node, which
# submits the processing job with the Inv* settings in the IMI config and
# waits for it.
#--------------------------------------------------------------------------
max_parallel_inversions=220

# Partition for satellite diagnostics and the inversion stage. Its nodes
# must have more memory than InversionMemory in the IMI config.
inversion_partition=ondem16

# Delete Jacobian run output already used by satellite diagnostics and the
# inversion stage (true or false). Keep false until pruning has been checked.
prune_output=true


#--------------------------------------------------------------------------
# Storage on FSx, for run_segments.sh. Within a segment, faces are started
# only while the space they will need fits under fsx_max_use_percent of FSx.
#--------------------------------------------------------------------------
fsx_max_use_percent=90         # Ultimately want to use 90
gb_per_element_day=0.055       # hourly Jacobian run output
gb_per_element_restart=0.0048  # one restart file
storage_safety_factor=1.1

# Release files from FSx that are safely in S3, to free space. They are
# loaded back from S3 if anything reads them. Needs sudo without a password.
release_met=true         # met fields every face is done with
release_products=false   # overpass and .pkl files, once exported to S3
release_done_faces=true  # restarts, hourly output and products of faces done with a segment
met_dir="${HOME}/fsx/dzhang-imi-gchp/ExtData/GEOS_0.25x0.3125/GEOS_FP"

#--------------------------------------------------------------------------
# Directories and files
#--------------------------------------------------------------------------
# Directory with the target face directories. 'output' must be an export DRA
# with new/changed/delete enabled
test_dir="${HOME}/fsx/dzhang-imi-gchp-output/output"

util_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
imi_dir="$(cd "${util_dir}/../.." && pwd)"

template_config="${imi_dir}/config-gchp-imi.yml"   # base config, edit by hand

# Nodes and cores for each Jacobian run, from the base IMI config
run_nodes=$(grep "^NUM_NODES:" "${template_config}" | awk '{print $2}')
run_cores=$(grep "^TOTAL_CORES:" "${template_config}" | awk '{print $2}')

# Face configs, in this directory relative to imi_dir. Configs made by
# setup_gchp_imi_jacobian.sh are named config_<face name>.yml ("full").
# Use "short" for names like config_T001.yml.
config_dir="configs_${test_name}"
config_name_style=full

mapping_file="${util_dir}/target_coords_num-sv.csv"

# Files linked into each template_run
env_file="${imi_dir}/gchp.env"
cleanup_restarts_script="${util_dir}/removeLocalRestarts.sh"

# Inputs downloaded from S3
dandan_s3_out="s3://dzhang-imi-gchp-output/output"
prior_emis_pattern="GEOSChem.Emissions.202501*"   # must cover StartDate to EndDate

# Repository branches to use (you must manually check these out!)
imi_branch=feature/gchp_14.7.1+bug_fixes_and_improvements
gchp_branch=feature/improved_14.7.1_for_gchp_imi

#--------------------------------------------------------------------------
# Helper functions
#--------------------------------------------------------------------------

# Return the directory name for a target face number, e.g. 7 -> ..._T007
face_name_for() {
    printf "%s_T%03d" "${target_prefix}" "$1"
}

# Return the number of state vector elements of a face, from the mapping file
# Argument: face number
num_sv_for() {
    sed -n "$((1 + $1))p" "${mapping_file}" | cut -d, -f3
}

# Return the config file for a face, relative to imi_dir
# Argument: face name, e.g. Global_1yr_2025_C36S10_T001
config_for() {
    if [[ "${config_name_style}" == short ]]; then
        echo "${config_dir}/config_${1##*_}.yml"
    else
        echo "${config_dir}/config_${1}.yml"
    fi
}

# Set one key in an IMI config file, e.g. change DoJacobian to true.
# Stops with an error if the key is not in the file.
# Arguments: config file, key, value
set_config() {
    if ! grep -q "^$2:" "$1"; then
        echo "ERROR: key $2 not found in $1"
        return 1
    fi
    sed -i "s|^$2:.*|$2: $3|" "$1"
}

# List the directories of all of this user's jobs in the Slurm queue
# (pending, running or requeued), one per line. Returns 1 if squeue fails,
# so that a failed query is never mistaken for an empty queue.
active_run_dirs() {
    local queue
    queue=$(squeue -h -u "${USER}" -o "%.400Z") || return 1
    awk '{print $1}' <<< "${queue}"
}

# Print the status of one Jacobian run directory, read from file job_status
#   ACTIVE          a job for this directory is in the Slurm queue
#   NO_CAP_RESTART  file cap_restart is missing
#   DONE            cap_restart has reached the end date
#   FAILED          the last job wrote FAILED to job_status
#   READY           can be submitted
# Arguments: run directory, end date YYYYMMDD, output of active_run_dirs
run_status() {
    local run_dir=$1
    local end_date=$2
    local active=$3
    local start_date last_result

    # Slurm stores the physical path, so compare with that
    if grep -qxF "$(cd -P "${run_dir}" && pwd)" <<< "${active}"; then
        echo ACTIVE
        return
    fi
    if [[ ! -f "${run_dir}/cap_restart" ]]; then
        echo NO_CAP_RESTART
        return
    fi
    start_date=$(cut -c 1-8 "${run_dir}/cap_restart")
    if [[ ! "${start_date}" < "${end_date}" ]]; then
        echo DONE
        return
    fi
    last_result=$(grep -E "^(DONE|FAILED)" "${run_dir}/job_status" 2>/dev/null | tail -n 1)
    if [[ "${last_result}" == FAILED* ]]; then
        echo FAILED
    else
        echo READY
    fi
}

# Submit one Jacobian run to an end date, from its own directory, and
# print the sbatch message
# Arguments: run directory, end date YYYYMMDD, maximum days
submit_run() {
    local run_dir=$1
    local end_date=$2
    local max_days=$3
    local run_name
    run_name=$(basename "${run_dir}")
    (cd "${run_dir}" && sbatch -p "${run_partition}" -t "${run_time}" \
        -N "${run_nodes}" -n "${run_cores}" \
        -J "${run_name}" "${run_name}.run" "${end_date}" "${max_days}")
}

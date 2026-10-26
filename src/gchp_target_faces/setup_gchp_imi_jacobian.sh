#!/bin/bash

# This script sets up IMI directories for all target faces. Steps are:
#   1. Makes one IMI config per face from the IMI config file in the IMI repo
#   2. Sets up the first face only. This compiles GCHP to get executables for varying
#      number of Jacobian tracres.
#   3. Checks that every GCHP executable exists, up to the max number of tracers
#   4. Submits setup jobs for the remaining faces
#
# Configurable settings, such as directory paths, are in seperate file gchp_imi_common.sh.
# Always on / always off IMI switches (DoSpinup, DisableRun0000, etc.) must be set
# by hand in the IMI config file prior to running this script. Make sure you check
# those settings in the IMI config file before proceeding, and this script will prompt
# you to do that, printing out the current values for you to review
#
# IMPORTANT: Before running this script: (1) clone the IMI, (2) checkout the correct branch,
# (3) clone GCHP within the IMI, (4) checkout its branches, and (5) update submodules.
# The commands for these steps are as follows:
#   cd ~
#   git clone https://github.com/geoschem/integrated_methane_inversion.git
#   cd integrated_methane_inversion
#   git checkout feature/gchp_14.7.1+bug_fixes_and_improvements
#   git clone -b feature/improved_14.7.1_for_gchp_imi https://github.com/geoschem/GCHP.git
#   cd GCHP
#   git submodule update --init --recursive
#
# Because GCHP is compiled locally in this script, you should run it in a tmux session.
# Here are steps to do that:
#  1. Download tmux if not done already (sudo dnf install tmux -y)
#  2. Open a new session called setup (tmux new -s setup) 
#  3. Run the script (./setup_gchp_imi_jacobian.sh)
#  4. Review the settings and continue, or quit to change config files and rerun
#  5. A message will print with text to tail follow the log. Copy this.
#  6. Exit the tmux session (Control-b d)
#  7. Copy the text to tail follow the log
#
# To see tmux sessions in progress: tmux ls
# To get back into the session    : tmux attach -t setup
# To delete the tmux session      : tmux kill-session -t setup
#
# Note that if you close your ssh connection to the cluster the tmux session
# will keep going and you can go back into it later.


set -e

# If starting from scratch is true, the directory with all target faces will be
# completely replaced from scratch
start_from_scratch=true

# Make sure the other utility scripts are located where this script is run
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

# Never delete directories that were not made for this test, e.g. runs
# and configs set up elsewhere.
if [[ "${start_from_scratch}" == "true" ]] && \
   [[ "$(basename "${test_dir}")" != "${test_name}" || "${config_dir}" != "configs_${test_name}" ]]; then
    echo "ERROR: start_from_scratch would delete ${test_dir} and ${imi_dir}/${config_dir}, which are not named after test ${test_name}. Stopping."
    exit 1
fi

if [[ -z "${setup_partition}" ]]; then
    echo "ERROR: Set slurm partition name in file gchp_imi_common.sh. This will be used to setup directories for all target faces except the first, which will be done on the head node."
    exit 1
fi

#--------------------------------------------------------------------------
# Check that the IMI is on the branch specified in gchp_imi_common.sh
#--------------------------------------------------------------------------
imi_current_branch=$(git -C "${imi_dir}" rev-parse --abbrev-ref HEAD)
if [[ "${imi_current_branch}" != "${imi_branch}" ]]; then
    echo "ERROR: ${imi_dir} is on branch '${imi_current_branch}', not ${imi_branch}."
    exit 1
fi

#--------------------------------------------------------------------------
# Check that GCHP exists and is on branch specified in gchp_imi_common.sh
#--------------------------------------------------------------------------
gchp_dir="${imi_dir}/GCHP"
if [[ ! -d "${gchp_dir}" ]]; then
    echo "ERROR: ${gchp_dir} does not exist. Clone it first:"
    echo "   cd ${imi_dir}"
    echo "   git clone -b ${gchp_branch} https://github.com/geoschem/GCHP.git"
    echo "   cd GCHP && git submodule update --init --recursive"
    exit 1
fi
gchp_current_branch=$(git -C "${gchp_dir}" rev-parse --abbrev-ref HEAD)
if [[ "${gchp_current_branch}" != "${gchp_branch}" ]]; then
    echo "ERROR: ${gchp_dir} is on branch '${gchp_current_branch}', not ${gchp_branch}."
    exit 1
fi

#--------------------------------------------------------------------------
# Function to create the IMI config file for one target face via
# copy and modification of the IMI config file in the IMI repository.
# Arguments: face name, target longitude, target latitude
#--------------------------------------------------------------------------
make_face_config() {
    local face_name=$1
    local lon=$2
    local lat=$3
    local face_config="${imi_dir}/$(config_for "${face_name}")"

    cp "${template_config}" "${face_config}"
    set_config "${face_config}" RunName            "\"${face_name}\""
    set_config "${face_config}" UseSpotRunScript   true
    set_config "${face_config}" TARGET_LON         "${lon}"
    set_config "${face_config}" TARGET_LAT         "${lat}"
    set_config "${face_config}" OutputPath         "\"${test_dir}\""
    set_config "${face_config}" NumJacobianTracers "${max_num_tracers}"
    set_config "${face_config}" CreateAutomaticRectilinearStateVectorFile false
    set_config "${face_config}" StateVectorFile    "\"${test_dir}/${face_name}/StateVector.nc\""
}

#--------------------------------------------------------------------------
# Function to print each key's value from the base config, or NOT FOUND if missing
#--------------------------------------------------------------------------
show_config_keys() {
    local key line
    for key in "$@"; do
        line=$(grep -m 1 "^${key}:" "${template_config}" || true)
        if [[ -n "${line}" ]]; then
            printf "  %-28s%s\n" "${key}:" "${line#*:}"
        else
            printf "  %-28s %s\n" "${key}:" "NOT FOUND"
        fi
    done
}

#--------------------------------------------------------------------------
# Show the user what is in the IMI config and ask them to review it
# Allow them to stop if it is not what they want.
#--------------------------------------------------------------------------
echo "=============================================================================="
echo "Review settings in base IMI config: ${template_config}"
echo "=============================================================================="
echo "Simulation:"
show_config_keys StartDate EndDate Checkpoint_Freq Met \
                 CS_RES STRETCH_GRID STRETCH_FACTOR
echo "Compute resources for GCHP runs:"
show_config_keys TOTAL_CORES NUM_NODES NUM_CORES_PER_NODE
echo "Input data:"
show_config_keys DataPath DataPathTROPOMI BlendedTROPOMI \
                 RestartFilePrefix BCpath BCversion
echo "Environment:"
show_config_keys CondaEnv CondaFile GEOSChemEnv
echo "Switches (stage switches are set by setup_face.sh):"
show_config_keys RunSetup DisableRun0000 KalmanMode \
                 PruneOutputDir SetupJacobianRuns \
                 SetupSpinupRun SetupInversion SetupPosteriorRun \
                 DoHemcoPriorEmis RegridHemcoPriorEmis DoSpinup \
                 DoJacobian ReDoJacobian DoSatelliteDiagnostics \
                 DoInversion DoPosterior DoPreview
echo ""
echo "Set per face by this script (from gchp_imi_common.sh):"
echo "  RunName, TARGET_LON, TARGET_LAT  from ${mapping_file}"
echo "  OutputPath:                 ${test_dir}"
echo "  NumJacobianTracers:         ${max_num_tracers}"
echo "  UseSpotRunScript:           true"
echo "  PruneOutputDir (inversion): ${prune_output}  (prune_output in gchp_imi_common.sh)"
echo "  StateVectorFile:            ${test_dir}/<face>/StateVector.nc"
echo "  Target faces:               ${start_face} to ${num_faces}"
echo ""
echo "Downloads and Slurm:"
echo "  S3 source:                  ${dandan_s3_out}"
echo "  Prior emissions pattern:    ${prior_emis_pattern}"
echo "  Setup partition:            ${setup_partition}"
echo "  IMI branch:                 ${imi_current_branch} ($(git -C "${imi_dir}" rev-parse --short HEAD))"
echo "  GCHP branch:                ${gchp_current_branch} ($(git -C "${gchp_dir}" rev-parse --short HEAD))"
echo "=============================================================================="
if [[ "${start_from_scratch}" == "true" ]]; then
    echo "WARNING: This will DELETE ${test_dir} and its face configs."
    echo "=============================================================================="
fi
read -p "Please review the settings in the IMI config file listed above. Continue? [y/N] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted by user."
    exit 1
fi

echo " "
if [[ "${start_from_scratch}" == "true" ]]; then
    rm -rf "${test_dir}" "${imi_dir}/${config_dir}"
fi
mkdir -p "${test_dir}/logs" "${imi_dir}/${config_dir}"

#--------------------------------------------------------------------------
# Make one config file per target face
#--------------------------------------------------------------------------
for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    file_row=$((1 + i_face))
    csv_line=$(sed -n "${file_row}p" "${mapping_file}")
    IFS=',' read -r target_lon target_lat num_sv <<< "${csv_line}"
    make_face_config "$(face_name_for ${i_face})" "${target_lon}" "${target_lat}"
done
echo "setup_gchp_imi_jacobian.sh: Created face configs in ${imi_dir}/${config_dir}"

#--------------------------------------------------------------------------
# Set up the first target face directory; this step includes one-time
# compilation of GCHP which will take a while
#--------------------------------------------------------------------------
first_face=$(face_name_for ${start_face})
first_log="${test_dir}/logs/setup_${first_face}.log"
echo "setup_gchp_imi_jacobian.sh: Setting up first face ${first_face} on this node (builds GCHP if needed)..."
echo "   Follow progress with: tail -f ${first_log}"
if ! nice bash setup_face.sh "${first_face}" > "${first_log}" 2>&1; then
    echo "ERROR: Setup of ${first_face} failed. See ${first_log}"
    exit 1
fi

#--------------------------------------------------------------------------
# Make sure all of the expected executables are there
#--------------------------------------------------------------------------
bin_dir="${imi_dir}/GCHP/build/bin"
for exe in gchp.default $(seq -f "gchp.%g" 10 10 ${max_num_tracers}); do
    if [[ ! -f "${bin_dir}/${exe}" ]]; then
        echo "ERROR: ${bin_dir}/${exe} is missing. See ${first_log} and ${imi_dir}/GCHP/build/build_geoschem.log"
        exit 1
    fi
done
echo "setup_gchp_imi_jacobian.sh: All GCHP executables found. Submitting the remaining faces."

#--------------------------------------------------------------------------
# Now that GCHP is built, setup directories for the remaining target faces
# by submitting one batch job per target face
#--------------------------------------------------------------------------
for ((i_face=start_face+1; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    sbatch -p "${setup_partition}" \
	   -c "${face_cores}" \
           -t "${face_time}" \
           -o "${test_dir}/logs/setup_${face_name}.log" \
           setup_face.sh "${face_name}"
done
echo "setup_gchp_imi_jacobian.sh: Done. Check progress with: tail -n 1 ${test_dir}/*/setup_status"

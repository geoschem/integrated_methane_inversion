#!/bin/bash
#SBATCH -N 1

# Set up the IMI directories for one target face:
#   1. Download the state vector mapping file from S3
#   2. Create the template and prior emissions run directories
#   3. Download the pre-generated prior emissions from S3
#   4. Create the Jacobian run directories
#
# This script is submitted within script setup_gchp_imi_jacobian.sh. Usage:
#   sbatch setup_face.sh <face_name>
#
# Setup progress is written to file <face dir>/setup_status. The last line
# is DONE on success, or FAILED with the line number on error.

set -e

#--------------------------------------------------------------------------
# Step 0: Preliminary setup
#--------------------------------------------------------------------------

# Slurm runs a copy of this script, so find the common file through the
# directory it was submitted from.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/gchp_imi_common.sh"

face_name=$1
face_dir="${test_dir}/${face_name}"
config=$(config_for "${face_name}")   # relative to imi_dir
status_file="${face_dir}/setup_status"

mkdir -p "${face_dir}"
echo "STARTED $(date) job ${SLURM_JOB_ID} on $(hostname)" > "${status_file}"
trap 'echo "FAILED $(date) at line ${LINENO} of setup_face.sh" >> "${status_file}"' ERR

# run_imi.sh must be run from the IMI directory, so go there now
cd "${imi_dir}"

#--------------------------------------------------------------------------
# Step 1: Download the state vector file
#--------------------------------------------------------------------------
echo "setup_face.sh: Downloading state vector for ${face_name}"
aws s3 cp "${dandan_s3_out}/${face_name}/StateVector.nc" "${face_dir}/StateVector.nc"

#--------------------------------------------------------------------------
# Step 2: Create template and prior emissions run directories
#--------------------------------------------------------------------------
echo "setup_face.sh: Creating template and prior emissions run directories"
set_config "${config}" SetupTemplateRundir       true
set_config "${config}" SetupHemcoPriorEmisRunDir true
set_config "${config}" SetupJacobianRuns         false
./run_imi.sh "${config}"

# Link utility scripts into the template run directory so the Jacobian
# run directories copied from it in step 4 get them too
cd "${face_dir}/template_run"
ln -nsf "${env_file}"
ln -nsf "${cleanup_restarts_script}"
cd "${imi_dir}"

#--------------------------------------------------------------------------
# Step 3: Download prior emissions (needed to create Jacobian run dirs)
#--------------------------------------------------------------------------
echo "setup_face.sh: Downloading prior emissions"
aws s3 sync "${dandan_s3_out}/${face_name}/hemco_prior_emis/OutputDir/" \
    "${face_dir}/hemco_prior_emis/OutputDir/" \
    --exclude "*" --include "${prior_emis_pattern}"
if ! ls "${face_dir}/hemco_prior_emis/OutputDir/"GEOSChem.Emissions.* >/dev/null 2>&1; then
    echo "ERROR: No prior emissions downloaded matching ${prior_emis_pattern}"
    echo "FAILED $(date): no prior emissions downloaded" >> "${status_file}"
    exit 1
fi

#--------------------------------------------------------------------------
# Step 4: Create Jacobian run directories
#--------------------------------------------------------------------------
echo "setup_face.sh: Creating Jacobian run directories"
set_config "${config}" SetupTemplateRundir       false
set_config "${config}" SetupHemcoPriorEmisRunDir false
set_config "${config}" SetupJacobianRuns         true
./run_imi.sh "${config}"

echo "DONE $(date)" >> "${status_file}"

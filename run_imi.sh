#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=2000
#SBATCH --mail-type=END
#SBATCH -o "imi_output.log"
#SBATCH --open-mode=append

## Uncomment to use PBS
###PBS -l nodes=1,ncpus=1
###PBS -o "imi_output.log"

# This script will run the Integrated Methane Inversion (IMI) with GEOS-Chem.
# For documentation, see https://imi.readthedocs.io.
#
# Authors: Daniel Varon, Melissa Sulprizio, Lucas Estrada, Will Downs


## Config log rotation - preserves history across re-launches with the same config
if [[ $# == 1 ]]; then
    _imi_cfg="$1"
else
    _imi_cfg="config.yml"
fi
_imi_tag=$(basename "$_imi_cfg" .yml)
_imi_log="imi_output.log"
_imi_sentinel="imi_output.log.tag"

if [[ -s "$_imi_log" && -f "$_imi_sentinel" ]]; then
    _imi_prev=$(<"$_imi_sentinel")
    if [[ "$_imi_prev" != "$_imi_tag" ]]; then
        _imi_ts=$(date +%Y%m%d_%H%M%S)
        cp "$_imi_log" "imi_output.${_imi_prev}.${_imi_ts}.bak.log"
        # Truncate instead of delete, because slurm has the file open
        : > "$_imi_log"
    fi
fi
# Always copy the current tag into the sentinel file
echo "$_imi_tag" > "$_imi_sentinel"

printf "\n================================================================================\n"
printf "IMI invocation: %s  job=%s  config=%s\n" \
    "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "${SLURM_JOB_ID:-no-slurm}" "$_imi_tag"
printf "================================================================================\n"

unset _imi_cfg _imi_tag _imi_log _imi_sentinel _imi_prev _imi_ts

##=======================================================================
## Import Shell functions
##=======================================================================
source src/utilities/common.sh
source src/components/setup_component/setup.sh
source src/components/template_component/template.sh
source src/components/template_component/compile.sh
source src/components/statevector_component/statevector.sh
source src/components/hemco_prior_emis_component/hemco_prior_emis.sh
source src/components/preview_component/preview.sh
source src/components/spinup_component/spinup.sh
source src/components/jacobian_component/jacobian.sh
source src/components/inversion_component/inversion.sh
source src/components/posterior_component/posterior.sh
source src/components/kalman_component/kalman.sh
source src/components/osse_component/osse_sim.sh

# trap and exit on errors
trap 'imi_failed $LINENO' ERR

start_time=$(date)
setup_start=$(date +%s)

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE (run_imi.sh) ===\n"

# Check if user has specified a configuration file
if [[ $# == 1 ]]; then
    ConfigFile=$1
else
    ConfigFile="config.yml"
fi

# Get the conda environment name and source file
# These variables are sourced manually because
# we need the python environment to parse the yaml file
PythonEnv=$(grep '^PythonEnv:' ${ConfigFile} |
    sed 's/PythonEnv://' |
    sed 's/#.*//' |
    sed 's/^[[:space:]]*//' |
    tr -d '"')

# Load conda/mamba/micromamba and append the current directory to PYTHONPATH
source $PythonEnv
export PYTHONPATH=${PYTHONPATH}:$(pwd -P)

# Parsing the config file
eval $(python src/utilities/parse_yaml.py ${ConfigFile})

if [[ -z "$GEOSChemEnv" ]]; then
    printf "\nWarning: GEOS-Chem environment not specified in config file.\n"
    printf "GEOS-Chem dependencies are assumed to be preloaded\n"
else
    # Load environment for compiling and running GEOS-Chem
    if [ ! -f "${GEOSChemEnv}" ]; then
        printf "\nGEOS-Chem environment file ${GEOSChemEnv} does not exist!"
        printf "\nIMI $RunName Aborted\n"
        exit 1
    else
        printf "\nLoading GEOS-Chem environment: ${GEOSChemEnv}\n"
            source ${GEOSChemEnv}
    fi

    # If scheduler is PBS, get the list of needed sites
    if [[ "$SchedulerType" = "PBS" ]]; then
        convert_sbatch_to_pbs
        sed -i -e "s/SLURM_ARRAY_TASK_ID/PBS_ARRAY_INDEX/g" ${OutputPath}/src/geoschem_run_scripts/run_jacobian_simulations.sh
        sed -i '/^export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK/s/^/# /' ${OutputPath}/src/geoschem_run_scripts/run.template
    fi
fi

# Check all necessary config variables are present
python src/utilities/sanitize_input_yaml.py $ConfigFile || imi_failed

# Set path to IMI runs
RunDirs="${OutputPath}/${RunName}"

##=======================================================================
## Standard settings
##=======================================================================

# In safe mode check whether selected options will overwrite existing files
if "$SafeMode"; then

    # Check if directories exist before creating them
    if ([ -d "${RunDirs}/spinup_run" ] && "$SetupSpinupRun") ||
        ([ -d "${RunDirs}/jacobian_runs" ] && "$SetupJacobianRuns") ||
        ([ -d "${RunDirs}/inversion" ] && "$SetupInversion") ||
        ([ -d "${RunDirs}/posterior_run" ] && "$SetupPosteriorRun"); then

        printf "\nERROR: Run directories in ${RunDirs}/"
        printf "\n   already exist. Please change RunName or change the"
        printf "\n   Setup* options to false in the IMI config file.\n"
        printf "\n  To proceed, and overwrite existing run directories, set"
        printf "\n  SafeMode in the config file to false.\n"
        printf "\nIMI $RunName Aborted\n"
        exit 1
    fi

    # Check if output from previous runs exists
    if ([ -d "${RunDirs}/spinup_run" ] && "$DoSpinup") ||
        ([ -d "${RunDirs}/jacobian_runs" ] && "$DoJacobian") ||
        ([ -d "${RunDirs}/inversion" ] && "$DoInversion") ||
        ([ -d "${RunDirs}/posterior_run/OutputDir/" ] && "$DoPosterior"); then
        printf "\nWARNING: Output files in ${RunDirs}/"
        printf "\n  may be overwritten. Please change RunName in the IMI"
        printf "\n  config file to avoid overwriting files.\n"
        printf "\n  To proceed, and overwrite existing output files, set"
        printf "\n  SafeMode in the config file to false.\n"
        printf "\nIMI $RunName Aborted\n"
        exit 1
    fi
fi

# Path to inversion setup
InversionPath=$(pwd -P)
ConfigPath=${InversionPath}/${ConfigFile}

# add inversion path to python path
export PYTHONPATH=${PYTHONPATH}:${InversionPath}

# Make run directory
mkdir -p -v ${RunDirs}

# Set/Collect information about the GEOS-Chem version, IMI version,
# and TROPOMI processor version
GEOSCHEM_VERSION=14.7.1
IMI_VERSION=$(git describe --tags)
TROPOMI_PROCESSOR_VERSION=$(grep 'VALID_TROPOMI_PROCESSOR_VERSIONS =' src/utilities/download_TROPOMI.py |
    sed 's/VALID_TROPOMI_PROCESSOR_VERSIONS = //' |
    tr -d '"')

# copy config file to run directory and add some run information to it for reference
cp $ConfigFile "${RunDirs}/config_${RunName}.yml"
echo ""
echo "## ================== IMI run information ==================" >>"${RunDirs}/config_${RunName}.yml"
echo "# Run with IMI version: ${IMI_VERSION}" >>"${RunDirs}/config_${RunName}.yml"
echo "# GEOS-Chem version: ${GEOSCHEM_VERSION}" >>"${RunDirs}/config_${RunName}.yml"
echo "# TROPOMI/blended processor version(s): ${TROPOMI_PROCESSOR_VERSION}" >>"${RunDirs}/config_${RunName}.yml"

##=======================================================================
##  Download the TROPOMI data
##=======================================================================
# Download TROPOMI or blended dataset from AWS
satelliteCache=${RunDirs}/satellite_data

if [[ -z "$DataPathObs" ]]; then
    mkdir -p -v $satelliteCache

    if [[ "$SatelliteProduct" == "BlendedTROPOMI" ]]; then
        downloadScript=src/utilities/download_blended_TROPOMI.py
    elif [[ "$SatelliteProduct" == "TROPOMI" ]]; then
        downloadScript=src/utilities/download_TROPOMI.py
    else
        printf "$SatelliteProduct is not currently supported for download"
    fi
    submit_job $SchedulerType true $RequestedMemory $RequestedCPUs $RequestedTime $downloadScript $StartDate $EndDate $satelliteCache
else
    # use existing tropomi data and create a symlink to it
    if [[ ! -L $satelliteCache ]]; then
        ln -s $DataPathObs $satelliteCache
    fi
fi

# Check to make sure there are no duplicate TROPOMI files (e.g., two files with the same orbit number but a different processor version)
python src/utilities/test_TROPOMI_dir.py $satelliteCache

##=======================================================================
##  Run the setup script
##=======================================================================
if "$RunSetup"; then
    echo 'calling setup_imi'
    setup_imi
    echo 'setup_imi done'
else
    echo 'skipping setup_imi'
fi
setup_end=$(date +%s)

##=======================================================================
##  Submit spinup simulation
##=======================================================================
if "$DoSpinup"; then
    echo 'calling run_spinup'
    run_spinup
    echo 'run_spinup done'
else
    echo 'skipping run_spinup'
fi

if ("$DoOSSE" && "$EnableOSSE"); then
    setup_osse
    run_osse
fi

##=======================================================================
##  Run Kalman Filter Mode
##=======================================================================
if "$KalmanMode"; then
    echo 'calling setup_kf'
    setup_kf
    echo 'setup_kf done'
    echo 'calling run_kf'
    run_kf
    echo 'run_kf done'
else
    echo 'skipping run_kf'
fi

##=======================================================================
##  Submit Jacobian simulation
##=======================================================================
if ("$DoJacobian" && ! "$KalmanMode"); then
    echo 'calling run_jacobian'
    run_jacobian
    echo 'run_jacobian done'
else
    echo 'skipping run_jacobian'
fi

##=======================================================================
##  Process data and run inversion
##=======================================================================
if ("$DoInversion" && ! "$KalmanMode"); then
    echo 'calling run_inversion'
    run_inversion
    echo 'run_inversion done'
else
    echo 'skipping run_inversion'
fi

##=======================================================================
##  Submit posterior simulation and process the output
##=======================================================================
if ("$DoPosterior" && ! "$KalmanMode"); then
    echo 'calling run_posterior'
    run_posterior
    echo 'run_posterior done'
else
    echo 'skipping run_posterior'
fi

printf "\n=== DONE RUNNING THE IMI ===\n"

# Run time
end_time=$(date)
printf "\nIMI started : %s" "$start_time"
printf "\nIMI ended   : %s" "$end_time"
printf "\n"

if ! "$KalmanMode"; then
    echo 'calling print_stats'
    print_stats
    echo 'print_stats done'
else
    echo 'skipping print_stats'
fi

# Copy output log to run directory for storage
if [[ -f ${InversionPath}/imi_output.log ]]; then
    cp "${InversionPath}/imi_output.log" "${RunDirs}/imi_output.log"
fi

# copy config file to run directory
cd $InversionPath
cp $ConfigFile "${RunDirs}/config_${RunName}.yml"

# Upload output to S3 if specified
if "$S3Upload"; then
    echo 'calling s3_upload.py to upload output to S3'
    python src/utilities/s3_upload.py $ConfigFile
    echo 's3 upload done'
else
    echo 'skipping s3 upload'
fi

exit 0

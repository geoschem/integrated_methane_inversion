#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=2000
#SBATCH -o "imi_output.log"

# This script will run the Integrated Methane Inversion (IMI) with GEOS-Chem.
# For documentation, see https://imi.readthedocs.io.
#
# Authors: Daniel Varon, Melissa Sulprizio, Lucas Estrada, Will Downs

##=======================================================================
## Import Shell functions
##=======================================================================
source src/utilities/common.sh
source src/components/setup_component/setup.sh
source src/components/template_component/template.sh
source src/components/statevector_component/statevector.sh
source src/components/hemco_prior_emis_component/hemco_prior_emis.sh
source src/components/preview_component/preview.sh
source src/components/spinup_component/spinup.sh
source src/components/jacobian_component/jacobian.sh
source src/components/inversion_component/inversion.sh
source src/components/posterior_component/posterior.sh
source src/components/kalman_component/kalman.sh

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
CondaEnv=$(grep '^CondaEnv:' ${ConfigFile} |
    sed 's/CondaEnv://' |
    sed 's/#.*//' |
    sed 's/^[[:space:]]*//' |
    tr -d '"')
CondaFile=$(eval echo $(grep '^CondaFile:' ${ConfigFile} |
    sed 's/CondaFile://' |
    sed 's/#.*//' |
    sed 's/^[[:space:]]*//' |
    tr -d '"'))

# Load conda/mamba/micromamba e.g. ~/.bashrc
source $CondaFile

# Activate Conda environment
printf "\nActivating conda environment: ${CondaEnv}\n"
conda activate ${CondaEnv}

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
GEOSCHEM_VERSION=14.6.2
IMI_VERSION=$(git describe --tags)
TROPOMI_PROCESSOR_VERSION=$(grep 'VALID_TROPOMI_PROCESSOR_VERSIONS =' src/utilities/download_TROPOMI.py |
    sed 's/VALID_TROPOMI_PROCESSOR_VERSIONS = //' |
    tr -d '"')

# copy config file to run directory and add some run information to it for reference
cp $ConfigFile "${RunDirs}/config_${RunName}.yml"
echo "## ================== IMI run information ==================" >>"${RunDirs}/config_${RunName}.yml"
echo "# Run with IMI version: ${IMI_VERSION}" >>"${RunDirs}/config_${RunName}.yml"
echo "# GEOS-Chem version: ${GEOSCHEM_VERSION}" >>"${RunDirs}/config_${RunName}.yml"
echo "# TROPOMI/blended processor version(s): ${TROPOMI_PROCESSOR_VERSION}" >>"${RunDirs}/config_${RunName}.yml"

##=======================================================================
##  Download the TROPOMI data
##=======================================================================

# Download TROPOMI or blended dataset from AWS
tropomiCache=${RunDirs}/satellite_data

if [[ -z "$DataPathTROPOMI" ]]; then
    mkdir -p -v $tropomiCache

    if "$BlendedTROPOMI"; then
        downloadScript=src/utilities/download_blended_TROPOMI.py
    else
        downloadScript=src/utilities/download_TROPOMI.py
    fi
    sbatch --mem $RequestedMemory \
        -c $RequestedCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -o imi_output.tmp \
        -W $downloadScript $StartDate $EndDate $tropomiCache
    wait
    cat imi_output.tmp >>${InversionPath}/imi_output.log
    rm imi_output.tmp
else
    # use existing tropomi data and create a symlink to it
    if [[ ! -L $tropomiCache ]]; then
        ln -s $DataPathTROPOMI $tropomiCache
    fi
fi

# Check to make sure there are no duplicate TROPOMI files (e.g., two files with the same orbit number but a different processor version)
python src/utilities/test_TROPOMI_dir.py $tropomiCache

##=======================================================================
##  Run the setup script
##=======================================================================
if "$RunSetup"; then
    setup_imi
fi
setup_end=$(date +%s)

##=======================================================================
##  Submit spinup simulation
##=======================================================================
if "$DoSpinup"; then
    run_spinup
fi

##=======================================================================
##  Run Kalman Filter Mode
##=======================================================================
if "$KalmanMode"; then
    setup_kf
    run_kf
fi

##=======================================================================
##  Submit Jacobian simulation
##=======================================================================
if ("$DoJacobian" && ! "$KalmanMode"); then
    run_jacobian
fi

##=======================================================================
##  Process data and run inversion
##=======================================================================
if ("$DoInversion" && ! "$KalmanMode"); then
    run_inversion
fi

##=======================================================================
##  Submit posterior simulation and process the output
##=======================================================================
if ("$DoPosterior" && ! "$KalmanMode"); then
    run_posterior
fi

printf "\n=== DONE RUNNING THE IMI ===\n"

# Run time
end_time=$(date)
printf "\nIMI started : %s" "$start_time"
printf "\nIMI ended   : %s" "$end_time"
printf "\n"

if ! "$KalmanMode"; then
    print_stats
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
    python src/utilities/s3_upload.py $ConfigFile
fi

exit 0

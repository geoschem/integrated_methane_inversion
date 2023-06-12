#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
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
source src/components/preview_component/preview.sh
source src/components/spinup_component/spinup.sh
source src/components/jacobian_component/jacobian.sh
source src/components/inversion_component/inversion.sh
source src/components/posterior_component/posterior.sh

# trap and exit on errors
trap 'imi_failed $LINENO' ERR

start_time=$(date)
setup_start=$(date +%s)

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE (run_imi.sh) ===\n"

# Check if user has specified a configuration file
if [[ $# == 1 ]] ; then
    ConfigFile=$1
else
    ConfigFile="config.yml"
fi

# Get configuration
source src/utilities/parse_yaml.sh
eval $(parse_yaml ${ConfigFile})

if ! "$isAWS"; then
    # Activate Conda environment
    printf "\nActivating conda environment: ${CondaEnv}\n"
    eval "$(conda shell.bash hook)"
    conda activate $CondaEnv
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
    if ([ -d "${RunDirs}/spinup_run" ] && "$SetupSpinupRun") || \
       ([ -d "${RunDirs}/jacobian_runs" ] && "$SetupJacobianRuns") || \
       ([ -d "${RunDirs}/inversion" ] && "$SetupInversion") || \
       ([ -d "${RunDirs}/posterior_run" ] && "$SetupPosteriorRun"); then
        
        printf "\nERROR: Run directories in ${RunDirs}/"
	printf "\n   already exist. Please change RunName or change the"
	printf "\n   Setup* options to false in the IMI config file.\n"
        printf "\nIMI $RunName Aborted\n"
        exit 1 
    fi

    # Check if output from previous runs exists
    if ([ -d "${RunDirs}/spinup_run" ] && "$DoSpinup") || \
       ([ -d "${RunDirs}/jacobian_runs" ] && "$DoJacobian") || \
       ([ -d "${RunDirs}/inversion" ] && "$DoInversion") || \
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

##=======================================================================
##  Download the TROPOMI data
##=======================================================================

# Download TROPOMI data from AWS. You will be charged if your ec2 instance is not in the eu-central-1 region.
mkdir -p -v ${RunDirs}
tropomiCache=${RunDirs}/data_TROPOMI
if "$isAWS"; then
    { # test if instance has access to TROPOMI bucket
        stdout=`aws s3 ls s3://meeo-s5p`
    } || { # catch 
        printf "\nError: Unable to connect to TROPOMI bucket. This is likely caused by misconfiguration of the ec2 instance iam role s3 permissions.\n"
        printf "IMI $RunName Aborted.\n"
        exit 1
    }
    mkdir -p -v $tropomiCache
    printf "Downloading TROPOMI data from S3\n"
    python src/utilities/download_TROPOMI.py $StartDate $EndDate $tropomiCache
    printf "\nFinished TROPOMI download\n"
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
if  "$DoSpinup"; then
    run_spinup
fi

##=======================================================================
##  Submit Jacobian simulation
##=======================================================================
if "$DoJacobian"; then
    run_jacobian
fi

##=======================================================================
##  Process data and run inversion
##=======================================================================
if "$DoInversion"; then
    run_inversion
fi

##=======================================================================
##  Submit posterior simulation and process the output
##=======================================================================
if "$DoPosterior"; then
    run_posterior
fi

printf "\n=== DONE RUNNING THE IMI ===\n"

# Run time
end_time=$(date)
printf "\nIMI started : %s" "$start_time"
printf "\nIMI ended   : %s" "$end_time"
printf "\n"
printf "\nRuntime statistics (s):"
printf "\n Setup     : $( [[ ! -z $setup_end ]] && echo $(( $setup_end - $setup_start )) || echo 0 )"
printf "\n Spinup     : $( [[ ! -z $spinup_end ]] && echo $(( $spinup_end - $spinup_start )) || echo 0 )"
printf "\n Jacobian     : $( [[ ! -z $jacobian_end ]] && echo $(( $jacobian_end - $jacobian_start )) || echo 0 )"
printf "\n Inversion     : $( [[ ! -z $inversion_end ]] && echo $(( $inversion_end - $inversion_start )) || echo 0 )"
printf "\n Posterior     : $( [[ ! -z $posterior_end ]] && echo $(( $posterior_end - $posterior_start )) || echo 0 )\n\n"

# copy output log to run directory for storage
if [[ -f ${InversionPath}/imi_output.log ]]; then
    cp "${InversionPath}/imi_output.log" "${RunDirs}/imi_output.log"
fi

# copy config file to run directory
cd $InversionPath
cp $ConfigFile "${RunDirs}/config_${RunName}.yml"

exit 0

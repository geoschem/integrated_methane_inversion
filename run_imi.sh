#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o "imi_output.log"

# This script will run the Integrated Methane Inversion (IMI) with GEOS-Chem.
# For documentation, see https://imi.readthedocs.io.
#
# Authors: Daniel Varon, Melissa Sulprizio, Lucas Estrada, Will Downs

# Error message for if the IMI fails
imi_failed() {
    printf "\nFATAL ERROR: IMI exiting."
    cp "${InversionPath}/imi_output.log" "${OutputPath}/${RunName}/imi_output.log"
    exit 1
}

start_time=$(date)
setup_start=$(date +%s)

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE ===\n"

# Check if user has specified a configuration file
if [[ $# == 1 ]] ; then
    ConfigFile=$1
else
    ConfigFile="config.yml"
fi

# Get configuration
source src/utilities/parse_yaml.sh
eval $(parse_yaml ${ConfigFile})

# Set path to IMI runs
RunDirs="${OutputPath}/${RunName}"

##=======================================================================
## Standard settings
##=======================================================================

# In safe mode check whether selected options will overwrite existing files
if "$SafeMode"; then
    if ([ -d "${RunDirs}/spinup_run" ] && "$DoSpinup") \
       || ([ -d "${RunDirs}/jacobian_runs" ] && "$DoJacobian") \
       || ([ -d "${RunDirs}/inversion" ] && "$DoInversion") \
       || ([ -d "${RunDirs}/posterior_run" ] && "$DoPosterior"); then
        
        printf "\nError: files in ${RunDirs}/ may be overwritten. Please change RunName in the IMI config file to avoid overwriting files."
        printf "To proceed, and overwrite existing files, set SafeMode in the config file to false." 
        printf "\nIMI $RunName Aborted"
        exit 1 
    fi
fi

# Path to inversion setup
InversionPath=$(pwd -P)

##=======================================================================
##  Download the TROPOMI data
##=======================================================================

# Download TROPOMI data from AWS. You will be charged if your ec2 instance is not in the eu-central-1 region.
if "$isAWS"; then
    { # test if instance has access to TROPOMI bucket
        stdout=`aws s3 ls s3://meeo-s5p`
    } || { # catch 
        printf "\nError: Unable to connect to TROPOMI bucket. This is likely caused by misconfiguration of the ec2 instance iam role s3 permissions.\n"
        printf "IMI $RunName Aborted."
        exit 1
    }
    tropomiCache=${RunDirs}/data_TROPOMI
    mkdir -p -v $tropomiCache
    python src/utilities/download_TROPOMI.py $StartDate $EndDate $tropomiCache
    printf "\nFinished TROPOMI download\n"
fi

##=======================================================================
##  Run the setup script
##=======================================================================

if "$RunSetup"; then

    printf "\n=== RUNNING SETUP SCRIPT ===\n"

    cd ${InversionPath}

    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Run the setup script
    ./setup_imi.sh ${ConfigFile}; wait;

    printf "\n=== DONE RUNNING SETUP SCRIPT ===\n"

fi
setup_end=$(date +%s)


##=======================================================================
##  Submit spinup simulation
##=======================================================================

spinup_start=$(date +%s)
if  "$DoSpinup"; then

    printf "\n=== SUBMITTING SPINUP SIMULATION ===\n"

    cd ${RunDirs}/spinup_run

    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Submit job to job scheduler
    sbatch -W ${RunName}_Spinup.run; wait;

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed

    printf "\n=== DONE SPINUP SIMULATION ===\n"
    
fi
spinup_end=$(date +%s)

##=======================================================================
##  Submit Jacobian simulation
##=======================================================================

jacobian_start=$(date +%s)
if "$DoJacobian"; then

    printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"

    cd ${RunDirs}/jacobian_runs

    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv} 
    fi

    # Submit job to job scheduler
    ./submit_jacobian_simulations_array.sh; wait;

    # check if any jacobians exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed

    printf "\n=== DONE JACOBIAN SIMULATIONS ===\n"

fi
jacobian_end=$(date +%s)

##=======================================================================
##  Process data and run inversion
##=======================================================================

inversion_start=$(date +%s)
if "$DoInversion"; then

    printf "\n=== RUNNING INVERSION ===\n"

    cd ${RunDirs}/inversion

    if ! "$isAWS"; then
        # Replace nCPUs, partitions

        # Activate Conda environment
        printf "\nActivating conda environment: ${CondaEnv}\n"
        conda activate $CondaEnv
    fi

    # Execute inversion driver script
    sbatch -W run_inversion.sh; wait;
        
    printf "\n=== DONE RUNNING INVERSION ===\n"

fi
inversion_end=$(date +%s)

##=======================================================================
##  Submit posterior simulation and process the output
##=======================================================================

posterior_start=$(date +%s)
if "$DoPosterior"; then

    cd ${RunDirs}/posterior_run
    
    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Submit job to job scheduler
    printf "\n=== SUBMITTING POSTERIOR SIMULATION ===\n"
    sbatch -W ${RunName}_Posterior.run; wait;
    printf "\n=== DONE POSTERIOR SIMULATION ===\n"

    cd ${RunDirs}/inversion

    # Fill missing data (first hour of simulation) in posterior output
    PosteriorRunDir="${RunDirs}/posterior_run"
    PrevDir="${RunDirs}/spinup_run"
    printf "\n=== Calling postproc_diags.py for posterior ===\n"
    python postproc_diags.py $RunName $PosteriorRunDir $PrevDir $StartDate; wait
    printf "\n=== DONE -- postproc_diags.py ===\n"

    # Build directory for hourly posterior GEOS-Chem output data
    mkdir -p data_converted_posterior
    mkdir -p data_geoschem_posterior
    GCsourcepth="${PosteriorRunDir}/OutputDir"
    GCDir="./data_geoschem_posterior"
    printf "\n=== Calling setup_gc_cache.py for posterior ===\n"
    python setup_gc_cache.py $StartDate $EndDate $GCsourcepth $GCDir; wait
    printf "\n=== DONE -- setup_gc_cache.py ===\n"

    # Sample GEOS-Chem atmosphere with TROPOMI
    function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
    function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
    LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
    LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
    LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
    LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)
    nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc)
    FetchTROPOMI="False"
    isPost="True"

    printf "\n=== Calling jacobian.py to sample posterior simulation (without jacobian sensitivity analysis) ===\n"
    python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $isPost; wait
    printf "\n=== DONE sampling the posterior simulation ===\n\n"

fi
posterior_end=$(date +%s)

# Remove temporary files
if "$isAWS"; then
    rm -f /home/ubuntu/foo.nc
fi

# Run time
end_time=$(date)
printf "\nIMI started: %s" "$start_time"
printf "\nIMI ended: %s\n\n" "$end_time"

printf "Statistics:"
printf "Setup runtime (s): $(( $setup_end - $setup_start ))"
printf "Spinup runtime (s): $(( $spinup_end - $spinup_start ))"
printf "Inversion runtime (s): $(( $inversion_end - $inversion_start ))"
printf "Jacobian runtime (s): $(( $jacobian_end - $jacobian_start ))"
printf "Posterior runtime (s): $(( $posterior_end - $posterior_start ))"

# copy output log to run directory for storage
cp "${InversionPath}/imi_output.log" "${RunDirs}/imi_output.log"

exit 0

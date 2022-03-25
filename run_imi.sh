#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o "imi_output.log"

# This script will run the Integrated Methane Inversion (IMI) with GEOS-Chem.
# For documentation, see https://imi.readthedocs.io.
#
# Authors: Daniel Varon, Melissa Sulprizio, Lucas Estrada, Will Downs

start_time=$(date)

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

##=======================================================================
## Standard settings
##=======================================================================

# In safe mode check whether selected options will overwrite existing files
if "$SafeMode"; then
    if ([ -d "${OutputPath}/${RunName}/spinup_run" ] && "$DoSpinup") \
       || ([ -d "${OutputPath}/${RunName}/jacobian_runs" ] && "$DoJacobian") \
       || ([ -d "${OutputPath}/${RunName}/inversion" ] && "$DoInversion") \
       || ([ -d "${OutputPath}/${RunName}/posterior_run" ] && "$DoPosterior"); then
        
        echo "Error: files in ${OutputPath}/${RunName}/ may be overwritten. Please change RunName in the IMI config file to avoid overwriting files."
        echo "To proceed, and overwrite existing files, set SafeMode in the config file to false." 
        echo "IMI $RunName Aborted"
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
        printf "Error: Unable to connect to TROPOMI bucket. This is likely caused by misconfiguration of the ec2 instance iam role s3 permissions."
        printf "IMI $RunName Aborted."
        exit 1
    }
    tropomiCache=${OutputPath}/${RunName}/data_TROPOMI
    mkdir -p -v $tropomiCache
    python src/utilities/download_TROPOMI.py $StartDate $EndDate $tropomiCache
    printf "Finished TROPOMI download"
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
    ./setup_imi.sh; wait;

    printf "\n=== DONE RUNNING SETUP SCRIPT ===\n"

fi

##=======================================================================
##  Submit spinup simulation
##=======================================================================

if  "$DoSpinup"; then

    printf "\n=== SUBMITTING SPINUP SIMULATION ===\n"

    cd ${OutputPath}/${RunName}/spinup_run

    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Submit job to job scheduler
    sbatch -W ${RunName}_Spinup.run; wait;

    printf "=== DONE SPINUP SIMULATION ===\n"
    
fi

##=======================================================================
##  Submit Jacobian simulation
##=======================================================================

if "$DoJacobian"; then

    printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"

    cd ${OutputPath}/${RunName}/jacobian_runs

    if ! "$isAWS"; then
        # Replace nCPUs, partitions

        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv} 
    fi

    # Submit job to job scheduler
    ./submit_jacobian_simulations_array.sh; wait;

    printf "=== DONE JACOBIAN SIMULATIONS ===\n"

fi

##=======================================================================
##  Process data and run inversion
##=======================================================================

if "$DoInversion"; then

    printf "\n=== RUNNING INVERSION ===\n"

    cd ${OutputPath}/${RunName}/inversion

    if ! "$isAWS"; then
        # Replace nCPUs, partitions

        # Activate Conda environment
        printf "Activating conda environment: ${CondaEnv}\n"
        conda activate $CondaEnv
    fi

    # Execute inversion driver script
    sbatch -W run_inversion.sh; wait;
        
    printf "=== DONE RUNNING INVERSION ===\n"

fi

##=======================================================================
##  Submit posterior simulation and process the output
##=======================================================================

if "$DoPosterior"; then

    cd ${OutputPath}/${RunName}/posterior_run
    
    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Submit job to job scheduler
    printf "\n=== SUBMITTING POSTERIOR SIMULATION ===\n"
    sbatch -W ${RunName}_Posterior.run; wait;
    printf "=== DONE POSTERIOR SIMULATION ===\n"

    cd ${OutputPath}/${RunName}/inversion

    # Fill missing data (first hour of simulation) in posterior output
    PosteriorRunDir="${OutputPath}/${RunName}/posterior_run"
    PrevDir="${OutputPath}/${RunName}/spinup_run"
    printf "\n=== Calling postproc_diags.py for posterior ===\n"
    python postproc_diags.py $RunName $PosteriorRunDir $PrevDir $StartDate; wait
    printf "=== DONE -- postproc_diags.py ===\n"

    # Build directory for hourly posterior GEOS-Chem output data
    mkdir -p data_converted_posterior
    mkdir -p data_geoschem_posterior
    GCsourcepth="${PosteriorRunDir}/OutputDir"
    GCDir="./data_geoschem_posterior"
    printf "\n=== Calling setup_gc_cache.py for posterior ===\n"
    python setup_gc_cache.py $StartDate $EndDate $GCsourcepth $GCDir; wait
    printf "=== DONE -- setup_gc_cache.py ===\n"

    # Sample GEOS-Chem atmosphere with TROPOMI
    function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
    function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
    LonMinInvDomain=$(ncmin lon ${OutputPath}/${RunName}/StateVector.nc)
    LonMaxInvDomain=$(ncmax lon ${OutputPath}/${RunName}/StateVector.nc)
    LatMinInvDomain=$(ncmin lat ${OutputPath}/${RunName}/StateVector.nc)
    LatMaxInvDomain=$(ncmax lat ${OutputPath}/${RunName}/StateVector.nc)
    nElements=$(ncmax StateVector ${OutputPath}/${RunName}/StateVector.nc)
    FetchTROPOMI="False"
    isPost="True"

    printf "\n=== Calling jacobian.py to sample posterior simulation (without jacobian sensitivity analysis) ===\n"
    python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $isPost; wait
    printf "=== DONE sampling the posterior simulation ===\n\n"

fi

# Remove temporary files
if "$isAWS"; then
    rm -f /home/ubuntu/foo.nc
fi

# Run time
end_time=$(date)
printf "\nIMI started: %s" "$start_time"
printf "\nIMI ended: %s\n\n" "$end_time"

exit 0

#!/bin/bash

# This script will run a CH4 analytical inversion with GEOS-Chem.
# (mps, 2/20/2021)
# (djv, 12/7/2021)

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE ===\n"

# Function to parse yaml files from shell script
# By Stefan Farestam via stackoverflow:
# https://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# Get configuration
eval $(parse_yaml config.yml)
# For reference, this defines the following environment variables:
# General: $isAWS, $RunName
# Period of interest: $StartDate, $EndDate, $SpinupMonths
# Region of interest: $LonMin, $LonMax, $LatMin, $LatMax
# Inversion: $PriorError, $ObsError, $Gamma
# Grid: $Res, $Met, $HalfPolar, $Levs, $NestedGrid, $REGION, $Buffer
# Setup modules: $CreateStateVectorFile, $SetupTemplateRundir, $SetupSpinupRun, $SetupJacobianRuns, $SetupInversion, $SetupPosteriorRun
# Run modules: $RunSetup, $DoSpinup, $DoJacobian, $DoInversion, $DoPosterior
# State vector: $BufferDeg, $nBufferClusters, $LandThreshold
# If custom state vec file: $StateVectorFile, $LonMinCustomStateVector, $LonMaxCustomStateVector, $LatMinCustomStateVector, $LatMaxCustomStateVector
# Harvard-Cannon: $nCPUs, $partition

# My path
if "$isAWS"; then
    MyPath="/home/ubuntu/CH4_Workflow"
    SetupPath="/home/ubuntu/setup_CH4"
else
    MyPath="/n/holyscratch01/jacob_lab/msulprizio/CH4"
    SetupPath="FILL"
fi

## ======================================================================
## Specific to Harvard's Cannon cluster
## ======================================================================

# Path to inversion setup
InversionPath=$(pwd -P)

# Environment files
NCOEnv="${InversionPath}/envs/Harvard-Cannon/gcc.ifort17_cannon.env"
GCCEnv="${InversionPath}/envs/Harvard-Cannon/gcc.gfortran10.2_cannon.env"
CondaEnv="ch4_inv"
FortranCompiler="~/env/envs/gcc_cmake.ifort17_openmpi_cannon.env"

##=======================================================================
##  Run the setup script
##=======================================================================
if "$RunSetup"; then

    printf "\n=== RUNNING SETUP SCRIPT ===\n"

    cd ${SetupPath}

    if ! "$isAWS"; then
        # Load fortran compiler
        source ${FortranCompiler}
    fi

    # Run the setup script
    ./setup_ch4_inversion.sh; wait;

fi

##=======================================================================
##  Submit spinup simulation
##=======================================================================
if  "$DoSpinup"; then

    printf "\n=== SUBMITTING SPINUP SIMULATION ===\n"

    cd ${MyPath}/${RunName}/spinup_run

    if ! "$isAWS"; then
        # Replace nCPUs, partitions
        
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GCCEnv}
    fi

    # Submit job to job scheduler
    sbatch ${RunName}_Spinup.run; wait;

    printf "\n=== DONE SPINUP SIMULATION ===\n"
    
fi

##=======================================================================
##  Submit Jacobian simulation
##=======================================================================
if "$DoJacobian"; then

    printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"

    cd ${MyPath}/${RunName}/jacobian_runs

    if ! "$isAWS"; then
        # Replace nCPUs, partitions

        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GCCEnv} 
    fi

    # Submit job to job scheduler
    ./submit_jacobian_simulations_array.sh; wait;

    printf "\n=== DONE JACOBIAN SIMULATIONS ===\n"

fi

##=======================================================================
##  Process data and run inversion
##=======================================================================
if "$DoInversion"; then

    printf "\n=== RUNNING INVERSION ===\n"

    cd ${MyPath}/${RunName}/inversion

    if ! "$isAWS"; then
        # Replace nCPUs, partitions

        # Activate Conda environment
        printf "Activating conda environment: ${CondaEnv}\n"
        source activate $CondaEnv
    fi

    # Execute inversion driver script
    sbatch run_inversion.sh; wait;
        
    printf "=== DONE RUNNING INVERSION ===\n"

fi

##=======================================================================
##  Submit posterior simulation and process output
##=======================================================================
if "$DoPosterior"; then

    cd ${MyPath}/${RunName}/posterior_run
    
    if ! "$isAWS"; then
        # Replace nCPUs, partitions

        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GCCEnv}
    fi

    # Submit job to job scheduler
    sbatch ${RunName}_Posterior.run; wait;

    cd ${MyPath}/${RunName}/inversion

    # Fill missing data (first hour of simulation) in posterior output
    PosteriorRunDir="${MyPath}/${RunName}/posterior_run"
    PrevDir="${MyPath}/${RunName}/spinup_run"
    printf "Calling postproc_diags.py for posterior\n"
    python postproc_diags.py $RunName $PosteriorRunDir $PrevDir $StartDate; wait
    printf "DONE -- postproc_diags.py\n\n"

    # Build directory for hourly posterior GEOS-Chem output data
    mkdir -p data_GC_posterior
    GCsourcepth="${PosteriorRunDir}/OutputDir"
    GCDir="./data_GC_posterior"
    printf "Calling setup_GCdatadir.py for posterior\n"
    python setup_GCdatadir.py $StartDate $EndDate $GCsourcepth $GCDir; wait
    printf "DONE -- setup_GCdatadir.py\n\n"

    useSensi="True"
    python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $FetchTROPOMI $useSensi; wait

fi

exit 0
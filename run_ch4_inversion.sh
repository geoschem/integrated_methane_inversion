#!/bin/bash

# This script will run a CH4 analytical inversion with GEOS-Chem.
# (mps, 2/20/2021)
# (djv, 12/7/2021)

start_time=$(date)

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE ===\n"

# Get configuration
source parse_yaml.sh
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
##  Download the TROPOMI data
##=======================================================================
if "$isAWS"; then
    Sat_datadir=${MyPath}/${RunName}/data_TROPOMI
    mkdir -p -v $Sat_datadir
    python PostprocessingScripts/CH4_TROPOMI_INV/download_TROPOMI.py $StartDate $EndDate $Sat_datadir
fi

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

    printf "\n=== DONE RUNNING SETUP SCRIPT ===\n"

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
    sbatch -W ${RunName}_Spinup.run; wait;

    printf "=== DONE SPINUP SIMULATION ===\n"
    
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

    printf "=== DONE JACOBIAN SIMULATIONS ===\n"

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
    sbatch -W run_inversion.sh; wait;
        
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
    printf "\n=== SUBMITTING POSTERIOR SIMULATION ===\n"
    sbatch -W ${RunName}_Posterior.run; wait;
    printf "=== DONE POSTERIOR SIMULATION ===\n"

    cd ${MyPath}/${RunName}/inversion

    # Fill missing data (first hour of simulation) in posterior output
    PosteriorRunDir="${MyPath}/${RunName}/posterior_run"
    PrevDir="${MyPath}/${RunName}/spinup_run"
    printf "\n=== Calling postproc_diags.py for posterior ===\n"
    python postproc_diags.py $RunName $PosteriorRunDir $PrevDir $StartDate; wait
    printf "=== DONE -- postproc_diags.py ===\n"

    # Build directory for hourly posterior GEOS-Chem output data
    mkdir -p data_converted_posterior
    mkdir -p data_GC_posterior
    GCsourcepth="${PosteriorRunDir}/OutputDir"
    GCDir="./data_GC_posterior"
    printf "\n=== Calling setup_GCdatadir.py for posterior ===\n"
    python setup_GCdatadir.py $StartDate $EndDate $GCsourcepth $GCDir; wait
    printf "=== DONE -- setup_GCdatadir.py ===\n"

    # Sample GEOS-Chem atmosphere with TROPOMI
    function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
    # Here we use bc to do floating point arithmetic
    LonMinInvDomain=`echo "$LonMin-$BufferDeg" | bc`
    LonMaxInvDomain=`echo "$LonMax+$BufferDeg" | bc`
    LatMinInvDomain=`echo "$LatMin-$BufferDeg" | bc`
    LatMaxInvDomain=`echo "$LatMax+$BufferDeg" | bc`
    StateVectorFilePath="${MyPath}/${RunName}/StateVector.nc"
    if ! "$CreateStateVectorFile"; then
        function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
        LonMinInvDomain=$(ncmin lon $StateVectorFile)
        LonMaxInvDomain=$(ncmax lon $StateVectorFile)
        LatMinInvDomain=$(ncmin lat $StateVectorFile)
        LatMaxInvDomain=$(ncmax lat $StateVectorFile)
        StateVectorFilePath=$StateVectorFile
    fi
    nElements=$(ncmax StateVector $StateVectorFilePath)
    FetchTROPOMI="False"
    isPost="True"

    printf "\n=== Calling jacobian.py to sample posterior simulation (without jacobian sensitivity analysis) ===\n"
    python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $FetchTROPOMI $isPost; wait
    printf "=== DONE sampling the posterior simulation ===\n\n"

fi

end_time=$(date)
printf "\nIMI started: %s" "$start_time"
printf "\nIMI ended: %s\n\n" "$end_time"

exit 0
#!/bin/bash

#SBATCH -o run_inversion_%j.out

##=======================================================================
## Parse config.yml file
##=======================================================================

send_error() {
    file=`basename "$0"`
    printf "\nInversion Error: on line number ${1} of ${file}: IMI exiting."
    echo "Error Status: 1" > .error_status_file.txt
    exit 1
}

# remove error status file if present
rm -f .error_status_file.txt

# trap and exit on errors
trap 'send_error $LINENO' ERR

printf "\n=== PARSING CONFIG FILE ===\n"

invPath={INVERSION_PATH}
configFile={CONFIG_FILE}

# Get configuration
#  This defines $StartDate, $EndDate, $nBufferClusters, $RunName
#  It also define $PriorError, $ObsError, $Gamma, $PrecomputedJacobian
#  Parsing the config file here facilitates generation of inversion ensembles
#  All that needs to be done is to edit the config file for $PriorError,
#   $ObsError, and $Gamma
#  Make sure $PrecomputedJacobian is true, and then re-run this script
#   (or run_imi.sh with only the $DoInversion module switched on in config.yml).
PythonEnv=$(grep '^PythonEnv:' ${invPath}/${configFile} |
    sed 's/PythonEnv://' |
    sed 's/#.*//' |
    sed 's/^[[:space:]]*//' |
    tr -d '"')
echo $PythonEnv
source ${invPath}/${PythonEnv}
eval $(python ${invPath}/src/utilities/parse_yaml.py ${invPath}/${configFile})

#=======================================================================
# Configuration (these settings generated on initial setup)
#=======================================================================
LonMinInvDomain={LON_MIN}
LonMaxInvDomain={LON_MAX}
LatMinInvDomain={LAT_MIN}
LatMaxInvDomain={LAT_MAX}
nElements={STATE_VECTOR_ELEMENTS}
nTracers={NUM_JACOBIAN_TRACERS}
OutputPath={OUTPUT_PATH}
Res={RES}
JacobianRunsDir="${OutputPath}/${RunName}/jacobian_runs"
PriorRunDir="${JacobianRunsDir}/${RunName}_0000"
BackgroundRunDir="${JacobianRunsDir}/${RunName}_background"
PosteriorRunDir="${OutputPath}/${RunName}/posterior_run"
StateVectorFile={STATE_VECTOR_PATH}
GCDir="${OutputPath}/${RunName}/inversion/data_geoschem"
GCVizDir="${OutputPath}/${RunName}/inversion/data_geoschem_prior"
JacobianDir="${OutputPath}/${RunName}/inversion/data_converted"
sensiCache="${OutputPath}/${RunName}/inversion/data_sensitivities"
satelliteCache="${OutputPath}/${RunName}/satellite_data"
period_i={PERIOD}

# For Kalman filter: assume first inversion period (( period_i = 1 )) by default
# Switch is flipped to false automatically if (( period_i > 1 ))
FirstSimSwitch=$1

printf "\n=== EXECUTING RUN_INVERSION.SH ===\n"
    
#=======================================================================
# Error checks
#=======================================================================

# Make sure specified paths exist
if [[ ! -d ${JacobianRunsDir} ]]; then
    printf "${JacobianRunsDir} does not exist. Please fix JacobianRunsDir in run_inversion.sh.\n"
    exit 1
fi
if [[ ! -f ${StateVectorFile} ]]; then
    printf "${StateVectorFile} does not exist. Please fix StateVectorFile in run_inversion.sh.\n"
    exit 1
fi

#=======================================================================
# Setup GC data directory in workdir
#=======================================================================

printf "Calling setup_gc_cache.py\n"
if "$LognormalErrors"; then
    # for lognormal errors we use the clean background run
    GCsourcepth="${BackgroundRunDir}/OutputDir"
    PriorOutputDir="${PriorRunDir}/OutputDir"
    # also need the prior cache so that we can visualize the prior simulation
    python ${OutputPath}/${RunName}/inversion/setup_gc_cache.py $StartDate $EndDate $PriorOutputDir $GCVizDir; wait
else
    # for normal errors we use the prior run
    GCsourcepth="${PriorRunDir}/OutputDir"
fi

export PYTHONPATH=${PYTHONPATH}:${OutputPath}
python ${OutputPath}/${RunName}/inversion/setup_gc_cache.py $StartDate $EndDate $GCsourcepth $GCDir; wait
printf "DONE -- setup_gc_cache.py\n\n"

#=======================================================================
# setup geoschem cache for pseudo observations if doing OSSE
#=======================================================================
if "$EnableOSSE"; then
    RunDirOSSE="${OutputPath}/${RunName}/osse_observations_run"
    GCDirOSSE="./data_geoschem_osse"
    mkdir -p  $GCDirOSSE
    # If simulating observations, we need to postprocess the observation data
    python postproc_diags.py $RunName $RunDirOSSE $PrevDir $StartDate $Res; wait
    python setup_gc_cache.py $StartDate $EndDate "${RunDirOSSE}/OutputDir" $GCDirOSSE; wait
fi

#=======================================================================
# Generate Jacobian matrix files 
#=======================================================================

printf "Calling jacobian.py\n"
isPost="False"
if ! "$PrecomputedJacobian"; then
    buildJacobian="True"
    jacobian_sf="None"
else
    buildJacobian="False"
    jacobian_sf=${OutputPath}/${RunName}/inversion/jacobian_scale_factors.npy
fi

python -u ${OutputPath}/${RunName}/inversion/jacobian.py ${OutputPath}/${RunName}/inversion ${invPath}/${configFile} $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $Species $satelliteCache $SatelliteProduct $UseWaterObs $isPost $period_i $buildJacobian False; wait
if "$LognormalErrors"; then
    # for lognormal error visualization of the prior we sample the prior run
    # without constructing the jacobian matrix
    python ${OutputPath}/${RunName}/inversion/jacobian.py ${OutputPath}/${RunName}/inversion ${invPath}/${configFile} $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $Species $satelliteCache $SatelliteProduct $UseWaterObs $isPost $period_i False True; wait
fi
printf " DONE -- jacobian.py\n\n"

#=======================================================================
# Do inversion
#=======================================================================
if "$LognormalErrors"; then
    # for lognormal errors we merge our y, y_bkgd and partial K matrices
    python ${OutputPath}/${RunName}/inversion/merge_partial_k.py $JacobianDir $StateVectorFile ${OutputPath}/${RunName}/config_${RunName}.yml $PrecomputedJacobian

    # then we run the inversion
    printf "Calling lognormal_invert.py\n"
    python ${OutputPath}/${RunName}/inversion/lognormal_invert.py ${invPath}/${configFile} $StateVectorFile $jacobian_sf
    printf "DONE -- lognormal_invert.py\n\n"
else
    posteriorSF="./inversion_result.nc"
    python_args=(${OutputPath}/${RunName}/inversion/invert.py ${OutputPath}/${RunName}/config_${RunName}.yml $nElements $JacobianDir $posteriorSF $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $Res $jacobian_sf $StateVectorFile)

    printf "Calling invert.py\n"
    python "${python_args[@]}"; wait
    printf "DONE -- invert.py\n\n"
    #=======================================================================
    # Create gridded posterior scaling factor netcdf file
    #=======================================================================
    GriddedPosterior="./gridded_posterior.nc"

    printf "Calling make_gridded_posterior.py\n"
    python ${OutputPath}/${RunName}/inversion/make_gridded_posterior.py $posteriorSF $StateVectorFile $GriddedPosterior; wait
    printf "DONE -- make_gridded_posterior.py\n\n"
fi

printf "Exiting run_inversion.sh"

exit 0

#!/bin/bash

#SBATCH -N 1
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
#  This defines $StartDate, $EndDate, $nBufferClusters, $RunName, $isAWS
#  It also define $PriorError, $ObsError, $Gamma, $PrecomputedJacobian
#  Parsing the config file here facilitates generation of inversion ensembles
#  All that needs to be done is to edit the config file for $PriorError,
#   $ObsError, and $Gamma
#  Make sure $PrecomputedJacobian is true, and then re-run this script
#   (or run_imi.sh with only the $DoInversion module switched on in config.yml).

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
SpinupDir="${OutputPath}/${RunName}/spinup_run"
JacobianRunsDir="${OutputPath}/${RunName}/jacobian_runs"
PriorRunDir="${JacobianRunsDir}/${RunName}_0000"
BackgroundRunDir="${JacobianRunsDir}/${RunName}_background"
PosteriorRunDir="${OutputPath}/${RunName}/posterior_run"
StateVectorFile={STATE_VECTOR_PATH}
GCDir="./data_geoschem"
GCVizDir="./data_geoschem_prior"
JacobianDir="./data_converted"
sensiCache="./data_sensitivities"
tropomiCache="${OutputPath}/${RunName}/satellite_data"
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
# Postprocess the SpeciesConc and LevelEdgeDiags files from GEOS-Chem
#=======================================================================

printf "Calling postproc_diags.py, FSS=$FirstSimSwitch\n"
if "$FirstSimSwitch"; then
    if [[ ! -d ${SpinupDir} ]]; then
	printf "${SpinupDir} does not exist. Please fix SpinupDir or set FirstSimSwitch to False in run_inversion.sh.\n"
	exit 1
    fi
    PrevDir=$SpinupDir
else
    PrevDir=$PosteriorRunDir
    if [[ ! -d ${PosteriorRunDir} ]]; then
	printf "${PosteriorRunDir} does not exist. Please fix PosteriorRunDir in run_inversion.sh.\n"
	exit 1
    fi
fi
printf "  - Hour 0 for ${StartDate} will be obtained from ${PrevDir}\n"

if ! "$PrecomputedJacobian"; then

    # Postprocess all the Jacobian simulations
    python postproc_diags.py $RunName $JacobianRunsDir $PrevDir $StartDate $Res; wait

else

    # Only postprocess the Prior simulation
    python postproc_diags.py $RunName $PriorRunDir $PrevDir $StartDate $Res; wait
    if "$LognormalErrors"; then
        # for lognormal errors we need to postprocess the background run too
        python postproc_diags.py $RunName $BackgroundRunDir $PrevDir $StartDate $Res; wait
    fi
fi
printf "DONE -- postproc_diags.py\n\n"

#=======================================================================
# Calculate GEOS-Chem sensitivities and save to sensitivities directory
#=======================================================================

if ! "$PrecomputedJacobian"; then
    # add an argument to calc_sensi.py if optimizing BCs and/or OH
    if "$OptimizeBCs"; then
        pertBCs=$PerturbValueBCs
    else
	pertBCs=0.0
    fi
    if "$OptimizeOH"; then
        pertOH=$PerturbValueOH
    else
	pertOH=0.0
    fi
    PerturbPath=${OutputPath}/${RunName}/archive_perturbation_sfs/flat_pert_sf_${period_i}.npy
    python_args=(calc_sensi.py $nElements $nTracers $PerturbPath $StartDate $EndDate $JacobianRunsDir $RunName $sensiCache $pertBCs $pertOH)
    printf "Calling calc_sensi.py\n"
    python "${python_args[@]}"; wait
    printf "DONE -- calc_sensi.py\n\n"
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
    python setup_gc_cache.py $StartDate $EndDate $PriorOutputDir $GCVizDir; wait
else
    # for normal errors we use the prior run
    GCsourcepth="${PriorRunDir}/OutputDir"
fi

python setup_gc_cache.py $StartDate $EndDate $GCsourcepth $GCDir; wait
printf "DONE -- setup_gc_cache.py\n\n"

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
    jacobian_sf=./jacobian_scale_factors.npy

fi

python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $BlendedTROPOMI $isPost $buildJacobian False; wait
if "$LognormalErrors"; then
    # for lognormal error visualization of the prior we sample the prior run
    # without constructing the jacobian matrix
    python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $BlendedTROPOMI $isPost False True; wait
fi
printf " DONE -- jacobian.py\n\n"

#=======================================================================
# Do inversion
#=======================================================================
if "$OptimizeBCs"; then
    ErrorBCs=$PriorErrorBCs
else
    ErrorBCs=0.0
fi
if "$OptimizeOH"; then
    ErrorOH=$PriorErrorOH
else
    ErrorOH=0.0
fi

if "$LognormalErrors"; then
    # for lognormal errors we merge our y, y_bkgd and partial K matrices
    python merge_partial_k.py $JacobianDir $StateVectorFile $ObsError $PrecomputedJacobian

    # then we run the inversion
    printf "Calling lognormal_invert.py\n"
    python lognormal_invert.py ${invPath}/${configFile} $StateVectorFile $jacobian_sf
    printf "DONE -- lognormal_invert.py\n\n"
else
    posteriorSF="./inversion_result.nc"
    python_args=(invert.py $nElements $JacobianDir $posteriorSF $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $PriorError $ObsError $Gamma $Res $jacobian_sf $ErrorBCs $ErrorOH)

    printf "Calling invert.py\n"
    python "${python_args[@]}"; wait
    printf "DONE -- invert.py\n\n"
    #=======================================================================
    # Create gridded posterior scaling factor netcdf file
    #=======================================================================
    GriddedPosterior="./gridded_posterior.nc"

    printf "Calling make_gridded_posterior.py\n"
    python make_gridded_posterior.py $posteriorSF $StateVectorFile $GriddedPosterior; wait
    printf "DONE -- make_gridded_posterior.py\n\n"
fi

exit 0

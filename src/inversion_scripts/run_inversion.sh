#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH -t 0-03:00
#SBATCH --mem 4000
#SBATCH -o run_inversion_%j.out
#SBATCH -e run_inversion_%j.err

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE ===\n"

# Get configuration (relative paths assume script is run from the inversion directory)
source ../../../integrated_methane_inversion/src/utilities/parse_yaml.sh
eval $(parse_yaml ../../../integrated_methane_inversion/config.yml)
# This defines $StartDate, $EndDate, $nBufferClusters, $RunName, $isAWS
# It also define $PriorError, $ObsError, $Gamma, $PrecomputedJacobian
# Parsing the config file here facilitates generation of inversion ensembles
# All that needs to be done is to edit the config file for $PriorError, $ObsError, and $Gamma, 
# make sure $PrecomputedJacobian is true, and then re-run this script (or run_imi.sh
# with only the $DoInversion module switched on in config.yml).

#=======================================================================
# Configuration (these settings generated on initial setup)
#=======================================================================
LonMinInvDomain={LON_MIN}
LonMaxInvDomain={LON_MAX}
LatMinInvDomain={LAT_MIN}
LatMaxInvDomain={LAT_MAX}
nElements={STATE_VECTOR_ELEMENTS}
MyPath={MY_PATH}
Res={RES}
SpinupDir="${MyPath}/${RunName}/spinup_run"
JacobianRunsDir="${MyPath}/${RunName}/jacobian_runs"
PosteriorRunDir="${MyPath}/${RunName}/posterior_run"
StateVectorFile={STATE_VECTOR_PATH}
GCDir="./data_GC"
JacobianDir="./data_converted"
sensiCache="./data_sensitivities"
tropomiCache="${MyPath}/${RunName}/data_TROPOMI"

# Only matters for Kalman filter inversions, to be implemented in a future version of the IMI
FirstSimSwitch=true

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

if ! "$PrecomputedJacobian"; then

    printf "Calling postproc_diags.py, FSS=$FirstSimSwitch\n"
    if "$FirstSimSwitch"; then
        if [[ ! -d ${SpinupDir} ]]; then
        printf "${SpinupDir} does not exist. Please fix SpinupDir or set FirstSimSwitch to False in run_inversion.sh.\n"
        exit 1
        fi
        PrevDir=$SpinupDir
    else
        PrevDir=%$PosteriorRunDir
        if [[ ! -d ${PosteriorRunDir} ]]; then
        printf "${PosteriorRunDir} does not exist. Please fix PosteriorRunDir in run_inversion.sh.\n"
        exit 1
        fi
    fi
    printf "  - Hour 0 for ${StartDate} will be obtained from ${PrevDir}\n"

    python postproc_diags.py $RunName $JacobianRunsDir $PrevDir $StartDate; wait
    printf "DONE -- postproc_diags.py\n\n"

fi

#=======================================================================
# Calculate GEOS-Chem sensitivities and save to sensitivities directory
#=======================================================================

if ! "$PrecomputedJacobian"; then

    # 50% perturbation
    Perturbation=0.5

    printf "Calling calc_sensi.py\n"
    python calc_sensi.py $nElements $Perturbation $StartDate $EndDate $JacobianRunsDir $RunName $sensiCache; wait
    printf "DONE -- calc_sensi.py\n\n"

fi

#=======================================================================
# Setup GC data directory in workdir
#=======================================================================

if ! "$PrecomputedJacobian"; then

    GCsourcepth="${JacobianRunsDir}/${RunName}_0000/OutputDir"

    printf "Calling setup_gc_cache.py\n"
    python setup_gc_cache.py $StartDate $EndDate $GCsourcepth $GCDir; wait
    printf "DONE -- setup_gc_cache.py\n\n"

fi

#=======================================================================
# Generate Jacobian matrix files 
#=======================================================================

if ! "$PrecomputedJacobian"; then

    printf "Calling jacobian.py\n"
    isPost="False"
    python jacobian.py $StartDate $EndDate $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $isPost; wait
    printf " DONE -- jacobian.py\n\n"

fi

#=======================================================================
# Do inversion
#=======================================================================

posteriorSF="./inversion_result.nc"

printf "Calling invert.py\n"
python invert.py $nElements $JacobianDir $posteriorSF $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $PriorError $ObsError $Gamma $Res; wait
printf "DONE -- invert.py\n\n"

#=======================================================================
# Create gridded posterior scaling factor netcdf file
#=======================================================================
GriddedPosterior="./gridded_posterior.nc"

printf "Calling make_gridded_posterior.py\n"
python make_gridded_posterior.py $posteriorSF $StateVectorFile $GriddedPosterior; wait
printf "DONE -- make_gridded_posterior.py\n\n"

exit 0

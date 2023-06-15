#!/bin/bash

# Functions available in this file include:
#   - setup_inversion 
#   - run_inversion 

# Description: Setup inversion run directory
# Usage:
#   setup_inversion
setup_inversion() {
    printf "\n=== SETTING UP INVERSION DIRECTORY ===\n"
    
    cd ${OutputPath}/$RunName
    mkdir -p -v inversion
    mkdir -p inversion/data_converted
    mkdir -p inversion/data_geoschem
    mkdir -p inversion/data_sensitivities
    mkdir -p inversion/data_visualization
    mkdir -p inversion/operators
    
    cp ${InversionPath}/src/inversion_scripts/calc_sensi.py inversion/
    cp ${InversionPath}/src/inversion_scripts/invert.py inversion/
    cp ${InversionPath}/src/inversion_scripts/jacobian.py inversion/
    cp ${InversionPath}/src/inversion_scripts/operators/* inversion/operators/
    cp ${InversionPath}/src/inversion_scripts/make_gridded_posterior.py inversion/
    cp ${InversionPath}/src/inversion_scripts/postproc_diags.py inversion/
    cp ${InversionPath}/src/inversion_scripts/setup_gc_cache.py inversion/
    cp ${InversionPath}/src/inversion_scripts/utils.py inversion/
    cp ${InversionPath}/src/inversion_scripts/run_inversion.sh inversion/
    cp ${InversionPath}/src/notebooks/visualization_notebook.ipynb inversion/
    sed -i -e "s:{INVERSION_PATH}:${InversionPath}:g" \
           -e "s:{CONFIG_FILE}:${ConfigFile}:g" \
           -e "s:{STATE_VECTOR_ELEMENTS}:${nElements}:g" \
           -e "s:{OUTPUT_PATH}:${OutputPath}:g" \
           -e "s:{STATE_VECTOR_PATH}:../StateVector.nc:g" \
           -e "s:{LON_MIN}:${LonMinInvDomain}:g" \
           -e "s:{LON_MAX}:${LonMaxInvDomain}:g" \
           -e "s:{LAT_MIN}:${LatMinInvDomain}:g" \
           -e "s:{LAT_MAX}:${LatMaxInvDomain}:g" \
           -e "s:{RES}:${gridResLong}:g" inversion/run_inversion.sh

    printf "\n=== DONE SETTING UP INVERSION DIRECTORY ===\n"
}

# Description: Run inversion
# Usage:
#   run_inversion
run_inversion() {
    inversion_start=$(date +%s)
    printf "\n=== RUNNING INVERSION ===\n"

    cd ${RunDirs}/inversion

    if ! "$isAWS"; then
        # Activate Conda environment
        printf "\nActivating conda environment: ${CondaEnv}\n"
        eval "$(conda shell.bash hook)"
        conda activate $CondaEnv
    fi

    # Execute inversion driver script
    sbatch --mem $SimulationMemory \
           -c $SimulationCPUs \
           -t $RequestedTime \
           -p $SchedulerPartition \
           -W run_inversion.sh; wait;

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO
        
    printf "\n=== DONE RUNNING INVERSION ===\n"
    inversion_end=$(date +%s)
}
#!/bin/bash

# Functions available in this file include:
#   - setup_inversion 
#   - run_inversion 
#   - run_notebooks 

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
           -e "s:{RES}:${Res}:g" inversion/run_inversion.sh

    if "$KalmanMode"; then
        # Rename inversion directory as template directory
        mv ${RunDirs}/inversion ${RunDirs}/inversion_template
    fi

    printf "\n=== DONE SETTING UP INVERSION DIRECTORY ===\n"
}

# Description: Run inversion
# Usage:
#   run_inversion
run_inversion() {
    inversion_start=$(date +%s)
    printf "\n=== RUNNING INVERSION ===\n"
    FirstSimSwitch=true
    if "$KalmanMode"; then
        cd ${RunDirs}/kf_inversions/period${period_i}
        # Modify inversion driver script to reflect current inversion period
        sed -i "s|data_TROPOMI\"|data_TROPOMI\"\n\n# Defined via run_kf.sh:\nStartDate=${StartDate_i}\nEndDate=${EndDate_i}|g" run_inversion.sh
        if (( period_i > 1 )); then
            FirstSimSwitch=false
        fi
    else
        cd ${RunDirs}/inversion
    fi

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
           -W run_inversion.sh $FirstSimSwitch; wait;

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO
        
    printf "\n=== DONE RUNNING INVERSION ===\n"
    inversion_end=$(date +%s)
}

# Description: Run visualization notebooks and export to html
# Usage:
#   run_notebooks
run_notebooks() {
    config_path=${InversionPath}/${ConfigFile}
    printf "\n=== RUNNING VISUALIZATION NOTEBOOKS ===\n"
    if "$KalmanMode"; then
        cd ${RunDirs}/kf_inversions/period${period_i}
    else
        cd ${RunDirs}/inversion
    fi
    # replace config file path in viz notebook
    sed -i 's|\/home\/ubuntu\/integrated_methane_inversion\/config.yml|'$config_path'|g' visualization_notebook.ipynb
    jupyter nbconvert --execute --to html visualization_notebook.ipynb
    printf "\n=== DONE RUNNING NOTEBOOKS ===\n"
}

#!/bin/bash

# Functions available in this file include:
#   - setup_posterior
#   - run_posterior

# Description: Setup posterior GCClassic run directory
# Usage:
#   setup_posterior
setup_posterior() {
    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml"
        exit 9999
    fi

    printf "\n=== CREATING POSTERIOR RUN DIRECTORY ===\n"

    cd ${RunDirs}

    # Define the run directory name
    PosteriorName="${RunName}_Posterior"

    # Make the directory
    runDir="posterior_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/* ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable
    ln -s ../GEOSChem_build/gcclassic .

    # Link to restart file
    RestartFileFromSpinup=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        ln -s $RestartFileFromSpinup Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
    else
        RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
        ln -s $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
        if "$UseBCsForRestart"; then
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n"
        fi
    fi

    # Update settings in HEMCO_Config.rc
    if "$LognormalErrors"; then
        gridded_posterior_filename="gridded_posterior_ln.nc"
    else
        gridded_posterior_filename="gridded_posterior.nc"
    fi

    # Apply posterior scale factors to emissions
    # But exclude soil absorption from the application of the scale factors
    # Read in soil absorption separately with MeMo_SOIL_ABSORPTION
    sed -i -e "s|\.\./\.\.|\.\.|g" \
        -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" \
        -e "s|--> Emis_PosteriorSF       :       false|--> Emis_PosteriorSF       :       true|g" \
        -e "s|--> UseTotalPriorEmis      :       false|--> UseTotalPriorEmis      :       true|g" \
        -e "/(((MeMo_SOIL_ABSORPTION/i ))).not.UseTotalPriorEmis" \
        -e "/)))MeMo_SOIL_ABSORPTION/a (((.not.UseTotalPriorEmis" \
        -e "s|gridded_posterior.nc|${RunDirs}/inversion/${gridded_posterior_filename}|g" \
        -e "s|GFED                   : on|GFED                   : off|g" HEMCO_Config.rc

    # Turn on LevelEdgeDiags output
    # Output daily restarts to avoid trouble at month boundaries
    if "$HourlyCH4"; then
        sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
            -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
            -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
            -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' \
            -e 's/Restart.frequency:          '\''End'\''/Restart.frequency:          00000001 000000/g' \
            -e 's/Restart.duration:           '\''End'\''/Restart.duration:           00000001 000000/g' HISTORY.rc
    fi

    ### Turn on observation operators if requested, for posterior run
    activate_observations

    # Create run script from template
    sed -e "s:namename:${PosteriorName}:g" \
        -e "s:##:#:g" ch4_run.template >${PosteriorName}.run
    chmod 755 ${PosteriorName}.run
    rm -f ch4_run.template

    ### Perform dry run if requested
    if "$PosteriorDryRun"; then
        printf "\nExecuting dry-run for posterior run...\n"
        ./gcclassic --dryrun &>log.dryrun
        # prevent restart file from getting downloaded since
        # we don't want to overwrite the one we link to above
        sed -i '/GEOSChem.Restart/d' log.dryrun
        ./download_data.py log.dryrun aws
    fi

    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING POSTERIOR RUN DIRECTORY ===\n"
}

# Description: Run posterior simulation and process output
# Usage:
#   run_posterior
run_posterior() {
    posterior_start=$(date +%s)
    cd ${RunDirs}/posterior_run

    if $LognormalErrors; then
        inversion_result_filename="inversion_result_ln.nc"
    else
        inversion_result_filename="inversion_result.nc"
    fi

    if "$OptimizeBCs"; then
        if "$KalmanMode"; then
            inv_result_path="${RunDirs}/kf_inversions/period${period_i}/${inversion_result_filename}"
        else
            inv_result_path="${RunDirs}/inversion/${inversion_result_filename}"
        fi
        # set BC optimal delta values
        PerturbBCValues=$(generate_optimized_BC_values $inv_result_path)
        # add BC optimization delta to boundary condition edges
        sed -i -e "s|CH4_boundary_condition_ppb_increase_NSEW:.*|CH4_boundary_condition_ppb_increase_NSEW: ${PerturbBCValues}|g" \
            -e "s|perturb_CH4_boundary_conditions: false|perturb_CH4_boundary_conditions: true|g" geoschem_config.yml

        printf "\n=== BC OPTIMIZATION: BC optimized perturbation values for NSEW set to: ${PerturbBCValues} ===\n"
    fi

    if "$OptimizeOH"; then
        if "$KalmanMode"; then
            inv_result_path="${RunDirs}/kf_inversions/period${period_i}/${inversion_result_filename}"
        else
            inv_result_path="${RunDirs}/inversion/${inversion_result_filename}"
        fi
        # set OH optimal delta values
        PerturbOHValue=$(generate_optimized_OH_value $inv_result_path)
        # add OH optimization delta to boundary condition edges
        sed -i -e "s| OH_pert_factor  1.0| OH_pert_factor  ${PerturbOHValue}|g" HEMCO_Config.rc
        printf "\n=== OH OPTIMIZATION: OH optimized perturbation value set to: ${PerturbOHValue} ===\n"
    fi

    # Submit job to job scheduler
    printf "\n=== SUBMITTING POSTERIOR SIMULATION ===\n"
    sbatch --mem $RequestedMemory \
        -c $RequestedCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -W ${RunName}_Posterior.run
    wait

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    printf "\n=== DONE POSTERIOR SIMULATION ===\n"
    if "$KalmanMode"; then
        cd ${RunDirs}/kf_inversions/period${period_i}
        if ((period_i == 1)); then
            PrevDir="${RunDirs}/spinup_run"
        else
            PrevDir="${RunDirs}/posterior_run"
        fi
    else
        StartDate_i=$StartDate
        EndDate_i=$EndDate
        cd ${RunDirs}/inversion
        PrevDir="${RunDirs}/spinup_run"
    fi

    # Fill missing data (first hour of simulation) in posterior output
    PosteriorRunDir="${RunDirs}/posterior_run"
    printf "\n=== Calling postproc_diags.py for posterior ===\n"
    python ${InversionPath}/src/inversion_scripts/postproc_diags.py $RunName $PosteriorRunDir $PrevDir $StartDate_i $Res
    wait
    printf "\n=== DONE -- postproc_diags.py ===\n"

    # Build directory for hourly posterior GEOS-Chem output data
    mkdir -p data_converted_posterior
    mkdir -p data_visualization_posterior
    mkdir -p data_geoschem_posterior
    GCsourcepth="${PosteriorRunDir}/OutputDir"
    GCDir="./data_geoschem_posterior"
    printf "\n=== Calling setup_gc_cache.py for posterior ===\n"
    python ${InversionPath}/src/inversion_scripts/setup_gc_cache.py $StartDate_i $EndDate_i $GCsourcepth $GCDir
    wait
    printf "\n=== DONE -- setup_gc_cache.py ===\n"

    # Sample GEOS-Chem atmosphere with TROPOMI
    LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
    LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
    LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
    LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)
    nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc)
    if "$OptimizeBCs"; then
        nElements=$((nElements + 4))
    fi
    if "$OptimizeOH"; then
        nElements=$((nElements + 1))
    fi
    FetchTROPOMI="False"
    isPost="True"
    buildJacobian="False"

    printf "\n=== Calling jacobian.py to sample posterior simulation (without jacobian sensitivity analysis) ===\n"
    python ${InversionPath}/src/inversion_scripts/jacobian.py $StartDate_i $EndDate_i $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $BlendedTROPOMI $isPost $buildJacobian False
    wait
    printf "\n=== DONE sampling the posterior simulation ===\n\n"
    posterior_end=$(date +%s)

    # convert vizualization notebooks to html
    run_notebooks
}

# Description: Generates the updated NSEW perturbation to apply to domain edge BCs
# Usage:
#   generate_optimized_BC_values <path-to-inversion-result> <bc-pert-value>
generate_optimized_BC_values() {
    if $OptimizeOH; then
        python -c "import sys; import xarray;\
       xhat = xarray.load_dataset(sys.argv[1])['xhat'].values[-5:-1];\
       print(xhat.tolist())" $1
    else
        python -c "import sys; import xarray;\
       xhat = xarray.load_dataset(sys.argv[1])['xhat'].values[-4:];\
       print(xhat.tolist())" $1
    fi
}

# Description: Generates the updated perturbation to apply to OH
# Usage:
#   generate_optimized_OH_values <path-to-inversion-result> <oh-pert-value>
generate_optimized_OH_value() {
    python -c "import sys; import xarray;\
    xhat = xarray.load_dataset(sys.argv[1])['xhat'].values[-1:];\
    print(xhat.tolist()[0])" $1
}

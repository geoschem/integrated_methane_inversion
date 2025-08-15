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
    if "$UseGCHP"; then
        RunDuration=$(get_run_duration "$StartDate" "$EndDate")
        sed -i -e "s/Run_Duration=\"[0-9]\{8\} 000000\"/Run_Duration=\"${RunDuration} 000000\"/" \
            -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
            -e "s/^TOTAL_CORES=.*/TOTAL_CORES=${TOTAL_CORES}/" \
            -e "s/^NUM_NODES=.*/NUM_NODES=${NUM_NODES}/" \
            -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=${NUM_CORES_PER_NODE}/" \
            setCommonRunSettings.sh
        ln -nsf ../GEOSChem_build/gchp .
    else
        ln -nsf ../GEOSChem_build/gcclassic .
    fi

    # Link to restart file
    if "$UseGCHP"; then
        RestartFileFromSpinup=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.c${CS_RES}.nc4
    else
        RestartFileFromSpinup=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    fi
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        RestartFile=$RestartFileFromSpinup
    else
        if "$UseBCsForRestart"; then
            if "$UseGCHP"; then
                # regrid restart file to GCHP resolution
                TROPOMIBC="${RestartFilePrefix}${StartDate}_0000z.nc4"
                TemplatePrefix="${RunDirs}/${runDir}/Restarts/GEOSChem.Restart.20190101_0000z"
                FilePrefix="GEOSChem.Restart.${StartDate}_0000z"
                cd "${RunDirs}/CS_grids"
                TROPOMIBC72="temp_tropomi-bc.nc4"
                python ${InversionPath}/src/utilities/regrid_vertgrid_47-to-72.py $TROPOMIBC $TROPOMIBC72
                regrid_tropomi-BC-restart_gcc2gchp ${TROPOMIBC72} ${TemplatePrefix} ${FilePrefix} ${CS_RES} ${STRETCH_GRID} ${STRETCH_FACTOR} ${TARGET_LAT} ${TARGET_LON}
                RestartFile="${RunDirs}/CS_grids/${FilePrefix}.c${CS_RES}.nc4"
                cd $runDir
            else
                RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
                sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
                printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n"
            fi
        fi
    fi
    if "$UseGCHP"; then
        ln -nsf $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.c${CS_RES}.nc4
    else
        ln -nsf $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
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
        -e "s|--> Emis_PosteriorSF       :       false|--> Emis_PosteriorSF       :       true|g" \
        -e "s|--> UseTotalPriorEmis      :       false|--> UseTotalPriorEmis      :       true|g" \
        -e "s|gridded_posterior.nc|${RunDirs}/inversion/${gridded_posterior_filename}|g" \
        -e "s|GFED                   : on|GFED                   : off|g" HEMCO_Config.rc

    if "$UseGCHP"; then
        sed -i -e "s|^#EMIS_SF|EMIS_SF|g" ExtData.rc
        sed -i -e "s|\.\./\.\.|\.\.|g" \
            -e "s|gridded_posterior.nc|${RunDirs}/inversion/${gridded_posterior_filename}|g" ExtData.rc
    fi

    if [ "$OptimizeSoil" != true ]; then
        sed -i -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" \
            -e "/(((MeMo_SOIL_ABSORPTION/i ))).not.UseTotalPriorEmis" \
            -e "/)))MeMo_SOIL_ABSORPTION/a (((.not.UseTotalPriorEmis" \
            HEMCO_Config.rc
        if "$UseGCHP"; then
            sed -i -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" \
                ExtData.rc
        fi
    fi
    # Turn on LevelEdgeDiags output
    # Output daily restarts to avoid trouble at month boundaries
    if "$HourlyCH4"; then
        if "$UseGCHP"; then
            sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                -e 's/LevelEdgeDiags.frequency:.*/LevelEdgeDiags.frequency:      010000/g' \
                -e 's/LevelEdgeDiags.duration:.*/LevelEdgeDiags.duration:       240000/g' \
                HISTORY.rc
            # turn on monthly checkpoint
            sed -i -e 's/^Midrun_Checkpoint=OFF/Midrun_Checkpoint=ON/' \
                -e 's/^Midrun_Checkpoint=.*/Midrun_Checkpoint=monthly/' \
                setCommonRunSettings.sh

        else
            sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
                -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
                -e 's/Restart.frequency:          '\''End'\''/Restart.frequency:          00000001 000000/g' \
                -e 's/Restart.duration:           '\''End'\''/Restart.duration:           00000001 000000/g' HISTORY.rc
        fi
    fi

    ### Turn on observation operators if requested, for posterior run
    activate_observations

    # Create run script from template
    if "$UseGCHP"; then
        sed -e "s:namename:${PosteriorName}:g" \
            -e "s:##:#:g" gchp_ch4_run.template >${PosteriorName}.run
        chmod 755 ${PosteriorName}.run
        rm -f gchp_ch4_run.template
    else
        sed -e "s:namename:${PosteriorName}:g" \
            -e "s:##:#:g" ch4_run.template >${PosteriorName}.run
        chmod 755 ${PosteriorName}.run
        rm -f ch4_run.template
    fi

    ### Perform dry run if requested
    if [ "$UseGCHP" != "true" ]; then
        if "$PosteriorDryRun"; then
            printf "\nExecuting dry-run for posterior run...\n"
            ./gcclassic --dryrun &>log.dryrun
            # prevent restart file from getting downloaded since
            # we don't want to overwrite the one we link to above
            sed -i '/GEOSChem.Restart/d' log.dryrun
            python download_gc_data.py log.dryrun aws
        fi
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

    printf "\n=== SETTING UP POSTERIOR OPTIMIZATION ===\n"

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

        printf "\n--- BC OPTIMIZATION ---\n"
        printf "BC optimized perturbation values for NSEW set to: ${PerturbBCValues}\n"
    fi

    if "$OptimizeOH"; then
        if "$KalmanMode"; then
            inv_result_path="${RunDirs}/kf_inversions/period${period_i}/${inversion_result_filename}"
        else
            inv_result_path="${RunDirs}/inversion/${inversion_result_filename}"
        fi

        # set OH optimal delta values
        PerturbOHValue=$(generate_optimized_OH_value $inv_result_path)

        printf "\n=== OH OPTIMIZATION ===\n"
        if "$isRegional"; then
            # Apply single OH scale factor to entire region
            sed -i -e "s| OH_pert_factor  1.0| OH_pert_factor  ${PerturbOHValue}|g" HEMCO_Config.rc
            printf "OH optimized perturbation value set to: ${PerturbOHValue}\n"
        else
            # Apply hemispheric OH perturbation values using mask file
            Output_fpath="./gridded_posterior_oh_scale.nc"
            oh_sfs=($PerturbOHValue)
            Hemis_mask_fpath="${DataPath}/HEMCO/MASKS/v2024-08/hemisphere_mask.01x01.nc"
            OptimizeNorth='True'
            OptimizeSouth='True'
            gridded_optimized_OH ${oh_sfs[0]} ${oh_sfs[1]} $Hemis_mask_fpath $Output_fpath $OptimizeNorth $OptimizeSouth $STRETCH_GRID $STRETCH_FACTOR $TARGET_LAT $TARGET_LON
            
            # Modify OH scale factor in HEMCO config
            sed -i -e "s| OH_pert_factor  1.0 - - - xy 1 1| OH_pert_factor ${Output_fpath} oh_scale 2000\/1\/1\/0 C xy 1 1|g" HEMCO_Config.rc
            
            if "$UseGCHP"; then
                # add entry in ExtData.rc for GCHP
                sed -i -e "s|^#OH_pert_factor.*|OH_pert_factor 1 N Y - none none oh_scale ${Output_fpath}|" ExtData.rc
            fi

            printf "OH optimized perturbation values set to:\n"
            printf " ${oh_sfs[0]} for Northern Hemisphere\n"
            printf " ${oh_sfs[1]} for Southern Hemisphere\n"
        fi

    fi

    # Submit job to job scheduler
    printf "\n=== SUBMITTING POSTERIOR SIMULATION ===\n"
    if "$UseGCHP"; then
        sbatch --mem $RequestedMemory \
            -N $NUM_NODES \
            -n $TOTAL_CORES \
            -t $RequestedTime \
            -p $SchedulerPartition \
            -W ${RunName}_Posterior.run
    else
        sbatch --mem $RequestedMemory \
            -c $RequestedCPUs \
            -t $RequestedTime \
            -p $SchedulerPartition \
            -W ${RunName}_Posterior.run
    fi
    wait

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    printf "\n=== DONE POSTERIOR SIMULATION ===\n"

    if "$KalmanMode"; then
        cd ${RunDirs}/kf_inversions/period${period_i}
    else
        StartDate_i=$StartDate
        EndDate_i=$EndDate
        cd ${RunDirs}/inversion
    fi

    PosteriorRunDir="${RunDirs}/posterior_run"
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
    if "$UseGCHP"; then
        LonMinInvDomain=-180
        LonMaxInvDomain=180
        LatMinInvDomain=-90
        LatMaxInvDomain=90
    else
        LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
        LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
        LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
        LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)
    fi
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
    # fill kf_period with dummy number here
    kf_period=1

    printf "\n=== Calling jacobian.py to sample posterior simulation (without jacobian sensitivity analysis) ===\n"
    python ${InversionPath}/src/inversion_scripts/jacobian.py ${ConfigPath} $StartDate_i $EndDate_i $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $BlendedTROPOMI $UseWaterObs $isPost $kf_period $buildJacobian False
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
        if $isRegional; then
            python -c "import sys; import xarray;\
            xhat = xarray.load_dataset(sys.argv[1])['xhat'].values[-5:-1];\
            print(xhat.tolist())" $1
        else
            python -c "import sys; import xarray;\
            xhat = xarray.load_dataset(sys.argv[1])['xhat'].values[-6:-2];\
            print(xhat.tolist())" $1
        fi
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
    if $isRegional; then
        python -c "import sys; import xarray;\
        xhat = xarray.load_dataset(sys.argv[1])['xhat'].values[-1:];\
        print(xhat.tolist()[0])" $1
    else
        python -c "import sys; import xarray;\
        xhat = xarray.load_dataset(sys.argv[1])['xhat'].values[-2:];\
        print(xhat.tolist()[0], ' ', xhat.tolist()[1])" $1
    fi
}

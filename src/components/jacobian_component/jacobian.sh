#!/bin/bash

# Functions available in this file include:
#   - setup_jacobian
#   - run_jacobian
#   - create_simulation_dir
#   - generate_BC_perturb_values

# Description: Setup jacobian run directory
# Usage:
#   setup_jacobian
setup_jacobian() {
    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml"
        exit 9999
    fi
    printf "\n=== CREATING JACOBIAN RUN DIRECTORIES ===\n"

    cd ${RunDirs}

    # make dir for jacobian ics/bcs
    mkdir -p jacobian_lowbg_ics_bcs/Restarts
    if $isRegional; then
        mkdir -p jacobian_lowbg_ics_bcs/BCs
        OrigBCFile=${fullBCpath}/GEOSChem.BoundaryConditions.${StartDate}_0000z.nc4
        python ${InversionPath}/src/components/jacobian_component/make_jacobian_icbc.py $OrigBCFile ${RunDirs}/jacobian_lowbg_ics_bcs/BCs $StartDate $Species
    fi

    if [[ $PerturbationType = "eigenvector" ]] || [[ $PerturbationType = "grid" ]]; then
        # create lowbg restart file (normally this is done when running the
        # Jacobian for Kalman filter reasons, but when using eigenvectors,
        # the prior run is run first before the Jacobian directories.)
        OrigRestartFile=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
        if test -f "$OrigRestartFile" || "$DoSpinup"; then
            OrigRestartFile=$OrigRestartFile
        else
            OrigRestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
        fi
        python ${InversionPath}/src/components/jacobian_component/make_jacobian_icbc.py $OrigRestartFile ${RunDirs}/jacobian_lowbg_ics_bcs/Restarts $StartDate $Species
        cd jacobian_lowbg_ics_bcs/Restarts/
        if [ -f GEOSChem.BoundaryConditions.lowbg.${StartDate}_0000z.nc4 ]; then
            mv GEOSChem.BoundaryConditions.lowbg.${StartDate}_0000z.nc4 GEOSChem.Restart.lowbg.${StartDate}_0000z.nc4
        fi
        cd ../..
    fi

    # Create directory that will contain all Jacobian run directories
    mkdir -p -v jacobian_runs_${PerturbationType}

    if [ $NumJacobianTracers -gt 1 ]; then
        nRuns=$(calculate_num_jacobian_runs $NumJacobianTracers $nElements $OptimizeBCs $OptimizeOH)

        # Determine approx. number of CH4 tracers per Jacobian run
        printf "\nCombining Jacobian runs: Generating $nRuns run directories with approx. $NumJacobianTracers CH4 tracers (representing state vector elements) per run\n"
    else
        nRuns=$nElements
    fi

    # Copy run scripts
    cp ${InversionPath}/src/geoschem_run_scripts/submit_jacobian_simulations_array.sh jacobian_runs_${PerturbationType}/
    sed -i -e "s:{START}:0:g" \
        -e "s:{END}:${nRuns}:g" \
        -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs_${PerturbationType}/submit_jacobian_simulations_array.sh


    if [ $MaxSimultaneousRuns -gt 0 ]; then
        # Error check
        if [ $MaxSimultaneousRuns -gt $nRuns ]; then
            printf "\MaxSimultaneousRuns=${MaxSimultaneousRuns} is greater than the total runs=${nRuns}. Please modify MaxSimultenaousRuns in config.yml"
            exit 9999
        fi
        sed -i -e "s:{JOBS}:%${MaxSimultaneousRuns}:g" jacobian_runs_${PerturbationType}/submit_jacobian_simulations_array.sh
    else
        sed -i -e "s:{JOBS}::g" jacobian_runs_${PerturbationType}/submit_jacobian_simulations_array.sh
    fi


    cp ${InversionPath}/src/geoschem_run_scripts/run_prior_simulation.sh jacobian_runs_${PerturbationType}/
    sed -i -e "s:{RunName}:${RunName}:g" \
        -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs_${PerturbationType}/run_prior_simulation.sh
    cp ${InversionPath}/src/geoschem_run_scripts/run_bkgd_simulation.sh jacobian_runs_${PerturbationType}/
    sed -i -e "s:{RunName}:${RunName}:g" \
        -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs_${PerturbationType}/run_bkgd_simulation.sh

    # Initialize (x=0 is base run, i.e. no perturbation; x=1 is state vector element=1; etc.)
    x=0

    # Create jacobian run directories
    while [ $x -le $nRuns ]; do

        # Current state vector element
        xUSE=$x

        # Add zeros to string name
        if [ $x -lt 10 ]; then
            xstr="000${x}"
        elif [ $x -lt 100 ]; then
            xstr="00${x}"
        elif [ $x -lt 1000 ]; then
            xstr="0${x}"
        else
            xstr="${x}"
        fi
        create_simulation_dir

        # Increment
        x=$(($x + 1))
    done

    if "$LognormalErrors"; then
        x="background"
        xstr=$x
        create_simulation_dir
    fi

    # HN: We don't do this because I don't use multiplicative scale factors
    # if [[ $Species = "CO2" ]]; then
    #     set -e
    #     # update perturbation values before running jacobian simulations
    #     printf "\n=== UPDATING PERTURBATION SFs ===\n"
    #     python ${InversionPath}/src/components/jacobian_component/make_perturbation_sf.py $ConfigPath 1 $PerturbValue
    # fi

    printf "\n=== DONE CREATING JACOBIAN RUN DIRECTORIES ===\n"
}

# Description: Create simulation directory for defined xstr
# Usage:
#   create_simulation_dir
create_simulation_dir() {
    # Define the run directory name
    name="${RunName}_${xstr}"

    # Make the directory
    runDir="./jacobian_runs_${PerturbationType}/${name}"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/* ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each rundir
    ln -s ../../GEOSChem_build/gcclassic .

    # If it's the prior run, save out the satellite data 
    working_dir=$(pwd)
    in="MODEL_LEVEL_EDGE_DIR: '/nobackupnfs1/hnesser/CO2_inversion/jacobian_runs_grid//OutputDir'"
    out="MODEL_LEVEL_EDGE_DIR: '/nobackupnfs1/hnesser/CO2_inversion/jacobian_runs_${PerturbationType}/CO2_inversion_0000/OutputDir'"
    sed -i -e "s|${in}|${out}|g" config_satellite_operator.yaml
    # in="MODEL_CONCENTRATION_DIR: 'OutputDir'"
    # out="MODEL_CONCENTRATION_DIR: '${working_dir}/OutputDir'"
    # sed -i -e "s|${in}|${out}|g" config_satellite_operator.yaml
    if [[ $x -eq 0 ]]; then
        sed -i -e "s|SAVE_SATELLITE_DATA: 'False'|SAVE_SATELLITE_DATA: 'True'|g" \
            -e "s|SAVE_INTERPOLATION: 'False'|SAVE_INTERPOLATION: 'True'|g" \
        config_satellite_operator.yaml
    fi

    # link to restart file
    RestartFileFromSpinup=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        RestartFile=$RestartFileFromSpinup
    else
        RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
        if "$UseBCsForRestart"; then
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
        fi
    fi
    ln -s $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4

    # Modify HEMCO_Config.rc to turn off individual emission inventories
    # and use total emissions (without soil absorption) saved out from prior
    # emissions simulation instead. For the prior and OH sims we add soil
    # absorption back in below
    # If CO2, don't use total prior emissions because of a HEMCO bug. Also,
    # we turn off the AnalyticalInversion even though really this should be 
    # a flag that depends on whether eigenvectors are being used.
    if [[ $Species = "CO2" ]]; then
        BoolFlag='false'
    else
        printf "\nTurning on use of total prior emissions in HEMCO_Config.rc.\n"
        BoolFlag='true'
    fi
    sed -i -e "s|UseTotalPriorEmis      :       false|UseTotalPriorEmis      :       ${BoolFlag}|g" \
        -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" \
        -e "s|GFED                   : on|GFED                   : off|g" HEMCO_Config.rc
    # HN: Removing this because we removed state vector specification in 
    # HEMCO_Config one way or another, and otherwise the AnalyticalInversion
    # flag isn't a problem
    # if $PerturbEigenvectors; then
    #     sed -i -e "s|AnalyticalInversion    :       false|AnalyticalInversion    :       ${BoolFlag}|g" HEMCO_Config.rc
    # fi

    # We no longer need this, either
    # if [[ $x -gt 0 ]]; then
    #     sed -i -e "s|UseEigenvectors        :       false|UseEigenvectors        :       ${PerturbEigenvectors}|g" \
    #     HEMCO_Config.rc
    # fi

    # Add in Obspack
    if $ObsPack; then
        sed -i -e "s|\./obspack_co2_1_OCO2MIP_2018-11-28\.YYYYMMDD\.nc|$DataPathObsPack|g" geoschem_config.yml
        sed -i '/obspack:/!b;n;s/activate: false/activate: true/' geoschem_config.yml
    fi

    # Determine which elements are BC perturbations
    BC_elem=false
    bcThreshold=$nElements
    if "$OptimizeBCs"; then
        if "$OptimizeOH"; then
            bcThreshold=$(($nElements - 5))
        else
            bcThreshold=$(($nElements - 4))
        fi
    fi

    # Determine which element (if any) is an OH perturbation
    OH_elem=false
    ohThreshold=$nElements
    if "$OptimizeOH"; then
        ohThreshold=$(($nElements - 1))
    fi

    # Update settings in HISTORY.rc
    # Only save out hourly pressure fields to daily files for base run
    if [[ $x -eq 0 ]] || [[ "$x" = "background" ]]; then
        if "$HourlySpecies"; then
            sed -i -e 's/'\''Restart/#'\''Restart/g' \
                -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
                -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
                -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
        fi
    # For all other runs, just disable Restarts
    else
        if "$HourlySpecies"; then
            sed -i -e 's/'\''Restart/#'\''Restart/g' HISTORY.rc
        fi
    fi

    # for background simulation, disable the emissions
    # needed for lognormal error inversion
    if [ "$x" = "background" ]; then
        sed -i -e 's/EMISSIONS              :       true/EMISSIONS              :       false/g' \
            -e 's/GFED                   : on    CH4/GFED                   : off    CH4/g' HEMCO_Config.rc
    fi

    # Create run script from template
    sed -e "s:namename:${name}:g" run.template >${name}.run
    if [[ $x -eq 0 ]]; then
        sed -i -e "s|cleanup=true|cleanup=false|g" \
            -e "s|##|#|g" ${name}.run
        sed -i "/grid_perturbation_outputs/d" ${name}.run
    else
        sed -i "/tccon_operator/d" ${name}.run
    fi

    rm -f run.template
    chmod 755 ${name}.run

    ### Turn on observation operators if requested, only for base run
    if [[ $x -eq 0 ]] || [[ "$x" = "background" ]]; then
        activate_observations
    fi

    # Turn off sectoral emissions diagnostics since total emissions are
    # read in for jacobian runs
    sed -i -e "s:EmisCH4:#EmisCH4:g" HEMCO_Diagn.rc
    sed -i -e "s:#EmisCH4_Total:EmisCH4_Total:g" HEMCO_Diagn.rc
    sed -i -e "s:#EmisCH4_SoilAbsorb:EmisCH4_SoilAbsorb:g" HEMCO_Diagn.rc

    if is_number "$x"; then
        ### Perform dry run if requested, only for base run
        if [[ $x -eq 0 ]]; then
            if "$ProductionDryRun"; then
                printf "\nExecuting dry-run for production runs...\n"
                ./gcclassic --dryrun &>log.dryrun
                # prevent restart file from getting downloaded since
                # we don't want to overwrite the one we link to above
                sed -i '/GEOSChem.Restart/d' log.dryrun
                ./download_data.py log.dryrun aws
            fi
        fi

        # Determine start and end element numbers for this run directory
        if [[ $x -eq 0 ]]; then
            # if using 1 tracer per simulation. Or is the prior simulation.
            start_element=$x
            end_element=$x
        else
            start_element=$((end_element + 1))
            # calculate tracer end based on the number of tracers and bc/oh thresholds
            # Note: the prior simulation, BC simulations, and OH simulation get their
            # own dedicated simulation, so end_element is the same as start_element
            end_element=$(calculate_tracer_end $start_element $nElements $NumJacobianTracers $bcThreshold $ohThreshold)
        fi

        # Perturb OH if this is the OH perturbations simulation
        if [ $start_element -gt $ohThreshold ]; then
            OH_elem=true
            sed -i -e "s| OH_pert_factor  1.0| OH_pert_factor  ${PerturbValueOH}|g" HEMCO_Config.rc
        fi

        # If the current state vector element is one of the BC state vector elements, then
        # turn on BC optimization for the corresponding edge
        if [[ $start_element -gt $bcThreshold ]] && [[ "$OH_elem" = false ]]; then
            BC_elem=true
            PerturbBCValues=$(generate_BC_perturb_values $bcThreshold $start_element $PerturbValueBCs)
            sed -i -e "s|CH4_boundary_condition_ppb_increase_NSEW:.*|CH4_boundary_condition_ppb_increase_NSEW: ${PerturbBCValues}|g" \
                -e "s|perturb_CH4_boundary_conditions: false|perturb_CH4_boundary_conditions: true|g" geoschem_config.yml
        fi
    fi

    # the prior, OH perturbation, and background simulations need to have soil absorption
    # and, in the case, of kalman mode the prior is scaled by the nudged scale factors
    if [[ $x -eq 0 ]] || [[ "$x" = "background" ]] || [[ $OH_elem = true ]]; then
        # Use MeMo soil absorption for the prior simulation
        sed -i -e "/(((MeMo_SOIL_ABSORPTION/i ))).not.UseTotalPriorEmis" \
            -e "/)))MeMo_SOIL_ABSORPTION/a (((.not.UseTotalPriorEmis" HEMCO_Config.rc
        if "$KalmanMode"; then
            # Use nudged scale factors for the prior simulation and OH simulation for kalman mode
            sed -i -e "s|--> Emis_PosteriorSF       :       false|--> Emis_PosteriorSF       :       true|g" \
                -e "s|--> UseTotalPriorEmis      :       false|--> UseTotalPriorEmis      :       true|g" \
                -e "s|gridded_posterior.nc|${RunDirs}/ScaleFactors.nc|g" HEMCO_Config.rc
        fi

    else
        # set lowbg boundary conditions and restarts for all other perturbation simulations
        # Note that we use the timecycle flag C to avoid having to make additional files
        RestartFile=${RunDirs}/jacobian_lowbg_ics_bcs/Restarts/GEOSChem.Restart.lowbg.${StartDate}_0000z.nc4
        BCFilelowbg=${RunDirs}/jacobian_lowbg_ics_bcs/BCs/GEOSChem.BoundaryConditions.lowbg.${StartDate}_0000z.nc4
        BCSettingslowbg="SpeciesBC_${Species}  1980-2024/1-12/1-31/* C xyz 1 ${Species} - 1 1"
        sed -i -e "s|.*GEOSChem\.BoundaryConditions.*|\* BC_${Species} ${BCFilelowbg} ${BCSettingslowbg}|g" HEMCO_Config.rc
        # create symlink to lowbg restart file
        ln -sf $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
        # Also, set emissions to zero for default tracer by applying new ZERO scale factor
        sed -i -e "/1 NEGATIVE -1.0 - - - xy 1 1/a 5 ZERO 0.0 - - - xy 1 1" \
            -e "s|CH4 - 1 500|CH4 5 1 500|g" \
            -e "s|CO2 - \([1-3]\) 500|CO2 5 \1 500|g" HEMCO_Config.rc 
        # HON this works for totalprioremis and for HEMCO defined emissions
    fi

    # Modify restart and BC entries in HEMCO_Config.rc to look for CH4/CO2 only
    # instead of all advected species
    sed -i -e "s/SPC_/SPC_${Species}/g" -e "s/?ALL?/${Species}/g" -e "s/EFYO xyz 1 \*/EFYO xyz 1 ${Species}/g" HEMCO_Config.rc
    sed -i -e "s/BC_ /BC_${Species} /g" -e "s/?ADV?/${Species}/g" -e "s/EFY xyz 1 \*/EFY xyz 1 ${Species}/g" HEMCO_Config.rc
    
    # Initialize previous lines to search
    GcPrevLine='- '${Species}
    HcoPrevLine1='EFYO xyz 1 '${Species}' - 1 '
    HcoPrevLine2='(((CMS_FLUX'
    # # if ! $PerturbEigenvectors; then
    # if [[ $PerturbationType = "grid" ]]; then
    #     HcoPrevLine2='0 FOSSIL '
    # elif [[ $PerturbationType = "eigenvector" ]]; then
    #     HcoPrevLine2='${Species} 5 1 500'
    # else
    #     printf "PerturbationType $PerturbationType not supported"
    #     exit 1
    # fi
    # HcoPrevLine3='Perturbations.txt - - - xy count 1'
    HcoPrevLine4='\* BC_'${Species}
    PertPrevLine='DEFAULT    0     0.0'

    # Loop over state vector element numbers for this run and add each element
    # as a $Species tracer in the configuraton files
    if is_number "$x"; then
        if [ $x -gt 0 ] && [ "$BC_elem" = false ] && [ "$OH_elem" = false ]; then
            i=0
            add_new_tracer
            for i in $(seq $start_element $end_element); do
                add_new_tracer
            done
        fi
        sed -i -e "s|CO2_Emis_Prior_0000|CO2_Emis_Prior|g" HEMCO_Config.rc
    fi

    # Navigate back to top-level directory
    cd ../..
}

# Description: Add new tracers to a simulation
# Usage: add_new_tracer
add_new_tracer() {
    if [ $i -lt 10 ]; then
        istr="000${i}"
    elif [ $i -lt 100 ]; then
        istr="00${i}"
    elif [ $i -lt 1000 ]; then
        istr="0${i}"
    else
        istr="${i}"
    fi

    # Start HEMCO scale factor ID at 2000 to avoid conflicts with
    # preexisting scale factors/masks
    # if ! ${PerturbEigenvectors}; then
    if [[ $PerturbationType = "grid" ]]; then
        SFnum=$((2000 + i))
    else
        SFnum='-'
    fi

    # Add lines to geoschem_config.yml
    # Spacing in GcNewLine is intentional
    GcNewLine='\
      - '${Species}'_'$istr
    sed -i -e "/$GcPrevLine/a $GcNewLine" geoschem_config.yml
    GcPrevLine='- '${Species}'_'$istr

    # Add lines to species_database.yml
    if [ ${Species} = "CH4" ]; then
        SpcNextLine='CHBr3:'
    elif [ ${Species} = "CO2" ]; then
        SpcNextLine='CO_PROP'
    else
        printf "${Species} is not supported."
    fi

    if [[ $Species = "CH4" ]]; then
        bg_vv="1.8e-6"
        fullname="Methane"
    elif [[ $Species = "CO2" ]]; then
        bg_vv="3.55e-4"
        fullname="Carbon dioxide"
    fi
    SpcNewLines=${Species}'_'$istr':\n  << : *'${Species}'properties\n  Background_VV: '${bg_vv}'\n  FullName: '${fullname}''
    sed -i -e "s|$SpcNextLine|$SpcNewLines\n$SpcNextLine|g" species_database.yml

    HcoNewLine1='\
* SPC_'${Species}'_'$istr' - - - - - - '${Species}'_'$istr' - 1 1'
    sed -i -e "/$HcoPrevLine1/a $HcoNewLine1" HEMCO_Config.rc
    HcoPrevLine1='SPC_'${Species}'_'$istr

    # Add perturbation lines to HEMCO_Config.yml
    if [[ $PerturbationType = "grid" ]]; then
        pert_str='perturbations'
    elif [[ $PerturbationType = "eigenvector" ]]; then
        pert_str='eigenvectors0'
    fi

    HcoNewLine2='\
0 '${Species}'_Emis_Prior_'$istr' '${RunDirs}'/'${pert_str}'/perturbation_'$istr'.nc pert 1990-2025/1-12/1/0 C xy kg/m2/s '${Species}'_'$istr' - 1 500'
    sed -i "/$HcoPrevLine2/a $HcoNewLine2" HEMCO_Config.rc
    HcoPrevLine2=${Species}'_'$istr' - 1 500'

    if $isRegional; then
        HcoNewLine4='\
    * BC_'${Species}'_'$istr' - - - - - - '${Species}'_'$istr' - 1 1'
        sed -i -e "/$HcoPrevLine4/a $HcoNewLine4" HEMCO_Config.rc
        HcoPrevLine4='BC_'${Species}'_'$istr
    fi
}

# Description: Run jacobian simulations
# Usage:
#   run_jacobian
run_jacobian() {

    pushd ${RunDirs}

    # Copy run scripts
    # need to re-copy since config vars are
    # hardcoded and redojacobian might have changed
    cp ${InversionPath}/src/geoschem_run_scripts/run_jacobian_simulations.sh jacobian_runs_${Perturbation_Type}/
    sed -i -e "s:{RunName}:${RunName}:g" \
        -e "s:{InversionPath}:${InversionPath}:g" \
        -e "s:{KalmanMode}:${KalmanMode}:g" \
        -e "s:{EndDate}:${EndDate}:g" \
        -e "s:{ReDoJacobian}:${ReDoJacobian}:g" jacobian_runs_${Perturbation_Type}/run_jacobian_simulations.sh

    popd

    if ! "$PrecomputedJacobian"; then
        jacobian_start=$(date +%s)
        if "$KalmanMode"; then
            jacobian_period=${period_i}
        else
            jacobian_period=1
        fi

        set -e
        # update perturbation values before running jacobian simulations
        printf "\n=== UPDATING PERTURBATION SFs ===\n"
        python ${InversionPath}/src/components/jacobian_component/make_perturbation_sf.py $ConfigPath $jacobian_period $PerturbValue

        cd ${RunDirs}/jacobian_runs_${Perturbation_Type}

        # # create lowbg restart file
        OrigRestartFile=$(readlink ${RunName}_0000/Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4)
        python ${InversionPath}/src/components/jacobian_component/make_jacobian_icbc.py $OrigRestartFile ${RunDirs}/jacobian_lowbg_ics_bcs/Restarts $StartDate $Species
        cd ${RunDirs}/jacobian_lowbg_ics_bcs/Restarts/
        if [ -f GEOSChem.BoundaryConditions.lowbg.${StartDate}_0000z.nc4 ]; then
            mv GEOSChem.BoundaryConditions.lowbg.${StartDate}_0000z.nc4 GEOSChem.Restart.lowbg.${StartDate}_0000z.nc4
        fi
        cd ${RunDirs}/jacobian_runs_${Perturbation_Type}
        set +e

        printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"
        # Submit job to job scheduler
        source submit_jacobian_simulations_array.sh

        if "$LognormalErrors"; then
            submit_job $SchedulerType false run_bkgd_simulation.sh
            wait
        fi

        # check if any jacobians exited with non-zero exit code
        [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

        printf "\n=== DONE JACOBIAN SIMULATIONS ===\n"
        jacobian_end=$(date +%s)
    else
        # Add symlink pointing to jacobian matrix files from the reference
        # inversion w/ precomputed Jacobian
        if "$KalmanMode"; then
            cd ${RunDirs}/kf_inversions/period${period_i}
            precomputedJacobianCachePrefix=${ReferenceRunDir}/kf_inversions/period${period_i}
        else
            cd ${RunDirs}/inversion
            precomputedJacobianCachePrefix=${ReferenceRunDir}/inversion
        fi

        precomputedJacobianCache=${precomputedJacobianCachePrefix}/data_converted
        ln -s $precomputedJacobianCache data_converted_reference

        # Run the prior simulation
        cd ${JacobianRunsDir}

        # Submit prior simulation to job scheduler
        printf "\n=== SUBMITTING PRIOR SIMULATION ===\n"
        submit_job $SchedulerType true run_prior_simulation.sh
        printf "=== DONE PRIOR SIMULATION ===\n"

        # Run the background simulation if lognormal errors enabled
        if "$LognormalErrors"; then
            printf "\n=== SUBMITTING BACKGROUND SIMULATION ===\n"
            submit_job $SchedulerType false run_bkgd_simulation.sh
            printf "=== DONE BACKGROUND SIMULATION ===\n"
        fi

        # Get Jacobian scale factors
        python ${InversionPath}/src/inversion_scripts/get_jacobian_scalefactors.py $period_i $RunDirs $ReferenceRunDir $ConfigPath
        wait
        printf "Got Jacobian scale factors\n"
    fi
}

# Description: Print perturbation string for BC optimization
#   based on the current state vector element
#   Returns [float, float, float, float]
# Usage:
#   generate_BC_perturb_values <bcThreshold> <element-number> <pert-value>
generate_BC_perturb_values() {
    python -c "import sys;\
    bc_perturb = [0.0, 0.0, 0.0, 0.0];\
    bcThreshold = int(sys.argv[1]) + 1;\
    element = int(sys.argv[2]);\
    pert_index = element % bcThreshold;\
    bc_perturb[pert_index] = float(sys.argv[3]);\
    print(bc_perturb)" $1 $2 $3
}

# Description: Print end element for multitracer perturbation runs
#   based on the current starting element, number of tracers, and whether
#   it is an OH or BC perturbation run
#   Returns int
# Usage:
#   calculate_tracer_end <start-element> <n-elements> <number-tracers> <bcThreshold> <ohThreshold>
calculate_tracer_end() {
    python -c "
import sys
start_elem = int(sys.argv[1])
n_elems = int(sys.argv[2])
nTracers = int(sys.argv[3])
bcThreshold = int(sys.argv[4])
ohThreshold = int(sys.argv[5])
end_elem = start_elem + nTracers - 1
# Ensure end element is within bounds
if end_elem > n_elems:
    end_elem = n_elems
# If this is a BC or OH perturbation run, only perturb the current element
if start_elem > bcThreshold or start_elem > ohThreshold:
    end_elem = start_elem
else:
    while end_elem > bcThreshold or end_elem > ohThreshold:
        end_elem -= 1
print(end_elem)
" $1 $2 $3 $4 $5
}

# Description: Print number of jacobian runs for multitracer perturbation runs
#   based on the number of targeted tracers per simulation, number of state
#   vector elements, and whether OH and BC are optimized. Returns an int.
# Usage:
#   calculate_num_jacobian_runs <num-tracers> <number-elements> <bc-optimized> <oh-optimized>
calculate_num_jacobian_runs() {
    python -c "
import sys
import math
nTracers = int(sys.argv[1])
nElements = int(sys.argv[2])
bcOptimized = sys.argv[3].lower() == 'true'
ohOptimized = sys.argv[4].lower() == 'true'
numStandaloneRuns = 0
if bcOptimized:
    numStandaloneRuns += 4
if ohOptimized:
    numStandaloneRuns += 1
nRuns = math.ceil((nElements - numStandaloneRuns) / nTracers)
nRuns += numStandaloneRuns
print(nRuns)
" $1 $2 $3 $4
}

is_number() {
    local s="$1"
    [[ $s =~ ^[0-9]+$ ]]
}
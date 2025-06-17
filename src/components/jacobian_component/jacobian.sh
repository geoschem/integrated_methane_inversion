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
    mkdir -p jacobian_1ppb_ics_bcs/Restarts
    if "$isRegional"; then
        mkdir -p jacobian_1ppb_ics_bcs/BCs
        OrigBCFile="${fullBCpath}/GEOSChem.BoundaryConditions.${StartDate}_0000z.nc4"
        python ${InversionPath}/src/components/jacobian_component/make_jacobian_icbc.py $ConfigPath $OrigBCFile ${RunDirs}/jacobian_1ppb_ics_bcs/BCs $StartDate
    fi

    # create 1ppb restart file
    # link to restart file
    if "$UseGCHP"; then
        OrigRestartFile="${RunDirs}/CS_grids/GEOSChem.Restart.${StartDate}_0000z.c${CS_RES}.nc4"
        if [ ! -f "$OrigRestartFile" ]; then
            # regrid restart file to GCHP resolution
            TROPOMIBC=${RestartFilePrefix}${StartDate}_0000z.nc4
            Template="${RunDirs}/${runDir}/Restarts/GEOSChem.Restart.20190101_0000z.c${CS_RES}.nc4"
            FilePrefix="GEOSChem.Restart.${StartDate}_0000z"
            cd ${RunDirs}/CS_grids
            ./regrid_tropomi-BC-restart_gcc2gchp.sh ${TROPOMIBC} ${Template} ${FilePrefix} ${CS_RES}
            RestartFile="${RunDirs}/CS_grids/${FilePrefix}.c${CS_RES}.nc4"
            cd ${RunDirs}
        fi
    else
        OrigRestartFile="${RestartFilePrefix}${StartDate}_0000z.nc4"
    fi
    python ${InversionPath}/src/components/jacobian_component/make_jacobian_icbc.py $ConfigPath $OrigRestartFile ${RunDirs}/jacobian_1ppb_ics_bcs/Restarts $StartDate
    cd ${RunDirs}/jacobian_1ppb_ics_bcs/Restarts/
    if [ "$UseGCHP" != "true" ]; then
        if [ -f GEOSChem.BoundaryConditions.1ppb.${StartDate}_0000z.nc4 ]; then
            mv GEOSChem.BoundaryConditions.1ppb.${StartDate}_0000z.nc4 GEOSChem.Restart.1ppb.${StartDate}_0000z.nc4
            ncrename -v SpeciesBC_CH4,SpeciesRst_CH4 GEOSChem.Restart.1ppb.${StartDate}_0000z.nc4
        fi
    fi
    cd ${RunDirs}

    # Create directory that will contain all Jacobian run directories
    mkdir -p -v jacobian_runs

    if ! "$PrecomputedJacobian"; then
        if [ $NumJacobianTracers -gt 1 ]; then
            nRuns=$(calculate_num_jacobian_runs $NumJacobianTracers $nElements $OptimizeBCs $OptimizeOH $isRegional)

            # Determine approx. number of CH4 tracers per Jacobian run
            printf "\nCombining Jacobian runs: Generating $nRuns run directories with approx. $NumJacobianTracers CH4 tracers (representing state vector elements) per run\n"
        else
            nRuns=$nElements
        fi
    else 
        # only need to set up the prior run directory
        nRuns=0
    fi

    # Copy run scripts
    cp ${InversionPath}/src/geoschem_run_scripts/submit_jacobian_simulations_array.sh jacobian_runs/
    sed -i -e "s:{START}:0:g" \
        -e "s:{END}:${nRuns}:g" \
        -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/submit_jacobian_simulations_array.sh
    if [ $MaxSimultaneousRuns -gt 0 ]; then
        # Error check
        if [ $MaxSimultaneousRuns -gt $nRuns ]; then
            printf "\MaxSimultaneousRuns=${MaxSimultaneousRuns} is greater than the total runs=${nRuns}. Please modify MaxSimultenaousRuns in config.yml"
            exit 9999
        fi
        sed -i -e "s:{JOBS}:%${MaxSimultaneousRuns}:g" jacobian_runs/submit_jacobian_simulations_array.sh
    else
        sed -i -e "s:{JOBS}::g" jacobian_runs/submit_jacobian_simulations_array.sh
    fi
    cp ${InversionPath}/src/geoschem_run_scripts/run_prior_simulation.sh jacobian_runs/
    sed -i -e "s:{RunName}:${RunName}:g" \
        -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/run_prior_simulation.sh
    cp ${InversionPath}/src/geoschem_run_scripts/run_bkgd_simulation.sh jacobian_runs/
    sed -i -e "s:{RunName}:${RunName}:g" \
        -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/run_bkgd_simulation.sh

    if "$KalmanMode"; then
        jacobian_period=${period_i}
    else
        jacobian_period=1
    fi

    set -e
    # generate gridded perturbation values for all state vector elements
    printf "\n=== GENERATE GRIDDED PERTURBATION SFs ===\n"
    python ${InversionPath}/src/components/jacobian_component/make_perturbation_sf.py $ConfigPath $jacobian_period $PerturbValue
    printf "\n=== DONE GENERATE GRIDDED PERTURBATION SFs ===\n"

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

    printf "\n=== DONE CREATING JACOBIAN RUN DIRECTORIES ===\n"
}
# Description: Create simulation directory for defined xstr
# Usage:
#   create_simulation_dir
create_simulation_dir() {
    # Define the run directory name
    name="${RunName}_${xstr}"

    # Make the directory
    runDir="./jacobian_runs/${name}"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/* ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each rundir
    if "$UseGCHP"; then
        sed -i -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
            -e "s/^TOTAL_CORES=.*/TOTAL_CORES=${TOTAL_CORES}/" \
            -e "s/^NUM_NODES=.*/NUM_NODES=${NUM_NODES}/" \
            -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=${NUM_CORES_PER_NODE}/" \
            setCommonRunSettings.sh
        ln -nsf ../../GEOSChem_build/gchp .
    else
        ln -nsf ../../GEOSChem_build/gcclassic .
    fi

    # link to restart file
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
                Template="${RunDirs}/${runDir}/Restarts/GEOSChem.Restart.20190101_0000z.c${CS_RES}.nc4"
                FilePrefix="GEOSChem.Restart.${StartDate}_0000z"
                cd "${RunDirs}/CS_grids"
                ./regrid_tropomi-BC-restart_gcc2gchp.sh ${TROPOMIBC} ${Template} ${FilePrefix} ${CS_RES}
                RestartFile="${RunDirs}/CS_grids/${FilePrefix}.c${CS_RES}.nc4"
                cd "${RunDirs}/jacobian_runs/${name}"
            else
                RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
                sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            fi
        fi
    fi

    if "$UseGCHP"; then
        ln -nsf $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.c${CS_RES}.nc4
    else
        ln -nsf $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
    fi

    # Modify HEMCO_Config.rc to turn off individual emission inventories
    # and use total emissions (without soil absorption) saved out from prior
    # emissions simulation instead. For the prior and OH sims we add soil
    # absorption back in below
    sed -i -e "s|UseTotalPriorEmis      :       false|UseTotalPriorEmis      :       true|g" \
        -e "s|AnalyticalInversion    :       false|AnalyticalInversion    :       true|g" \
        -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" \
        -e "s|GFED                   : on|GFED                   : off|g" HEMCO_Config.rc

    if "$UseGCHP"; then
        sed -i -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" ExtData.rc
    fi
    # Determine which elements are BC perturbations
    BC_elem=false
    bcThreshold=$nElements
    if "$OptimizeBCs"; then
        if "$OptimizeOH"; then
            if "$isRegional"; then
                bcThreshold=$(($nElements - 5))
            else
                bcThreshold=$(($nElements - 6))
            fi
        else
            bcThreshold=$(($nElements - 4))
        fi
    fi

    # Determine which element (if any) is an OH perturbation
    OH_elem=false
    ohThreshold=$nElements
    if "$OptimizeOH"; then
        if "$isRegional"; then
            ohThreshold=$(($nElements - 1))
        else
            ohThreshold=$(($nElements - 2))
        fi
    fi

    # Update settings in HISTORY.rc
    # Only save out hourly pressure fields to daily files for base run
    if [[ $x -eq 0 ]] || [[ "$x" = "background" ]]; then
        if "$UseGCHP"; then
            sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                -e 's/LevelEdgeDiags.frequency:.*/LevelEdgeDiags.frequency:      010000/g' \
                -e 's/LevelEdgeDiags.duration:.*/LevelEdgeDiags.duration:       240000/g' HISTORY.rc
        else
            sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
                -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' HISTORY.rc
        fi
    fi
    # disable Restart for all runs
    if [ "$UseGCHP" != "true" ]; then
        if "$HourlyCH4"; then
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
    if "$UseGCHP"; then
        sed -e "s:namename:${name}:g" gchp_ch4_run.template >${name}.run
        rm -f gchp_ch4_run.template
    else
        sed -e "s:namename:${name}:g" ch4_run.template >${name}.run
        rm -f ch4_run.template
    fi
    chmod 755 ${name}.run

    ### Turn on observation operators if requested, only for base run
    if [[ $x -eq 0 ]] || [[ "$x" = "background" ]]; then
        activate_observations
    fi

    # Turn off emissions diagnostics to save disk space
    # These should remain unchanged from hemco_prior_emis
    sed -i -e "s:EmisCH4:#EmisCH4:g" HEMCO_Diagn.rc

    if is_number "$x"; then
        ### Perform dry run if requested, only for base run
        if [[ $x -eq 0 ]]; then
            if [ "$UseGCHP" != "true" ]; then
                if "$ProductionDryRun"; then
                    printf "\nExecuting dry-run for production runs...\n"
                    ./gcclassic --dryrun &>log.dryrun
                    # prevent restart file from getting downloaded since
                    # we don't want to overwrite the one we link to above
                    sed -i '/GEOSChem.Restart/d' log.dryrun
                    python download_gc_data.py log.dryrun aws
                fi
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

        # Perturb OH if this is an OH state vector element
        if [ $start_element -gt $ohThreshold ]; then
            OH_elem=true
            if "$isRegional"; then
                # Perturb OH everywhere if this is a reginoal simulation
                sed -i -e "s| OH_pert_factor  1.0| OH_pert_factor ${PerturbValueOH}|g" HEMCO_Config.rc
            else
                # Perturb OH by hemisphere if this is a global simulation
                # Apply hemispheric OH perturbation values using mask file
                Output_fpath="./gridded_perturbation_oh_scale.nc"
                Hemis_mask_fpath="${DataPath}/HEMCO/MASKS/v2024-08/hemisphere_mask.01x01.nc"
                OptimizeNorth='False'
                OptimizeSouth='False'
                if [ $start_element -eq $((ohThreshold + 1)) ]; then
                    OptimizeNorth='True'
                else
                    OptimizeSouth='True'
                fi
                gridded_optimized_OH $PerturbValueOH $PerturbValueOH $Hemis_mask_fpath $Output_fpath $OptimizeNorth $OptimizeSouth
                
                # Modify OH scale factor in HEMCO config
                sed -i -e "s| OH_pert_factor  1.0 - - - xy 1 1| OH_pert_factor ${Output_fpath} oh_scale 2000\/1\/1\/0 C xy 1 1|g" HEMCO_Config.rc
                
                if "$UseGCHP"; then
                    # add entry in ExtData.rc for GCHP
                    sed -i -e "s|^#OH_pert_factor.*|OH_pert_factor 1 N Y 2000-01-01-T00:00:00 none none oh_scale ${Output_fpath}|" ExtData.rc
                fi
            fi
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

        # create a break in EMISSIONS logic block for MeMo in background simulation
        if [[ "$x" = "background" ]]; then
            sed -i -e "/(((MeMo_SOIL_ABSORPTION/i )))EMISSIONS" HEMCO_Config.rc
            sed -i -e "/)))MeMo_SOIL_ABSORPTION/a (((EMISSIONS" HEMCO_Config.rc
        fi

        if "$KalmanMode"; then
            # Use nudged scale factors for the prior simulation and OH simulation for kalman mode
            sed -i -e "s|--> Emis_PosteriorSF       :       false|--> Emis_PosteriorSF       :       true|g" \
                -e "s|--> UseTotalPriorEmis      :       false|--> UseTotalPriorEmis      :       true|g" \
                -e "s|gridded_posterior.nc|${RunDirs}/ScaleFactors.nc|g" HEMCO_Config.rc
            if "$UseGCHP"; then
                sed -i -e "s|gridded_posterior.nc|./RunDirs/ScaleFactors.nc|g" ExtData.rc
            fi
        fi

    else
        # set 1ppb CH4 boundary conditions and restarts for all other perturbation simulations
        # Note that we use the timecycle flag C to avoid having to make additional files
        if "$isRegional"; then
            BCFile1ppb=${RunDirs}/jacobian_1ppb_ics_bcs/BCs/GEOSChem.BoundaryConditions.1ppb.${StartDate}_0000z.nc4
            BCSettings1ppb="SpeciesBC_CH4  1980-2021/1-12/1-31/* C xyz 1 CH4 - 1 1"
            sed -i -e "s|.*GEOSChem\.BoundaryConditions.*|\* BC_CH4 ${BCFile1ppb} ${BCSettings1ppb}|g" HEMCO_Config.rc
        fi
        # create symlink to 1ppb restart file
        if "$UseGCHP"; then
            RestartFile1ppb=${RunDirs}/jacobian_1ppb_ics_bcs/Restarts/GEOSChem.Restart.1ppb.${StartDate}_0000z.c${CS_RES}.nc4
            restart_fpath="Restarts/GEOSChem.Restart.1ppb.${StartDate}_0000z.c${CS_RES}.nc4"
            # Conditionally add tracers if missing
            all_tracers_exist $start_element $end_element "$restart_fpath" || {
                add_jacobian_tracers_restart_for_gchp "$start_element" "$end_element" "$RestartFile1ppb" "$restart_fpath"
            }
            cd Restarts
            ln -nsf "GEOSChem.Restart.1ppb.${StartDate}_0000z.c${CS_RES}.nc4" "GEOSChem.Restart.${StartDate}_0000z.c${CS_RES}.nc4"
            cd ..
        else
            RestartFile1ppb=${RunDirs}/jacobian_1ppb_ics_bcs/Restarts/GEOSChem.Restart.1ppb.${StartDate}_0000z.nc4
            RestartFile=$RestartFile1ppb
            ln -nsf $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
        fi
        # Also, set emissions to zero for default CH4 tracer by applying ZERO scale factor (id 5)
        sed -i -e "s|CH4 - 1 500|CH4 5 1 500|g" HEMCO_Config.rc
    fi

    # Modify restart and BC entries in HEMCO_Config.rc to look for CH4 only
    # instead of all advected species
    if [ "$UseGCHP" != "true" ]; then
        sed -i -e "s/SPC_/SPC_CH4/g" -e "s/?ALL?/CH4/g" -e "s/EFYO xyz 1 \*/EFYO xyz 1 CH4/g" HEMCO_Config.rc
        if "$isRegional"; then
            sed -i -e "s/BC_ /BC_CH4 /g" -e "s/?ADV?/CH4/g" -e "s/EFY xyz 1 \*/EFY xyz 1 CH4/g" HEMCO_Config.rc
        fi
    fi

    # Initialize previous lines to search
    GcPrevLine='- CH4'
    HcoPrevLine1='EFYO xyz 1 CH4 - 1 '
    HcoPrevLine2='1 500'
    HcoPrevLine3="#200N SCALE_ELEM_000N ${RunDirs}/StateVector.nc StateVector 2000/1/1/0 C xy 1 1 N"
    HcoPrevLine4='\* BC_CH4'
    ExtPrevLine3="#SCALE_ELEM_000N  1 N Y 2000-01-01T00:00:00 none none StateVector ./RunDirs/StateVector.nc"
    
    # Loop over state vector element numbers for this run and add each element
    # as a CH4 tracer in the configuraton files
    if is_number "$x"; then
        if [ $x -gt 0 ] && [ "$BC_elem" = false ] && [ "$OH_elem" = false ]; then
            for i in $(seq $start_element $end_element); do
                add_new_tracer
            done
        fi
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
    SFnum=$((2000 + i))

    # Add lines to geoschem_config.yml
    # Spacing in GcNewLine is intentional
    GcNewLine='\
      - CH4_'$istr
    sed -i -e "/$GcPrevLine/a $GcNewLine" geoschem_config.yml
    GcPrevLine='- CH4_'$istr

    # Add lines to species_database.yml
    SpcNextLine='CHBr3:'
    SpcNewLines='CH4_'$istr':\n  << : *CH4properties\n  Background_VV: 1.8e-6\n  FullName: Methane'
    sed -i -e "s|$SpcNextLine|$SpcNewLines\n$SpcNextLine|g" species_database.yml

    # Add lines for new tracers to HEMCO_Config.rc
    HcoNewLine2='0 CH4_Emis_Prior_'$istr' - - - - - - CH4_'$istr' '4/$SFnum' 1 500'
    sed -i -e "\|$HcoPrevLine2|a $HcoNewLine2" HEMCO_Config.rc
    HcoPrevLine2=$HcoNewLine2

    HcoNewLine3="$SFnum SCALE_ELEM_$istr ${RunDirs}/StateVector.nc StateVector 2000/1/1/0 C xy 1 1 $i"
    sed -i -e "\|$HcoPrevLine3|a $HcoNewLine3" HEMCO_Config.rc
    HcoPrevLine3=$HcoNewLine3

    if [ "$UseGCHP" != "true" ]; then
        # Add lines for restarts of new tracers to HEMCO_Config.rc
        HcoNewLine1='* SPC_CH4_'$istr' - - - - - - CH4_'$istr' - 1 1'
        sed -i -e "/$HcoPrevLine1/a $HcoNewLine1" HEMCO_Config.rc
        HcoPrevLine1='SPC_CH4_'$istr
        if "$isRegional"; then
            HcoNewLine4='* BC_CH4_'$istr' - - - - - - CH4_'$istr' - 1 1'
            sed -i -e "/$HcoPrevLine4/a $HcoNewLine4" HEMCO_Config.rc
            HcoPrevLine4='BC_CH4_'$istr
        fi
    else
        # Add lines for new tracers to ExtData.rc
        ExtNewLine3="SCALE_ELEM_$istr  1 N Y 2000-01-01T00:00:00 none none StateVector ./RunDirs/StateVector.nc"
        sed -i -e "\|$ExtPrevLine3|a $ExtNewLine3" ExtData.rc
        ExtPrevLine3=$ExtNewLine3
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
    cp ${InversionPath}/src/geoschem_run_scripts/run_jacobian_simulations.sh jacobian_runs/
    sed -i -e "s:{RunName}:${RunName}:g" \
        -e "s:{InversionPath}:${InversionPath}:g" \
        -e "s:{KalmanMode}:${KalmanMode}:g" \
        -e "s:{EndDate}:${EndDate}:g" \
        -e "s:{ReDoJacobian}:${ReDoJacobian}:g" jacobian_runs/run_jacobian_simulations.sh

    popd

    if "$KalmanMode"; then
        jacobian_period=${period_i}
    else
        jacobian_period=1
    fi

    if ! "$PrecomputedJacobian"; then

        cd ${RunDirs}/jacobian_runs
        jacobian_start=$(date +%s)

        set +e

        printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"
        # Submit job to job scheduler
        if [ "$UseGCHP" = "true" ]; then
            # Replace the line with -c $RequestedCPUs \ by three lines for multiple nodes
            sed -i -E "s|-c[[:space:]]+\\\$RequestedCPUs[[:space:]]*\\\\|-N \$NUM_NODES \\\
-n \$TOTAL_CORES \\\
-c 1 \\\
|" submit_jacobian_simulations_array.sh
        else
            # Just append -N 1 \ after the line with -c $RequestedCPUs \
            sed -i -E "/-c[[:space:]]+\\\$RequestedCPUs[[:space:]]*\\\\/a -N 1 \\" submit_jacobian_simulations_array.sh
        fi

        if "$LognormalErrors"; then
            if "$UseGCHP"; then
                sbatch --mem $RequestedMemory \
                    -c 1 \
                    -N $NUM_NODES \
                    -n $TOTAL_CORES \
                    -t $RequestedTime \
                    -p $SchedulerPartition \
                    -W run_bkgd_simulation.sh
            else
                sbatch --mem $RequestedMemory \
                    -c $RequestedCPUs \
                    -N 1 \
                    -t $RequestedTime \
                    -p $SchedulerPartition \
                    -W run_bkgd_simulation.sh
            fi
            wait
        fi

        # check if any jacobians exited with non-zero exit code
        [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

        printf "\n=== DONE JACOBIAN SIMULATIONS ===\n"
        jacobian_end=$(date +%s)
    else
        set +e
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
        JacobianRunsDir=${RunDirs}/jacobian_runs
        cd ${JacobianRunsDir}

        # Submit prior simulation to job scheduler
        printf "\n=== SUBMITTING PRIOR SIMULATION ===\n"
        if "$UseGCHP"; then
            sbatch --mem $RequestedMemory \
                -N $NUM_NODES \
                -n $TOTAL_CORES \
                -c 1 \
                -t $RequestedTime \
                -o imi_output.tmp \
                -p $SchedulerPartition \
                -W run_prior_simulation.sh
        else
            sbatch --mem $RequestedMemory \
                -c $RequestedCPUs \
                -N 1 \
                -t $RequestedTime \
                -o imi_output.tmp \
                -p $SchedulerPartition \
                -W run_prior_simulation.sh
        fi
        wait
        cat imi_output.tmp >>${InversionPath}/imi_output.log
        rm imi_output.tmp
        # check if prior simulation exited with non-zero exit code
        [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

        printf "=== DONE PRIOR SIMULATION ===\n"

        # Run the background simulation if lognormal errors enabled
        if "$LognormalErrors"; then
            printf "\n=== SUBMITTING BACKGROUND SIMULATION ===\n"
            if "$UseGCHP"; then
                sbatch --mem $RequestedMemory \
                    -c 1 \
                    -N $NUM_NODES \
                    -n $TOTAL_CORES \
                    -t $RequestedTime \
                    -p $SchedulerPartition \
                    -W run_bkgd_simulation.sh
            else
                sbatch --mem $RequestedMemory \
                    -c $RequestedCPUs \
                    -N 1 \
                    -t $RequestedTime \
                    -p $SchedulerPartition \
                    -W run_bkgd_simulation.sh
            fi
            wait
            # check if background simulation exited with non-zero exit code
            [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO
            printf "=== DONE BACKGROUND SIMULATION ===\n"
        fi

        # Get Jacobian scale factors
        python ${InversionPath}/src/inversion_scripts/get_jacobian_scalefactors.py $jacobian_period $RunDirs $ReferenceRunDir $KalmanMode
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
isRegional = sys.argv[5].lower() == 'true'
numStandaloneRuns = 0
if bcOptimized:
    numStandaloneRuns += 4
if ohOptimized:
    if isRegional:
        numStandaloneRuns += 1
    else:
        numStandaloneRuns += 2
nRuns = math.ceil((nElements - numStandaloneRuns) / nTracers)
nRuns += numStandaloneRuns
print(nRuns)
" $1 $2 $3 $4 $5
}

is_number() {
    local s="$1"
    [[ $s =~ ^[0-9]+$ ]]
}

# Description: add restarts for all Jacobian tracers for GCHP
# Usage:
#   add_jacobian_tracers_restart_for_gchp <start_element> <end_element> <org_restart_fpath> <new_restart_fpath>
add_jacobian_tracers_restart_for_gchp(){
    local start=$1
    local end=$2
    local org_restart_fpath=$3
    local new_restart_fpath=$4

    # Make a copy of the original file to work on
    cp "$org_restart_fpath" "$new_restart_fpath"

    for ((i=start; i<=end; i++)); do
        tracer=$(printf "SPC_CH4_%04d" "$i")
        
        # Duplicate variable using ncap2
        ncap2 -O -s "${tracer}=SPC_CH4;" "$new_restart_fpath" "$new_restart_fpath"
        
        # Optionally update long_name attribute
        ncatted -O -a long_name,"$tracer",o,c,"Dry mixing ratio of species CH4_$(printf "%04d" "$i")" "$new_restart_fpath"
    done
}

# Function to check if all SPC_CH4_#### variables exist
all_tracers_exist() {
    local start=$1
    local end=$2
    local file=$3

    if [ ! -f "$file" ]; then
        echo "File $file does not exist."
        return 1
    fi

    # Extract variable names from file header
    local vars
    vars=$(ncdump -h "$file" | awk '/(float|double|int|char|byte)[[:space:]]+SPC_CH4_[0-9]{4}/ {gsub(/\(.*/, "", $2); print $2}')

    for ((i=start; i<=end; i++)); do
        tracer=$(printf "SPC_CH4_%04d" "$i")
        if ! echo "$vars" | grep -q "^$tracer$"; then
            echo "Missing tracer: $tracer"
            return 1
        fi
    done

    return 0
}
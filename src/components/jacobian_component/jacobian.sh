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

    # Create directory that will contain all Jacobian run directories
    mkdir -p -v jacobian_runs

    if [ $NumJacobianRuns -gt 0 ]; then
	nRuns=$((NumJacobianRuns-1))

	# Determine approx. number of CH4 tracers per Jacobian run
	nTracers=$((nElements/NumJacobianRuns))
	printf "\nCombining Jacbian runs: Generating $NumJacobianRuns run directories with approx. $nTracers CH4 tracers (representing state vector elements) per run\n"
    else
	nRuns=$nElements
	nTracers=1
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

    # Initialize (x=0 is base run, i.e. no perturbation; x=1 is state vector element=1; etc.)
    x=0
    # Link to restart file
    RestartFileFromSpinup=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        RestartFile=$RestartFileFromSpinup
    else
        RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
        if "$UseBCsForRestart"; then
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
        fi
    fi
    PertRestartFile=$(python ${InversionPath}/src/components/jacobian_component/make_jacobian_icbc.py $RestartFile ${RunTemplate}/Restarts $StartDate)
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
    ln -s ../../GEOSChem_build/gcclassic .

    # TODO change to $PertRestartFile
    # link to restart file
    ln -s $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4

    # Modify HEMCO_Config.rc to turn off individual emission inventories
    # and use total emissions saved out from prior emissions simulation
    # instead
    printf "\nTurning on use of total prior emissions in HEMCO_Config.rc.\n"
    sed -i -e "s|UseTotalPriorEmis      :       false|UseTotalPriorEmis      :       true|g" \
           -e "s|AnalyticalInversion    :       false|AnalyticalInversion    :       true|g" \
           -e "s|GFED                   : on|GFED                   : off|g" HEMCO_Config.rc

    # Apply perturbations to total emissions with soil absorption removed
    # In this case, we still need to read soil absorption for overall CH4 flux
    #  so remove from the UseTotalPriorEmis brackets
    sed -i -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" HEMCO_Config.rc

    # TODO -- ask melissa what to do about soil absorption -- do we need this?
    # sed -i -e "/(((MeMo_SOIL_ABSORPTION/a )))UseTotalPriorEmis" \
    #        -e "/)))MeMo_SOIL_ABSORPTION/a (((UseTotalPriorEmis" HEMCO_Config.rc

    # Update settings in HISTORY.rc
    # Only save out hourly pressure fields to daily files for base run
    if [[ $x -eq 0 ]] || [[ "$x" = "background" ]]; then
        if "$HourlyCH4"; then
            sed -i -e 's/'\''Restart/#'\''Restart/g' \
                -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
                -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
                -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
        fi
    # For all other runs, just disable Restarts
    else
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
    sed -e "s:namename:${name}:g" ch4_run.template >${name}.run
    rm -f ch4_run.template
    chmod 755 ${name}.run

    ### Turn on observation operators if requested, only for base run
    if [[ $x -eq 0 ]] || [[ "$x" = "background" ]]; then
        activate_observations
    fi

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

    # BC optimization setup
    if "$OptimizeBCs"; then
        if "$OptimizeOH"; then
            bcThreshold=$(($nElements - 5))
        else
            bcThreshold=$(($nElements - 4))
        fi
        # The last four state vector elements are reserved for BC optimization of NSEW
        # domain edges. If the current state vector element is one of these, then
        # turn on BC optimization for the corresponding edge and revert emission perturbation
        if [[ $x -gt $bcThreshold ]]; then
            PerturbBCValues=$(generate_BC_perturb_values $bcThreshold $x $PerturbValueBCs)
            sed -i -e "s|CH4_boundary_condition_ppb_increase_NSEW:.*|CH4_boundary_condition_ppb_increase_NSEW: ${PerturbBCValues}|g" \
                -e "s|perturb_CH4_boundary_conditions: false|perturb_CH4_boundary_conditions: true|g" geoschem_config.yml
        fi
    fi

    if "$OptimizeOH"; then
        # The last state vector element is reserved for OH optimization.
        # If this is the current state vector element, then modify the OH
        # perturb value in HEMCO_Config.rc and revert emission perturbation.
        OHthreshold=$(($nElements - 1))
        if [ $x -gt $OHthreshold ]; then
            sed -i -e "s| OH_pert_factor  1.0| OH_pert_factor  ${PerturbValueOH}|g" HEMCO_Config.rc
        fi
    fi

    # Turn off sectoral emissions diagnostics since total emissions are
    # read in for jacobian runs
    sed -i -e "s:EmisCH4:#EmisCH4:g" HEMCO_Diagn.rc
    sed -i -e "s:#EmisCH4_Total:EmisCH4_Total:g" HEMCO_Diagn.rc
    
    # Determine start and end element numbers for this run directory
    if [ $NumJacobianRuns -lt 0 ]; then
	start_element=$x
	end_element=$x
    else
	if [ $x -eq 0 ]; then
	    start_element=0
	else
	    start_element=$(( end_element + 1 ))
	fi
	if [ $x -eq $nRuns ]; then
	    end_element=$nElements
	else
	    end_element=$(( start_element + nTracers ))
	fi
    fi

    # Modify restart and BC entries in HEMCO_Config.rc to look for CH4 only
    # instead of all advected species
    sed -i -e "s/SPC_/SPC_CH4/g"  -e "s/?ALL?/CH4/g" -e "s/EFYO xyz 1 \*/EFYO xyz 1 CH4/g" HEMCO_Config.rc
    sed -i -e "s/BC_ /BC_CH4 /g"  -e "s/?ADV?/CH4/g" -e "s/EFY xyz 1 \*/EFY xyz 1 CH4/g" HEMCO_Config.rc

    # Initialize previous lines to search
    GcPrevLine='- CH4'
    HcoPrevLine1='EFYO xyz 1 CH4 - 1 '
    HcoPrevLine2='CH4 - 1 500'
    HcoPrevLine3='Perturbations.txt - - - xy count 1'
    HcoPrevLine4='SpeciesBC_CH4'
    PertPrevLine='DEFAULT    0     0.0'

    # TODO: figure out what to do about the prior simulation
    # by default remove all emissions
    sed -i -e "s/DEFAULT    0     1.0/$PertPrevLine/g" Perturbations.txt

		
    # Loop over state vector element numbers for this run and add each element
    # as a CH4 tracer in the configuraton files
    for i in $(seq $start_element $end_element); do

	if [ $i -lt 10 ]; then
	    istr="000${i}"
	elif [ $x -lt 100 ]; then
	    istr="00${i}"
	elif [ $x -lt 1000 ]; then
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

	# Add lines to HEMCO_Config.yml
	HcoNewLine1='\
* SPC_CH4_'$istr' - - - - - - CH4_'$istr' - 1 1'
	sed -i -e "/$HcoPrevLine1/a $HcoNewLine1" HEMCO_Config.rc
	HcoPrevLine1='SPC_CH4_'$istr

	HcoNewLine2='\
0 CH4_Emis_Prior_'$istr' - - - - - - CH4_'$istr' '$SFnum' 1 500'
	sed -i "/$HcoPrevLine2/a $HcoNewLine2" HEMCO_Config.rc
	HcoPrevLine2='CH4_'$istr' '$SFnum' 1 500'

	HcoNewLine3='\
'$SFnum' SCALE_ELEM_'$istr' Perturbations_'$istr'.txt - - - xy count 1'
	sed -i "/$HcoPrevLine3/a $HcoNewLine3" HEMCO_Config.rc
	HcoPrevLine3='SCALE_ELEM_'$istr' Perturbations_'$istr'.txt - - - xy count 1'

	HcoNewLine4='\
* BC_CH4_'$istr' - - - - - - CH4_'$istr' - 1 1'
	sed -i -e "/$HcoPrevLine4/a $HcoNewLine4" HEMCO_Config.rc
	HcoPrevLine4='BC_CH4_'$istr

	# Add new Perturbations.txt and update
	cp Perturbations.txt Perturbations_${istr}.txt
	PertNewLine='\
ELEM_'$istr'  '$i'     '1.0''
	sed -i "/$PertPrevLine/a $PertNewLine" Perturbations_${istr}.txt

    done
   
    # Navigate back to top-level directory
    cd ../..
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

    if ! "$PrecomputedJacobian"; then
        jacobian_start=$(date +%s)
        if "$KalmanMode"; then
            jacobian_period=${period_i}
        else
            jacobian_period=1
        fi
        
        # update perturbation values before running jacobian simulations
        printf "\n=== UPDATING PERTURBATION SFs ===\n"
        python ${InversionPath}/src/components/jacobian_component/make_perturbation_sf.py $ConfigPath $jacobian_period 

        cd ${RunDirs}/jacobian_runs

        printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"
        # Submit job to job scheduler
        source submit_jacobian_simulations_array.sh

        if "$LognormalErrors"; then
            sbatch --mem $RequestedMemory \
                -c $RequestedCPUs \
                -t $RequestedTime \
                -p $SchedulerPartition \
                -W run_bkgd_simulation.sh
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
        sbatch --mem $RequestedMemory \
            -c $RequestedCPUs \
            -t $RequestedTime \
            -p $SchedulerPartition \
            -W run_prior_simulation.sh
        wait
        cat imi_output.tmp >>${InversionPath}/imi_output.log
        rm imi_output.tmp
        printf "=== DONE PRIOR SIMULATION ===\n"

        # Run the background simulation if lognormal errors enabled
        if "$LognormalErrors"; then
            printf "\n=== SUBMITTING BACKGROUND SIMULATION ===\n"
            sbatch --mem $RequestedMemory \
                -c $RequestedCPUs \
                -t $RequestedTime \
                -p $SchedulerPartition \
                -W run_bkgd_simulation.sh
            wait
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


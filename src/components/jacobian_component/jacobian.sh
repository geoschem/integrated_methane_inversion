#!/bin/bash

# Functions available in this file include:
#   - setup_jacobian 
#   - run_jacobian 
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

    # Copy run scripts
    cp ${InversionPath}/src/geoschem_run_scripts/run_jacobian_simulations.sh jacobian_runs/
    sed -i -e "s:{RunName}:${RunName}:g" \
           -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/run_jacobian_simulations.sh
    cp ${InversionPath}/src/geoschem_run_scripts/submit_jacobian_simulations_array.sh jacobian_runs/
    sed -i -e "s:{START}:0:g" \
           -e "s:{END}:${nElements}:g" \
           -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/submit_jacobian_simulations_array.sh
    cp ${InversionPath}/src/geoschem_run_scripts/run_prior_simulation.sh jacobian_runs/
    sed -i -e "s:{RunName}:${RunName}:g" \
           -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/run_prior_simulation.sh
    cp ${InversionPath}/src/geoschem_run_scripts/run_bkgd_simulation.sh jacobian_runs/
    sed -i -e "s:{RunName}:${RunName}:g" \
           -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/run_bkgd_simulation.sh

    # Initialize (x=0 is base run, i.e. no perturbation; x=1 is state vector element=1; etc.)
    x=0

    # Create run directory for each state vector element so we can
    # apply the perturbation to each
    while [ $x -le $nElements ]; do

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
	    x=$[$x+1]
    done

    if "$LognormalErrors"; then
        x="background"
        xstr=$x
        create_simulation_dir
    fi

    printf "\n=== DONE CREATING JACOBIAN RUN DIRECTORIES ===\n"
}

create_simulation_dir() {
    # Define the run directory name
	name="${RunName}_${xstr}"

	# Make the directory
	runDir="./jacobian_runs/${name}"
	mkdir -p -v ${runDir}

	# Copy run directory files
	cp -r ${RunTemplate}/*  ${runDir}
	cd $runDir

	# Link to GEOS-Chem executable instead of having a copy in each rundir
	rm -rf gcclassic
	ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFileFromSpinup=${RunDirs}/spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        ln -s $RestartFileFromSpinup Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
	else
	    RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
	    ln -s $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
	    if "$UseBCsForRestart"; then
		    sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
        fi
	fi
   
	# Update settings in geoschem_config.yml except for the base run
    if [ $x -ne 0 ] && [ "$x" != "background" ]; then
	    sed -i -e "s|emission_perturbation_factor: 1.0|emission_perturbation_factor: ${PerturbValue}|g" \
	           -e "s|state_vector_element_number: 0|state_vector_element_number: ${xUSE}|g" geoschem_config.yml
    fi

    # BC optimization setup
    if "$OptimizeBCs"; then
        bcThreshold=$(($nElements - 4))
        # The last four state vector elements are reserved for BC optimization of NSEW
        # domain edges. If the current state vector element is one of these, then
        # turn on BC optimization for the corresponding edge and revert emission perturbation
        if [ $x -gt $bcThreshold ]; then
            PerturbBCValues=$(generate_BC_perturb_values $bcThreshold $x $PerturbValueBCs)
            sed -i -e "s|CH4_boundary_condition_ppb_increase_NSEW:.*|CH4_boundary_condition_ppb_increase_NSEW: ${PerturbBCValues}|g" \
                -e "s|perturb_CH4_boundary_conditions: false|perturb_CH4_boundary_conditions: true|g" \
                -e "s|emission_perturbation_factor: ${PerturbValue}|emission_perturbation_factor: 1.0|g" \
                -e "s|state_vector_element_number: ${xUSE}|state_vector_element_number: 0|g" geoschem_config.yml
        fi
    fi 

	# Update settings in HISTORY.rc
	# Only save out hourly pressure fields to daily files for base run
	if [ $x -eq 0 ] || [ "$x" = "background" ]; then
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
        # TODO: OH_pert_factor sed will be unnecessary once geoschem Emissions logical is fixed
        # until  then we remove the OH_pert_factor line and add it back in before at line 877
        # sed -i '/2 OH_pert_factor  1.0 - - - xy 1 1/d' HEMCO_Config.rc
        sed -i -e 's/EMISSIONS              :       true/EMISSIONS              :       false/g' \
               -e 's/GFED                   : on    CH4/GFED                   : off    CH4/g'  HEMCO_Config.rc
        
        # TODO: remove this hack once geoschem Emissions logical is fixed
        # echo '2 OH_pert_factor  1.0 - - - xy 1 1' > add.txt
        # sed -i -e '877r add.txt' HEMCO_Config.rc
        # rm add.txt
    fi

	# Create run script from template
	sed -e "s:namename:${name}:g" ch4_run.template > ${name}.run
	rm -f ch4_run.template
	chmod 755 ${name}.run

    ### Turn on observation operators if requested, only for base run
    if [ $x -eq 0 ] || [ "$x" = "background" ]; then
    	activate_observations
    fi

    ### Perform dry run if requested, only for base run
    if [ $x -eq 0 ]; then
        if "$ProductionDryRun"; then
            printf "\nExecuting dry-run for production runs...\n"
            ./gcclassic --dryrun &> log.dryrun
            ./download_data.py log.dryrun aws
        fi
    fi

    # Navigate back to top-level directory
	cd ../..
}

# Description: Run jacobian simulations
# Usage:
#   run_jacobian
run_jacobian() {
    if ! "$PrecomputedJacobian"; then
        jacobian_start=$(date +%s)
        printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"

        cd ${RunDirs}/jacobian_runs

        if ! "$isAWS"; then
            # Load environment with modules for compiling GEOS-Chem Classic
            source ${GEOSChemEnv} 
        fi

        # Submit job to job scheduler
        source submit_jacobian_simulations_array.sh

        if "$LognormalErrors"; then
            sbatch --mem $SimulationMemory \
                -c $SimulationCPUs \
                -t $RequestedTime \
                -p $SchedulerPartition \
                -W run_bkgd_simulation.sh; wait;
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

        precomputedJacobianCacheNormal=${precomputedJacobianCachePrefix}/data_converted
        precomputedJacobianCacheLognormal=${precomputedJacobianCachePrefix}/data_converted_background
        
        if [ -d $precomputedJacobianCacheNormal ]; then
            ln -s $precomputedJacobianCacheNormal data_converted_reference
        fi

        if [ -d $precomputedJacobianCacheLognormal ]; then
            ln -s $precomputedJacobianCacheLognormal data_converted_reference_background
        fi


        # Run the prior simulation
        cd ${JacobianRunsDir}
            
        if ! "$isAWS"; then
            # Load environment with modules for compiling GEOS-Chem Classic
            source ${GEOSChemEnv}
        fi

        # Submit prior simulation to job scheduler
        printf "\n=== SUBMITTING PRIOR SIMULATION ===\n"
        sbatch --mem $SimulationMemory \
                -c $SimulationCPUs \
                -t $RequestedTime \
                -p $SchedulerPartition \
                -W run_prior_simulation.sh; wait;
        printf "=== DONE PRIOR SIMULATION ===\n"

        # Get Jacobian scale factors
        python ${InversionPath}/src/inversion_scripts/get_jacobian_scalefactors.py $period_i $RunDirs $ReferenceRunDir; wait
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


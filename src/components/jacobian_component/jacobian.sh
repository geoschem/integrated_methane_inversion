#!/bin/bash

# Functions available in this file include:
#   - setup_jacobian 
#   - run_jacobian 

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
   
	# Update settings in geoschem_config.yml
	sed -i -e "s|emission_perturbation: 1.0|emission_perturbation: ${PerturbValue}|g" \
	       -e "s|state_vector_element_number: 0|state_vector_element_number: ${xUSE}|g" geoschem_config.yml

	# Update settings in HISTORY.rc
	# Only save out hourly pressure fields to daily files for base run
	if [ $x -eq 0 ]; then
	    if "$HourlyCH4"; then
                sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                       -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
                       -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
                       -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
	    fi
	fi

	# Create run script from template
	sed -e "s:namename:${name}:g" ch4_run.template > ${name}.run
	rm -f ch4_run.template
	chmod 755 ${name}.run

    ### Turn on observation operators if requested, only for base run
    if [ $x -eq 0 ]; then
    	if "$GOSAT"; then
		OLD="GOSAT: false"
		NEW="GOSAT: true"
		sed -i "s/$OLD/$NEW/g" geoschem_config.yml
    	fi
    	if "$TCCON"; then
		OLD="TCCON: false"
		NEW="TCCON: true"
		sed -i "s/$OLD/$NEW/g" geoschem_config.yml
    	fi
    	if "$AIRS"; then
		OLD="AIR: false"
		NEW="AIR: true"
		sed -i "s/$OLD/$NEW/g" geoschem_config.yml
	fi
	if "$PLANEFLIGHT"; then
		mkdir -p Plane_Logs
		sed -i "/planeflight/{N;s/activate: false/activate: true/}" geoschem_config.yml
	
		OLD="flight_track_file: Planeflight.dat.YYYYMMDD"
		NEW="flight_track_file: Planeflights\/Planeflight.dat.YYYYMMDD"
		sed -i "s/$OLD/$NEW/g" geoschem_config.yml
		OLD="output_file: plane.log.YYYYMMDD"
		NEW="output_file: Plane_Logs\/plane.log.YYYYMMDD"
		sed -i "s/$OLD/$NEW/g" geoschem_config.yml
    	fi
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

	# Increment
	x=$[$x+1]

    done

    printf "\n=== DONE CREATING JACOBIAN RUN DIRECTORIES ===\n"
}

# Description: Run jacobian simulations
# Usage:
#   run_jacobian
run_jacobian() {
    jacobian_start=$(date +%s)
    printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"

    cd ${RunDirs}/jacobian_runs

    if ! "$isAWS"; then
        # Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv} 
    fi

    # Submit job to job scheduler
    source submit_jacobian_simulations_array.sh

    # check if any jacobians exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO

    printf "\n=== DONE JACOBIAN SIMULATIONS ===\n"
    jacobian_end=$(date +%s)
}

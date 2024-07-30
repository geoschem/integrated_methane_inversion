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
    mkdir -p jacobian_1ppb_ics_bcs/BCs
    OrigBCFile=${fullBCpath}/GEOSChem.BoundaryConditions.${StartDate}_0000z.nc4
    python ${InversionPath}/src/components/jacobian_component/make_jacobian_icbc.py $OrigBCFile ${RunDirs}/jacobian_1ppb_ics_bcs/BCs $StartDate

    # Create directory that will contain all Jacobian run directories
    mkdir -p -v jacobian_runs

    if [ $NumJacobianTracers -gt 1 ]; then
        nRuns=$(calculate_num_jacobian_runs $NumJacobianTracers $nElements $OptimizeBCs $OptimizeOH)

        # Determine approx. number of CH4 tracers per Jacobian run
        printf "\nCombining Jacobian runs: Generating $nRuns run directories with approx. $NumJacobianTracers CH4 tracers (representing state vector elements) per run\n"
    else
        nRuns=$nElements
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
    printf "\nTurning on use of total prior emissions in HEMCO_Config.rc.\n"
    sed -i -e "s|UseTotalPriorEmis      :       false|UseTotalPriorEmis      :       true|g" \
        -e "s|AnalyticalInversion    :       false|AnalyticalInversion    :       true|g" \
        -e "s|EmisCH4_Total|EmisCH4_Total_ExclSoilAbs|g" \
        -e "s|GFED                   : on|GFED                   : off|g" HEMCO_Config.rc

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
        # set 1ppb CH4 boundary conditions and restarts for all other perturbation simulations
        # Note that we use the timecycle flag C to avoid having to make additional files
        RestartFile=${RunDirs}/jacobian_1ppb_ics_bcs/Restarts/GEOSChem.Restart.1ppb.${StartDate}_0000z.nc4
        BCFile1ppb=${RunDirs}/jacobian_1ppb_ics_bcs/BCs/GEOSChem.BoundaryConditions.1ppb.${StartDate}_0000z.nc4
        BCSettings1ppb="SpeciesBC_CH4  1980-2021/1-12/1-31/* C xyz 1 CH4 - 1 1"
        sed -i -e "s|.*GEOSChem\.BoundaryConditions.*|\* BC_CH4 ${BCFile1ppb} ${BCSettings1ppb}|g" HEMCO_Config.rc
        # create symlink to 1ppb restart file
        ln -sf $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
        # Also, set emissions to zero for default CH4 tracer by applying new ZERO scale factor
        sed -i -e "/1 NEGATIVE       -1.0 - - - xy 1 1/a 5 ZERO            0.0 - - - xy 1 1" \
            -e "s|CH4 - 1 500|CH4 5 1 500|g" HEMCO_Config.rc
    fi

    # Modify restart and BC entries in HEMCO_Config.rc to look for CH4 only
    # instead of all advected species
    sed -i -e "s/SPC_/SPC_CH4/g" -e "s/?ALL?/CH4/g" -e "s/EFYO xyz 1 \*/EFYO xyz 1 CH4/g" HEMCO_Config.rc
    sed -i -e "s/BC_ /BC_CH4 /g" -e "s/?ADV?/CH4/g" -e "s/EFY xyz 1 \*/EFY xyz 1 CH4/g" HEMCO_Config.rc

    # Initialize previous lines to search
    GcPrevLine='- CH4'
    HcoPrevLine1='EFYO xyz 1 CH4 - 1 '
    HcoPrevLine2='CH4 5 1 500'
    HcoPrevLine3='Perturbations.txt - - - xy count 1'
    HcoPrevLine4='\* BC_CH4'
    PertPrevLine='DEFAULT    0     0.0'

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

    # by default remove all emissions except for in the prior simulation
    # and the OH perturbation simulation
    if [ $x -gt 0 ]; then
        sed -i -e "s/DEFAULT    0     1.0/$PertPrevLine/g" Perturbations.txt
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

    # Add new Perturbations.txt and update for non prior runs
    cp Perturbations.txt Perturbations_${istr}.txt
    if [ $x -gt 0 ]; then
        PertNewLine='\
ELEM_'$istr'  '$i'     '0.0''
        sed -i "/$PertPrevLine/a $PertNewLine" Perturbations_${istr}.txt
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

        cd ${RunDirs}/jacobian_runs

        # create 1ppb restart file
        OrigRestartFile=$(readlink ${RunName}_0000/Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4)
        python ${InversionPath}/src/components/jacobian_component/make_jacobian_icbc.py $OrigRestartFile ${RunDirs}/jacobian_1ppb_ics_bcs/Restarts $StartDate
        cd ${RunDirs}/jacobian_1ppb_ics_bcs/Restarts/
        if [ -f GEOSChem.BoundaryConditions.1ppb.${StartDate}_0000z.nc4 ]; then
            mv GEOSChem.BoundaryConditions.1ppb.${StartDate}_0000z.nc4 GEOSChem.Restart.1ppb.${StartDate}_0000z.nc4
        fi
        cd ${RunDirs}/jacobian_runs
        set +e

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

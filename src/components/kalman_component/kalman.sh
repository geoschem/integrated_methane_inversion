#!/bin/bash

# Functions available in this file include:
#   - setup_kf
#   - run_kf

# Description: Setup Kalman filter prereqiuisites
# Usage:
#   setup_kf
setup_kf() {
    StateVectorFile="${RunDirs}/StateVector.nc"

    # Create a parent directory for the Kalman filter inversions
    # Include a link to the state vector file for use with run_inversion.sh
    mkdir -p ${RunDirs}/kf_inversions
    ln -sf $StateVectorFile $RunDirs/kf_inversions/StateVector.nc

    # Define Kalman filter update periods
    python ${InversionPath}/src/kf_scripts/make_periods_csv.py $StartDate $EndDate $UpdateFreqDays $RunDirs; wait

    # Create unit scale factor file
    python ${InversionPath}/src/kf_scripts/make_unit_sf.py $StateVectorFile $RunDirs; wait

    # Create directory to archive prior scale factors for each inversion period
    mkdir -p ${RunDirs}/archive_sf

    # Number of state vector elements
    nElements=$(ncmax StateVector ${StateVectorFile})
}

run_kf() {

    if ("$DoJacobian" && "$DoInversion" && "$DoPosterior"); then
        # Key directories
        JacobianRunsDir="${RunDirs}/jacobian_runs"
        PriorRunDir="${JacobianRunsDir}/${RunName}_0000"
        PosteriorRunDir="${RunDirs}/posterior_run"
        InversionDir="${RunDirs}/inversion_template"

        PeriodsFile="${RunDirs}/periods.csv"
        nPeriods=$(($(wc -l < ${PeriodsFile}) - 1))


        for ((i=FirstPeriod;i<=nPeriods;i++)); do

            ##=======================================================================
            ##  Setup (dates, emission scale factors)
            ##=======================================================================

            # Print current period
            echo -e "\nPeriod ${i}"

            # Create inversion directory for the period
            cp -r ${RunDirs}/inversion_template/. ${RunDirs}/kf_inversions/period${i}

            # Get Start/End dates of current period from periods.csv
            ithLine=$(sed "$((i+1))q;d" $PeriodsFile)
            ithDates=(${ithLine//,/ })
            StartDate_i=${ithDates[0]}
            EndDate_i=${ithDates[1]}
            echo "Start, End: $StartDate_i, $EndDate_i"

            # Set dates in geoschem_config.yml for prior, perturbation, and posterior runs
            python ${InversionPath}/src/kf_scripts/change_dates.py $StartDate_i $EndDate_i $JacobianRunsDir; wait
            python ${InversionPath}/src/kf_scripts/change_dates.py $StartDate_i $EndDate_i $PosteriorRunDir; wait
            echo "Edited Start/End dates in geoschem_config.yml for prior/perturbed/posterior simulations: $StartDate_i to $EndDate_i"

            # Prepare initial (prior) emission scale factors for the current period
            ConfigPath=${InversionPath}/${ConfigFile}
            python ${InversionPath}/src/kf_scripts/prepare_sf.py $ConfigPath $i ${RunDirs} $NudgeFactor; wait

            ##=======================================================================
            ##  Submit all Jacobian simulations OR submit only the Prior simulation
            ##=======================================================================
            
            # run jacobian simulation for the given period
            run_jacobian

            # run inversion for the given period
            run_inversion

            # Update ScaleFactor.nc with the new posterior scale factors before running the posterior simulation
            # NOTE: This also creates the posterior_sf_period{i}.nc file in archive_sf/
            python ${InversionPath}/src/kf_scripts/multiply_posteriors.py $i ${RunDirs}; wait
            echo "Multiplied posterior scale factors over record"

            # Print total posterior emissions
            python ${InversionPath}/src/kf_scripts/print_posterior_emissions.py $ConfigPath $i ${RunDirs}; wait

            run_posterior

            # Make a copy of the posterior output/diags files for postproc_diags.py
            copydir="${PosteriorRunDir}/OutputDir"
            cp ${copydir}/GEOSChem.SpeciesConc.${EndDate_i}_0000z.nc4 ${copydir}/GEOSChem.SpeciesConc.Copy.${EndDate_i}_0000z.nc4
            cp ${copydir}/GEOSChem.LevelEdgeDiags.${EndDate_i}_0000z.nc4 ${copydir}/GEOSChem.LevelEdgeDiags.Copy.${EndDate_i}_0000z.nc4
            echo "Made a copy of the final posterior SpeciesConc and LevelEdgeDiags files"

            # Make link to restart file from posterior run directory in each Jacobian run directory
            for ((x=0;x<=nElements;x++)); do
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
                ln -sf ${PosteriorRunDir}/Restarts/GEOSChem.Restart.${EndDate_i}_0000z.nc4 ${JacobianRunsDir}/${RunName}_${xstr}/Restarts/.
            done
            echo "Copied posterior restart to $((x-1)) Jacobian run directories for next iteration"

            cd ${InversionPath}

            # Delete unneeded daily restart files from Jacobian and posterior directories
            # python ${InversionPath}/src/kf_scripts/cull_restarts.py $JacobianRunsDir $PosteriorRunDir $StartDate_i $EndDate_i

            # Move to next time step
            echo -e "Moving to next iteration\n"


        done
    else
        printf "Exiting: DoJacobian, DoInversion, and DoPosterior must all be set to true to continue with KalmanMode\n"
    fi

}

# run_period() {


# }

#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o "imi_output.log"

# This script will run the Integrated Methane Inversion (IMI) with GEOS-Chem.
# For documentation, see https://imi.readthedocs.io.
#
# Authors: Daniel Varon, Melissa Sulprizio, Lucas Estrada, Will Downs

##=======================================================================
## Import Shell functions
##=======================================================================
source src/utilities/common.sh
source src/components/setup_component/setup.sh
source src/components/template_component/template.sh
source src/components/statevector_component/statevector.sh
source src/components/preview_component/preview.sh
source src/components/spinup_component/spinup.sh
source src/components/jacobian_component/jacobian.sh
source src/components/inversion_component/inversion.sh
source src/components/posterior_component/posterior.sh

# trap and exit on errors
trap 'imi_failed $LINENO' ERR

start_time=$(date)
setup_start=$(date +%s)

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE (run_imi.sh) ===\n"

# Check if user has specified a configuration file
if [[ $# == 1 ]] ; then
    ConfigFile=$1
else
    ConfigFile="config.yml"
fi

# Get configuration
source src/utilities/parse_yaml.sh
eval $(parse_yaml ${ConfigFile})

if ! "$isAWS"; then
    # Activate Conda environment
    printf "\nActivating conda environment: ${CondaEnv}\n"
    eval "$(conda shell.bash hook)"
    conda activate $CondaEnv
fi

# Check all necessary config variables are present
python src/utilities/sanitize_input_yaml.py $ConfigFile || imi_failed

# Set path to IMI runs
RunDirs="${OutputPath}/${RunName}"

##=======================================================================
## Standard settings
##=======================================================================

# In safe mode check whether selected options will overwrite existing files
if "$SafeMode"; then

    # Check if directories exist before creating them
    if ([ -d "${RunDirs}/spinup_run" ] && "$SetupSpinupRun") || \
       ([ -d "${RunDirs}/jacobian_runs" ] && "$SetupJacobianRuns") || \
       ([ -d "${RunDirs}/inversion" ] && "$SetupInversion") || \
       ([ -d "${RunDirs}/posterior_run" ] && "$SetupPosteriorRun"); then
        
        printf "\nERROR: Run directories in ${RunDirs}/"
        printf "\n   already exist. Please change RunName or change the"
        printf "\n   Setup* options to false in the IMI config file.\n"
        printf "\nIMI $RunName Aborted\n"
        exit 1 
    fi

    # Check if output from previous runs exists
    if ([ -d "${RunDirs}/spinup_run" ] && "$DoSpinup") || \
       ([ -d "${RunDirs}/jacobian_runs" ] && "$DoJacobian") || \
       ([ -d "${RunDirs}/inversion" ] && "$DoInversion") || \
       ([ -d "${RunDirs}/posterior_run/OutputDir/" ] && "$DoPosterior"); then
        printf "\nWARNING: Output files in ${RunDirs}/" 
        printf "\n  may be overwritten. Please change RunName in the IMI"
        printf "\n  config file to avoid overwriting files.\n"
        printf "\n  To proceed, and overwrite existing output files, set"
        printf "\n  SafeMode in the config file to false.\n" 
        printf "\nIMI $RunName Aborted\n"
        exit 1 
    fi
fi

# Path to inversion setup
InversionPath=$(pwd -P)

##=======================================================================
##  Download the TROPOMI data
##=======================================================================

# Download TROPOMI data from AWS. You will be charged if your ec2 instance is not in the eu-central-1 region.
mkdir -p -v ${RunDirs}
tropomiCache=${RunDirs}/data_TROPOMI
if "$isAWS"; then
    { # test if instance has access to TROPOMI bucket
        stdout=`aws s3 ls s3://meeo-s5p`
    } || { # catch 
        printf "\nError: Unable to connect to TROPOMI bucket. This is likely caused by misconfiguration of the ec2 instance iam role s3 permissions.\n"
        printf "IMI $RunName Aborted.\n"
        exit 1
    }
    mkdir -p -v $tropomiCache
    printf "Downloading TROPOMI data from S3\n"
    python src/utilities/download_TROPOMI.py $StartDate $EndDate $tropomiCache
    printf "\nFinished TROPOMI download\n"
else
    # use existing tropomi data and create a symlink to it
    if [[ ! -L $tropomiCache ]]; then
	ln -s $DataPathTROPOMI $tropomiCache
    fi
fi

# Check to make sure there are no duplicate TROPOMI files (e.g., two files with the same orbit number but a different processor version)
python src/utilities/test_TROPOMI_dir.py $tropomiCache

##=======================================================================
##  Run the setup script
##=======================================================================
if "$RunSetup"; then
    setup_imi
fi
setup_end=$(date +%s)

##=======================================================================
##  Submit spinup simulation
##=======================================================================
if  "$DoSpinup"; then
    run_spinup
fi

##=======================================================================
##  Initiate Kalman filter
##=======================================================================

if ("$DoJacobian" && "$DoInversion" && "$DoPosterior"); then

    # Activate Conda environment
    printf "\nActivating conda environment: ${CondaEnv}\n"
    eval "$(conda shell.bash hook)"
    if "$isAWS"; then
        source $CondaFile
    fi
    conda activate $CondaEnv

    # Key files and directories
    JacobianRunsDir="${RunDirs}/jacobian_runs"
    PriorRunDir="${JacobianRunsDir}/${RunName}_0000"
    PosteriorRunDir="${RunDirs}/posterior_run"
    StateVectorFile="${RunDirs}/StateVector.nc"
    InversionDir="${RunDirs}/inversion_template"

    # Create a parent directory for the Kalman filter inversions
    # Include a link to the state vector file for use with run_inversion.sh
    mkdir -p ${RunDirs}/kf_inversions
    ln -sf $StateVectorFile $RunDirs/kf_inversions/StateVector.nc

    # Define Kalman filter update periods
    python ${InversionPath}/src/components/kalman_component/make_periods_csv.py $StartDate $EndDate $UpdateFreqDays $RunDirs; wait
    PeriodsFile="${RunDirs}/periods.csv"
    nPeriods=$(($(wc -l < ${PeriodsFile}) - 1))

    # Create unit scale factor file
    python ${InversionPath}/src/components/kalman_component/make_unit_sf.py $StateVectorFile $RunDirs; wait

    # Create directory to archive prior scale factors for each inversion period
    mkdir -p ${RunDirs}/archive_sf

    # Number of state vector elements
    if ! "$isAWS"; then
        # Load environment with NCO
        source ${NCOEnv}
    fi
    function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
    nElements=$(ncmax StateVector ${StateVectorFile})

    # Kalman filter loop
    # NOTE: FirstPeriod is defined in config.yml
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
        python ${InversionPath}/src/components/kalman_component/change_dates.py $StartDate_i $EndDate_i $JacobianRunsDir; wait
        python ${InversionPath}/src/components/kalman_component/change_dates.py $StartDate_i $EndDate_i $PosteriorRunDir; wait
        echo "Edited Start/End dates in geoschem_config.yml for prior/perturbed/posterior simulations: $StartDate_i to $EndDate_i"

        # Prepare initial (prior) emission scale factors for the current period
        ConfigPath=${InversionPath}/${ConfigFile}
        python ${InversionPath}/src/components/kalman_component/prepare_sf.py $ConfigPath $i ${RunDirs} $NudgeFactor; wait

        ##=======================================================================
        ##  Submit all Jacobian simulations OR submit only the Prior simulation
        ##=======================================================================

        jacobian_start=$(date +%s)
        if "$DoJacobian"; then

            if ! "$PrecomputedJacobian"; then

                printf "\n=== SUBMITTING JACOBIAN SIMULATIONS ===\n"

                cd ${JacobianRunsDir}

                if ! "$isAWS"; then
                    # Load environment with modules for compiling GEOS-Chem Classic
                    source ${GEOSChemEnv} 
                fi

                # Submit job to job scheduler
                source submit_jacobian_simulations_array.sh; wait;

                # check if any jacobians exited with non-zero exit code
                [ ! -f ".error_status_file.txt" ] || imi_failed

                printf "=== DONE JACOBIAN SIMULATIONS ===\n"

            else

                # Replace the (empty) data_sensitivities folder with a symlink to the
                # sensitivities from the reference inversion w/ precomputed Jacobian.
                cd ${RunDirs}/kf_inversions/period${i}

                if "$KalmanMode"; then
                    precomputedSensiCache=${ReferenceRunDir}/kf_inversions/period${i}/data_sensitivities
                else
                    # Need this alternative when run_imi.sh and run_kf.sh get combined... TODO
                    precomputedSensiCache=${ReferenceRunDir}/inversion/data_sensitivities
                fi
                # mv rather than rm, to prevent accidental deletion of original data_sensitivities/ ?
                mv data_sensitivities temp_dir
                ln -s $precomputedSensiCache data_sensitivities

                # Run the prior simulation
                cd ${JacobianRunsDir}
            
                if ! "$isAWS"; then
                    # Load environment with modules for compiling GEOS-Chem Classic
                    source ${GEOSChemEnv}
                fi

                # Submit prior simulation to job scheduler
                printf "\n=== SUBMITTING PRIOR SIMULATION ===\n"
                sbatch -W run_prior_simulation.sh; wait;
                printf "=== DONE PRIOR SIMULATION ===\n"

                # Get Jacobian scale factors
                python ${InversionPath}/src/inversion_scripts/get_jacobian_scalefactors.py $i $RunDirs $ReferenceRunDir; wait
                printf "Got Jacobian scale factors\n"

            fi

        fi
        jacobian_end=$(date +%s)

        ##=======================================================================
        ##  Process data and run inversion
        ##=======================================================================

        inversion_start=$(date +%s)
        if "$DoInversion"; then

            printf "\n=== RUNNING INVERSION ===\n"

            cd ${RunDirs}/kf_inversions/period${i}

            if ! "$isAWS"; then
                # Activate Conda environment
                printf "\nActivating conda environment: ${CondaEnv}\n"
                eval "$(conda shell.bash hook)"
                conda activate $CondaEnv
            fi

            # Modify inversion driver script to reflect current inversion period
            sed -i "s|data_TROPOMI\"|data_TROPOMI\"\n\n# Defined via run_kf.sh:\nStartDate=${StartDate_i}\nEndDate=${EndDate_i}|g" run_inversion.sh
            if (( i > 1 )); then
                sed -i "s,FirstSimSwitch=true,FirstSimSwitch=false,g" run_inversion.sh
            fi
            
            # Execute inversion driver script
            sbatch -W run_inversion.sh; wait;
                
            printf "=== DONE RUNNING INVERSION ===\n\n"

        fi
        inversion_end=$(date +%s)

        # Update ScaleFactor.nc with the new posterior scale factors before running the posterior simulation
        # NOTE: This also creates the posterior_sf_period{i}.nc file in archive_sf/
        python ${InversionPath}/src/components/kalman_component/multiply_posteriors.py $i ${RunDirs}; wait
        echo "Multiplied posterior scale factors over record"

        # Print total posterior emissions
        python ${InversionPath}/src/components/kalman_component/print_posterior_emissions.py $ConfigPath $i ${RunDirs}; wait

        ##=======================================================================
        ##  Submit posterior simulation and process the output
        ##=======================================================================

        posterior_start=$(date +%s)
        if "$DoPosterior"; then

            cd ${PosteriorRunDir}
            
            if ! "$isAWS"; then
                # Load environment with modules for compiling GEOS-Chem Classic
                source ${GEOSChemEnv}
            fi

            # Submit job to job scheduler
            printf "\n=== SUBMITTING POSTERIOR SIMULATION ===\n"
            sbatch -W ${RunName}_Posterior.run; wait;
            printf "=== DONE POSTERIOR SIMULATION ===\n"

            cd ${RunDirs}/kf_inversions/period${i}

            # Fill missing data (first hour of simulation) in posterior output
            if (( i == 1 )); then
                PrevDir="${RunDirs}/spinup_run"
            else
                PrevDir="${RunDirs}/posterior_run"
            fi
            printf "\n=== Calling postproc_diags.py for posterior ===\n"
            python ${InversionPath}/src/inversion_scripts/postproc_diags.py $RunName $PosteriorRunDir $PrevDir $StartDate_i; wait
            printf "=== DONE -- postproc_diags.py ===\n"

            # Build directory for hourly posterior GEOS-Chem output data
            mkdir -p data_converted_posterior
            mkdir -p data_geoschem_posterior
            mkdir -p data_visualization_posterior
            GCsourcepth="${PosteriorRunDir}/OutputDir"
            GCDir="./data_geoschem_posterior"
            printf "\n=== Calling setup_gc_cache.py for posterior ===\n"
            python ${InversionPath}/src/inversion_scripts/setup_gc_cache.py $StartDate_i $EndDate_i $GCsourcepth $GCDir; wait
            printf "=== DONE -- setup_gc_cache.py ===\n"

            if ! "$isAWS"; then
                # Load environment with NCO
                source ${NCOEnv}
            fi

            # Sample GEOS-Chem atmosphere with TROPOMI
            function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
            function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
            LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
            LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
            LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
            LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)
            nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc)
            rm ~/foo.nc
            FetchTROPOMI="False"
            isPost="True"

            printf "\n=== Calling jacobian.py to sample posterior simulation (without jacobian sensitivity analysis) ===\n"
            python ${InversionPath}/src/inversion_scripts/jacobian.py $StartDate_i $EndDate_i $LonMinInvDomain $LonMaxInvDomain $LatMinInvDomain $LatMaxInvDomain $nElements $tropomiCache $isPost; wait
            printf "=== DONE sampling the posterior simulation ===\n\n"

        fi
        posterior_end=$(date +%s)

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

fi

exit 0

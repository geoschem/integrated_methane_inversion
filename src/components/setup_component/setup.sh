#!/bin/bash

# Functions available in this file include:
#   - setup_imi
#   - activate_observations

# Description:
#   This script will set up an Integrated Methane Inversion (IMI) with GEOS-Chem.
#   For documentation, see https://imi.readthedocs.io.
# Usage:
#   setup_imi
setup_imi() {
    printf "\n=== RUNNING SETUP ===\n"

    cd ${InversionPath}

    if ! "$isAWS"; then
		if [ ! -f "${GEOSChemEnv}" ]; then
			printf "\nGEOS-Chem environment file does not exist!"
			printf "\nIMI $RunName Aborted\n"
			exit 1
		else
	        # Load environment with modules for compiling GEOS-Chem Classic
    	    source ${GEOSChemEnv}
    	fi
    fi

    ##=======================================================================
    ## Standard settings
    ##=======================================================================
    # Start and end date for the spinup simulation
    SpinupStart=$(date --date="${StartDate} -${SpinupMonths} month" +%Y%m%d)
    SpinupEnd=${StartDate}

    # Use global boundary condition files for initial conditions
    UseBCsForRestart=true

    printf "\nActivating conda environment: ${CondaEnv}\n"
    eval "$(conda shell.bash hook)"
    if "$isAWS"; then
        # Get max process count for spinup, production, and run_inversion scripts
        output=$(echo $(slurmd -C))
        array=($output)
        cpu_str=$(echo ${array[1]})
        cpu_count=$(echo ${cpu_str:5})

        # With sbatch reduce cpu_count by 1 to account for parent sbatch process 
        # using 1 core 
        if "$UseSlurm"; then 
            cpu_count="$((cpu_count-1))"
        fi

        # Source Conda environment file
        source $CondaFile

    fi

    # Activate Conda environment
    conda activate $CondaEnv

    ##=======================================================================
    ## Download Boundary Conditions files if requested
    ##=======================================================================
    fullBCpath="${BCpath}/${BCversion}"
    if "$BCdryrun"; then

        mkdir -p ${fullBCpath}

        if "$DoSpinup"; then
            START=${SpinupStart}
        else
            START=${StartDate}
        fi
        printf "\nDownloading boundary condition data for $START to $EndDate\n"
        python src/utilities/download_bc.py ${START} ${EndDate} ${fullBCpath} ${BCversion}

    fi

    ##=======================================================================
    ## Download initial restart file if requested
    ##=======================================================================
    if "$RestartDownload"; then
        RestartFile=${RestartFilePrefix}${SpinupStart}_0000z.nc4
        if [ ! -f "$RestartFile" ]; then
            aws s3 cp --request-payer=requester s3://imi-boundary-conditions/GEOSChem.BoundaryConditions.${SpinupStart}_0000z.nc4 $RestartFile
        fi
        RestartFilePreview=${RestartFilePreviewPrefix}${StartDate}_0000z.nc4
        if [ ! -f "$RestartFilePreview" ]; then
            aws s3 cp --request-payer=requester s3://imi-boundary-conditions/GEOSChem.BoundaryConditions.${StartDate}_0000z.nc4 $RestartFilePreview
        fi
    fi

    ##=======================================================================
    ## Define met and grid fields for HEMCO_Config.rc
    ##=======================================================================
    if [ "$Met" == "geosfp" ]; then
        metUC="GEOSFP"
        metDir="GEOS_FP"
        native="0.25x0.3125"
        constYr="2011"
        LandCoverFileExtension="nc"
    elif [ "$Met" == "merra2" ]; then
        metUC="MERRA2"
        metDir="MERRA2"
        native="0.5x0.625"
        constYr="2015"
        LandCoverFileExtension="nc4"
    fi

    if [ "$Res" = "4x5" ]; then
        gridRes="${Res}"
        gridResLong="4.0x5.0"
    elif [ "$Res" == "2x2.5" ]; then
        gridRes="2x25"
        gridResLong="2.0x2.5"
    elif [ "$Res" == "0.5x0.625" ]; then
        gridRes="05x0625"
        gridResLong="${Res}"
    elif [ "$Res" == "0.25x0.3125" ]; then
        gridRes="025x03125"
        gridResLong="${Res}"
    fi

    if [ -z "$NestedRegion" ]; then
        gridDir="$Res"
    else
        gridDir="${Res}_${NestedRegion}"
    fi

    # Define path to GEOS-Chem run directory files
    GCClassicPath="${InversionPath}/GCClassic"
    RunFilesPath="${GCClassicPath}/run"

    # Create working directory if it doesn't exist yet
    RunDirs="${OutputPath}/${RunName}"
    if [ ! -d "${RunDirs}" ]; then
        mkdir -p -v ${RunDirs}
    fi

    ##=======================================================================
    ## Create state vector file
    ##=======================================================================

    if "$CreateAutomaticRectilinearStateVectorFile"; then
        create_statevector
    else
        # Copy custom state vector to $RunDirs directory for later use
        cp $StateVectorFile ${RunDirs}/StateVector.nc
    fi

    # Determine number of elements in state vector file
    nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc) 
    printf "\nNumber of state vector elements in this inversion = ${nElements}\n\n"

    # Define inversion domain lat/lon bounds
    LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
    LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
    LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
    LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)
    Lons="${LonMinInvDomain}, ${LonMaxInvDomain}"
    Lats="${LatMinInvDomain}, ${LatMaxInvDomain}"

    ##=======================================================================
    ## Set up template run directory
    ##=======================================================================
    runDir="template_run"
    RunTemplate="${RunDirs}/${runDir}"
    if "$SetupTemplateRundir"; then
        setup_template
    fi

    ##=======================================================================
    ## Reduce state vector dimension
    ##=======================================================================
    if "$ReducedDimensionStateVector"; then
        reduce_dimension
    fi

    ##=======================================================================
    ##  Set up IMI preview run directory
    ##=======================================================================
    preview_start=$(date +%s)
    if  "$DoPreview"; then
        run_preview
    fi 
    preview_end=$(date +%s)

    ##=======================================================================
    ##  Set up spinup run directory
    ##=======================================================================
    if  "$SetupSpinupRun"; then
        setup_spinup
    fi

    ##=======================================================================
    ##  Set up posterior run directory
    ##=======================================================================
    if  "$SetupPosteriorRun"; then
        setup_posterior
    fi

    ##=======================================================================
    ##  Set up Jacobian run directories
    ##=======================================================================
    if "$SetupJacobianRuns"; then
        setup_jacobian
    fi

    ##=======================================================================
    ##  Setup inversion directory
    ##=======================================================================
    if "$SetupInversion"; then
        setup_inversion
    fi

    # Remove temporary files
    if "$isAWS"; then
        rm -f /home/ubuntu/foo.nc
    fi

    # Run time
    echo "Statistics (setup):"
    echo "Preview runtime (s): $(( $preview_end - $preview_start ))"
    echo "Note: this is part of the Setup runtime reported by run_imi.sh"
    printf "\n=== DONE RUNNING SETUP SCRIPT ===\n"
}

# Description: Turn on switches for extra observation operators
#   Works on geoschem_config.yml file in the current directory  
# Usage:
#   activate_observations
activate_observations() {
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

}


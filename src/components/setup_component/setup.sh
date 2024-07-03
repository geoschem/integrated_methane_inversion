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

    ##=======================================================================
    ## Standard settings
    ##=======================================================================
    # Start and end date for the spinup simulation
    SpinupStart=$(date --date="${StartDate} -${SpinupMonths} month" +%Y%m%d)
    SpinupEnd=${StartDate}

    # Use global boundary condition files for initial conditions
    UseBCsForRestart=true    

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
            aws s3 cp --no-sign-request s3://imi-boundary-conditions/GEOSChem.BoundaryConditions.${SpinupStart}_0000z.nc4 $RestartFile
        fi
        RestartFilePreview=${RestartFilePreviewPrefix}${StartDate}_0000z.nc4
        if [ ! -f "$RestartFilePreview" ]; then
            aws s3 cp --no-sign-request s3://imi-boundary-conditions/GEOSChem.BoundaryConditions.${StartDate}_0000z.nc4 $RestartFilePreview
        fi
    fi

    ##=======================================================================
    ## Define met and grid fields for HEMCO_Config.rc
    ##=======================================================================
    if [[ "$Met" == "GEOSFP" || "$Met" == "GEOS-FP" || "$Met" == "geosfp" ]]; then
        metDir="GEOS_FP"
        native="0.25x0.3125"
        constYr="2011"
        LandCoverFileExtension="nc"
    elif [[ "$Met" == "MERRA2" || "$Met" == "MERRA-2" || "$Met" == "merra2" ]]; then
        metDir="MERRA2"
        native="0.5x0.625"
        constYr="2015"
        LandCoverFileExtension="nc4"
    else
	printf "\nERROR: Meteorology field ${Met} is not supported by the IMI. "
	printf "\n Options are GEOSFP or MERRA2.\n"
	exit 1
    fi

    if [ "$Res" == "0.25x0.3125" ]; then
        gridDir="${Res}"
        gridFile="025x03125"
    elif [ "$Res" == "0.5x0.625" ]; then
        gridDir="${Res}" 
        gridFile="05x0625"
    elif [ "$Res" == "2.0x2.5" ]; then
        gridDir="2x2.5"
        gridFile="2x25"
    elif [ "$Res" = "4.0x5.0" ]; then
        gridDir="4x5"
        gridFile="4x5"
    else
	printf "\nERROR: Grid resolution ${Res} is not supported by the IMI. "
	printf "\n Options are 0.25x0.3125, 0.5x0.625, 2.0x2.5, or 4.0x5.0.\n"
	exit 1
    fi
    # Use cropped met for regional simulations instead of using global met
    if "$isRegional"; then
        gridDir="${gridDir}_${RegionID}"
    fi

    # Clone defined version of GCClassic
    # Define path to GEOS-Chem run directory files
    cd "${InversionPath}"
    if [ ! -d "GCClassic" ]; then
        git clone https://github.com/geoschem/GCClassic.git
        cd GCClassic
        git checkout ${GEOSCHEM_VERSION}
        git submodule update --init --recursive
        cd ..
    else
        cd GCClassic
        if grep -Fq "VERSION ${GEOSCHEM_VERSION}" CMakeLists.txt; then
            printf "\nGCClassic already exists and is the correct version ${GEOSCHEM_VERSION}.\n"
        else
            printf "\nERROR: GCClassic already exists but is not version ${GEOSCHEM_VERSION}.\n"
            exit 1
        fi
        cd ..
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
    ## Create or copy state vector file
    ##=======================================================================

    if "$CreateAutomaticRectilinearStateVectorFile"; then
        create_statevector
    else
        # Copy custom state vector to $RunDirs directory for later use
        printf "\nCopying state vector file\n"
        cp -v $StateVectorFile ${RunDirs}/StateVector.nc
    fi

    # Determine number of elements in state vector file
    nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc)
    if "$OptimizeBCs"; then
	nElements=$((nElements+4))
    fi
    if "$OptimizeOH";then
	nElements=$((nElements+1))
    fi
    printf "\nNumber of state vector elements in this inversion = ${nElements}\n\n"

    # Define inversion domain lat/lon bounds
    LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
    LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
    LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
    LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)

    ##=======================================================================
    ## Set up template run directory
    ##=======================================================================
    runDir="template_run"
    RunTemplate="${RunDirs}/${runDir}"
    if "$SetupTemplateRundir"; then
        setup_template
    fi

    ##=======================================================================
    ## Generate Prior Emissions
    ##=======================================================================
    if "$DoPriorEmis"; then
       run_prior
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

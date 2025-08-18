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
            python src/utilities/download_aws_file.py s3://imi-boundary-conditions/${BCversion}/GEOSChem.BoundaryConditions.${SpinupStart}_0000z.nc4 $RestartFile
        fi
    fi

    ##=======================================================================
    ## Define met and grid fields for HEMCO_Config.rc
    ##=======================================================================
    if [[ "$Met" == "GEOSFP" || "$Met" == "GEOS-FP" || "$Met" == "geosfp" ]]; then
        metDir="GEOS_FP"
        native="0.25x0.3125"
        metsuffix="025x03125"
        constYr="2011"
        LandCoverFileExtension="nc"
    elif [[ "$Met" == "MERRA2" || "$Met" == "MERRA-2" || "$Met" == "merra2" ]]; then
        metDir="MERRA2"
        native="0.5x0.625"
        metsuffix="05x0625"
        constYr="2015"
        LandCoverFileExtension="nc4"
    else
        printf "\nERROR: Meteorology field ${Met} is not supported by the IMI. "
        printf "\n Options are GEOSFP or MERRA2.\n"
        exit 1
    fi

    if [ "$UseGCHP" != "true" ]; then
        if [ "$Res" == "0.125x0.15625" ]; then
            gridDir="${Res}"
            gridFile="0125x015625" 
        elif [ "$Res" == "0.25x0.3125" ]; then
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
            printf "\n Options are 0.125x0.15625, 0.25x0.3125, 0.5x0.625, 2.0x2.5, or 4.0x5.0.\n"
            exit 1
        fi
    else
        gridDir="${native}"
        gridFile="${metsuffix}"
    fi
    # Use cropped met for regional simulations instead of using global met
    if [ "$RegionID" != "" ]; then
        gridDir="${gridDir}_${RegionID}"
    fi

    # Clone defined version of GCHP or GCClassic
    # Define path to GEOS-Chem run directory files
    cd "${InversionPath}"
    if "$UseGCHP"; then
        if [ ! -d "GCHP" ]; then
            git clone https://github.com/geoschem/GCHP.git
            cd GCHP
            git checkout ${GEOSCHEM_VERSION}
            git submodule update --init --recursive
            cd ..
        else
            cd GCHP
            if grep -Fq "VERSION ${GEOSCHEM_VERSION}" CMakeLists.txt; then
                printf "\nGCHP already exists and is the correct version ${GEOSCHEM_VERSION}.\n"
            else
                printf "\nERROR: GCHP already exists but is not version ${GEOSCHEM_VERSION}.\n"
                exit 1
            fi
            cd ..
        fi
    else
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
    fi

    # Define path to GEOS-Chem run directory files
    GCHPPath="${InversionPath}/GCHP"
    GCClassicPath="${InversionPath}/GCClassic"
    
    # Create working directory if it doesn't exist yet
    RunDirs="${OutputPath}/${RunName}"
    if [ ! -d "${RunDirs}" ]; then
        mkdir -p -v ${RunDirs}
    fi

    ##=======================================================================
    ## Create regridding weights to CS grid
    ##=======================================================================
    if "$UseGCHP"; then
        CSgridDir="${RunDirs}/CS_grids"
        if [ ! -d "$CSgridDir" ]; then
            mkdir -p -v "$CSgridDir"
        fi
        cd "$CSgridDir"

        prefix=$(get_GridSpec_prefix "$CS_RES" "$STRETCH_GRID" "$STRETCH_FACTOR" "$TARGET_LAT" "$TARGET_LON")
        gridspec_fname="${CSgridDir}/${prefix}_gridspec.nc"
        if [ -f "$gridspec_fname" ]; then
            echo "GridSpec file of already exists: $gridspec_fname"
        else
            echo "Generating grid spec: $gridspec_fname"
            if "$STRETCH_GRID"; then
                gridspec-create sgcs -s "${STRETCH_FACTOR}" -t "${TARGET_LAT}" "${TARGET_LON}" "$CS_RES" > /dev/null 2>&1
            else
                gridspec-create gcs "$CS_RES" > /dev/null 2>&1
            fi
        fi

        # Create CS grid file only if it doesn't exist
        gridfpath="${CSgridDir}/grids.c${CS_RES}.nc"
        if [ -f "$gridfpath" ]; then
            echo "CS grid file already exists: $gridfpath"
        else
            echo "Creating CS grid file: $gridfpath"
            generate_grid_from_GridSpec "$CS_RES" "$gridfpath" "$STRETCH_GRID" "$STRETCH_FACTOR" "$TARGET_LAT" "$TARGET_LON"
        fi

        # Generate regridding weights only if not already present
        dst_grid=$gridspec_fname
        if "$STRETCH_GRID"; then
            weightsfile="regrid_weights_${native}_to_c${CS_RES}_s${STRETCH_FACTOR}_${TARGET_LAT}N_${TARGET_LON}E_conserve.nc"
        else
            weightsfile="regrid_weights_${native}_to_c${CS_RES}_conserve.nc"
        fi

        regridding_method="conserve"

        if [ -f "$weightsfile" ]; then
            echo "Regridding weights file already exists: $weightsfile"
        else
            if [ "$metDir" = "GEOS_FP" ]; then
                gridspec-create latlon -b -180 -90 180 90 -pc -hp -dc -o . 721 1152 > /dev/null 2>&1
                src_grid="regular_lat_lon_721x1152.nc"
            elif [ "$metDir" = "MERRA2" ]; then
                gridspec-create latlon -b -180 -90 180 90 -pc -hp -dc -o . 361 576 > /dev/null 2>&1
                src_grid="regular_lat_lon_361x576.nc"
            else
                echo "Unsupported metDir: $metDir"
                exit 1
            fi

            echo "Generating regridding weights file: $weightsfile"
            ESMF_RegridWeightGen -s "$src_grid" -d "$dst_grid" -m "$regridding_method" -w "$weightsfile" > /dev/null 2>&1
            rm -f PET*
        fi
        cd "${InversionPath}"
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
        nElements=$((nElements + 4))
    fi
    if "$OptimizeOH"; then
        if "$isRegional"; then
            nElements=$((nElements + 1))
        else
            nElements=$((nElements + 2))
        fi
    fi
    printf "\nNumber of state vector elements in this inversion = ${nElements}\n\n"

    # Define inversion domain lat/lon bounds
    if [ "$UseGCHP" != "true" ]; then
        LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
        LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
        LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
        LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)
    else
        LonMinInvDomain=-180
        LonMaxInvDomain=180
        LatMinInvDomain=-90
        LatMaxInvDomain=90
    fi

    # Define custom Kalman filter periods
    if "$KalmanMode"; then
        if ! "$MakePeriodsCSV"; then
            printf "Copying custom periods.csv to the run directory.\n"
            cp $CustomPeriodsCSV ${RunDirs}/
        fi
    fi

    ##=======================================================================
    ## Set up template run directory
    ##=======================================================================
    runDir="template_run"
    RunTemplate="${RunDirs}/${runDir}"
    if "$SetupTemplateRundir"; then
        setup_template
    fi

    ##=======================================================================
    ## Generate Prior Emissions using a HEMCO standalone run or GCHP prior run
    ##=======================================================================
    if "$DoHemcoPriorEmis"; then
        if [ "$UseGCHP" != "true" ]; then
            run_hemco_prior_emis
        else
            setup_prior_gchp
            run_prior_gchp $StartDate $EndDate
        fi
    fi

    ##=======================================================================
    ## Reduce state vector dimension
    ##=======================================================================
    if "$ReducedDimensionStateVector"; then
        reduce_dimension # to do: adapt to GCHP
    fi

    ##=======================================================================
    ##  Run the IMI preview
    ##=======================================================================
    preview_start=$(date +%s)
    if "$DoPreview"; then
        run_preview
    fi
    preview_end=$(date +%s)

    ##=======================================================================
    ##  Set up spinup run directory
    ##=======================================================================
    if "$SetupSpinupRun"; then
        setup_spinup
    fi

    ##=======================================================================
    ##  Set up posterior run directory
    ##=======================================================================
    if "$SetupPosteriorRun"; then
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

    # Run time
    echo "Statistics (setup):"
    echo "Preview runtime (s): $(($preview_end - $preview_start))"
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
    if "$DoObsPack"; then
        sed -i "/obspack/{N;s/activate: false/activate: true/}" geoschem_config.yml
        ln -s ${DataPath}/Observations/ObsPack ObsPack
        OLD="input_file: .\/obspack_co2_1_OCO2MIP_2018-11-28.YYYYMMDD.nc"
        NEW="input_file: .\/ObsPack\/CH4\/YYYY\/obspack_ch4.YYYYMMDD.nc"
        sed -i "s|$OLD|$NEW|g" geoschem_config.yml
    fi 

}

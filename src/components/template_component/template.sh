#!/bin/bash

# Functions available in this file include:
#   - setup_template

# Description: Setup template GCClassic run directory
# Usage:
#   setup_template
setup_template() {
    printf "\n=== CREATING TEMPLATE RUN DIRECTORY ===\n"

    if "$UseGCHP"; then
        cd ${GCHPPath}/run
    else
        cd ${GCClassicPath}/run
    fi

    # The createRunDir.sh script assumes the file ~/.geoschem/config exists
    # and contains the path to GEOS-Chem input data
    export GC_USER_REGISTERED=true
    if [[ ! -f ${HOME}/.geoschem/config ]]; then
        mkdir -p ${HOME}/.geoschem
        echo "export GC_DATA_ROOT=${DataPath}" >>${HOME}/.geoschem/config
        source ${HOME}/.geoschem/config
    fi

    if [[ -d ${RunTemplate} ]]; then
        printf "\nERROR: ${RunTemplate} already exists. Please remove or set 'SetupTemplateRunDir: false' in config.yml.\n"
        exit 9999
    fi

    # Commands to feed to createRunDir.sh
    if [[ "$Met" == "MERRA2" || "$Met" == "MERRA-2" || "$Met" == "merra2" ]]; then
        metNum="1"
    elif [[ "$Met" == "GEOSFP" || "$Met" == "GEOS-FP" || "$Met" == "geosfp" ]]; then
	# Add y to metNum to skip prompt about GEOS-FP discontinuity
        metNum="2\ny"
    else
        printf "\nERROR: Meteorology field ${Met} is not supported by the IMI. "
        printf "\n Options are GEOSFP or MERRA2.\n"
        exit 1
    fi

    if "$UseGCHP"; then
        if [ "${metNum}" == "1" ]; then
            cmd="5\n2\n${metNum}\n${RunDirs}\n${runDir}\nn\n"
        else
            # GEOSFP: Use daily files pre-processed for GEOS-Chem
            cmd="5\n2\n${metNum}\ny\n1\n1\n${RunDirs}\n${runDir}\nn\n"
        fi
    else
        if [ "$Res" = "4.0x5.0" ]; then
            cmd="9\n${metNum}\n1\n2\n${RunDirs}\n${runDir}\nn\n"
        elif [ "$Res" == "2.0x2.5" ]; then
            cmd="9\n${metNum}\n2\n2\n${RunDirs}\n${runDir}\nn\n"
        elif [ "$Res" == "0.5x0.625" ]; then
            if "$isRegional"; then
                # Use NA domain by default and adjust lat/lon below
                cmd="9\n${metNum}\n3\n4\n2\n${RunDirs}\n${runDir}\nn\n"
            else
                cmd="9\n${metNum}\n3\n1\n2\n${RunDirs}\n${runDir}\nn\n"
            fi
        elif [ "$Res" == "0.25x0.3125" ]; then
            if "$isRegional"; then
                # Use NA domain by default and adjust lat/lon below
                cmd="9\n${metNum}\n4\n4\n2\n${RunDirs}\n${runDir}\nn\n"
            else
                cmd="9\n${metNum}\n4\n1\n2\n${RunDirs}\n${runDir}\nn\n"
            fi
        elif [ "$Res" == "0.125x0.15625" ]; then
            if "$isRegional"; then
                # Use NA domain by default and adjust lat/lon below
                cmd="9\n${metNum}\n5\n4\n2\n${RunDirs}\n${runDir}\nn\n" #regional run
            else
                cmd="9\n${metNum}\n5\n1\n2\n${RunDirs}\n${runDir}\nn\n"
            fi
        else
            printf "\nERROR: Grid resolution ${Res} is not supported by the IMI. "
            printf "\n Options are 0.125x0.15625, 0.25x0.3125, 0.5x0.625, 2.0x2.5, or 4.0x5.0.\n"
            exit 1
        fi
    fi

    # Create run directory
    printf ${cmd} | ./createRunDir.sh >>createRunDir.log 2>&1
    rm -f createRunDir.log
    printf "\nCreated ${RunTemplate}\n"

    cd ${RunTemplate}

    # Copy download script to run directory
    cp ${InversionPath}/src/utilities/download_gc_data.py download_gc_data.py

    if "$UseGCHP"; then
        RunDuration=$(get_run_duration "$StartDate" "$EndDate")
        echo "$StartDate 000000" > cap_restart
        sed -i -e "s/Run_Duration=\"[0-9]\{8\} 000000\"/Run_Duration=\"${RunDuration} 000000\"/" \
            -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
            -e "s/^TOTAL_CORES=.*/TOTAL_CORES=${TOTAL_CORES}/" \
            -e "s/^NUM_NODES=.*/NUM_NODES=${NUM_NODES}/" \
            -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=${NUM_CORES_PER_NODE}/" \
            -e 's/^AutoUpdate_Diagnostics=.*$/AutoUpdate_Diagnostics=OFF/' \
            setCommonRunSettings.sh
        sed -i -e "s/monthly:[[:space:]]*1/monthly:        0/g" HISTORY.rc
    else
        # Modify geoschem_config.yml based on settings in config.yml
        sed -i -e "s:20190101:${StartDate}:g" \
            -e "s:20190201:${EndDate}:g" geoschem_config.yml
    fi

    if "$isRegional"; then
        # Adjust lat/lon bounds because GEOS-Chem defines the domain
        # based on grid cell edges (not centers) for the lat/lon bounds
        Lons=$(calculate_geoschem_domain lon ${RunDirs}/StateVector.nc ${LonMinInvDomain} ${LonMaxInvDomain})
        Lats=$(calculate_geoschem_domain lat ${RunDirs}/StateVector.nc ${LatMinInvDomain} ${LatMaxInvDomain})
        sed -i -e "s:-130.0,  -60.0:${Lons}:g" \
            -e "s:9.75,  60.0:${Lats}:g" geoschem_config.yml
    fi

    # Update time cycling flags to use most recent year
    sed -i "s/RF xy/C xy/g" HEMCO_Config.rc

    # Modify path to state vector file in HEMCO_Config.rc
    OLD=" ./StateVector.nc"
    NEW=" ${RunDirs}/StateVector.nc"
    sed -i -e "s@$OLD@$NEW@g" HEMCO_Config.rc

    if "$KalmanMode"; then
        jacobian_period=${period_i}
    else
        jacobian_period=1
    fi
    scale_OLD=" ./gridded_pert_scale_1.nc"
    scale_NEW=" ${RunDirs}/archive_perturbation_sfs/gridded_pert_scale_${jacobian_period}.nc"
    sed -i -e "s@$scale_OLD@$scale_NEW@g" HEMCO_Config.rc

    if [ "$UseGCHP" = true ]; then
        # Too long file name could be problematic in GCHP
        ln -nsf "${RunDirs}" RunDirs
        NEW_Ext=" ./RunDirs/StateVector.nc"
        sed -i -e "s@$OLD@$NEW_Ext@g" ExtData.rc
        scale_NEW_Ext=" ./RunDirs/archive_perturbation_sfs/gridded_pert_scale_${jacobian_period}.nc"
        sed -i -e "s@$scale_OLD@$scale_NEW_Ext@g" ExtData.rc
    fi

    # Modify HEMCO_Config.rc if running Kalman filter
    if "$KalmanMode"; then
        sed -i -e "s|gridded_posterior.nc|${RunDirs}/ScaleFactors.nc|g" HEMCO_Config.rc
        if "$UseGCHP"; then
            sed -i -e "s|gridded_posterior.nc|${RunDirs}/ScaleFactors.nc|g" ExtData.rc
        fi
    fi

    # Modify HEMCO_Config.rc based on settings in config.yml
    # Use cropped met fields (add the region to both METDIR and the met files)
    if [ "$RegionID" != "" ]; then
	if [ "$Res" != "0.125x0.15625" ]; then
           sed -i -e "s:GEOS_${Res}_NA:GEOS_${Res}_${RegionID}:g" HEMCO_Config.rc.gmao_metfields
           sed -i -e "s:\$RES.NA:\$RES.${RegionID}:g" HEMCO_Config.rc.gmao_metfields
        # Modify the METDIR for 0.125x0.15625 simulation
        elif [ "$Res" = "0.125x0.15625" ]; then
           sed -i -e "s:GEOS_0.25x0.3125_NA\/GEOS_FP:GEOS_0.25x0.3125_${RegionID}\/GEOS_FP:g" HEMCO_Config.rc.gmao_metfields_0125
           OLD="GEOS_0.125x0.15625_NA/GEOS_FP"
           NEW="GEOS_0.125x0.15625_${RegionID}/GEOS_FP_DerivedWinds"
           sed -i "s|$OLD|$NEW|g" HEMCO_Config.rc.gmao_metfields_0125
           sed -i '/METDIR/d' HEMCO_Config.rc
	fi
   fi

    # By default, only output emissions at the end of the simulation
    sed -i -e "s|DiagnFreq:                   Monthly|DiagnFreq:                   End|g" HEMCO_Config.rc
    # do not output Emissions collection
    if "$UseGCHP"; then
        sed -i "s/'Emissions',/#'Emissions',/" HISTORY.rc
    fi
    # Add a new ZERO scale factor for use in jacobian simulations
    sed -i -E '/^1[[:space:]]+NEGATIVE[[:space:]]+-1\.0([[:space:]]+-){3}[[:space:]]+xy[[:space:]]+1[[:space:]]+1/a 5 ZERO      0.0 - - - xy 1 1' HEMCO_Config.rc

    # Modify path to BC files
    sed -i -e "s:\$ROOT/SAMPLE_BCs/v2021-07/CH4:${fullBCpath}:g" HEMCO_Config.rc

    # If reading total prior emissions (as in the jacobian and posterior), read a new file each month
    sed -i -e "s|EmisCH4_Total \$YYYY/\$MM/\$DD/0|EmisCH4_Total 1900-2050/1-12/1-31/0|g" HEMCO_Config.rc

    # Temporary fix: Modify path to HEMCO prior emissions (the path is currently
    # hardcoded in the template HEMCO config file in GEOS-Chem)
    sed -i -e "s|prior_run|hemco_prior_emis|g" HEMCO_Config.rc
    if "$UseGCHP"; then
        sed -i -e "s|prior_run|hemco_prior_emis|g" ExtData.rc
    fi
    
    # Modify HISTORY.rc - comment out diagnostics that aren't needed
    sed -i -e "s:'CH4':#'CH4':g" \
        -e "s:'Metrics:#'Metrics:g" \
        -e "s:'StateMet:#'StateMet:g" HISTORY.rc

    # If turned on, save out hourly CH4 concentrations to daily files
    # use time-average mode
    if "$HourlyCH4"; then
        if "$UseGCHP"; then
            sed -i -e 's/SpeciesConc.frequency:.*/SpeciesConc.frequency:      010000/g' \
                -e 's/SpeciesConc.duration:.*/SpeciesConc.duration:       240000/g' HISTORY.rc
        else
            sed -i -e 's/SpeciesConc.frequency:      00000100 000000/SpeciesConc.frequency:      00000000 010000/g' \
                -e 's/SpeciesConc.duration:       00000100 000000/SpeciesConc.duration:       00000001 000000/g' HISTORY.rc
        fi
    fi

    # Remove sample restart file; GCHP restarts are just soft links and are needed later as template
    if [ "$UseGCHP" != true ]; then
        rm -f Restarts/GEOSChem.Restart.20190101_0000z.nc4
    fi

    # Copy template run script
    if "$UseGCHP"; then
        cp ${InversionPath}/src/geoschem_run_scripts/gchp_ch4_run.template .
    else
        cp ${InversionPath}/src/geoschem_run_scripts/ch4_run.template .
    fi

    # Compile GEOS-Chem and store executable in GEOSChem_build directory
    if [[ -f "../GEOSChem_build/gchp" || -f "../GEOSChem_build/gcclassic" ]]; then
        printf "\nGEOS-Chem executable is already built and stored in GEOSChem_build\n"
        rm -rf build
    else
        cd build
        if "$UseGCHP"; then
            printf "\nCompiling GCHP...\n"
            cmake ${InversionPath}/GCHP >>build_geoschem.log 2>&1
        else
            printf "\nCompiling GEOS-Chem...\n"
            cmake ${InversionPath}/GCClassic >>build_geoschem.log 2>&1
        fi
        cmake . -DRUNDIR=.. -DMECH=carbon >>build_geoschem.log 2>&1
        make -j install >>build_geoschem.log 2>&1
        cd ..
        if "$UseGCHP"; then
            if [[ -f gchp ]]; then
                mkdir ../GEOSChem_build
                mv -v gchp ../GEOSChem_build/
                mv build/CMakeCache.txt ../GEOSChem_build
                mv build/build_geoschem.log ../GEOSChem_build
                rm -rf build
            else
                printf "\nGCHP build failed! \n\nSee ${RunDirs}/GEOSChem_build/build_geoschem.log for details\n"
                exit 999
            fi
        else
            if [[ -f gcclassic ]]; then
                mv build_info ../GEOSChem_build
                mv -v gcclassic ../GEOSChem_build/
                mv build/build_geoschem.log ../GEOSChem_build
                rm -rf build
            else
                printf "\nGEOS-Chem build failed! \n\nSee ${RunDirs}/GEOSChem_build/build_geoschem.log for details\n"
                exit 999
            fi
        fi
        printf "\nDone compiling GEOS-Chem \n\nSee ${RunDirs}/GEOSChem_build for details\n\n"
    fi
    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING TEMPLATE RUN DIRECTORY ===\n"
}

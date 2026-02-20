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
            cmd="3\n2\n${metNum}\n${RunDirs}\n${runDir}\nn\n"
        else
            # GEOSFP: Use daily files pre-processed for GEOS-Chem
            cmd="3\n2\n${metNum}\n1\n1\n${RunDirs}\n${runDir}\nn\n"
        fi
    else
        if [ "$Res" = "4.0x5.0" ]; then
            cmd="3\n2\n${metNum}\n1\n2\n${RunDirs}\n${runDir}\nn\n"
        elif [ "$Res" == "2.0x2.5" ]; then
            cmd="3\n2\n${metNum}\n2\n2\n${RunDirs}\n${runDir}\nn\n"
        elif [ "$Res" == "0.5x0.625" ]; then
            if "$isRegional"; then
                # Use NA domain by default and adjust lat/lon below
                cmd="3\n2\n${metNum}\n3\n4\n2\n${RunDirs}\n${runDir}\nn\n"
            else
                cmd="3\n2\n${metNum}\n3\n1\n2\n${RunDirs}\n${runDir}\nn\n"
            fi
        elif [ "$Res" == "0.25x0.3125" ]; then
            if "$isRegional"; then
                # Use NA domain by default and adjust lat/lon below
                cmd="3\n2\n${metNum}\n4\n4\n2\n${RunDirs}\n${runDir}\nn\n"
            else
                cmd="3\n2\n${metNum}\n4\n1\n2\n${RunDirs}\n${runDir}\nn\n"
            fi
        elif [ "$Res" == "0.125x0.15625" ]; then
            if "$isRegional"; then
                # Use NA domain by default and adjust lat/lon below
                cmd="3\n2\n${metNum}\n5\n4\n2\n${RunDirs}\n${runDir}\nn\n" #regional run
            else
                cmd="3\n2\n${metNum}\n5\n1\n2\n${RunDirs}\n${runDir}\nn\n"
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
        sed -i -e "s/^BEG_DATE:.*/BEG_DATE:     ${StartDate} 000000/" \
            -e "s/^END_DATE:.*/END_DATE:     ${EndDate} 000000/" CAP.rc
        echo "$StartDate 000000" > cap_restart
        sed -i -e "s/Run_Duration=\"[0-9]\{8\} 000000\"/Run_Duration=\"${RunDuration} 000000\"/" \
            -e "s/^CS_RES=.*/CS_RES=${CS_RES}/" \
            -e "s/^TOTAL_CORES=.*/TOTAL_CORES=${TOTAL_CORES}/" \
            -e "s/^NUM_NODES=.*/NUM_NODES=${NUM_NODES}/" \
            -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=${NUM_CORES_PER_NODE}/" \
            -e 's/^AutoUpdate_Diagnostics=.*$/AutoUpdate_Diagnostics=OFF/' \
            setCommonRunSettings.sh
        # turn on monthly checkpoint
        sed -i -e 's/^Midrun_Checkpoint=OFF/Midrun_Checkpoint=ON/' \
            -e 's/^Checkpoint_Freq=.*/Checkpoint_Freq=monthly/' \
            setCommonRunSettings.sh
        sed -i -e "s/monthly:[[:space:]]*1/monthly:        0/g" HISTORY.rc
        if "$STRETCH_GRID"; then
            sed -i -e "s/^STRETCH_GRID=.*/STRETCH_GRID=ON/" \
                -e "s/^STRETCH_FACTOR=.*/STRETCH_FACTOR=${STRETCH_FACTOR}/" \
                -e "s/^TARGET_LAT=.*/TARGET_LAT=${TARGET_LAT}/" \
                -e "s/^TARGET_LON=.*/TARGET_LON=${TARGET_LON}/" \
                setCommonRunSettings.sh
        fi
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

    # Too long file name could be problematic
    ln -nsf ${RunDirs} RunDirs
    # Modify path to state vector file in HEMCO_Config.rc
    OLD=" ./StateVector.nc"
    NEW=" ./RunDirs/StateVector.nc"
    sed -i -e "s@$OLD@$NEW@g" HEMCO_Config.rc

    if "$KalmanMode"; then
        jacobian_period=${period_i}
    else
        jacobian_period=1
    fi

    scale_OLD=" ./gridded_pert_scale_1.nc"
    scale_NEW=" ./RunDirs/archive_perturbation_sfs/gridded_pert_scale_${jacobian_period}.nc"
    sed -i -e "s@$scale_OLD@$scale_NEW@g" HEMCO_Config.rc

    if [ "$UseGCHP" = true ]; then
        sed -i -e "s@^#SCALE_PERT@SCALE_PERT@g" ExtData.rc
        NEW_Ext=" ./RunDirs/StateVector.nc"
        sed -i -e "s@$OLD@$NEW_Ext@g" ExtData.rc
        scale_NEW_Ext=" ./RunDirs/archive_perturbation_sfs/gridded_pert_scale_${jacobian_period}.nc"
        sed -i -e "s@$scale_OLD@$scale_NEW_Ext@g" ExtData.rc
    fi

    # Modify HEMCO_Config.rc if running Kalman filter
    if "$KalmanMode"; then
        sed -i -e "s|gridded_posterior.nc|./RunDirs/ScaleFactors.nc|g" HEMCO_Config.rc
        if "$UseGCHP"; then
            sed -i -e "s|gridded_posterior.nc|./RunDirs/ScaleFactors.nc|g" ExtData.rc
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
        sed -i -e "s|^#CH4_Emis_Prior|CH4_Emis_Prior|g" ExtData.rc
        sed -i -e "s|prior_run|hemco_prior_emis|g" ExtData.rc
    fi
    
    # Modify HISTORY.rc - comment out diagnostics that aren't needed
    sed -i -e "s:'Carbon':#'Carbon':g" \
        -e "s:'Metrics':#'Metrics':g" \
        -e "s:'StateMet':#'StateMet':g" HISTORY.rc

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

    # Compile GEOS-Chem/GCHP and store in <src>/build/bin directory
    if "$UseGCHP"; then
        printf "\nCompiling GCHP...\n"
        execname="gchp"
        build_dir=${GCHPPath}/build
        KPP_dir=${GCHPPath}/src/GCHP_GridComp/GEOSChem_GridComp/geos-chem/KPP
    else
        printf "\nCompiling GEOS-Chem...\n"
        execname="gcclassic"
        build_dir=${GCClassicPath}/build
        KPP_dir=${GCClassicPath}/src/GEOS-Chem/KPP
    fi

    mkdir -p ${build_dir}
    cd ${build_dir}
    # first build a default exexutable without Jacobian tracers
    if [ -f "bin/${execname}.default" ]; then
        echo "Executable bin/${execname}.default already exists — skipping rebuild."
    else
        echo "Building ${execname}.default ..."

        # remove CMakeCache.txt once to initialize JACOBIAN cmake option
        rm -f CMakeCache.txt
        cd "${KPP_dir}/carbon"
        if [ ! -f carbon.eqn.default ]; then
            mv carbon.eqn carbon.eqn.default
        fi
        ln -nsf carbon.eqn.default carbon.eqn
        # initialize KPP mechanism to be the default carbon equations
        cd ${KPP_dir}
        ./build_mechanism.sh carbon >> "${build_dir}/build_geoschem.log" 2>&1 \
            || { echo "ERROR: build_mechanism.sh carbon failed." >&2; exit 1; }

        cd "${build_dir}"
        cmake .. >> build_geoschem.log 2>&1
        cmake . -DMECH=carbon -DJACOBIAN=n >> build_geoschem.log 2>&1
        make -j >> build_geoschem.log 2>&1

        mv "bin/${execname}" "bin/${execname}.default"
    fi

    # then expand a series of carbon equations for Jacobian tracers
    # sanity check on NumJacobianTracers
    if ! [[ "$NumJacobianTracers" =~ ^[0-9]+$ ]] || [ "$NumJacobianTracers" -eq 0 ]; then
        echo "ERROR: NumJacobianTracers must be a positive integer > 0." >&2
        exit 1
    fi
    # determine number of executables to be built based on NumJacobianTracers
    # get the ceiling number rounded to 10s
    upper=$(( (NumJacobianTracers + 9) / 10 * 10 ))

    cd "${build_dir}"
    cmake . -DMECH=carbon -DJACOBIAN=y >> build_geoschem.log 2>&1
    for n in $(seq 10 10 $upper); do
        if [ -f "bin/${execname}.${n}" ]; then
            echo "Executable bin/${execname}.${n} already exists — skipping rebuild."
        else
            echo "Building ${execname}.${n} ..."

            cd ${KPP_dir}/carbon
            ${InversionPath}/src/utilities/expand_carbon_eqn.py \
                carbon.eqn.default ${n} > carbon.eqn.${n}
            ln -nsf carbon.eqn.${n} carbon.eqn
            # generate KPP carbon mechanism
            cd ${KPP_dir}
            ./build_mechanism.sh carbon >> "${build_dir}/build_geoschem.log" 2>&1 \
                || { echo "ERROR: build_mechanism.sh carbon failed." >&2; exit 1; }
            # build GCC/GCHP with expanded carbon mechanism for Jacobian tracers
            cd ${build_dir}
            make -j >>build_geoschem.log 2>&1
            mv bin/${execname} bin/${execname}.${n}
        fi
    done

    printf "\nDone compiling GEOS-Chem \n\nSee ${build_dir}/build_geoschem.log for details\n\n"

    # Navigate back to working directory
    cd ${RunDirs}

    printf "\n=== DONE CREATING TEMPLATE RUN DIRECTORY ===\n"
}

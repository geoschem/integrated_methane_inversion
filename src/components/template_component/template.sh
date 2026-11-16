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

    #---------------------------------------------------------------------
    # Create GEOS-Chem run directory 
    #---------------------------------------------------------------------
    simNum=3  # Carbon simulation

    # Species
    if [[ "$Species" == "CH4" ]]; then
    	spcNum=2
    elif [[ "$Species" == "CO2" ]]; then
    	spcNum=3
    elif [[ "$Species" == "CH4_CO2" ]]; then
    	spcNum=1
    else
    	printf "\nERROR: Species ${Species} is not supported by the IMI. "
    	printf "\n Options are CH4, CO2, or CH4_CO2"
    	exit 1
    fi

    # Meteorology field
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
	# GCHP: Hardcode command to feed to createrundir.sh
        if [ "${metNum}" == "1" ]; then
            cmd="5\n2\n${metNum}\n${RunDirs}\n${runDir}\nn\n"
        else
            # GEOSFP: Use daily files pre-processed for GEOS-Chem
            cmd="5\n2\n${metNum}\ny\n1\n1\n${RunDirs}\n${runDir}\nn\n"
        fi
    else
	# GCClassic: Allow for grid resolution and regional options
	if [[ "$Res" == "4.0x5.0" ]]; then
            resNum=1
	elif [[ "$Res" == "2.0x2.5" ]]; then
            resNum=2
	elif [[ "$Res" == "0.5x0.625" ]]; then
	    resNum=3
	elif [[ "$Res" == "0.25x0.3125" ]]; then
	    resNum=4
	elif [[ "$Res" == "0.125x0.15625" ]]; then
	    resNum=5
        else
            printf "\nERROR: Grid resolution ${Res} is not supported by the IMI. "
            printf "\n Options are 0.125x0.15625, 0.25x0.3125, 0.5x0.625, 2.0x2.5, or 4.0x5.0.\n"
            exit 1
        fi
	if "$isRegional"; then
	    regionNum=4  # Use NA domain by default and adjust lat/lon below
	else
	    regionNum=1
	fi

	# GCClassic: Command to feed to createRunDir.sh
	if [[ "$Res" == "4.0x5.0" || "$Res" == "2.0x2.5" ]]; then
	    cmd="${simNum}\n${spcNum}\n${metNum}\n${resNum}\n2\n${RunDirs}\n${runDir}\nn\n"
	else
	    cmd="${simNum}\n${spcNum}\n${metNum}\n${resNum}\n${regionNum}\n2\n${RunDirs}\n${runDir}\nn\n"
	fi
    fi
    
    # Create run directory
    printf ${cmd} | ./createRunDir.sh >>createRunDir.log 2>&1
    rm -f createRunDir.log
    printf "\nCreated ${RunTemplate}\n"

    #---------------------------------------------------------------------
    # Modify default settings in run directory
    #---------------------------------------------------------------------
    cd ${RunTemplate}

    # Copy download script to run directory
    cp ${InversionPath}/src/utilities/download_gc_data.py download_gc_data.py

    if "$UseGCHP"; then
	# GCHP: Modify run directory files based on settings in config.yml
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

	# Set stretched grid paramaters from IMI config
        if "$STRETCH_GRID"; then
            sed -i -e "s/^STRETCH_GRID=.*/STRETCH_GRID=ON/" \
                -e "s/^STRETCH_FACTOR=.*/STRETCH_FACTOR=${STRETCH_FACTOR}/" \
                -e "s/^TARGET_LAT=.*/TARGET_LAT=${TARGET_LAT}/" \
                -e "s/^TARGET_LON=.*/TARGET_LON=${TARGET_LON}/" \
                setCommonRunSettings.sh
        fi
    else
	# GCClassic: modify geoschem_config.yml based on settings in config.yml
	if [ "$Res" == "0.125x0.15625" ]; then
	    sed -i -e "s:20230101:${StartDate}:g" \
                   -e "s:20230201:${EndDate}:g" geoschem_config.yml
	else
	    sed -i -e "s:20190101:${StartDate}:g" \
                   -e "s:20190201:${EndDate}:g" geoschem_config.yml
	fi
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
    scale_NEW=" ${RunDirs}/gridded_perturbation_sf.nc"
    sed -i -e "s@$scale_OLD@$scale_NEW@g" HEMCO_Config.rc

    if [ "$UseGCHP" = true ]; then
        sed -i -e "s@^#SCALE_PERT@SCALE_PERT@g" ExtData.rc
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
           sed -i -e "s:GEOS_0.25x0.3125_NA:GEOS_0.25x0.3125_${RegionID}:g" HEMCO_Config.rc.gmao_metfields_0125
           sed -i -e "s:GEOS_0.125x0.15625_NA:GEOS_0.125x0.15625_${RegionID}:g" HEMCO_Config.rc.gmao_metfields_0125
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
    # use a space beforehand to avoid adding multiple #
    sed -i -e "s: 'CH4': #'CH4':g" \
        -e "s: 'ProdLoss: #'ProdLoss:g" \
        -e "s: 'StateMet: #'StateMet:g" \
        -e "s: 'SpeciesConcMND: #'SpeciesConcMND:g" \
        -e "s: 'Met_PEDGEDRY: #'Met_PEDGEDRY:g" \
        -e "s: 'Met_PFICU: #'Met_PFICU:g" \
        -e "s: 'Met_PFILSAN: #'Met_PFILSAN:g" \
        -e "s: 'Met_PFLCU: #'Met_PFLCU:g" \
        -e "s: 'Met_PFLLSAN: #'Met_PFLLSAN:g" HISTORY.rc

    # If turned on, save out hourly CH4 concentrations to daily files
    # use time-average mode
    if "$HourlySpecies"; then
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
        cp ${InversionPath}/src/geoschem_run_scripts/run.template .
    fi

    # Navigate back to working directory
    cd ${RunDirs}

    printf "\n=== DONE CREATING TEMPLATE RUN DIRECTORY ===\n"
}

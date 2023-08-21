#!/bin/bash

# Functions available in this file include:
#   - setup_template 

# Description: Setup template GCClassic run directory
# Usage:
#   setup_template
setup_template() {
    printf "\n=== CREATING TEMPLATE RUN DIRECTORY ===\n"

    cd ${GCClassicPath}/run

    # The createRunDir.sh script assumes the file ~/.geoschem/config exists
    # and contains the path to GEOS-Chem input data
	export GC_USER_REGISTERED=true
    if [[ ! -f ${HOME}/.geoschem/config ]]; then
	mkdir -p ${HOME}/.geoschem
	echo "export GC_DATA_ROOT=${DataPath}" >> ${HOME}/.geoschem/config
	source ${HOME}/.geoschem/config
    fi

    if [[ -d ${RunTemplate} ]]; then
	printf "\nERROR: ${RunTemplate} already exists. Please remove or set 'SetupTemplateRunDir: false' in config.yml.\n"
	exit 9999
    fi

    # Commands to feed to createRunDir.sh
    # Create a GEOS-FP 0.25x0.3125 nested NA CH4 run directory by default
    # (Grid and meteorology fields will be replaced below by the settings
    #  in config.yml)
    cmd="3\n2\n4\n4\n2\n${RunDirs}\n${runDir}\nn\n"

    # Create run directory
    printf ${cmd} | ./createRunDir.sh >> createRunDir.log 2>&1
    rm -f createRunDir.log
    printf "\nCreated ${RunTemplate}\n"

    cd ${RunTemplate}

    if "$isAWS"; then
	# Update GC data download to silence output from aws commands
	sed -i "s/command: 'aws s3 cp --request-payer=requester '/command: 'aws s3 cp --request-payer=requester --only-show-errors '/" download_data.yml
    fi

    # Modify geoschem_config.yml based on settings in config.yml
    sed -i -e "s:20190101:${StartDate}:g" \
           -e "s:20190201:${EndDate}:g" \
           -e "s:GEOSFP:${metUC}:g" \
           -e "s:0.25x0.3125:${gridResLong}:g" \
           -e "s:-130.0,  -60.0:${Lons}:g" \
           -e "s:9.75,  60.0:${Lats}:g" geoschem_config.yml

    # For CH4 inversions always turn analytical inversion on
    sed -i "/analytical_inversion/{N;s/activate: false/activate: true/}" geoschem_config.yml

    # Also turn on analytical inversion option in HEMCO_Config.rc
    OLD="--> AnalyticalInv          :       false"
    NEW="--> AnalyticalInv          :       true "
    sed -i "s/$OLD/$NEW/g" HEMCO_Config.rc

    # Modify path to state vector file in HEMCO_Config.rc
    OLD=" StateVector.nc"
    NEW=" ${RunDirs}/StateVector.nc"
    sed -i -e "s@$OLD@$NEW@g" HEMCO_Config.rc

    # Turn other options on/off according to settings above
    if "$UseEmisSF"; then
	OLD="use_emission_scale_factor: false"
	NEW="use_emission_scale_factor: true"
	sed -i "s/$OLD/$NEW/g" geoschem_config.yml
    fi
    if "$UseOHSF"; then
	OLD="use_OH_scale_factors: false"
	NEW="use_OH_scale_factors: true"
	sed -i "s/$OLD/$NEW/g" geoschem_config.yml
    fi

    # Modify HEMCO_Config.rc based on settings in config.yml
    # Use cropped met fields (add the region to both METDIR and the met files)
    if [ ! -z "$NestedRegion" ]; then
	sed -i -e "s:GEOS_0.25x0.3125\/GEOS_FP:GEOS_${native}_${NestedRegion}\/${metDir}:g" HEMCO_Config.rc
	sed -i -e "s:GEOS_0.25x0.3125\/GEOS_FP:GEOS_${native}_${NestedRegion}\/${metDir}:g" HEMCO_Config.rc.gmao_metfields
        sed -i -e "s:\$RES:\$RES.${NestedRegion}:g" HEMCO_Config.rc.gmao_metfields
    fi

    # Determine length of inversion period in days
    InvPeriodLength=$(( ( $(date -d ${EndDate} "+%s") - $(date -d ${StartDate} "+%s") ) / 86400))

    # If inversion period is < 32 days, use End diagnostic output frequency
    if (( ${InvPeriodLength} < 32 )); then
        sed -i -e "s|DiagnFreq:                   Monthly|DiagnFreq:                   End|g" HEMCO_Config.rc
    fi

    # Modify path to BC files
    sed -i -e "s:\$ROOT/SAMPLE_BCs/v2021-07/CH4:${fullBCpath}:g" HEMCO_Config.rc
    if [ "$NestedGrid" == "false" ]; then
        OLD="--> GC_BCs                 :       true "
        NEW="--> GC_BCs                 :       false"
        sed -i "s/$OLD/$NEW/g" HEMCO_Config.rc
    fi

    # Modify HISTORY.rc
    sed -i -e "s:'CH4':#'CH4':g" \
           -e "s:'Metrics:#'Metrics:g" HISTORY.rc
    
    # If turned on, save out hourly CH4 concentrations to daily files
    if "$HourlyCH4"; then
        sed -i -e 's/SpeciesConc.frequency:      00000100 000000/SpeciesConc.frequency:      00000000 010000/g' \
    	       -e 's/SpeciesConc.duration:       00000100 000000/SpeciesConc.duration:       00000001 000000/g' \
               -e 's/SpeciesConc.mode:           '\''time-averaged/SpeciesConc.mode:           '\''instantaneous/g' HISTORY.rc
    fi

    if ! "$isAWS"; then
	# Load environment with modules for compiling GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Remove sample restart file
    rm -f Restarts/GEOSChem.Restart.20190101_0000z.nc4

    # Copy template run script
    cp ${InversionPath}/src/geoschem_run_scripts/ch4_run.template .
    
    # Compile GEOS-Chem and store executable in template run directory
    printf "\nCompiling GEOS-Chem...\n"
    cd build
    cmake ${InversionPath}/GCClassic >> build_geoschem.log 2>&1
    cmake . -DRUNDIR=..  >> build_geoschem.log 2>&1 
    make -j install >> build_geoschem.log 2>&1
    cd ..
    if [[ -f gcclassic ]]; then
        rm -rf build
        mv build_info ../GEOSChem_build_info
    else
        printf "\nGEOS-Chem build failed! \n\nSee ${RunTemplate}/build/build_geoschem.log for details\n"
        exit 999
    fi
    printf "\nDone compiling GEOS-Chem \n\nSee ${RunDirs}/GEOSChem_build_info for details\n\n"
    
    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING TEMPLATE RUN DIRECTORY ===\n"
}

#!/bin/bash

# This script will set up an Integrated Methane Inversion (IMI) with GEOS-Chem.
# For documentation, see https://imi.readthedocs.io.
#
# Calling Sequence:
#  ./setup_imi.sh [CONFIG-FILE]
#
# Authors: Daniel Varon, Melissa Sulprizio, Lucas Estrada, Will Downs

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE (setup_imi.sh) ===\n"

# Check if user has specified a configuration file
if [[ $# == 1 ]] ; then
    ConfigFile=$1
else
    ConfigFile="config.yml"
fi

# Get configuration
source src/utilities/parse_yaml.sh
eval $(parse_yaml ${ConfigFile})

##=======================================================================
## Standard settings
##=======================================================================

# Path to inversion setup
InversionPath=$(pwd -P)

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
if "$BCdryrun"; then

    mkdir -p ${BCpath}

    if "$DoSpinup"; then
        START=${SpinupStart}
    else
        START=${StartDate}
    fi
    printf "\nDownloading boundary condition data for $START to $EndDate\n"
    python src/utilities/download_bc.py ${START} ${EndDate} ${BCpath}

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

    printf "\n=== CREATING RECTANGULAR STATE VECTOR FILE ===\n"
    
    # Use GEOS-FP or MERRA-2 CN file to determine ocean/land grid boxes
    LandCoverFile="${DataPath}/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.${NestedRegion}.${LandCoverFileExtension}"

    if "$isAWS"; then
	# Download land cover file
	s3_lc_path="s3://gcgrid/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.${NestedRegion}.${LandCoverFileExtension}"
	aws s3 cp --request-payer=requester ${s3_lc_path} ${LandCoverFile}
    fi

    # Output path and filename for state vector file
    StateVectorFName="StateVector.nc"

    # Create state vector file
    cd ${RunDirs}

    # Copy state vector creation script to working directory
    cp ${InversionPath}/src/utilities/make_state_vector_file.py .
    chmod 755 make_state_vector_file.py

    printf "\nCalling make_state_vector_file.py\n"
    python make_state_vector_file.py $LandCoverFile $StateVectorFName $LatMin $LatMax $LonMin $LonMax $BufferDeg $LandThreshold $nBufferClusters

    printf "\n=== DONE CREATING RECTANGULAR STATE VECTOR FILE ===\n"

else

    # Copy custom state vector to $RunDirs directory for later use
    cp $StateVectorFile ${RunDirs}/StateVector.nc

fi

if ! "$isAWS"; then
    # Load environment with NCO
    source ${NCOEnv}
fi

# Determine number of elements in state vector file
function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc)
rm ~/foo.nc
printf "\nNumber of state vector elements in this inversion = ${nElements}\n\n"

# Define inversion domain lat/lon bounds
function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
LonMinInvDomain=$(ncmin lon ${RunDirs}/StateVector.nc)
LonMaxInvDomain=$(ncmax lon ${RunDirs}/StateVector.nc)
LatMinInvDomain=$(ncmin lat ${RunDirs}/StateVector.nc)
LatMaxInvDomain=$(ncmax lat ${RunDirs}/StateVector.nc)
rm ~/foo.nc
Lons="${LonMinInvDomain}, ${LonMaxInvDomain}"
Lats="${LatMinInvDomain}, ${LatMaxInvDomain}"

##=======================================================================
## Set up template run directory
##=======================================================================

runDir="template_run"
RunTemplate="${RunDirs}/${runDir}"

if "$SetupTemplateRundir"; then

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
           -e "s:geosfp:${Met}:g" \
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

    # Modify HEMCO_Config.rc if running Kalman filter
    if "$KalmanMode"; then
    sed -i -e "s|use_emission_scale_factor: false|use_emission_scale_factor: true|g" geoschem_config.yml
    sed -i -e "s|--> Emis_ScaleFactor       :       false|--> Emis_ScaleFactor       :       true|g" \
           -e "s|gridded_posterior.nc|${RunDirs}/ScaleFactors.nc|g" HEMCO_Config.rc
    fi

    # Turn other options on/off according to settings above
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

    # Modify HEMCO_Config.rc based on settings in config.yml
    # Use cropped met fields (add the region to both METDIR and the met files)
    if [ ! -z "$NestedRegion" ]; then
	sed -i -e "s:GEOS_${native}:GEOS_${native}_${NestedRegion}:g" HEMCO_Config.rc
	sed -i -e "s:GEOS_${native}:GEOS_${native}_${NestedRegion}:g" HEMCO_Config.rc.gmao_metfields
        sed -i -e "s:\$RES:\$RES.${NestedRegion}:g" HEMCO_Config.rc.gmao_metfields
    fi

    # Determine length of inversion period in days
    InvPeriodLength=$(( ( $(date -d ${EndDate} "+%s") - $(date -d ${StartDate} "+%s") ) / 86400))

    # If inversion period is < 32 days, use End diagnostic output frequency
    if (( ${InvPeriodLength} < 32 )); then
        sed -i -e "s|DiagnFreq:                   Monthly|DiagnFreq:                   End|g" HEMCO_Config.rc
    fi

    # Modify path to BC files
    sed -i -e "s:\$ROOT/SAMPLE_BCs/v2021-07/CH4:${BCpath}:g" HEMCO_Config.rc
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

fi # SetupTemplateRunDir

##=======================================================================
##  Set up IMI preview run directory
##=======================================================================

preview_start=$(date +%s)
if  "$DoPreview"; then

    if ! "$isAWS"; then
	# Load environment with modules for running GEOS-Chem Classic
        source ${GEOSChemEnv}
    fi

    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml\n" 
        exit 9999
    fi

    printf "\n=== CREATING IMI PREVIEW RUN DIRECTORY ===\n"

    cd ${RunDirs}
    
    # Define the preview run name
    PreviewName="${RunName}_Preview"

    # Make the directory
    runDir="preview_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/*  ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFilePreview=${RestartFilePreviewPrefix}${StartDate}_0000z.nc4
    ln -s $RestartFilePreview Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
    if "$UseBCsForRestart"; then
        sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
    fi

    # End date for the preview simulation
    PreviewEnd=$(date --date="${StartDate} +1 day" +%Y%m%d)

    # Update settings in geoschem_config.yml
    sed -i -e "s|${EndDate}|${PreviewEnd}|g" geoschem_config.yml
    sed -i "/analytical_inversion/{N;s/activate: true/activate: false/}" geoschem_config.yml
    sed -i -e "s|use_emission_scale_factor: true|use_emission_scale_factor: false|g" geoschem_config.yml

    # Update settings in HEMCO_Config.rc
    sed -i -e "s|DiagnFreq:                   Monthly|DiagnFreq:                   End|g" HEMCO_Config.rc
    sed -i -e "s|--> Emis_ScaleFactor       :       true|--> Emis_ScaleFactor       :       false|g" HEMCO_Config.rc

    # Create run script from template
    sed -e "s:namename:${PreviewName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PreviewName}.run
    chmod 755 ${PreviewName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${PreviewName}.run
    fi

    ### Perform dry run if requested
    if "$PreviewDryRun"; then
        printf "\nExecuting dry-run for preview run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
    fi

    printf "\n=== DONE CREATING PREVIEW RUN DIRECTORY ===\n"

    ##===============##
    ##  Run preview  ##
    ##===============##

    printf "\n=== RUNNING IMI PREVIEW ===\n"

    # Submit preview GEOS-Chem job to job scheduler
    if "$UseSlurm"; then
        sbatch -W ${RunName}_Preview.run; wait;
    else
        ./${RunName}_Preview.run
    fi
    # Run preview script
    config_path=${InversionPath}/${ConfigFile}
    state_vector_path=${RunDirs}/StateVector.nc
    preview_dir=${RunDirs}/${runDir}
    tropomi_cache=${RunDirs}/data_TROPOMI
    preview_file=${InversionPath}/src/inversion_scripts/imi_preview.py

    # if running end to end script with sbatch then use
    # sbatch to take advantage of multiple cores 
    if "$UseSlurm"; then
        # set number of cores to run preview with
        if "$isAWS"; then
            sed -i -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${InversionPath}/src/inversion_scripts/imi_preview.py
        fi
        export PYTHONPATH=${PYTHONPATH}:${InversionPath}/src/inversion_scripts/
        chmod +x $preview_file
        sbatch -W $preview_file $InversionPath $config_path $state_vector_path $preview_dir $tropomi_cache; wait;
    else
        python $preview_file $InversionPath $config_path $state_vector_path $preview_dir $tropomi_cache
    fi
    printf "\n=== DONE RUNNING IMI PREVIEW ===\n"

    # Escape condition for DOFS threshold? Read diagnostics file for expectedDOFS variable
    eval $(parse_yaml ${preview_dir}/preview_diagnostics.txt) 
    if [ 1 -eq "$(echo "${expectedDOFS} < ${DOFSThreshold}" | bc)" ]; then  
        printf "\nExpected DOFS = ${expectedDOFS} are less than DOFSThreshold = ${DOFSThreshold}. Exiting.\n"
        printf "Consider increasing the inversion period, increasing the prior error, or using another prior inventory.\n"
        exit 0
    fi

    # Navigate back to top-level directory
    cd ..

fi 
preview_end=$(date +%s)

##=======================================================================
##  Set up spinup run directory
##=======================================================================

if  "$SetupSpinupRun"; then

    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml\n" 
        exit 9999
    fi

    printf "\n=== CREATING SPINUP RUN DIRECTORY ===\n"
    
    cd ${RunDirs}
    
    # Define the run directory name
    SpinupName="${RunName}_Spinup"

    # Make the directory
    runDir="spinup_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/*  ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFile=${RestartFilePrefix}${SpinupStart}_0000z.nc4
    ln -s $RestartFile Restarts/GEOSChem.Restart.${SpinupStart}_0000z.nc4
    if "$UseBCsForRestart"; then
        sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
	printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n" 
    fi
    
    # Update settings in geoschem_config.yml
    sed -i -e "s|${StartDate}|${SpinupStart}|g" \
           -e "s|${EndDate}|${SpinupEnd}|g" geoschem_config.yml
    sed -i "/analytical_inversion/{N;s/activate: true/activate: false/}" geoschem_config.yml
    sed -i -e "s|use_emission_scale_factor: true|use_emission_scale_factor: false|g" geoschem_config.yml

    # Update settings in HEMCO_Config.rc
    sed -i -e "s|--> Emis_ScaleFactor       :       true|--> Emis_ScaleFactor       :       false|g" HEMCO_Config.rc

    # Turn on LevelEdgeDiags output
    if "$HourlyCH4"; then
        sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
               -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
               -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
               -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
    fi

    # Create run script from template
    sed -e "s:namename:${SpinupName}:g" \
        -e "s:##:#:g" ch4_run.template > ${SpinupName}.run
    chmod 755 ${SpinupName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${SpinupName}.run
    fi

    ### Perform dry run if requested
    if "$SpinupDryrun"; then
        printf "\nExecuting dry-run for spinup run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
    fi
    
    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING SPINUP RUN DIRECTORY ===\n"

fi # SetupSpinupRun

##=======================================================================
##  Set up posterior run directory
##=======================================================================

if  "$SetupPosteriorRun"; then

    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml" 
        exit 9999
    fi

    printf "\n=== CREATING POSTERIOR RUN DIRECTORY ===\n"
    
    cd ${RunDirs}
    
    # Define the run directory name
    PosteriorName="${RunName}_Posterior"

    # Make the directory
    runDir="posterior_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp -r ${RunTemplate}/*  ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFileFromSpinup=../../spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        ln -s $RestartFileFromSpinup Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
    else
        RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
        ln -s $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
        if "$UseBCsForRestart"; then
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n" 
        fi
    fi
    
    # Update settings in geoschem_config.yml
    sed -i "/analytical_inversion/{N;s/activate: true/activate: false/}" geoschem_config.yml
    sed -i "s/use_emission_scale_factor: false/use_emission_scale_factor: true/g" geoschem_config.yml
    
    # Update settings in HEMCO_Config.rc
    sed -i -e "s|--> Emis_ScaleFactor       :       false|--> Emis_ScaleFactor       :       true|g" \
           -e "s|gridded_posterior.nc|${RunDirs}/inversion/gridded_posterior.nc|g" HEMCO_Config.rc

    # Turn on LevelEdgeDiags output
    if "$HourlyCH4"; then
        sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
               -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
               -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
               -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
    fi

    # Create run script from template
    sed -e "s:namename:${PosteriorName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PosteriorName}.run
    chmod 755 ${PosteriorName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${PosteriorName}.run
    fi

    ### Perform dry run if requested
    if "$PosteriorDryRun"; then
        printf "\nExecuting dry-run for posterior run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
    fi
    
    # Navigate back to top-level directory
    cd ..

    printf "\n=== DONE CREATING POSTERIOR RUN DIRECTORY ===\n"
    
fi # SetupPosteriorRun

##=======================================================================
##  Set up Jacobian run directories
##=======================================================================

if "$SetupJacobianRuns"; then

    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml" 
        exit 9999
    fi

    printf "\n=== CREATING JACOBIAN RUN DIRECTORIES ===\n"
    
    cd ${RunDirs}

    # Create directory that will contain all Jacobian run directories
    mkdir -p -v jacobian_runs

    # Copy run scripts
    cp ${InversionPath}/src/geoschem_run_scripts/run_jacobian_simulations.sh jacobian_runs/
    sed -i -e "s:{RunName}:${RunName}:g" \
           -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/run_jacobian_simulations.sh
    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" jacobian_runs/run_jacobian_simulations.sh
    fi
    cp ${InversionPath}/src/geoschem_run_scripts/submit_jacobian_simulations_array.sh jacobian_runs/
    sed -i -e "s:{START}:0:g" \
           -e "s:{END}:${nElements}:g" \
           -e "s:{InversionPath}:${InversionPath}:g" jacobian_runs/submit_jacobian_simulations_array.sh

    # Initialize (x=0 is base run, i.e. no perturbation; x=1 is state vector element=1; etc.)
    x=0

    # Create run directory for each state vector element so we can
    # apply the perturbation to each
    while [ $x -le $nElements ]; do

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

	# Define the run directory name
	name="${RunName}_${xstr}"

	# Make the directory
	runDir="./jacobian_runs/${name}"
	mkdir -p -v ${runDir}

	# Copy run directory files
	cp -r ${RunTemplate}/*  ${runDir}
	cd $runDir

	# Link to GEOS-Chem executable instead of having a copy in each rundir
	rm -rf gcclassic
	ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFileFromSpinup=../../../spinup_run/Restarts/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        ln -s $RestartFileFromSpinup Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
	else
	    RestartFile=${RestartFilePrefix}${StartDate}_0000z.nc4
	    ln -s $RestartFile Restarts/GEOSChem.Restart.${StartDate}_0000z.nc4
	    if "$UseBCsForRestart"; then
		sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            fi
	fi
   
	# Update settings in geoschem_config.yml
	sed -i -e "s|emission_perturbation: 1.0|emission_perturbation: ${PerturbValue}|g" \
	       -e "s|state_vector_element_number: 0|state_vector_element_number: ${xUSE}|g" geoschem_config.yml

	# Update settings in HISTORY.rc
	# Only save out hourly pressure fields to daily files for base run
	if [ $x -eq 0 ]; then
	    if "$HourlyCH4"; then
                sed -i -e 's/#'\''LevelEdgeDiags/'\''LevelEdgeDiags/g' \
                       -e 's/LevelEdgeDiags.frequency:   00000100 000000/LevelEdgeDiags.frequency:   00000000 010000/g' \
                       -e 's/LevelEdgeDiags.duration:    00000100 000000/LevelEdgeDiags.duration:    00000001 000000/g' \
                       -e 's/LevelEdgeDiags.mode:        '\''time-averaged/LevelEdgeDiags.mode:        '\''instantaneous/g' HISTORY.rc
	    fi
	fi

	# Create run script from template
	sed -e "s:namename:${name}:g" ch4_run.template > ${name}.run
	rm -f ch4_run.template
	chmod 755 ${name}.run

    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c 1:g" ${name}.run

        sed -i -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c 1:g" ../run_jacobian_simulations.sh
    fi

    ### Perform dry run if requested, only for base run
    if [ $x -eq 0 ]; then
        if "$ProductionDryRun"; then
            printf "\nExecuting dry-run for production runs...\n"
            ./gcclassic --dryrun &> log.dryrun
            ./download_data.py log.dryrun aws
        fi
    fi

	# Navigate back to top-level directory
	cd ../..

	# Increment
	x=$[$x+1]

    done

    printf "\n=== DONE CREATING JACOBIAN RUN DIRECTORIES ===\n"

fi  # SetupJacobianRuns

##=======================================================================
##  Setup inversion directory
##=======================================================================

if "$SetupInversion"; then

    printf "\n=== SETTING UP INVERSION DIRECTORY ===\n"
    
    cd ${OutputPath}/$RunName
    mkdir -p -v inversion
    mkdir -p inversion/data_converted
    mkdir -p inversion/data_geoschem
    mkdir -p inversion/data_sensitivities
    mkdir -p inversion/data_visualization
    mkdir -p inversion/operators
    
    cp ${InversionPath}/src/inversion_scripts/calc_sensi.py inversion/
    cp ${InversionPath}/src/inversion_scripts/invert.py inversion/
    cp ${InversionPath}/src/inversion_scripts/jacobian.py inversion/
    cp ${InversionPath}/src/inversion_scripts/operators/* inversion/operators/
    cp ${InversionPath}/src/inversion_scripts/make_gridded_posterior.py inversion/
    cp ${InversionPath}/src/inversion_scripts/postproc_diags.py inversion/
    cp ${InversionPath}/src/inversion_scripts/setup_gc_cache.py inversion/
    cp ${InversionPath}/src/inversion_scripts/utils.py inversion/
    cp ${InversionPath}/src/inversion_scripts/run_inversion.sh inversion/
    cp ${InversionPath}/src/notebooks/visualization_notebook.ipynb inversion/
    sed -i -e "s:{INVERSION_PATH}:${InversionPath}:g" \
           -e "s:{CONFIG_FILE}:${ConfigFile}:g" \
           -e "s:{STATE_VECTOR_ELEMENTS}:${nElements}:g" \
           -e "s:{OUTPUT_PATH}:${OutputPath}:g" \
           -e "s:{STATE_VECTOR_PATH}:../StateVector.nc:g" \
           -e "s:{LON_MIN}:${LonMinInvDomain}:g" \
           -e "s:{LON_MAX}:${LonMaxInvDomain}:g" \
           -e "s:{LAT_MIN}:${LatMinInvDomain}:g" \
           -e "s:{LAT_MAX}:${LatMaxInvDomain}:g" \
           -e "s:{RES}:${gridResLong}:g" inversion/run_inversion.sh

    if "$isAWS"; then
        sed -i -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -n 1:#SBATCH -n ${cpu_count}:g" inversion/run_inversion.sh
    fi
    
    printf "\n=== DONE SETTING UP INVERSION DIRECTORY ===\n"

fi #SetupInversion

# Remove temporary files
if "$isAWS"; then
    rm -f /home/ubuntu/foo.nc
fi

# Run time
echo "Statistics (setup):"
echo "Preview runtime (s): $(( $preview_end - $preview_start ))"
echo "Note: this is part of the Setup runtime reported by run_imi.sh"
exit 0

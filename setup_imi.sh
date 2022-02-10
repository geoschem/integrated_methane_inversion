#!/bin/bash

# This script will set up an Integrated Methane Inversion (IMI) with GEOS-Chem.
# For documentation, see https://integrated-methane-inversion.readthedocs.io.
#
# Authors: Daniel Varon, Melissa Sulprizio, Lucas Estrada, Will Downs

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE ===\n"

# Get configuration
source src/utilities/parse_yaml.sh
eval $(parse_yaml config.yml)
# For reference, this defines the following environment variables:
# General: $isAWS, $RunName, $UseSlurm
# Period of interest: $StartDate, $EndDate, $SpinupMonths
# Region of interest: $LonMin, $LonMax, $LatMin, $LatMax
# Inversion: $PriorError, $ObsError, $Gamma, $PrecomputedJacobian
# Grid: $Res, $Met, $HalfPolar, $Levs, $NestedGrid, $REGION, $Buffer
# Setup modules: $CreateStateVectorFile, $SetupTemplateRundir, $SetupSpinupRun, $SetupJacobianRuns, $SetupInversion, $SetupPosteriorRun
# Run modules: $RunSetup, $DoSpinup, $DoJacobian, $DoInversion, $DoPosterior
# State vector: $BufferDeg, $nBufferClusters, $LandThreshold
# If custom state vec file: $StateVectorFile, $LonMinCustomStateVector, $LonMaxCustomStateVector, $LatMinCustomStateVector, $LatMaxCustomStateVector
# Harvard-Cannon: $nCPUs, $partition

##=======================================================================
## Standard settings
##=======================================================================

# Path to inversion setup
InversionPath=$(pwd -P)

# AWS only: Download missing GEOS-Chem input data from S3
# You will be charged if your ec2 instance is not in the us-east-1 region.
if "$isAWS"; then
    PreviewDryRun=true     # Met fields/emissions for preview run
    SpinupDryrun=true      # Met fields/emissions for spinup run
    ProductionDryRun=true  # Met fields/emissions for production runs
    PosteriorDryRun=true   # Met fields/emissions for posterior run
    BCdryrun=true          # Boundary condition files
else
    PreviewDryRun=false
    SpinupDryrun=false
    ProductionDryRun=false
    PosteriorDryRun=false
    BCdryrun=false
fi

# Start and end date for the spinup simulation
SpinupStart=$(date --date="${StartDate} -${SpinupMonths} month" +%Y%m%d)
SpinupEnd=${StartDate}

# Path where you want to set up CH4 inversion code and run directories
if "$isAWS"; then
    MyPath="/home/ubuntu/imi_output_dir"
    CondaEnv="geo"
else
    MyPath="/n/holyscratch01/jacob_lab/$USER"

    # Environment files (specific to Harvard's Cannon cluster)
    NCOEnv="${InversionPath}/envs/Harvard-Cannon/gcc.ifort17_cannon.env"
    GCCEnv="${InversionPath}/envs/Harvard-Cannon/gcc.gfortran10.2_cannon.env"
    CondaEnv="ch4_inv" # See envs/README to create this environment
fi

# Path to find non-emissions input data
if "$isAWS"; then
    DataPath="/home/ubuntu/ExtData"
else
    DataPath="/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData"
fi

# Path to initial restart file
UseBCsForRestart=true
if "$isAWS"; then
    RestartDownload=true # automatically download restart file
    RestartFile="${DataPath}/BoundaryConditions/GEOSChem.BoundaryConditions.${SpinupStart}_0000z.nc4"
else
    RestartDownload=false
    RestartFile="/n/seasasfs02/CH4_inversion/InputData/BoundaryConditions/OutputDir_bias_corrected_dk_2/GEOSChem.BoundaryConditions.${SpinupStart}_0000z.nc4"
fi
    
# Path to boundary condition files (for nested grid simulations)
# Must put backslash before $ in $YYYY$MM$DD to properly work in sed command
if "$isAWS"; then
    BCfiles="${DataPath}/BoundaryConditions"
else
    BCfiles="/n/seasasfs02/CH4_inversion/InputData/BoundaryConditions/OutputDir_bias_corrected_dk_2/GEOSChem.BoundaryConditions.\$YYYY\$MM\$DD_0000z.nc4"
fi

# Activate Conda environment
printf "Activating conda environment: ${CondaEnv}\n"
if "$isAWS"; then
    source /home/ubuntu/miniconda/etc/profile.d/conda.sh
    conda activate $CondaEnv
else
    source activate $CondaEnv
fi

# Get max process count for spinup, production, and run_inversion scripts
if "$isAWS"; then
    output=$(echo $(slurmd -C))
    array=($output)
    cpu_str=$(echo ${array[1]})
    cpu_count=$(echo ${cpu_str:5})
    
    # with sbatch reduce cpu_count by 1 to account for parent sbatch process 
    # using 1 core 
    if "$UseSlurm"; then 
        cpu_count="$((cpu_count-1))"
    fi
fi

## Jacobian settings
PerturbValue=1.5

# Apply scale factors from a previous inversion?
UseEmisSF=false
UseSeparateWetlandSF=false
UseOHSF=false

# Turn on planeflight diagnostic in GEOS-Chem?
PLANEFLIGHT=false

# Save out hourly diagnostics from GEOS-Chem?
# For use in satellite operators via post-processing -- required for TROPOMI
# inversions
HourlyCH4=true

# Turn on old observation operators in GEOS-Chem?
# These will save out text files comparing GEOS-Chem to observations, but have
# to be manually incorporated into the IMI
GOSAT=false
TCCON=false
AIRS=false

##=======================================================================
## Download Boundary Conditions files if requested
##=======================================================================

if "$BCdryrun"; then

    mkdir -p ${BCfiles}

    if "$DoSpinup"; then
        START=${SpinupStart}
    else
        START=${StartDate}
    fi
    echo "Downloading boundary condition data for $START to $EndDate"
    python src/utilities/download_bc.py ${START} ${EndDate} ${BCfiles}

fi

##=======================================================================
## Download initial restart file if requested
##=======================================================================

if "$RestartDownload"; then
    if [ ! -f "$RestartFile" ]; then
        aws s3 cp --request-payer=requester s3://imi-boundary-conditions/GEOSChem.BoundaryConditions.${SpinupStart}_0000z.nc4 $RestartFile
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

if [ -z "$REGION" ]; then
    gridDir="$Res"
else
    gridDir="${Res}_${REGION}"
fi

# Define path to GEOS-Chem run directory files
GCClassicPath="${InversionPath}/GCClassic"
RunFilesPath="${GCClassicPath}/run"

# Create working directory if it doesn't exist yet
mkdir -p -v ${MyPath}/${RunName}

##=======================================================================
## Create state vector file
##=======================================================================

if "$CreateStateVectorFile"; then

    printf "\n=== CREATING RECTANGULAR STATE VECTOR FILE ===\n"
    
    # Use GEOS-FP or MERRA-2 CN file to determine ocean/land grid boxes
    LandCoverFile="${DataPath}/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.${REGION}.${LandCoverFileExtension}"

    # Download land cover file
    s3_lc_path="s3://gcgrid/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.${REGION}.${LandCoverFileExtension}"
    aws s3 cp --request-payer=requester ${s3_lc_path} ${LandCoverFile}

    # Output path and filename for state vector file
    StateVectorFName="StateVector.nc"

    # Create state vector file
    cd ${MyPath}/$RunName

    # Copy state vector creation script to working directory
    cp ${InversionPath}/src/inversion_scripts/make_state_vector_file.py .
    chmod 755 make_state_vector_file.py

    printf "Calling make_state_vector_file.py\n"
    python make_state_vector_file.py $LandCoverFile $StateVectorFName $LatMin $LatMax $LonMin $LonMax $BufferDeg $LandThreshold $nBufferClusters

    printf "=== DONE CREATING RECTANGULAR STATE VECTOR FILE ===\n"

else

    # Copy custom state vector to $MyPath/$RunName directory for later use
    cp $StateVectorFile ${MyPath}/${RunName}/StateVector.nc

fi

# Load environment with NCO
if ! "$isAWS"; then
    source ${NCOEnv}
fi

# Determine number of elements in state vector file
function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
nElements=$(ncmax StateVector ${MyPath}/${RunName}/StateVector.nc)
printf "\n Number of state vector elements in this inversion= ${nElements}\n"

# Define inversion domain lat/lon bounds
function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
LonMinInvDomain=$(ncmin lon ${MyPath}/${RunName}/StateVector.nc)
LonMaxInvDomain=$(ncmax lon ${MyPath}/${RunName}/StateVector.nc)
LatMinInvDomain=$(ncmin lat ${MyPath}/${RunName}/StateVector.nc)
LatMaxInvDomain=$(ncmax lat ${MyPath}/${RunName}/StateVector.nc)
Lons="${LonMinInvDomain} ${LonMaxInvDomain}"
Lats="${LatMinInvDomain} ${LatMaxInvDomain}"

# Purge software modules if not on AWS
if ! "$isAWS"; then
    module purge
fi

##=======================================================================
## Set up template run directory
##=======================================================================

RunTemplate="${MyPath}/${RunName}/template_run"

if "$SetupTemplateRundir"; then

    printf "\n=== CREATING TEMPLATE RUN DIRECTORY ===\n"

    cd ${MyPath}/${RunName}

    # Create template run directory
    mkdir -p -v ${RunTemplate}

    # Copy run directory files directly from GEOS-Chem repository
    cp -RLv ${RunFilesPath}/input.geos.templates/input.geos.CH4 ${RunTemplate}/input.geos
    cp -RLv ${RunFilesPath}/HISTORY.rc.templates/HISTORY.rc.CH4 ${RunTemplate}/HISTORY.rc
    cp -RLv ${RunFilesPath}/runScriptSamples/ch4_run.template ${RunTemplate}
    cp -RLv ${RunFilesPath}/getRunInfo ${RunTemplate}/
    cp -RLv ${RunFilesPath}/Makefile ${RunTemplate}/
    cp -RLv ${RunFilesPath}/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.CH4 ${RunTemplate}/HEMCO_Diagn.rc
    cp -RLv ${RunFilesPath}/HEMCO_Config.rc.templates/HEMCO_Config.rc.CH4 ${RunTemplate}/HEMCO_Config.rc
    cp -RLv ${RunFilesPath}/../shared/download_data.py ${RunTemplate}/
    cp -RLv ${RunFilesPath}/../shared/download_data.yml ${RunTemplate}/
    cp -RLv ${GCClassicPath}/src/GEOS-Chem/run/shared/species_database.yml ${RunTemplate}/

    cd $RunTemplate

    # Update settings in input.geos
    sed -i -e "s:{DATE1}:${StartDate}:g" \
           -e "s:{DATE2}:${EndDate}:g" \
           -e "s:{TIME1}:000000:g" \
           -e "s:{TIME2}:000000:g" \
           -e "s:{MET}:${Met}:g" \
           -e "s:{DATA_ROOT}:${DataPath}:g" \
           -e "s:{SIM}:CH4:g" \
           -e "s:{RES}:${gridResLong}:g" \
           -e "s:{LON_RANGE}:${Lons}:g" \
           -e "s:{LAT_RANGE}:${Lats}:g" \
           -e "s:{CENTER_180}:T:g" \
           -e "s:{HALF_POLAR}:${HalfPolar}:g" \
           -e "s:{NLEV}:${Levs}:g" \
           -e "s:{NESTED_SIM}:${NestedGrid}:g" \
           -e "s:{BUFFER_ZONE}:${Buffer}:g" input.geos
    if [ "$NestedGrid" == "T" ]; then
	sed -i -e "s|timestep \[sec\]: 600|timestep \[sec\]: 300|g" \
           -e "s|timestep \[sec\]: 1200|timestep \[sec\]: 600|g" input.geos
    fi

    # For CH4 inversions always turn analytical inversion on
    OLD="Do analytical inversion?: F"
    NEW="Do analytical inversion?: T"
    sed -i "s/$OLD/$NEW/g" input.geos

    # Turn on analytical inversion option in HEMCO_Config.rc also
    OLD="--> AnalyticalInv          :       false"
    NEW="--> AnalyticalInv          :       true "
    sed -i "s/$OLD/$NEW/g" HEMCO_Config.rc

    # Modify path to state vector file in HEMCO_Config.rc
    OLD=" StateVector.nc"
    NEW=" ${MyPath}/${RunName}/StateVector.nc"
    echo $NEW
    sed -i -e "s@$OLD@$NEW@g" HEMCO_Config.rc

    # Turn other options on/off according to settings above
    if "$GOSAT"; then
	OLD="Use GOSAT obs operator? : F"
	NEW="Use GOSAT obs operator? : T"
	sed -i "s/$OLD/$NEW/g" input.geos
    fi
    if "$TCCON"; then
	OLD="Use TCCON obs operator? : F"
	NEW="Use TCCON obs operator? : T"
	sed -i "s/$OLD/$NEW/g" input.geos
    fi
    if "$AIRS"; then
	OLD="Use AIRS obs operator?  : F"
	NEW="Use AIRS obs operator?  : T"
	sed -i "s/$OLD/$NEW/g" input.geos
    fi
    if "$UseEmisSF"; then
	OLD=" => Use emis scale factr: F"
	NEW=" => Use emis scale factr: T"
	sed -i "s/$OLD/$NEW/g" input.geos
    fi
    if "$UseSeparateWetlandSF"; then
	OLD=" => Use sep. wetland SFs: F"
	NEW=" => Use sep. wetland SFs: T"
	sed -i "s/$OLD/$NEW/g" input.geos
    fi
    if "$UseOHSF"; then
	OLD=" => Use OH scale factors: F"
	NEW=" => Use OH scale factors: T"
	sed -i "s/$OLD/$NEW/g" input.geos
    fi
    if "$PLANEFLIGHT"; then
	mkdir -p Plane_Logs
	OLD="Turn on plane flt diag? : F"
	NEW="Turn on plane flt diag? : T"
	sed -i "s/$OLD/$NEW/g" input.geos
	OLD="Flight track info file  : Planeflight.dat.YYYYMMDD"
	NEW="Flight track info file  : Planeflights\/Planeflight.dat.YYYYMMDD"
	sed -i "s/$OLD/$NEW/g" input.geos
	OLD="Output file name        : plane.log.YYYYMMDD"
	NEW="Output file name        : Plane_Logs\/plane.log.YYYYMMDD"
	sed -i "s/$OLD/$NEW/g" input.geos
    fi

    # Set up HEMCO_Config.rc
    sed -i -e "/### BEGIN SECTION SETTINGS/r ${GCClassicPath}/run/HEMCO_Config.rc.templates/header.gmao"                     HEMCO_Config.rc
    sed -i -e "/# --- Meteorology fields for FlexGrid ---/r ${GCClassicPath}/run/HEMCO_Config.rc.templates/met_fields.gmao"  HEMCO_Config.rc

    # Determine length of inversion period in days
    InvPeriodLength=$(( ( $(date -d ${EndDate} "+%s") - $(date -d ${StartDate} "+%s") ) / 86400))

    # If inversion period is < 32 days, use End diagnostic output frequency
    if (( ${InvPeriodLength} < 32 )); then
        sed -i -e "s|DiagnFreq:                   Monthly|DiagnFreq:                   End|g" HEMCO_Config.rc
    fi
    # Otherwise the Monthly diagnostic output frequency is used
    sed -i -e "s:{VERBOSE}:0:g" \
           -e "s:{WARNINGS}:1:g" \
           -e "s:{DATA_ROOT}:${DataPath}:g" \
           -e "s:{GRID_DIR}:${gridDir}:g" \
           -e "s:{MET_DIR}:${metDir}:g" \
           -e "s:{NATIVE_RES}:${native}:g" \
           -e "s:\$ROOT/SAMPLE_BCs/v2019-05/CH4/GEOSChem.BoundaryConditions.\$YYYY\$MM\$DD_\$HH\$MNz.nc4:${BCfiles}:g" HEMCO_Config.rc

    if [ ! -z "$REGION" ]; then
        sed -i -e "s:\$RES:\$RES.${REGION}:g" HEMCO_Config.rc
    fi
    if [ "$NestedGrid" == "T" ]; then
        OLD="--> GC_BCs                 :       false "
        NEW="--> GC_BCs                 :       true  "
        sed -i "s/$OLD/$NEW/g" HEMCO_Config.rc
    fi

    # Set up HISTORY.rc
    # Use monthly output by default
    sed -i -e "s:{FREQUENCY}:00000100 000000:g" \
           -e "s:{DURATION}:00000100 000000:g" \
           -e "s:'CH4':#'CH4':g" \
           -e "s:'Metrics:#'Metrics:g" HISTORY.rc
    
    # If turned on, save out hourly CH4 concentrations to daily files
    if "$HourlyCH4"; then
        sed -i -e 's/SpeciesConc.frequency:      00000100 000000/SpeciesConc.frequency:      00000000 010000/g' \
    	       -e 's/SpeciesConc.duration:       00000100 000000/SpeciesConc.duration:       00000001 000000/g' \
               -e 's/SpeciesConc.mode:           '\''time-averaged/SpeciesConc.mode:           '\''instantaneous/g' HISTORY.rc
    fi

    # Load environment with modules for compiling GEOS-Chem Classic
    if ! "$isAWS"; then
        source ${GCCEnv}
    fi
    
    # Compile GEOS-Chem and store executable in template run directory
    mkdir build; cd build
    cmake ${InversionPath}/GCClassic
    cmake . -DRUNDIR=..
    make -j install
    cd ..
    if [[ -f gcclassic ]]; then
        rm -rf build
        mv build_info ../GEOSChem_build_info
    else
        echo "GEOS-Chem build failed"
        exit 999
    fi

    # Purge software modules
    if ! "$isAWS"; then
        module purge
    fi
    
    # Navigate back to top-level directory
    cd ..

    printf "=== DONE CREATING TEMPLATE RUN DIRECTORY ===\n"

fi # SetupTemplateRunDir

##=======================================================================
##  Set up IMI preview run directory
##=======================================================================

if  "$DoPreview"; then

    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/input.geos ]]; then
        echo "Template run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml" 
        exit 9999
    fi

    printf "\n=== CREATING IMI PREVIEW RUN DIRECTORY ===\n"

    cd ${MyPath}/${RunName}
    
    # Define the preview run name
    PreviewName="${RunName}_Preview"

    # Make the directory
    runDir="preview_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp ${RunTemplate}/*  ${runDir}
    cd $runDir
    mkdir -p OutputDir
    mkdir -p Restarts

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    ln -s $RestartFile GEOSChem.Restart.${SpinupStart}_0000z.nc4
    if "$UseBCsForRestart"; then
        sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
    fi

    # End date for the preview simulation
    PreviewEnd=$(date --date="${SpinupStart} +1 day" +%Y%m%d)

    # Update settings in input.geos
    sed -i -e "s|${StartDate}|${SpinupStart}|g" \
           -e "s|${EndDate}|${PreviewEnd}|g" \
           -e "s|Do analytical inversion?: T|Do analytical inversion?: F|g" \
           -e "s|{PERTURBATION}|1.0|g" \
           -e "s|{ELEMENT}|0|g" input.geos

    # Update settings in HEMCO_Config.rc
    sed -i -e "s|DiagnFreq:                   Monthly|DiagnFreq:                   End|g" HEMCO_Config.rc

    # Create run script from template
    sed -e "s:namename:${PreviewName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PreviewName}.run
    chmod 755 ${PreviewName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
        sed -i -e "/#SBATCH -p huce_cascade/d" \
               -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${PreviewName}.run
    fi

    ### Perform dry run if requested
    if "$PreviewDryRun"; then
        printf "Executing dry-run for preview run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
    fi

    printf "=== DONE CREATING PREVIEW RUN DIRECTORY ===\n"

    ##===============##
    ##  Run preview  ##
    ##===============##

    printf "\n=== RUNNING IMI PREVIEW ===\n"

    # Submit preview GEOS-Chem job to job scheduler
    sbatch -W ${RunName}_Preview.run; wait;

    # Run preview script
    config_path=${InversionPath}/config.yml
    state_vector_path=${MyPath}/${RunName}/StateVector.nc
    preview_dir=${MyPath}/${RunName}/${runDir}
    tropomi_cache=${MyPath}/${RunName}/data_TROPOMI
    preview_file=${InversionPath}/src/inversion_scripts/imi_preview.py

    # if running end to end script with sbatch then use
    # sbatch to take advantage of multiple cores 
    if "$UseSlurm"; then
        # set number of cores to run preview with
        if "$isAWS"; then
            sed -i -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${InversionPath}/src/inversion_scripts/imi_preview.py
        fi
        cd ${InversionPath}/src/inversion_scripts/
        chmod +x $preview_file
        sbatch -W $preview_file $config_path $state_vector_path $preview_dir $tropomi_cache $cpu_count; wait;
        cd $runDir
    else
        python $preview_file $config_path $state_vector_path $preview_dir $tropomi_cache $cpu_count
    fi
    printf "=== DONE RUNNING IMI PREVIEW ===\n"

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

##=======================================================================
##  Set up spinup run directory
##=======================================================================

if  "$SetupSpinupRun"; then

    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/input.geos ]]; then
        echo "Template run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml" 
        exit 9999
    fi

    printf "\n=== CREATING SPINUP RUN DIRECTORY ===\n"
    
    cd ${MyPath}/${RunName}
    
    # Define the run directory name
    SpinupName="${RunName}_Spinup"

    # Make the directory
    runDir="spinup_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp ${RunTemplate}/*  ${runDir}
    cd $runDir
    mkdir -p OutputDir
    mkdir -p Restarts

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    ln -s $RestartFile GEOSChem.Restart.${SpinupStart}_0000z.nc4
    if "$UseBCsForRestart"; then
        sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
    fi
    
    # Update settings in input.geos
    sed -i -e "s|${StartDate}|${SpinupStart}|g" \
           -e "s|${EndDate}|${SpinupEnd}|g" \
           -e "s|Do analytical inversion?: T|Do analytical inversion?: F|g" \
           -e "s|{PERTURBATION}|1.0|g" \
           -e "s|{ELEMENT}|0|g" input.geos

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
        sed -i -e "/#SBATCH -p huce_cascade/d" \
               -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${SpinupName}.run
    fi

    ### Perform dry run if requested
    if "$SpinupDryrun"; then
        printf "Executing dry-run for spinup run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
    fi
    
    # Navigate back to top-level directory
    cd ..

    printf "=== DONE CREATING SPINUP RUN DIRECTORY ===\n"

fi # SetupSpinupRun

##=======================================================================
##  Set up posterior run directory
##=======================================================================

if  "$SetupPosteriorRun"; then

    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/input.geos ]]; then
        echo "Template run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml" 
        exit 9999
    fi

    printf "\n=== CREATING POSTERIOR RUN DIRECTORY ===\n"
    
    cd ${MyPath}/${RunName}
    
    # Define the run directory name
    PosteriorName="${RunName}_Posterior"

    # Make the directory
    runDir="posterior_run"
    mkdir -p -v ${runDir}

    # Copy run directory files
    cp ${RunTemplate}/*  ${runDir}
    cd $runDir
    mkdir -p OutputDir
    mkdir -p Restarts

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFileFromSpinup=../spinup_run/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        ln -s $RestartFileFromSpinup GEOSChem.Restart.${StartDate}_0000z.nc4
    else
        if "$UseBCsForRestart"; then
            RestartFile=${DataPath}/BoundaryConditions/GEOSChem.BoundaryConditions.${StartDate}_0000z.nc4
            ln -s $RestartFile GEOSChem.Restart.${StartDate}_0000z.nc4
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n" 
        fi
    fi
    
    # Update settings in input.geos
    sed -i -e "s|Do analytical inversion?: T|Do analytical inversion?: F|g" \
           -e "s|Use emis scale factr    : F|Use emis scale factr    : T|g" \
           -e "s|{PERTURBATION}|1.0|g" \
           -e "s|{ELEMENT}|0|g" input.geos

    # Update settings in HEMCO_Config.rc
    sed -i -e "s|--> Emis_ScaleFactor       :       false|--> Emis_ScaleFactor       :       true|g" \
           -e "s|gridded_posterior.nc|${MyPath}/${RunName}/inversion/gridded_posterior.nc|g" HEMCO_Config.rc

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
        sed -i -e "/#SBATCH -p huce_cascade/d" \
               -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${PosteriorName}.run
    fi

    ### Perform dry run if requested
    if "$PosteriorDryRun"; then
        printf "Executing dry-run for posterior run...\n"
        ./gcclassic --dryrun &> log.dryrun
        ./download_data.py log.dryrun aws
    fi
    
    # Navigate back to top-level directory
    cd ..

    printf "=== DONE CREATING POSTERIOR RUN DIRECTORY ===\n"
    
fi # SetupPosteriorRun

##=======================================================================
##  Set up Jacobian run directories
##=======================================================================

if "$SetupJacobianRuns"; then

    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/input.geos ]]; then
        echo "Template run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml" 
        exit 9999
    fi

    printf "\n=== CREATING JACOBIAN RUN DIRECTORIES ===\n"
    
    cd ${MyPath}/${RunName}

    # Create directory that will contain all Jacobian run directories
    mkdir -p -v jacobian_runs

    # Copy run scripts
    cp ${RunFilesPath}/runScriptSamples/run_jacobian_simulations.sh jacobian_runs/
    sed -i -e "s:{RunName}:${RunName}:g" jacobian_runs/run_jacobian_simulations.sh
    if "$isAWS"; then
        sed -i -e "/#SBATCH -p huce_cascade/d" \
               -e "/#SBATCH -t/d" jacobian_runs/run_jacobian_simulations.sh
    fi
    cp ${RunFilesPath}/runScriptSamples/submit_jacobian_simulations_array.sh jacobian_runs/
    sed -i -e "s:{START}:0:g" -e "s:{END}:${nElements}:g" jacobian_runs/submit_jacobian_simulations_array.sh

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
	cp ${RunTemplate}/*  ${runDir}
	cd $runDir
	mkdir -p OutputDir
	mkdir -p Restarts

	# Link to GEOS-Chem executable instead of having a copy in each rundir
	rm -rf gcclassic
	ln -s ${RunTemplate}/gcclassic .

    # Link to restart file
    RestartFileFromSpinup=../../spinup_run/GEOSChem.Restart.${SpinupEnd}_0000z.nc4
    if test -f "$RestartFileFromSpinup" || "$DoSpinup"; then
        ln -s $RestartFileFromSpinup GEOSChem.Restart.${StartDate}_0000z.nc4
	else
	    if "$UseBCsForRestart"; then
            RestartFile=${DataPath}/BoundaryConditions/GEOSChem.BoundaryConditions.${StartDate}_0000z.nc4
            ln -s $RestartFile GEOSChem.Restart.${StartDate}_0000z.nc4
            sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
            printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n" 
        fi
	fi
   
	# Update settings in input.geos
	sed -i -e "s:{PERTURBATION}:${PerturbValue}:g" \
	       -e "s:{ELEMENT}:${xUSE}:g" input.geos

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
        sed -i -e "/#SBATCH -p huce_cascade/d" \
               -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c 1:g" ${name}.run

        sed -i -e "/#SBATCH -p huce_cascade/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c 1:g" ../run_jacobian_simulations.sh
    fi

    ### Perform dry run if requested, only for base run
    if [ $x -eq 0 ]; then
        if "$ProductionDryRun"; then
            printf "Executing dry-run for production runs...\n"
            ./gcclassic --dryrun &> log.dryrun
            ./download_data.py log.dryrun aws
        fi
    fi

	# Navigate back to top-level directory
	cd ../..

	# Increment
	x=$[$x+1]

    done

    printf "=== DONE CREATING JACOBIAN RUN DIRECTORIES ===\n"

fi  # SetupJacobianRuns

##=======================================================================
##  Setup inversion directory
##=======================================================================

if "$SetupInversion"; then

    printf "\n=== SETTING UP INVERSION DIRECTORY ===\n"
    
    cd ${MyPath}/$RunName
    mkdir -p -v inversion
    mkdir -p inversion/data_converted
    mkdir -p inversion/data_GC
    mkdir -p inversion/Sensi
    if "$isAWS"; then
        cp -rfP /home/ubuntu/backup_files/input_data/ ${MyPath}/
    else
        mkdir -p inversion/data_TROPOMI
        ln -s /n/holylfs05/LABS/jacob_lab/lshen/CH4/TROPOMI/data inversion/data_TROPOMI
    fi
    cp ${InversionPath}/src/inversion_scripts/calc_sensi.py inversion/
    cp ${InversionPath}/src/inversion_scripts/invert.py inversion/
    cp ${InversionPath}/src/inversion_scripts/jacobian.py inversion/
    cp ${InversionPath}/src/inversion_scripts/make_gridded_posterior.py inversion/
    cp ${InversionPath}/src/inversion_scripts/postproc_diags.py inversion/
    cp ${InversionPath}/src/inversion_scripts/setup_gc_cache.py inversion/
    cp ${InversionPath}/src/inversion_scripts/utils.py inversion/
    cp ${InversionPath}/src/inversion_scripts/run_inversion.sh inversion/
    cp ${InversionPath}/src/notebooks/visualization_notebook.ipynb inversion/
    sed -i -e "s:{STATE_VECTOR_ELEMENTS}:${nElements}:g" \
           -e "s:{MY_PATH}:${MyPath}:g" \
           -e "s:{STATE_VECTOR_PATH}:../StateVector.nc:g" \
           -e "s:{LON_MIN}:${LonMinInvDomain}:g" \
           -e "s:{LON_MAX}:${LonMaxInvDomain}:g" \
           -e "s:{LAT_MIN}:${LatMinInvDomain}:g" \
           -e "s:{LAT_MAX}:${LatMaxInvDomain}:g" \
           -e "s:{RES}:${gridResLong}:g" inversion/run_inversion.sh

    if "$isAWS"; then
        sed -i -e "/#SBATCH -p huce_intel/d" \
               -e "/#SBATCH -t/d" \
               -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -n 1:#SBATCH -n ${cpu_count}:g" inversion/run_inversion.sh
    fi
    
    printf "=== DONE SETTING UP INVERSION DIRECTORY ===\n"

fi #SetupInversion

# Remove temporary files
if "$isAWS"; then
    rm -f /home/ubuntu/foo.nc
fi

exit 0

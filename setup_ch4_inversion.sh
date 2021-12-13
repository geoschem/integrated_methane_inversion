#!/bin/bash

# This script will set up CH4 analytical inversions with GEOS-Chem. See
# setup_ch4_inversion_instructions.txt for details (mps, 2/20/2020)

##=======================================================================
## Parse config.yml file
##=======================================================================

printf "\n=== PARSING CONFIG FILE ===\n"

# Function to parse yaml files from shell script
# By Stefan Farestam via stackoverflow:
# https://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# Get configuration
eval $(parse_yaml config.yml)
# For reference, this defines the following environment variables:
# General: $isAWS, $RunName
# Period of interest: $StartDate, $EndDate, $SpinupMonths
# Region of interest: $LonMin, $LonMax, $LatMin, $LatMax
# Inversion: $PriorError, $ObsError, $Gamma
# Grid: $Res, $Met, $HalfPolar, $Levs, $NestedGrid, $REGION, $Buffer
# Setup modules: $CreateStateVectorFile, $SetupTemplateRundir, $SetupSpinupRun, $SetupJacobianRuns, $SetupInversion, $SetupPosteriorRun
# Run modules: $RunSetup, $DoSpinup, $DoJacobian, $DoInversion, $DoPosterior
# State vector: $BufferDeg, $nBufferClusters, $LandThreshold, $StateVectorFile
# Harvard-Cannon: $nCPUs, $partition

##=======================================================================
## Standard settings
##=======================================================================

# Path to inversion setup
InversionPath=$(pwd -P)

# AWS only: Download missing GEOS-Chem input data from S3 (you will be charged)
if "$isAWS"; then
    SpinupDryrun=true      # Met fields/emissions for spinup run
    ProductionDryRun=true  # Met fields/emissions for production runs
    PosteriorDryRun=true   # Met fields/emissions for posterior run
    BCdryrun=true          # Boundary condition files
else
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
    MyPath="/home/ubuntu/CH4_Workflow"
    CondaEnv="geo"
else
    MyPath="/n/holyscratch01/jacob_lab/msulprizio/CH4"

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
    python download_bc.py ${START} ${EndDate} ${BCfiles}

fi

##=======================================================================
## Download initial restart file if requested
##=======================================================================
if "$RestartDownload"; then
    if [ ! -f "$RestartFile" ]; then
	aws s3 cp --request-payer=requester s3://umi-bc-test/${RestartFile} $RestartFile
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
elif [ "$Met" == "merra2" ]; then
    metUC="MERRA2"
    metDir="MERRA2"
    native="0.5x0.625"
    constYr="2015"
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
mkdir -p -v ${MyPath}/$RunName

##=======================================================================
## Create state vector file
##=======================================================================
if "$CreateStateVectorFile"; then

    printf "\n=== CREATING STATE VECTOR FILE ===\n"
    
    # Use GEOS-FP or MERRA-2 CN file to determine ocean/land grid boxes
    LandCoverFile="${DataPath}/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.${REGION}.nc"

    # Output path and filename for state vector file
    StateVectorFile="StateVector.nc"

    # Create state vector file
    cd ${MyPath}/$RunName
    mkdir -p -v StateVectorFile # djv: why is StateVectorFile a directory here, it's .nc ?
    cd StateVectorFile

    # Copy state vector creation script to working directory
    cp ${InversionPath}/PostprocessingScripts/CH4_TROPOMI_INV/make_state_vector_file.py .
    chmod 755 make_state_vector_file.py

    # Activate Conda environment
    printf "Activating conda environment: ${CondaEnv}\n"
    if "$isAWS"; then
        source /home/ubuntu/miniconda/etc/profile.d/conda.sh
        conda activate $CondaEnv
    else
        source activate $CondaEnv
    fi

    printf "Calling make_state_vector_file.py\n"
    python make_state_vector_file.py $LandCoverFile $StateVectorFile $LatMin $LatMax $LonMin $LonMax $BufferDeg $LandThreshold $nBufferClusters
    conda deactivate
    
    # Define inversion domain lat/lon bounds
    Lons="$(( LonMin-BufferDeg )) $(( LonMax+BufferDeg ))"
    Lats="$(( LatMin-BufferDeg )) $(( LatMax+BufferDeg ))"

    printf "=== DONE CREATING STATE VECTOR FILE ===\n"

else
    echo Need something to define Lons, Lats from the custom state vector file!
fi

# Load environment with NCO
if ! "$isAWS"; then
    source ${NCOEnv}
fi

# Determine number of elements in state vector file
function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
nElements=$(ncmax StateVector $StateVectorFile)
printf "\n Number of state vector elements in this inversion= ${nElements}\n"

# Purge software modules
if ! "$isAWS"; then
    module purge
fi

##=======================================================================
## Set up template run directory
##=======================================================================
if "$SetupTemplateRundir"; then

    printf "\n=== CREATING TEMPLATE RUN DIRECTORY ===\n"

    cd ${MyPath}/${RunName}

    # Create template run directory
    RunTemplate="template_run"
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
    cp -RLv ${GCClassicPath}/src/GEOS-Chem/run/shared/species_database.yml ${RunTemplate}/

    cd $RunTemplate
    mkdir -p OutputDir
    mkdir -p Restarts

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
    OLD=" StateVectors.nc"
    if "$CreateStateVectorFile"; then
	NEW=" ${MyPath}/${RunName}/StateVectorFile/${StateVectorFile}"
    else
	NEW=" ${StateVectorFile}"
    fi
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
    # Use monthly emissions diagnostic output by default
    sed -i -e "s:End:Monthly:g" \
           -e "s:{VERBOSE}:0:g" \
           -e "s:{WARNINGS}:1:g" \
           -e "s:{DATA_ROOT}:${DataPath}:g" \
           -e "s:{GRID_DIR}:${gridDir}:g" \
           -e "s:{MET_DIR}:${metDir}:g" \
           -e "s:{NATIVE_RES}:${native}:g" \
           -e "s:\$ROOT/SAMPLE_BCs/v2019-05/CH4/GEOSChem.BoundaryConditions.\$YYYY\$MM\$DD_\$HH\$MNz.nc4:${BCfiles}:g" HEMCO_Config.rc

    if [ ! -z "$REGION" ]; then
        sed -i -e "s:\$Res:\$Res.${REGION}:g" HEMCO_Config.rc
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
    rm -rf build

    # Purge software modules
    if ! "$isAWS"; then
        module purge
    fi
    
    # Navigate back to top-level directory
    cd ..

    printf "=== DONE CREATING TEMPLATE RUN DIRECTORY ===\n"

fi # SetupTemplateRunDir

##=======================================================================
##  Set up spinup run directory
##=======================================================================

if "$isAWS"; then
    # Get max process count for spinup, production, and run_inversion scripts
    output=$(echo $(slurmd -C))
    array=($output)
    cpu_str=$(echo ${array[1]})
    cpu_count=$(echo ${cpu_str:5})
fi

if  "$SetupSpinupRun"; then

    # Make sure template run directory exists
    if [ ! -d ${RunTemplate} ]; then
	echo "Template run directory does not exist. Please set 'SetupTemplateRundir=true' in setup_ch4_inversion.sh" 
	exit 9999
    fi

    printf "\n=== CREATING SPINUP RUN DIRECTORY ===\n"
    
    cd ${MyPath}/${RunName}
    
    # Define the run directory name
    SpinupName="${RunName}_Spinup"

    # Make the directory
    runDir="spinup_run"
    mkdir -p -v ${runDir}

    # Copy and point to the necessary data
    cp -r ${RunTemplate}/*  ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ../${RunTemplate}/gcclassic .

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

    # Create run script from template
    sed -e "s:namename:${SpinupName}:g" \
	-e "s:##:#:g" ch4_run.template > ${SpinupName}.run
    chmod 755 ${SpinupName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
	sed -i -e "/#SBATCH -p huce_intel/d" \
	       -e "/#SBATCH -t/d" \
	       -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${SpinupName}.run
    fi

    ### Perform dry run if requested
    if "$SpinupDryrun"; then
       printf "Executing dry-run for spinup run...\n"
       ./gcclassic --dryrun &> log.dryrun
       ./download_data.py --aws log.dryrun
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
    if [ ! -d ${RunTemplate} ]; then
	echo "Template run directory does not exist. Please set 'SetupTemplateRundir=true' in setup_ch4_inversion.sh" 
	exit 9999
    fi

    printf "\n=== CREATING POSTERIOR RUN DIRECTORY ===\n"
    
    cd ${MyPath}/${RunName}
    
    # Define the run directory name
    PosteriorName="${RunName}_Posterior"

    # Make the directory
    runDir="posterior_run"
    mkdir -p -v ${runDir}

    # Copy and point to the necessary data
    cp -r ${RunTemplate}/*  ${runDir}
    cd $runDir

    # Link to GEOS-Chem executable instead of having a copy in each run dir
    rm -rf gcclassic
    ln -s ../${RunTemplate}/gcclassic .

    # Link to restart file
    if "$DoSpinup"; then
       ln -s ../spinup_run/GEOSChem.Restart.${SpinupEnd}_0000z.nc4 GEOSChem.Restart.${StartDate}_0000z.nc4
    else
       ln -s $RestartFile GEOSChem.Restart.${StartDate}_0000z.nc4
       if "$UseBCsForRestart"; then
	   sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
	   printf "\nWARNING: Changing restart field entry in HEMCO_Config.rc to read the field from a boundary condition file. Please revert SpeciesBC_ back to SpeciesRst_ for subsequent runs.\n" 
       fi
    fi
    
    # Update settings in input.geos
    sed -i -e "s|Do analytical inversion?: T|Do analytical inversion?: F|g" \
           -e "s|{PERTURBATION}|1.0|g" \
           -e "s|{ELEMENT}|0|g" input.geos

    # Create run script from template
    sed -e "s:namename:${SpinupName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PosteriorName}.run
    chmod 755 ${PosteriorName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
	sed -i -e "/#SBATCH -p huce_intel/d" \
	       -e "/#SBATCH -t/d" \
	       -e "/#SBATCH --mem/d" \
	       -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${PosteriorName}.run
    fi
    
    # Print messages
    printf "\nNote: You will need to manually modify HEMCO_Config.rc to apply the appropriate scale factors.\n"

    ### Perform dry run if requested
    if "$PosteriorDryRun"; then
	   printf "Executing dry-run for posterior run...\n"
	   ./gcclassic --dryrun &> log.dryrun
	   ./download_data.py --aws log.dryrun
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
    if [ ! -d ${RunTemplate} ]; then
	echo "Template run directory does not exist. Please set 'SetupTemplateRundir=true' in setup_ch4_inversion.sh" 
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
	sed -i -e "/#SBATCH -p huce_intel/d" \
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

	# Copy and point to the necessary data
	cp -r ${RunTemplate}/*  ${runDir}
	cd $runDir

	# Link to GEOS-Chem executable instead of having a copy in each rundir
	rm -rf gcclassic
	ln -s ../../${RunTemplate}/gcclassic .

	# Link to restart file
	if "$DoSpinup"; then
	    ln -s ../../spinup_run/GEOSChem.Restart.${SpinupEnd}_0000z.nc4 GEOSChem.Restart.${StartDate}_0000z.nc4
	else
	    ln -s $RestartFile GEOSChem.Restart.${StartDate}_0000z.nc4
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

    ### Perform dry run if requested, only for base run
    if [ $x -eq 0 ]; then
        if "$ProductionDryRun"; then
            printf "Executing dry-run for production runs...\n"
            ./gcclassic --dryrun &> log.dryrun
            ./download_data.py --aws log.dryrun
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
	mkdir -p inversion/data_TROPOMI
	cp -rfP /home/ubuntu/backup_files/input_data/ ${MyPath}/
    else
	ln -s /n/holylfs05/LABS/jacob_lab/lshen/CH4/TROPOMI/data inversion/data_TROPOMI
    fi
    cp ${InversionPath}/PostprocessingScripts/CH4_TROPOMI_INV/*.py inversion/
    cp ${InversionPath}/PostprocessingScripts/CH4_TROPOMI_INV/run_inversion.sh inversion/
    sed -i -e "s:{START}:${StartDate}:g" \
           -e "s:{END}:${EndDate}:g" \
	   -e "s:{STATE_VECTOR_ELEMENTS}:${nElements}:g" \
	   -e "s:{BUFFER_CLUSTERS}:${nBufferClusters}:g" \
	   -e "s:{MY_PATH}:${MyPath}:g" \
	   -e "s:{RUN_NAME}:${RunName}:g" \
	   -e "s:{STATE_VECTOR_PATH}:${StateVectorFile}:g" \
	   -e "s:{LON_MIN}:${LonMin}:g" \
	   -e "s:{LON_MAX}:${LonMax}:g" \
	   -e "s:{LAT_MIN}:${LatMin}:g" \
	   -e "s:{LAT_MAX}:${LatMax}:g" \
	   -e "s:{PRIOR_ERROR}:${PriorError}:g" \
	   -e "s:{OBS_ERROR}:${ObsError}:g" \
	   -e "s:{GAMMA}:${Gamma}:g" \
	   -e "s:{IS_AWS}:${IsAWS}:g" inversion/run_inversion.sh

    if "$isAWS"; then
       sed -i -e "/#SBATCH -p huce_intel/d" \
	      -e "/#SBATCH -t/d" \
	      -e "/#SBATCH --mem/d" \
	      -e "s:#SBATCH -n 1:#SBATCH -n ${cpu_count}:g" inversion/run_inversion.sh
    fi
    
    printf "=== DONE SETTING UP INVERSION DIRECTORY ===\n"

fi #SetupInversion

# Copy sample cluster files (djv: remove cluster terminology)
if "$isAWS"; then
    cp -rfP /home/ubuntu/backup_files/cluster_files/* /home/ubuntu/ExtData/HEMCO/
fi

exit 0

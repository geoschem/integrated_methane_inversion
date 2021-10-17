#!/bin/bash

# This script will set up CH4 analytical inversions with GEOS-Chem. See
# setup_ch4_inversion_instructions.txt for details (mps, 2/20/2020)

##=======================================================================
## User settings **MODIFY AS NEEDED**
##=======================================================================

# Path to inversion setup
UMIpath=$(pwd -P)

# Run IMI on AWS? If false, a local cluster will be assumed
isAWS=false

# Turn on/off different steps. This will allow you to come back to this
# script and set up different stages later.
CreateClusterFile=true    # If false, set ClusterFile below
SetupTemplateRundir=true
SetupSpinupRun=true
SetupJacobianRuns=true
SetupInversion=true
SetupPosteriorRun=true

# AWS only: Download missing GEOS-Chem input data from S3 (you will be charged)
if "$isAWS"; then
    SpinupDryrun=true      # Met fields/emissions for spinup run
    ProductionDryRun=true  # Met fields/emissions for production runs
    BCdryrun=true          # Boundary condition files
else
    SpinupDryrun=false
    ProductionDryRun=false
    BCdryrun=false
fi

# Name for this run
RunName="Test_Permian"

# Start and end date for the spinup simulation
DoSpinup=true
SpinupStart=20180401
SpinupEnd=20180501

# Start and end date for the production simulations
StartDate=20180501
EndDate=20180508

# Path where you want to set up CH4 inversion code and run directories
if "$isAWS"; then
    MyPath="/home/ubuntu/CH4_Workflow"
else
    MyPath="/n/holyscratch01/jacob_lab/msulprizio/CH4"

    # Environment files (specific to Harvard's Cannon cluster)
    NCOEnv="${UMIpath}/envs/Harvard-Cannon/gcc.ifort17_cannon.env"
    GCCEnv="${UMIpath}/envs/Harvard-Cannon/gcc.gfortran10.2_cannon.env"
    CondaEnv="ch4_inv" # See envs/README to create this environment
fi

# Path to find non-emissions input data
if "$isAWS"; then
    DataPath="/home/ubuntu/ExtData"
else
    DataPath="/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData"
fi

# Path to cluster file
ClusterFile="Clusters.nc"

# Path to initial restart file
UseBCsForRestart=true
if "$isAWS"; then
    RestartDownload=true # automatically download restart file
    RestartFile="${DataPath}/BoundaryConditions/GEOSChem.BoundaryConditions.${SpinupStart}_0000z.nc4"
else
    RestartFile="/n/seasasfs02/CH4_inversion/InputData/BoundaryConditions/OutputDir_bias_corrected_dk_2/GEOSChem.BoundaryConditions.${SpinupStart}_0000z.nc4"
fi
    
# Path to boundary condition files (for nested grid simulations)
# Must put backslash before $ in $YYYY$MM$DD to properly work in sed command
if "$isAWS"; then
    BCfiles="${DataPath}/BoundaryConditions"
else
    BCfiles="/n/seasasfs02/CH4_inversion/InputData/BoundaryConditions/OutputDir_bias_corrected_dk_2/GEOSChem.BoundaryConditions.\$YYYY\$MM\$DD_0000z.nc4"
fi

# Grid settings (Permian Basin example)
Res="0.25x0.3125"
Met="geosfp"
LonMin=-111.0
LonMax=-95.0
Lons="${LonMin} ${LonMax}"
LatMin=24.0
LatMax=39.0
Lats="${LatMin} ${LatMax}"
HalfPolar="F"
Levs="47"
NestedGrid="T"
Region="NA"  # NA,AS,CH,EU
Buffer="3 3 3 3"

# Jacobian settings
PerturbValue="1.5"

# Inversion settings
PriorError=0.5
ObsError=15
Gamma=0.25

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
	Start=${SpinupStart}
    else
	START=${StartDate}
    fi
    echo "Downloading boundary condition data for $START to $EndDate"
    python download_bc.py ${Start} ${EndDate} ${BCfiles}
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
GCClassicPath="${UMIpath}/GCClassic"
RunFilesPath="${GCClassicPath}/run"

# Create working directory if it doesn't exist yet
mkdir -p -v ${MyPath}/$RunName

##=======================================================================
## Create cluster file
##=======================================================================
if "$CreateClusterFile"; then

    printf "\n=== CREATING CLUSTER FILE ===\n"
    
    # Use GEOS-FP or MERRA-2 CN file to determine ocean/land grid boxes
    LandCoverFile="${DataPath}/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.nc"
    LandThreshold=0.25

    # Output path and filename for cluster file
    ClusterFile="Clusters.nc"

    # Width of k-means buffer area in degrees (default=5, approx 500 km)
    BufferDeg=5

    # Number of clusters for k-means (default=8)
    kClusters=8

    # Create cluster file
    cd ${MyPath}/$RunName
    mkdir -p -v ClusterFile
    cd ClusterFile

    # Copy cluster creation script to working directory
    cp ${UMIpath}/PostprocessingScripts/CH4_TROPOMI_INV/make_cluster_file.py .
    chmod 755 make_cluster_file.py

    # Activate Conda environment
    printf "Activating conda environment: ${CondaEnv}\n"
    source activate $CondaEnv
    
    printf "Calling make_cluster_file.py\n"
    python make_cluster_file.py $LandCoverFile $ClusterFile $LatMin $LatMax $LonMin $LonMax $BufferDeg $LandThreshold $kClusters

    conda deactivate
    
    printf "=== DONE CREATING CLUSTER FILE ===\n"

fi

# Load environment with NCO
source ${NCOEnv}

# Determine number of clusters from file
function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/foo.nc | cut -f 3- -d ' ' ; }
nClusters=$(ncmax Clusters $ClusterFile)
printf "\n Number of clusters in this inversion= ${nClusters}\n"

# Purge software modules
module purge

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

    # Modify path to cluster file in HEMCO_Config.rc
    OLD=" Clusters.nc"
    if "$CreateClusterFile"; then
	NEW=" ${MyPath}/${RunName}/ClusterFile/${ClusterFile}"
    else
	NEW=" ${ClusterFile}"
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
    source ${GCCEnv}
    
    # Compile GEOS-Chem and store executable in template run directory
    mkdir build; cd build
    cmake ${UMIpath}/GCClassic
    cmake . -DRUNDIR=..
    make -j install
    cd ..
    rm -rf build

    # Purge software modules
    module purge
    
    # Navigate back to top-level directory
    cd ..

    printf "=== DONE CREATING TEMPLATE RUN DIRECTORY ===\n"

fi # SetupTemplateRunDir

##=======================================================================
##  Set up spinup run directory
##=======================================================================

# Get max process count for spinup, production, and run_inversion scripts
output=$(echo $(slurmd -C))
array=($output)
cpu_str=$(echo ${array[1]})
cpu_count=$(echo ${cpu_str:5})

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
	   -e "s|pertpert|1.0|g" \
           -e "s|clustnumclustnum|0|g" input.geos

    # Create run script from template
    sed -e "s:namename:${SpinupName}:g" \
	-e "s:##:#:g" ch4_run.template > ${SpinupName}.run
    chmod 755 ${SpinupName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
	sed -i -e "/#SBATCH -p huce_intel/d" \
	       -e "/#SBATCH -t/d" \
	       -e "/#SBATCH --mem/d" \
               -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${spinup_name}.run
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
	   -e "s|pertpert|1.0|g" \
           -e "s|clustnumclustnum|0|g" input.geos

    # Create run script from template
    sed -e "s:namename:${SpinupName}:g" \
	-e "s:##:#:g" ch4_run.template > ${PosteriorName}.run
    chmod 755 ${PosteriorName}.run
    rm -f ch4_run.template

    if "$isAWS"; then
	sed -i -e "/#SBATCH -p huce_intel/d" \
	       -e "/#SBATCH -t/d" \
	       -e "/#SBATCH --mem/d" \
	       -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${posterior_name}.run
    fi
    
    # Print messages
    printf "\nNote: You will need to manually modify HEMCO_Config.rc to apply the appropriate scale factors.\n"

    ### Perform dry run if requested
    if "$ProductionDryrun"; then
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
    sed -i -e "s:{START}:0:g" -e "s:{END}:${nClusters}:g" jacobian_runs/submit_jacobian_simulations_array.sh

    # Initialize (x=0 is base run, i.e. no perturbation; x=1 is cluster=1; etc.)
    x=0

    # Create run directory for each cluster so we can apply perturbation to each
    while [ $x -le $nClusters ];do

	# Current cluster
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
	sed -i -e "s:pertpert:${PerturbValue}:g" \
               -e "s:clustnumclustnum:${xUSE}:g" input.geos

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
    if "#isAWS"; then
	mkdir -p inversion/data_TROPOMI
	cp -rfP /home/ubuntu/backup_files/input_data/ ${MyPath}/
    else
	ln -s /n/holylfs/LABS/jacob_lab/lshen/CH4/TROPOMI/data inversion/data_TROPOMI ${MyPath}/data_TROPOMI
    fi
    cp ${UMIpath}/PostprocessingScripts/CH4_TROPOMI_INV/*.py inversion/
    cp ${UMIpath}/PostprocessingScripts/CH4_TROPOMI_INV/run_inversion.sh inversion/
    sed -i -e "s:{CLUSTERS}:${nClusters}:g" \
	   -e "s:{START}:${StartDate}:g" \
           -e "s:{END}:${EndDate}:g" \
	   -e "s:{MY_PATH}:${MyPath}:g" \
	   -e "s:{RUN_NAME}:${RunName}:g" \
	   -e "s:{CLUSTER_PATH}:${ClusterFile}:g" \
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

# Copy sample cluster files
if "$isAWS"; then
    cp -rfP /home/ubuntu/backup_files/cluster_files/* /home/ubuntu/ExtData/HEMCO/
fi

exit 0

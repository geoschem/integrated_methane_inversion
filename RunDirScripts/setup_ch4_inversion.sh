#!/bin/bash

# ------------------------------------------------------------------
### Set up GEOS-Chem for Jacobian run (mps, 2/20/2020)
# ------------------------------------------------------------------

##=======================================================================
## Set variables

GC_VERSION=12.7.1

# Path where you want to set up CH4 inversion code and run directories
MY_PATH="/n/holyscratch01/jacob_lab/msulprizio/CH4"

# Path to find non-emissions input data
DATA_PATH="/n/holylfs/EXTERNAL_REPOS/GEOS-CHEM/gcgrid/data/ExtData"

# Path to initial restart file
RESTART_FILE="/n/seasasfs02/CH4_inversion/InputData/Restarts/GEOSChem.Restart.20090101_0000z.nc4"

# Path to boundary condition files (for nested grid simulations)
BC_FILES="n/holyscratch01/jacob_lab/hnesser/GC_TROPOMI_bias_rundirs/Nested_NA/run_dirs/rundir_LSinv/BCs_nc/BC.$YYYY$MM$DD.nc"

# Start and end date fo the simulations
START_DATE=20180501
END_DATE=20180601

# Grid settings
RES="4x5"
MET="merra2"
LONS="-180.0 180.0"
LATS=" -90.0  90.0"
HPOLAR="T"
LEVS="47"
NEST="F"
BUFFER="0 0 0 0"
gridDir="${RES}" # for use METDIR in HEMCO_Config.rc

# Jacobian settings
START_I=0
END_I=0
pPERT="1.5"
RUN_NAME="CH4_Jacobian"
RUN_TEMPLATE="${RES}_template"

# Turn on observation operators and planeflight diagnostics?
GOSAT=false
TCCON=false
UseEmisSF=false
UseSeparateWetlandSF=false
UseOHSF=false
PLANEFLIGHT=false

### Number of Clusters
start=$START_I
stop=$END_I
x=$start

##=======================================================================
## Get source code and run directories

# Copy source code with CH4 analytical inversion updates to your space
# Make sure branch with latest CH4 inversion updates is checked out
cp -r /n/seasasfs02/CH4_inversion/Code.CH4_Inv .
cd Code.CH4_Inv
git checkout CH4_Analytical_Inversion
cd ..

# Copy Unit Tester to create run directory to your space
# Make sure branch with latest CH4 inversion updates is checked out
cp -r /n/seasasfs02/CH4_inversion/UnitTester.CH4_Inv .
cd UnitTester.CH4_Inv
git checkout CH4_Analytical_Inversion
cd ..

# Copy run directory with template files directly from unit tester
RUN_SCRIPTS="/n/seasasfs02/CH4_inversion/RunDirScripts"
mkdir -p $RUN_NAME
cd $RUN_NAME
mkdir -p run_dirs
cp ${RUN_SCRIPTS}/submit_array_jobs run_dirs/
sed -e "s:{RunName}:${RUN_NAME}:g" run_dirs/submit_array_jobs
cp ${RUN_SCRIPTS}/run_array_job run_dirs/
sed -e "s:{START}:${START_I}:g" -e "s:{END}:${END_I}:g" run_dirs/run_array_job
cp ${RUN_SCRIPTS}/rundir_check.sh run_dirs/
mkdir -p bin
if [ "$NEST" == "T" ]; then
  cp -rLv ${MY_PATH}/UnitTester.CH4_Inv/runs/${MET}_*_CH4_na $RUN_TEMPLATE
else
  cp -rLv ${MY_PATH}/UnitTester.CH4_Inv/runs/${RES}_CH4 $RUN_TEMPLATE
fi

# Set up template run directory
cd $RUN_TEMPLATE
cp ${MY_PATH}/UnitTester.CH4_Inv/runs/shared_inputs/Makefiles/Makefile .
cp ${MY_PATH}/UnitTester.CH4_Inv/perl/getRunInfo .
cp ${RUN_SCRIPTS}/run.template .
ln -s $RESTART_FILE .
mkdir OutputDir
cd ..

# Define met and grid fields for HEMCO_Config.rc
if [ "$MET" == "geosfp" ]; then
  metDir="GEOS_FP"
  native="0.25x0.3125"
elif [ "$MET" == "merra2" ]; then
  metDir="MERRA2"
  native="0.5x0.625"
fi

##=======================================================================
##  Create run directories

while [ $x -le $stop ];do

   ### Positive or negative perturbation
   if [ $x -eq -1 ]; then
      PERT="1.0"
      xUSE=$x
   else
      PERT=$pPERT
      xUSE=$x
   fi

   ### Add zeros to string name
   if [ $x -lt 10 ]; then
      xstr="000${x}"
   elif [ $x -lt 100 ]; then
      xstr="00${x}"
   elif [ $x -lt 1000 ]; then
      xstr="0${x}"
   else
      xstr="${x}"
   fi

   ### Define the run directory name
   name="${RUN_NAME}_${xstr}"

   ### Make the directory
   runDir="./run_dirs/${name}"
   mkdir -p ${runDir}
   mkdir -p ${runDir}/Plane_Logs
   mkdir -p ${runDir}/Restarts
   ln -s ${MY_PATH}/Code.CH4_Inv CodeDir

   ### Copy and point to the necessary data
   cp -r ${RUN_TEMPLATE}/*  ${runDir}
   cd $runDir

   ### Create input.geos file from template
   InputFile="input.geos.template"
   sed -e "s:{DATE1}:${START_DATE}:g" \
       -e "s:{DATE2}:${END_DATE}:g" \
       -e "s:{TIME1}:000000:g" \
       -e "s:{TIME2}:000000:g" \
       -e "s:{MET}:${MET}:g" \
       -e "s:{DATA_ROOT}:${DATA_PATH}:g" \
       -e "s:{SIM}:CH4:g" \
       -e "s:{RES}:${RES}:g" \
       -e "s:{LON_RANGE}:${LONS}:g" \
       -e "s:{LAT_RANGE}:${LATS}:g" \
       -e "s:{HALF_POLAR}:${HPOLAR}:g" \
       -e "s:{NLEV}:${LEVS}:g" \
       -e "s:{NESTED_SIM}:${NEST}:g" \
       -e "s:{BUFFER_ZONE}:${BUFFER}:g" \
       -e "s:pertpert:${PERT}:g" \
       -e "s:clustnumclustnum:${xUSE}:g" \
       $InputFile > input.geos.temp
   mv input.geos.temp input.geos
   rm input.geos.template
   
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
       OLD="Turn on plane flt diag? : F"
       NEW="Turn on plane flt diag? : T"
       sed -i "s/$OLD/$NEW/g" input.geos
       OLD="Flight track info file  : Planeflight.dat.YYYYMMDD"
       NEW="Flight track info file  : Planeflights/Planeflight.dat.YYYYMMDD"
       sed -i "s/$OLD/$NEW/g" input.geos
       OLD="Output file name        : plane.log.YYYYMMDD"
       NEW="Output file name        : Plane_Logs/plane.log.YYYYMMDD"
       sed -i "s/$OLD/$NEW/g" input.geos
   fi

   ### Set up HEMCO_Config.rc
   ### Use monthly emissions diagnostic output for now
   sed -e "s:End:Monthly:g" \
       -e "s:{VERBOSE}:0:g" \
       -e "s:{WARNINGS}:1:g" \
       -e "s:{DATA_ROOT}:${DATA_PATH}:g" \
       -e "s:{GRID_DIR}:${gridDir}:g" \
       -e "s:{MET_DIR}:${metDir}:g" \
       -e "s:{NATIVE_RES}:${native}:g" \
       -e "s:$ROOT/SAMPLE_BCs/v2019-05/CH4/GEOSChem.BoundaryConditions.$YYYY$MM$DD_$HH$MNz.nc4:${BC_FILES}:g" \
       HEMCO_Config.template > HEMCO_Config.rc
   rm HEMCO_Config.template

   ### Set up HISTORY.rc
   ### use monthly output for now
   sed -e "s:{FREQUENCY}:00000100 000000:g" \
       -e "s:{DURATION}:00000100 000000:g" \
       HISTORY.rc.template > HISTORY.rc
   rm HISTORY.rc.template
   
   ### Create run script from template
   sed -e "s:namename:${name}:g" run.template > ${name}.run
   chmod 755 ${name}.run

   ### Compile code when creating first run directory
   if [ $x -eq $START_I ]; then
       make realclean CODE_DIR=$MY_PATH/Code.CH4_Inv
       make -j4 build TIMERS=1 CODE_DIR=$MY_PATH/Code.CH4_Inv
       cp -av geos ../../bin/
   fi

   ### Navigate back to top-level directory
   cd ../..

   ### Increment
   x=$[$x+1]

   ### Print diagnostics
   echo "CREATED: ${runDir}"

done

echo "=== DONE CREATING JACOBIAN RUN DIRECTORIES ==="

exit 0

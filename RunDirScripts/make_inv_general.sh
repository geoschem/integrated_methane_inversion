#!/bin/bash

# ------------------------------------------------------------------
### Make run directories for Jacobian run
# ------------------------------------------------------------------

##=======================================================================
## 1. Set variables

START_I=0
END_I=0
pPERT="1.5"
BFname=run.template
nameB="J_emis"
resstr="4x5"
mcwd=`pwd`
geos="geos.mp"
fclust="Clusters_45_July.nc"

date1=20180701
date2=20180801
   
### Number of Clusters
start=$START_I
stop=$END_I
x=$start

### Define the base directory
baseDir="${mcwd}/${resstr}_template"

##=======================================================================
## 2. Create run directories

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

   ### Define the name
   name="${nameB}_${xstr}"

   ### Make the directory
   runDir="${mcwd}/run_dirs/${name}"
   mkdir -p ${runDir}
   mkdir -p ${runDir}/plane_logs
   mkdir -p ${runDir}/restarts
      
   ### Copy and point to the necessary data
   cp -r ${baseDir}/*         ${runDir}

   ### Create input.geos file from template
   if [ $x -eq 0 ]; then
     InputFile="input.geos.baserun.template"
   else
     InputFile="input.geos.template"
   fi

   cd $runDir
   sed -e "s:11111111:${date1}:g" \
       -e "s:22222222:${date2}:g" \
       -e "s:pertpert:${PERT}:g" \
       -e "s:clustnumclustnum:${xUSE}:g" \
       $InputFile > input.geos.temp
   mv input.geos.temp input.geos
   rm input.geos*.template

   ### Create run script from template
   sed -e "s:namename:${name}:g" \
       ${BFname} > ${name}.run
       chmod 755 ${runDir}/${name}.run

   ### Compile code when creating first run directory
   if [ $x -eq $START_I ]; then
       make realclean
       make -j${OMP_NUM_THREADS} mpbuild
       mv ${geos} ${mcwd}/bin/
   fi

   ### Create symbolic links to executable and cluster file
   ln -s ${mcwd}/bin/${geos} ${runDir}/
   ln -s ${mcwd}/$fclust ${runDir}/$fclust
   cd $mcwd

   ### Increment
   x=$[$x+1]

   ### Print diagnostics
   echo "CREATED: ./run_dirs/${name}/"

done

exit 0

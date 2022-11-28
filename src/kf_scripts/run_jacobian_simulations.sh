#!/bin/bash

#SBATCH -J CH4_Jacobian
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 2000
#SBATCH -t 0-02:00

### Run directory
RUNDIR="/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/production_permian/CH4_Jacobian/run_dirs"

### Get current task Id
x=${SLURM_ARRAY_TASK_ID}

### The task is going to run a series of 4 jobs. For each job:
for iter in {0..3}; do

  ### Select a cluster Id
  xi=$((4*x+iter))

  ### Add zeros to the cluster Id
  if [ $xi -lt 10 ]; then
     xstr="000${xi}"
  elif [ $xi -lt 100 ]; then
     xstr="00${xi}"
  elif [ $xi -lt 1000 ]; then
     xstr="0${xi}"
  else
     xstr="${xi}"
  fi
 
  ### Run GEOS-Chem in the directory corresponding to the cluster Id
  cd ${RUNDIR}/CH4_Jacobian_${xstr}
  ./CH4_Jacobian_${xstr}.run

done

exit 0

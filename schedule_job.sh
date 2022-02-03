#!/bin/bash

# Description: 
#   running `sbatch <file>; wait;` within a preexisting sbatch 
#   process on aws results in the child sbatch job getting stuck in pending 
#   because all available resources on the instance are utilized by the 
#   driving sbatch command. In this case the child sbatch process 
#   should be run directly. This function schedules a job either by using 
#   sbatch or direct call depending on whether the parent process is using slurm
# Usage:
#   scheduleJob <runscript-name>
#      runscript-name: script to schedule or run job for
scheduleJob() {

    if "$UseSlurm" && "$isAWS"; then
        ./$1
    else
        # Submit job to job scheduler
        sbatch -W $1; wait;
    fi
}
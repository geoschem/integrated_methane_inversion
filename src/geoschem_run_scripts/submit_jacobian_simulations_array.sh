#!/bin/bash
echo "running {END} jacobian simulations" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

if [[ $SchedulerType = "slurm" | $SchedulerType = "tmux" ]]; then
    sbatch --array={START}-{END}{JOBS} --mem $JacobianMemory \
        -c $JacobianCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -W run_jacobian_simulations.sh
elif [[ $SchedulerType = "PBS" ]]; then
    qsub -J {START}-{END}{JOBS} 
        -l nodes=1 \
        -l mem="$JacobianMemory" \
        -l ncpus=$JacobianCPUs \
        -l walltime=$RequestedTime \
        -l site=needed=$SitesNeeded \
        -l model=ivy \
        -sync y run_jacobian_simulations.sh; wait;
fi
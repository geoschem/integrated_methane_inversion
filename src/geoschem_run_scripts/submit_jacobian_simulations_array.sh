#!/bin/bash
echo "running {END} jacobian simulations" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

if [[ $SchedulerType = "slurm" || $SchedulerType = "tmux" ]]; then
    sbatch --array={START}-{END}{JOBS} --mem $RequestedMemory \
        -c $RequestedCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -o imi_output.tmp \
        --open-mode=append \
        -W run_jacobian_simulations.sh
elif [[ $SchedulerType = "PBS" ]]; then
    qsub -J {START}-{END}{JOBS} \
        -lselect=1:ncpus=$RequestedCPUs:mem="$RequestedMemory":model=ivy \
        -l walltime=$RequestedTime \
        -l site=needed=$SitesNeeded \
        -o imi_output.tmp \
        -Wblock=True run_jacobian_simulations.sh
fi

cat imi_output.tmp >> {InversionPath}/imi_output.log
rm imi_output.tmp

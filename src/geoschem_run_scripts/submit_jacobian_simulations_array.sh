#!/bin/bash
echo "running {END} jacobian simulations" >> {InversionPath}/ici_output.log

# remove error status file if present
rm -f .error_status_file.txt

if [[ $SchedulerType = "slurm" | $SchedulerType = "tmux" ]]; then
    sbatch --array={START}-{END}{JOBS} --mem $RequestedMemory \
        -c $RequestedCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -o ici_output.tmp \
        --open-mode=append \
        -W run_jacobian_simulations.sh
elif [[ $SchedulerType = "PBS" ]]; then
    qsub -J {START}-{END}{JOBS} 
        -l nodes=1 \
        -l mem="$RequestedMemory" \
        -l ncpus=$RequestedCPUs \
        -l walltime=$RequestedTime \
        -l site=needed=$SitesNeeded \
        -l model=ivy \
        -o ici_output.tmp \
        -j oe -k oe \  
        -W run_jacobian_simulations.sh
fi

cat ici_output.tmp >> {InversionPath}/ici_output.log
rm ici_output.tmp

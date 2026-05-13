#!/bin/bash
echo "running 9 jacobian simulations" >> /orcd/data/dvaron/001/vedrau/integrated_methane_inversion/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

if [[ $SchedulerType = "slurm" || $SchedulerType = "tmux" ]]; then
    sbatch --array=0-9 --mem $RequestedMemory \
        -c $RequestedCPUs \
        -N 1 \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -o imi_output.tmp \
        --open-mode=append \
        -W run_jacobian_simulations.sh
elif [[ $SchedulerType = "PBS" ]]; then
    qsub -J 0-9 \
        -lselect=1:ncpus=$RequestedCPUs:mem="$RequestedMemory":model=ivy \
        -l walltime=$RequestedTime \
        -l site=needed=$SitesNeeded \
        -o imi_output.tmp \
        -Wblock=True run_jacobian_simulations.sh
fi

cat imi_output.tmp >> /orcd/data/dvaron/001/vedrau/integrated_methane_inversion/imi_output.log
rm imi_output.tmp

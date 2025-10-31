#!/bin/bash
echo "running {END} jacobian simulations" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

if [ "$SchedulerType" == "slurm" ]; then
    sbatch --array={START}-{END}{JOBS} --mem $RequestedMemory \
    -c $RequestedCPUs \
    -t $RequestedTime \
    -p $SchedulerPartition \
    -o imi_output.tmp \
    --open-mode=append \
    -W run_jacobian_simulations.sh
else
    for i in $(seq {START} {END}); do
        export SLURM_ARRAY_TASK_ID=$i
        source run_jacobian_simulations.sh
    done
fi

cat imi_output.tmp >> {InversionPath}/imi_output.log
rm imi_output.tmp
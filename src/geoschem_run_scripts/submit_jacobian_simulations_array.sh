#!/bin/bash

# Determine which simulations of the job array to submit. This is every
# simulation unless EmulatedJacobian is enabled, in which case only a subset
# of the simulations is needed (see get_jacobian_array_indices in jacobian.sh)
JacobianArray=$(get_jacobian_array_indices {START} {END} \
    "${EmulatedJacobian:-false}" "$OptimizeBCs" "$OptimizeOH" "$isRegional")

echo "running jacobian simulations ${JacobianArray}" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

if [[ $SchedulerType = "slurm" || $SchedulerType = "tmux" ]]; then
    sbatch --array=${JacobianArray}{JOBS} --mem $RequestedMemory \
        -c $RequestedCPUs \
        -N 1 \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -o imi_output.tmp \
        --open-mode=append \
        -W run_jacobian_simulations.sh

    cat imi_output.tmp >> {InversionPath}/imi_output.log
    rm imi_output.tmp
elif [[ $SchedulerType = "PBS" ]]; then
    # PBS job arrays only accept a single index range, so submit one blocking
    # array job per contiguous group of indices (e.g. "0,1,13-16" -> 3 jobs)
    IFS=',' read -ra JacobianArrayGroups <<< "${JacobianArray}"
    JacobianArrayPids=()
    for i in "${!JacobianArrayGroups[@]}"; do
        group="${JacobianArrayGroups[$i]}"
        # PBS requires a range, so expand single indices (e.g. "0" -> "0-0")
        [[ "$group" == *-* ]] || group="${group}-${group}"
        qsub -J ${group}{JOBS} \
            -lselect=1:ncpus=$RequestedCPUs:mem="$RequestedMemory":model=ivy \
            -l walltime=$RequestedTime \
            -l site=needed=$SitesNeeded \
            -o imi_output_${i}.tmp \
            -Wblock=True run_jacobian_simulations.sh &
        JacobianArrayPids+=($!)
    done
    for pid in "${JacobianArrayPids[@]}"; do
        wait $pid
    done

    cat imi_output_*.tmp >> {InversionPath}/imi_output.log
    rm imi_output_*.tmp
fi

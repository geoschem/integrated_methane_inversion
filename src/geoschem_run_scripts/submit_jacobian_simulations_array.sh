#!/bin/bash
echo "running {END} jacobian simulations" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

if [[ $SchedulerType = "slurm" || $SchedulerType = "tmux" ]]; then
    if ! jacobian_submission=$(sbatch --parsable --array={START}-{END}{JOBS} --mem $RequestedMemory \
        -c $RequestedCPUs \
        -N 1 \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -o imi_output.tmp \
        --open-mode=append \
        run_jacobian_simulations.sh); then
        echo "ERROR: Failed to submit the Jacobian job array." >&2
        exit 1
    fi

    # --parsable returns job_id, or job_id;cluster for a remote cluster.
    jacobian_array_id=${jacobian_submission%%;*}
    if [[ ! $jacobian_array_id =~ ^[0-9]+$ ]]; then
        echo "ERROR: Invalid Jacobian array job ID: $jacobian_submission" >&2
        exit 1
    fi
    jacobian_cluster_args=()
    if [[ $jacobian_submission == *';'* ]]; then
        jacobian_cluster_args=(--clusters="${jacobian_submission#*;}")
    fi
    echo "Submitted Jacobian job array $jacobian_submission"

    # Avoid sbatch --wait: affected Slurm versions can return before all array
    # tasks finish. Wait for every queued state, including pending/completing.
    # Query array IDs rather than a specific job so a purged ID is simply absent.
    while true; do
        if ! jacobian_queued_ids=$(squeue \
            ${jacobian_cluster_args[@]+"${jacobian_cluster_args[@]}"} \
            --noheader --array --states=all --format='%F'); then
            echo "ERROR: Cannot check Jacobian array $jacobian_array_id; stopping IMI." >&2
            exit 1
        fi
        if ! grep -Fxq "$jacobian_array_id" <<< "$jacobian_queued_ids"; then
            break
        fi
        sleep 30
    done
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

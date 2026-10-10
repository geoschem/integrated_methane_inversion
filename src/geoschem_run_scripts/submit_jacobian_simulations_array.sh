#!/bin/bash
if ! jacobian_pending_ids=$(bash ./run_jacobian_simulations.sh --list-pending {START} {END}); then
    echo "ERROR: Failed to inspect Jacobian outputs before submission." >&2
    exit 1
fi

# remove error status file if present
rm -f .error_status_file.txt

if [[ -z $jacobian_pending_ids ]]; then
    echo "All Jacobian simulations are complete; no array submitted." | tee -a {InversionPath}/imi_output.log
    # This script is normally sourced by run_jacobian, but can also be executed.
    return 0 2>/dev/null || exit 0
fi
echo "Submitting Jacobian simulations: $jacobian_pending_ids" | tee -a {InversionPath}/imi_output.log

if [[ $SchedulerType = "slurm" || $SchedulerType = "tmux" ]]; then
    if ! jacobian_submission=$(sbatch --parsable --array="${jacobian_pending_ids}{JOBS}" --mem $RequestedMemory \
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
    qsub -J "${jacobian_pending_ids}{JOBS}" \
        -lselect=1:ncpus=$RequestedCPUs:mem="$RequestedMemory":model=ivy \
        -l walltime=$RequestedTime \
        -l site=needed=$SitesNeeded \
        -o imi_output.tmp \
        -Wblock=True run_jacobian_simulations.sh
fi

cat imi_output.tmp >> {InversionPath}/imi_output.log
rm imi_output.tmp

#!/bin/bash

# Determine which simulations of the job array to submit. This is every
# simulation unless EmulatedJacobian is enabled, in which case only a subset
# of the simulations is needed (see get_jacobian_array_indices in jacobian.sh)
JacobianArray=$(get_jacobian_array_indices {START} {END} \
    "${EmulatedJacobian:-false}" "$OptimizeBCs" "$OptimizeOH" "$isRegional")

echo "running jacobian simulations ${JacobianArray}" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

sbatch --array=${JacobianArray}{JOBS} --mem $RequestedMemory \
-c 1 \
-N $NUM_NODES \
-n $TOTAL_CORES \
-t $RequestedTime \
-p $SchedulerPartition \
-o imi_output.tmp \
--open-mode=append \
-W run_jacobian_simulations.sh

cat imi_output.tmp >> {InversionPath}/imi_output.log
rm imi_output.tmp

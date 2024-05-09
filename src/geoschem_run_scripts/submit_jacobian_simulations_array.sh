#!/bin/bash
echo "running {END} jacobian simulations" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

sbatch --array={START}-{END}{JOBS} --mem $JacobianMemory \
-c $JacobianCPUs \
-t $RequestedTime \
-p $JacobianSchedulerPartition \
-o imi_output.tmp \
--open-mode=append \
-W run_jacobian_simulations.sh

cat imi_output.tmp >> {InversionPath}/imi_output.log
rm imi_output.tmp
#!/bin/bash
echo "running {END} jacobian simulations" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

sbatch --array={START}-{END} --mem $JacobianMemory -c $JacobianCPUs -t $RequestedTime -W run_jacobian_simulations.sh

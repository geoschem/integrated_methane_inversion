#!/bin/bash
echo "running {END} jacobian simulations" >> {InversionPath}/imi_output.log

# remove error status file if present
rm -f .error_status_file.txt

sbatch --array={START}-{END} -W run_jacobian_simulations.sh

exit 0

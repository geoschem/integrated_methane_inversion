#!/bin/bash

sbatch --array={START}-{END} -W run_jacobian_simulations.sh

exit 0

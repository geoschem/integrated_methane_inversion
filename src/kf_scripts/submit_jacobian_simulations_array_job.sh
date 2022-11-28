#!/bin/bash

sbatch --array=0-60 -W run_jacobian_simulations.sh

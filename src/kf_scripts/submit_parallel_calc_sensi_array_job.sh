#!/bin/bash

sbatch --array=0-6 -W run_parallel_calc_sensi.sh $1 $2 $3 $4

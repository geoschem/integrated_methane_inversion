#!/bin/bash
#SBATCH -J test_job
#SBATCH -o test_job.out
#SBATCH -p debug
#SBATCH -N 1
#SBATCH --mem 100
#SBATCH --ntasks-per-node 1
#SBATCH -t 00:05:00

echo "Hello from Slurm job!"
sleep 3
echo "Job completed."

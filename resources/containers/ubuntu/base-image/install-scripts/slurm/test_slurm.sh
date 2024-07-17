#!/bin/bash
#SBATCH --job-name=test_job
#SBATCH --output=test_job.out
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --mem=100
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:05:00

echo "Hello from Slurm job!"
sleep 3
echo "Job completed."

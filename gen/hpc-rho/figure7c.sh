#!/bin/bash
#SBATCH --account=torch_pr_921_general
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=generate-det-jobs
#SBATCH --output=generate-det-rho-tau_%A_%a.out
#SBATCH --error=generate-det-spd-mu_%A_%a.err
#SBATCH --array=1-51

# Give Julia the same number of threads as allocated CPUs
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Print some info for debugging
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Julia threads: $JULIA_NUM_THREADS"

# Run the appropriate Julia script
julia generate-heat-rho-tau.jl $SLURM_ARRAY_TASK_ID

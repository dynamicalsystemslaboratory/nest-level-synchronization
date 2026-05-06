#!/bin/bash
#SBATCH --account=torch_pr_921_general
#SBATCH --mem=32G
#SBATCH --cpus-per-task=48
#SBATCH --job-name=generate-det-jobs
#SBATCH --output=generate-det-%A_%a.out
#SBATCH --error=generate-det-%A_%a.err
#SBATCH --array=1-17

# Give Julia the same number of threads as allocated CPUs
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Start from initial conditions?
SYSTEM_SIZE=$1
LOAD_STEPS=${2:-true}

# Print some info for debugging
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Julia threads: $JULIA_NUM_THREADS"

# Run the appropriate Julia script
julia generate-det-spd.jl $SYSTEM_SIZE $SLURM_ARRAY_TASK_ID $LOAD_STEPS

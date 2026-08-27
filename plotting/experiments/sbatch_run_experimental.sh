#!/bin/bash
#SBATCH --job-name=plot-exp
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --chdir=/storage/biology/projects/miller-lowry/beitner/binning-classification-visualization
#SBATCH --output=plots/sandbox/sbatch_%j.log
#SBATCH --error=plots/sandbox/sbatch_%j.err

set -e

# Activate visualization conda environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate visualization

echo "================================"
echo "Experimental Plot Generation"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "================================"

# Change to plotting directory
cd /storage/biology/projects/miller-lowry/beitner/binning-classification-visualization/plotting

# Run experimental plots
python experiments/run_experimental.py "$@"

echo "✓ Plotting job complete"

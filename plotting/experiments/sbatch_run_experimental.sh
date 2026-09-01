#!/bin/bash
#SBATCH --job-name=plot-exp
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --chdir=/storage/biology/projects/miller-lowry/beitner/binning-classification-visualization
#SBATCH --output=plots/sandbox/sbatch_%j.log
#SBATCH --error=plots/sandbox/sbatch_%j.err

set -e

# Change to plotting directory
cd /storage/biology/projects/miller-lowry/beitner/binning-classification-visualization/plotting

# Run with conda run (more reliable in SBATCH)
conda run -n visualization python experiments/run_experimental.py "$@"

echo "✓ Plotting job complete"

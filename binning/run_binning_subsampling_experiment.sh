#!/usr/bin/env bash
#SBATCH --job-name=binning_wrapper
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --output=logs/binning_wrapper_%j.out
#SBATCH --error=logs/binning_wrapper_%j.err
#SBATCH --chdir=/storage/biology/projects/miller-lowry/beitner/binning-classification-visualization

set -euo pipefail

# Use absolute path to avoid issues with SBATCH script copying
BINNING_DIR="/storage/biology/projects/miller-lowry/beitner/binning-classification-visualization/binning"
CONFIG_DIR="/storage/biology/projects/miller-lowry/beitner/binning-classification-visualization/configs"
DEFAULT_CONFIG="${CONFIG_DIR}/binning_subsampling_experiment.conf"
COMMON_RUNNER="${BINNING_DIR}/run_binning_experiment_common.sh"

if [[ "$#" -gt 0 && "$1" == "--config" ]]; then
    CONFIG_PATH="$2"
    shift 2
else
    CONFIG_PATH="$DEFAULT_CONFIG"
fi

bash "$COMMON_RUNNER" --config "$CONFIG_PATH" "$@"

if [[ "$#" -gt 0 && "$1" == "--config" ]]; then
    CONFIG_PATH="$2"
    shift 2
else
    CONFIG_PATH="$DEFAULT_CONFIG"
fi

bash "$COMMON_RUNNER" --config "$CONFIG_PATH" "$@"

#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/storage/biology/projects/miller-lowry/beitner"
WORK_DIR="${BASE_DIR}/binning-classification-wrapper"
DEFAULT_CONFIG="${WORK_DIR}/scripts/configs/binning_gold_standard_mag_library.conf"

if [[ "$#" -gt 0 && "$1" == "--config" ]]; then
    CONFIG_PATH="$2"
    shift 2
else
    CONFIG_PATH="$DEFAULT_CONFIG"
fi

bash "${WORK_DIR}/scripts/run_binning_experiment_common.sh" --config "$CONFIG_PATH" "$@"

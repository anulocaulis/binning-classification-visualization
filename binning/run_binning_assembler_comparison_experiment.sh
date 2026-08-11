#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/../configs"
DEFAULT_CONFIG="${CONFIG_DIR}/binning_assembler_comparison.conf"
COMMON_RUNNER="${SCRIPT_DIR}/run_binning_experiment_common.sh"

if [[ "$#" -gt 0 && "$1" == "--config" ]]; then
    CONFIG_PATH="$2"
    shift 2
else
    CONFIG_PATH="$DEFAULT_CONFIG"
fi

bash "$COMMON_RUNNER" --config "$CONFIG_PATH" "$@"

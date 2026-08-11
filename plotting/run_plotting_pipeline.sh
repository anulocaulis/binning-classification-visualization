#!/bin/bash
#SBATCH --job-name=plotting_pipeline
#SBATCH --output=logs/plotting_pipeline_%j.out
#SBATCH --error=logs/plotting_pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --partition=math-alderaan

# Master plotting pipeline wrapper.
#
# This SLURM job sources plotting_pipeline_config.sh and runs the core stats,
# filtering, and plotting steps in a fixed sequence with shared runtime and
# style configuration.
#
# Wrapped shell scripts:
#   1. filter_contigs_min_length.sh
#   2. filter_subsampling_contigs_min_length.sh
#   3. summary_stats.sh
#   4. ge1000_summary_stats.sh
#   5. subsample_assembly_summary_stats.sh
#   6. sequencing_run_summary_stats.sh
#
# Wrapped Python entry points:
#   7. plot_summary_stats.py
#   8. poster_figures.py
#   9. plot_subsample_sequencing_runs.py
#
# Execution order:
#   filter full assemblies -> filter subsample assemblies -> full stats ->
#   ge1000 stats -> subsample assembly stats -> sequencing-run stats ->
#   full summary plots -> filtered/ge1000 poster plots -> subsample plots
#
# The config file controls:
#   - input/output root paths
#   - expected sample/assembler counts for validation
#   - shared long/short/hybrid color palettes
#   - max rarefaction curve downsampling points
#   - step toggles and runtime/container settings

# Usage:
#   sbatch run_plotting_pipeline.sh [CONFIG_FILE]
# Example:
#   sbatch run_plotting_pipeline.sh ../configs/plotting_pipeline_config.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_PATH="${1:-$SCRIPT_DIR/../configs/plotting_pipeline_config.sh}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG_PATH"

export QC_TOOLS_CONTAINER="${QC_TOOLS_CONTAINER:-}"

mkdir -p "$SCRIPT_DIR/logs"
cd "$SCRIPT_DIR"

PLOT_IMPORT_CHECK='import pandas, matplotlib, seaborn'
RUN_MODE=""

if [[ -n "${CONDA_ENV_NAME:-}" ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "CONDA_ENV_NAME is set but 'conda' is unavailable in PATH." >&2
    exit 1
  fi
  if conda run -n "$CONDA_ENV_NAME" python -c "$PLOT_IMPORT_CHECK" >/dev/null 2>&1; then
    RUN_MODE="conda"
  else
    echo "Conda env '$CONDA_ENV_NAME' is missing pandas/matplotlib/seaborn." >&2
    exit 1
  fi
elif command -v "$PYTHON_BIN" >/dev/null 2>&1 && "$PYTHON_BIN" -c "$PLOT_IMPORT_CHECK" >/dev/null 2>&1; then
  RUN_MODE="host"
else
  echo "No plotting runtime found with pandas, matplotlib, and seaborn." >&2
  exit 1
fi

run_python() {
  local script_path="$1"
  shift

  if [[ "$RUN_MODE" == "conda" ]]; then
    conda run -n "$CONDA_ENV_NAME" python "$script_path" "$@"
  else
    "$PYTHON_BIN" "$script_path" "$@"
  fi
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    echo "[ERROR] Missing ${label}: $path" >&2
    exit 1
  fi
}

require_dir() {
  local path="$1"
  local label="$2"
  if [[ ! -d "$path" ]]; then
    echo "[ERROR] Missing ${label}: $path" >&2
    exit 1
  fi
}

validate_pipeline_setup() {
  require_file "$SCRIPT_DIR/filter_contigs_min_length.sh" "script"
  require_file "$SCRIPT_DIR/filter_subsampling_contigs_min_length.sh" "script"
  require_file "$SCRIPT_DIR/summary_stats.sh" "script"
  require_file "$SCRIPT_DIR/ge1000_summary_stats.sh" "script"
  require_file "$SCRIPT_DIR/subsample_assembly_summary_stats.sh" "script"
  require_file "$SCRIPT_DIR/sequencing_run_summary_stats.sh" "script"
  require_file "$SCRIPT_DIR/plot_summary_stats.py" "script"
  require_file "$SCRIPT_DIR/poster_figures.py" "script"
  require_file "$SCRIPT_DIR/plot_subsample_sequencing_runs.py" "script"

  if [[ -n "${QC_TOOLS_CONTAINER:-}" ]]; then
    require_file "$QC_TOOLS_CONTAINER" "QC container"
  fi

  if [[ "$RUN_FILTER_CONTIGS" -eq 1 || "$RUN_SUMMARY_STATS" -eq 1 || "$RUN_GE1000_SUMMARY_STATS" -eq 1 ]]; then
    require_dir "$FULL_ASSEMBLIES_DIR" "full assemblies directory"
    require_dir "$GE1000_ASSEMBLIES_DIR" "ge1000 assemblies directory"
  fi

  if [[ "$RUN_FILTER_SUBSAMPLING_CONTIGS" -eq 1 || "$RUN_SUBSAMPLE_ASSEMBLY_STATS" -eq 1 ]]; then
    require_dir "$SUBSAMPLING_DIR" "subsampling directory"
  fi

  if [[ "$RUN_SEQUENCING_RUN_STATS" -eq 1 ]]; then
    require_dir "$SEQUENCING_RUNS_DIR" "sequencing runs directory"
  fi

  if [[ "$RUN_POSTER_FIGURES" -eq 1 && ! -f "$GE1000_SUMMARY_LOG" && "$RUN_GE1000_SUMMARY_STATS" -eq 0 ]]; then
    echo "[ERROR] poster_figures.py needs $GE1000_SUMMARY_LOG or RUN_GE1000_SUMMARY_STATS=1" >&2
    exit 1
  fi

  if [[ "$RUN_PLOT_SUMMARY_STATS" -eq 1 && ! -f "$FULL_SUMMARY_LOG" && "$RUN_SUMMARY_STATS" -eq 0 ]]; then
    echo "[ERROR] plot_summary_stats.py needs $FULL_SUMMARY_LOG or RUN_SUMMARY_STATS=1" >&2
    exit 1
  fi

  if [[ "$RUN_PLOT_SUBSAMPLE_SEQUENCING" -eq 1 ]]; then
    if [[ ! -f "$SEQUENCING_CLEAN_TSV" && "$RUN_SEQUENCING_RUN_STATS" -eq 0 ]]; then
      echo "[ERROR] plot_subsample_sequencing_runs.py needs $SEQUENCING_CLEAN_TSV or RUN_SEQUENCING_RUN_STATS=1" >&2
      exit 1
    fi
    if [[ ! -f "$SUBSAMPLE_ASSEMBLY_LOG" && "$RUN_SUBSAMPLE_ASSEMBLY_STATS" -eq 0 ]]; then
      echo "[ERROR] plot_subsample_sequencing_runs.py needs $SUBSAMPLE_ASSEMBLY_LOG or RUN_SUBSAMPLE_ASSEMBLY_STATS=1" >&2
      exit 1
    fi
    if [[ ! -f "$FULL_SUMMARY_LOG" && "$RUN_SUMMARY_STATS" -eq 0 ]]; then
      echo "[ERROR] plot_subsample_sequencing_runs.py needs $FULL_SUMMARY_LOG or RUN_SUMMARY_STATS=1" >&2
      exit 1
    fi
  fi
}

validate_counts_from_log() {
  local log_path="$1"
  local expected_samples="$2"
  local expected_assemblers="$3"
  local sample_count assembler_count

  [[ -f "$log_path" ]] || return 0

  sample_count="$(awk -F': ' '/^Sample: /{print $2}' "$log_path" | sort -u | wc -l | tr -d ' ')"
  assembler_count="$(awk -F': ' '/^Assembler: /{print $2}' "$log_path" | sort -u | wc -l | tr -d ' ')"

  if [[ "$expected_samples" -gt 0 && "$sample_count" -ne "$expected_samples" ]]; then
    echo "[WARN] $log_path has $sample_count unique samples; expected $expected_samples" >&2
  fi
  if [[ "$expected_assemblers" -gt 0 && "$assembler_count" -ne "$expected_assemblers" ]]; then
    echo "[WARN] $log_path has $assembler_count unique assemblers; expected $expected_assemblers" >&2
  fi
}

echo "Using config: $CONFIG_PATH"
echo "Project root: $PROJECT_ROOT"
echo "Python run mode: $RUN_MODE"

validate_pipeline_setup

echo "Validation passed."

if [[ "$RUN_FILTER_CONTIGS" -eq 1 ]]; then
  echo "[STEP] filter_contigs_min_length.sh"
  bash "$SCRIPT_DIR/filter_contigs_min_length.sh"
fi

if [[ "$RUN_FILTER_SUBSAMPLING_CONTIGS" -eq 1 ]]; then
  echo "[STEP] filter_subsampling_contigs_min_length.sh"
  bash "$SCRIPT_DIR/filter_subsampling_contigs_min_length.sh" "$SUBSAMPLING_DIR"
fi

if [[ "$RUN_SUMMARY_STATS" -eq 1 ]]; then
  echo "[STEP] summary_stats.sh"
  bash "$SCRIPT_DIR/summary_stats.sh" "$FULL_ASSEMBLIES_DIR"
  validate_counts_from_log "$FULL_SUMMARY_LOG" "$EXPECTED_FULL_SAMPLE_COUNT" "$EXPECTED_ASSEMBLER_COUNT"
fi

if [[ "$RUN_GE1000_SUMMARY_STATS" -eq 1 ]]; then
  echo "[STEP] ge1000_summary_stats.sh"
  bash "$SCRIPT_DIR/ge1000_summary_stats.sh" "$GE1000_ASSEMBLIES_DIR" "$GE1000_SUMMARY_LOG"
fi

if [[ "$RUN_SUBSAMPLE_ASSEMBLY_STATS" -eq 1 ]]; then
  echo "[STEP] subsample_assembly_summary_stats.sh"
  bash "$SCRIPT_DIR/subsample_assembly_summary_stats.sh" "$SUBSAMPLING_DIR" "$SUBSAMPLE_ASSEMBLY_LOG"
  validate_counts_from_log "$SUBSAMPLE_ASSEMBLY_LOG" "$EXPECTED_SUBSAMPLE_SAMPLE_COUNT" "$EXPECTED_ASSEMBLER_COUNT"
fi

if [[ "$RUN_SEQUENCING_RUN_STATS" -eq 1 ]]; then
  echo "[STEP] sequencing_run_summary_stats.sh"
  bash "$SCRIPT_DIR/sequencing_run_summary_stats.sh" "$SEQUENCING_RUNS_DIR" "$SEQUENCING_LOG" "$SEQUENCING_CLEAN_TSV"
fi

if [[ "$RUN_PLOT_SUMMARY_STATS" -eq 1 ]]; then
  echo "[STEP] plot_summary_stats.py"
  mkdir -p "$FULL_OUT_ROOT/plots/summary_stats" "$FULL_OUT_ROOT/data"
  run_python "$SCRIPT_DIR/plot_summary_stats.py" \
    --log "$FULL_SUMMARY_LOG" \
    --outdir "$FULL_OUT_ROOT/plots/summary_stats" \
    --data-dir "$FULL_OUT_ROOT/data" \
    --plots-root "$FULL_OUT_ROOT/plots"
fi

if [[ "$RUN_POSTER_FIGURES" -eq 1 ]]; then
  echo "[STEP] poster_figures.py"
  mkdir -p "$FILTERED_OUT_ROOT/plots/ge1000_poster" "$FILTERED_OUT_ROOT/data"
  run_python "$SCRIPT_DIR/poster_figures.py" \
    --log "$GE1000_SUMMARY_LOG" \
    --input-csv "$FILTERED_OUT_ROOT/data/rarefaction_data.csv" \
    --summary-csv "$FILTERED_OUT_ROOT/data/parsed_summary_stats.csv" \
    --outdir "$FILTERED_OUT_ROOT/plots/ge1000_poster" \
    --data-dir "$FILTERED_OUT_ROOT/data" \
    --long-colors "$LONG_COLORS" \
    --short-colors "$SHORT_COLORS" \
    --hybrid-colors "$HYBRID_COLORS" \
    --max-tsv-points-per-curve "$MAX_TSV_POINTS_PER_CURVE"
fi

if [[ "$RUN_PLOT_SUBSAMPLE_SEQUENCING" -eq 1 ]]; then
  echo "[STEP] plot_subsample_sequencing_runs.py"
  mkdir -p "$SUBSAMPLE_OUT_ROOT/plots/sequencing_runs" "$SUBSAMPLE_OUT_ROOT/data"
  run_python "$SCRIPT_DIR/plot_subsample_sequencing_runs.py" \
    --input-tsv "$SEQUENCING_CLEAN_TSV" \
    --outdir "$SUBSAMPLE_OUT_ROOT/plots/sequencing_runs" \
    --data-dir "$SUBSAMPLE_OUT_ROOT/data" \
    --subsample-assembly-log "$SUBSAMPLE_ASSEMBLY_LOG" \
    --true-assembly-log "$FULL_SUMMARY_LOG" \
    --max-tsv-points-per-curve "$MAX_TSV_POINTS_PER_CURVE" \
    --long-colors "$LONG_COLORS" \
    --short-colors "$SHORT_COLORS" \
    --hybrid-colors "$HYBRID_COLORS"
fi

echo "Completed plotting pipeline."
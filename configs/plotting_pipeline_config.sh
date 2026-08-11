#!/usr/bin/env bash

# Central configuration for run_plotting_pipeline.sh.
# Adjust paths and counts here rather than editing the wrapper.

# shellcheck source=/dev/null
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$CONFIG_DIR/configs_master.conf"

PROJECT_ROOT="${PROJECT_PLOTTING_ROOT}"

# Container/runtime path used by stats/filtering scripts. Leave empty to let
# individual scripts search their built-in candidate locations.
QC_TOOLS_CONTAINER="${PROJECT_QC_TOOLS_CONTAINER}"

# Input roots
FULL_ASSEMBLIES_DIR="${PROJECT_FULL_ASSEMBLIES_DIR}"
GE1000_ASSEMBLIES_DIR="${PROJECT_GE1000_ASSEMBLIES_DIR}"
SUBSAMPLING_DIR="${PROJECT_SUBSAMPLING_DIR}"
SEQUENCING_RUNS_DIR="${PROJECT_SEQUENCING_RUNS_DIR}"

# Core generated logs/tables
FULL_SUMMARY_LOG="${PROJECT_ROOT}/summary_stats_log.txt"
GE1000_SUMMARY_LOG="${PROJECT_ROOT}/ge1000_summary_stats.txt"
SUBSAMPLE_ASSEMBLY_LOG="${PROJECT_ROOT}/subsample_assembly_summary_stats_log.txt"
SEQUENCING_LOG="${PROJECT_ROOT}/sequencing_run_summary_stats_log.txt"
SEQUENCING_CLEAN_TSV="${PROJECT_ROOT}/sequencing_run_summary_stats_log_clean.tsv"

# Output roots
FULL_OUT_ROOT="${PROJECT_ROOT}/full_assembly"
FILTERED_OUT_ROOT="${PROJECT_ROOT}/filtered_assembly"
SUBSAMPLE_OUT_ROOT="${PROJECT_ROOT}/subsampling"
PROJECT_REFINED_BIN_COUNTS_OUTDIR="${PROJECT_ROOT}/plots/refined_bin_counts"

# Project-level validation metadata. Set to 0 to disable a check.
EXPECTED_FULL_SAMPLE_COUNT="${PROJECT_EXPECTED_FULL_SAMPLE_COUNT}"
EXPECTED_SUBSAMPLE_SAMPLE_COUNT="${PROJECT_EXPECTED_SUBSAMPLE_SAMPLE_COUNT}"
EXPECTED_ASSEMBLER_COUNT="${PROJECT_EXPECTED_ASSEMBLER_COUNT}"

# Global plotting style
LONG_COLORS="#C6DBEF,#41B6C4,#08519C"
SHORT_COLORS="#A1D99B,#74C476,#238B45"
HYBRID_COLORS="#D0D1E6,#6A51A3"
MAX_TSV_POINTS_PER_CURVE=1000

# Runtime preferences
CONDA_ENV_NAME="${CONDA_ENV_NAME:-visualization}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

# Step toggles for the eight canonical scripts.
RUN_FILTER_CONTIGS=1
RUN_FILTER_SUBSAMPLING_CONTIGS=1
RUN_SUMMARY_STATS=1
RUN_GE1000_SUMMARY_STATS=1
RUN_SUBSAMPLE_ASSEMBLY_STATS=1
RUN_SEQUENCING_RUN_STATS=1
RUN_PLOT_SUMMARY_STATS=1
RUN_POSTER_FIGURES=1
RUN_PLOT_SUBSAMPLE_SEQUENCING=1
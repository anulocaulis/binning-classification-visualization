#!/bin/bash
#SBATCH --job-name=ge1000_stats
#SBATCH --output=logs/ge1000_stats_%j.out
#SBATCH --error=logs/ge1000_stats_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --partition=math-alderaan

# Usage:
#   sbatch ge1000_summary_stats.sh [ASSEMBLIES_DIR] [LOG_FILE]
# Example:
#   sbatch ge1000_summary_stats.sh /path/to/data/assemblies ge1000_summary_stats.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
# shellcheck source=/dev/null
if [[ -f "$SCRIPT_DIR/../configs/configs_master.conf" ]]; then
  source "$SCRIPT_DIR/../configs/configs_master.conf"
fi
ASSEMBLIES_DIR="${1:-$BASE_DIR/../data/assemblies}"
LOG_FILE="${2:-ge1000_summary_stats.txt}"
WORKDIR="$(pwd)"
CONTAINER=""

container_candidates=(
  "${QC_TOOLS_CONTAINER:-}"
  "${PROJECT_QC_TOOLS_CONTAINER:-}"
  "$BASE_DIR/containers/qc_tools_miniconda.sif"
  "$BASE_DIR/../containers/qc_tools_miniconda.sif"
  "$BASE_DIR/../Lowry-assemblies/containers/qc_tools_miniconda.sif"
)

for candidate in "${container_candidates[@]}"; do
  if [[ -n "$candidate" && -f "$candidate" ]]; then
    CONTAINER="$candidate"
    break
  fi
done

if [[ ! -f "$CONTAINER" ]]; then
  echo "Container not found. Set QC_TOOLS_CONTAINER or place qc_tools_miniconda.sif in a standard containers path." >&2
  exit 1
fi

if [[ ! -d "$ASSEMBLIES_DIR" ]]; then
  echo "Assemblies directory not found: $ASSEMBLIES_DIR" >&2
  exit 1
fi

mapfile -t assembly_files < <(
  {
    find "$ASSEMBLIES_DIR" -type f \( -name "*ge1000*.fa" -o -name "*ge1000*.fasta" \)
    # Include full Opera-MS assembly in the same summary log alongside ge1000 outputs.
    find "$ASSEMBLIES_DIR" -type f -path "*/assembly.opera_ms/opera_ms/assembly.fasta"
  } | sort -u
)

if [[ ${#assembly_files[@]} -eq 0 ]]; then
  echo "No ge1000 FASTA files (or full Opera-MS assembly FASTA) found under: $ASSEMBLIES_DIR" >&2
  exit 1
fi

# Recreate the log from scratch on each run to avoid duplicate appended blocks.
: > "$LOG_FILE"

processed_files=0

for assembly_file in "${assembly_files[@]}"; do
  rel_path="${assembly_file#"$ASSEMBLIES_DIR"/}"
  sample_dir="$(basename "$(dirname "$assembly_file")")"

  {
    echo "============================================================"
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Sample dir: $sample_dir"
    echo "Assembly: $assembly_file"
    echo "Relative path: $rel_path"
    echo "============================================================"

    singularity exec -B "$WORKDIR:$WORKDIR,$ASSEMBLIES_DIR:$ASSEMBLIES_DIR" "$CONTAINER" \
      stats.sh in="$assembly_file"

    echo
  } >> "$LOG_FILE" 2>&1

  processed_files=$((processed_files + 1))
done

echo "Appended stats for $processed_files FASTA files (ge1000 + full Opera-MS assembly) to $LOG_FILE"

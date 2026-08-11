#!/bin/bash
#SBATCH --job-name=subsample_asm_stats
#SBATCH --output=logs/subsample_asm_stats_%j.out
#SBATCH --error=logs/subsample_asm_stats_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --partition=math-alderaan

# Usage:
#   sbatch subsample_assembly_summary_stats.sh [SUBSAMPLING_DIR] [LOG_FILE]
# Example:
#   sbatch subsample_assembly_summary_stats.sh /path/to/data/subsampling subsample_assembly_summary_stats_log.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
# shellcheck source=/dev/null
if [[ -f "$SCRIPT_DIR/../configs/configs_master.conf" ]]; then
  source "$SCRIPT_DIR/../configs/configs_master.conf"
fi
SUBSAMPLING_DIR="${1:-$BASE_DIR/../data/subsampling}"
LOG_FILE="${2:-subsample_assembly_summary_stats_log.txt}"
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

if [[ ! -d "$SUBSAMPLING_DIR" ]]; then
  echo "Subsampling directory not found: $SUBSAMPLING_DIR" >&2
  exit 1
fi

mapfile -t sample_dirs < <(find "$SUBSAMPLING_DIR" -mindepth 1 -maxdepth 1 -type d -name "S*_subsample_*" | sort)
if [[ ${#sample_dirs[@]} -eq 0 ]]; then
  echo "No subsample directories found under: $SUBSAMPLING_DIR" >&2
  exit 1
fi

# Recreate the log from scratch on each run to avoid duplicate appended blocks.
: > "$LOG_FILE"

processed_files=0

extract_subsample_level() {
  local sample_name="$1"
  if [[ "$sample_name" =~ _subsample_([0-9]+)$ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "unknown"
  fi
}

infer_assembler_from_dir() {
  local assembly_dir="$1"
  local dirname
  dirname="$(basename "$assembly_dir")"
  if [[ "$dirname" == assembly.* ]]; then
    echo "${dirname#assembly.}"
    return
  fi
  echo "$dirname"
}

pick_primary_fasta() {
  local assembly_dir="$1"
  local candidate=""

  if [[ "$(basename "$assembly_dir")" == "opera_ms" ]]; then
    for name in \
      "intermediate_files/user_assembly/contigs.fasta" \
      "contigs.fasta" \
      "assembly.fasta" \
      "final.contigs.fa"; do
      candidate="$assembly_dir/$name"
      if [[ -f "$candidate" ]]; then
        echo "$candidate"
        return 0
      fi
    done
  fi

  for name in \
    "assembly.fasta" \
    "contigs.fasta" \
    "metamdbg.contigs.fasta" \
    "scaffolds.fasta" \
    "final.contigs.fa"; do
    candidate="$assembly_dir/$name"
    if [[ -f "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done

  return 1
}

for sample_dir in "${sample_dirs[@]}"; do
  sample="$(basename "$sample_dir")"
  subsample_level="$(extract_subsample_level "$sample")"

  mapfile -t assembly_dirs < <(
    find "$sample_dir" -mindepth 1 -maxdepth 1 -type d \
      \( -name "assembly.*" -o -name "opera_ms" \) \
      ! -name "*.failed_*" \
      | sort
  )

  if [[ ${#assembly_dirs[@]} -eq 0 ]]; then
    echo "Skipping $sample (no assembly directories found)" >&2
    continue
  fi

  for assembly_dir in "${assembly_dirs[@]}"; do
    assembler="$(infer_assembler_from_dir "$assembly_dir")"

    if ! assembly_file="$(pick_primary_fasta "$assembly_dir")"; then
      echo "Skipping $sample / $assembler (no canonical FASTA found in $assembly_dir)" >&2
      continue
    fi

    rel_path="${assembly_file#"$SUBSAMPLING_DIR"/}"

    {
      echo "============================================================"
      echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
      echo "Sample: $sample"
      echo "Subsample level: $subsample_level"
      echo "Assembler: $assembler"
      echo "Assembly: $assembly_file"
      echo "Relative path: $rel_path"
      echo "============================================================"

      singularity exec -B "$WORKDIR:$WORKDIR,$SUBSAMPLING_DIR:$SUBSAMPLING_DIR" "$CONTAINER" \
        stats.sh in="$assembly_file"

      echo
    } >> "$LOG_FILE" 2>&1

    processed_files=$((processed_files + 1))
  done
done

if [[ $processed_files -eq 0 ]]; then
  echo "No assembly FASTA files processed under: $SUBSAMPLING_DIR" >&2
  exit 1
fi

echo "Appended stats for $processed_files assembly FASTA files to $LOG_FILE"

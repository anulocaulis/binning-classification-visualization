#!/bin/bash
#SBATCH --job-name=sequencing_run_stats
#SBATCH --output=logs/sequencing_run_stats_%j.out
#SBATCH --error=logs/sequencing_run_stats_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --partition=math-alderaan

# Usage:
#   sbatch sequencing_run_summary_stats.sh [RUNS_DIR] [LOG_FILE] [CLEAN_TABLE_FILE]
# Example:
#   sbatch sequencing_run_summary_stats.sh /path/to/data sequencing_run_summary_stats_log.txt sequencing_run_summary_stats_clean.tsv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
# shellcheck source=/dev/null
if [[ -f "$SCRIPT_DIR/../configs/configs_master.conf" ]]; then
  source "$SCRIPT_DIR/../configs/configs_master.conf"
fi
RUNS_DIR="${1:-$BASE_DIR/../data}"
LOG_FILE="${2:-sequencing_run_summary_stats_log.txt}"
if [[ "$LOG_FILE" == *.* ]]; then
  CLEAN_TABLE_FILE_DEFAULT="${LOG_FILE%.*}_clean.tsv"
else
  CLEAN_TABLE_FILE_DEFAULT="${LOG_FILE}_clean.tsv"
fi
CLEAN_TABLE_FILE="${3:-$CLEAN_TABLE_FILE_DEFAULT}"
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

if [[ ! -d "$RUNS_DIR" ]]; then
  echo "Runs directory not found: $RUNS_DIR" >&2
  exit 1
fi

mapfile -t sample_dirs < <(find "$RUNS_DIR" -mindepth 1 -maxdepth 1 -type d -name "S*" | sort)
if [[ ${#sample_dirs[@]} -eq 0 ]]; then
  echo "No sample directories found under: $RUNS_DIR" >&2
  exit 1
fi

processed_files=0
processed_samples=0

if [[ ! -f "$CLEAN_TABLE_FILE" ]]; then
  echo -e "timestamp\tsample\tplatform\tinput\trelative_path\tfile_type\treads\ttotal_bases\ttotal_gbp\tstatus" > "$CLEAN_TABLE_FILE"
fi

infer_platform() {
  local run_file="$1"
  local fname
  fname="$(basename "$run_file" | tr '[:upper:]' '[:lower:]')"
  if [[ "$fname" == *.bam || "$fname" == *long* ]]; then
    echo "ONT"
    return
  fi
  if [[ "$fname" == *short* || "$fname" == *illumina* || "$fname" == *.fastq || "$fname" == *.fq || "$fname" == *.fastq.gz || "$fname" == *.fq.gz ]]; then
    echo "Illumina"
    return
  fi
  echo "unknown"
}

count_fastq_reads_bases() {
  local run_file="$1"
  if [[ "$run_file" == *.gz ]]; then
    gzip -cd "$run_file" | awk 'NR % 4 == 2 {reads += 1; bases += length($0)} END {printf "%d\t%d\n", reads + 0, bases + 0}'
  else
    awk 'NR % 4 == 2 {reads += 1; bases += length($0)} END {printf "%d\t%d\n", reads + 0, bases + 0}' "$run_file"
  fi
}

count_bam_reads_bases() {
  local run_file="$1"
  singularity exec -B "$WORKDIR:$WORKDIR,$RUNS_DIR:$RUNS_DIR" "$CONTAINER" \
    samtools view "$run_file" | awk '{reads += 1; bases += length($10)} END {printf "%d\t%d\n", reads + 0, bases + 0}'
}

for sample_dir in "${sample_dirs[@]}"; do
  sample="$(basename "$sample_dir")"

  mapfile -t run_files < <(
    find "$sample_dir" -type f \
      \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" -o -name "*.bam" \) \
      | sort
  )

  if [[ ${#run_files[@]} -eq 0 ]]; then
    echo "Skipping $sample (no sequencing run files found)" >&2
    continue
  fi

  processed_samples=$((processed_samples + 1))

  for run_file in "${run_files[@]}"; do
    rel_path="${run_file#"$RUNS_DIR"/}"
    file_type="${run_file##*.}"
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
    platform="$(infer_platform "$run_file")"
    reads=0
    total_bases=0
    total_gbp="0.000000"
    status="ok"

    {
      echo "============================================================"
      echo "Timestamp: $timestamp"
      echo "Sample: $sample"
      echo "Input: $run_file"
      echo "Relative path: $rel_path"
      echo "Platform: $platform"
      echo "============================================================"

      if [[ "$run_file" == *.bam ]]; then
        bam_counts="$(count_bam_reads_bases "$run_file")"
        reads="${bam_counts%%$'\t'*}"
        total_bases="${bam_counts##*$'\t'}"
        total_gbp="$(awk -v b="$total_bases" 'BEGIN {printf "%.6f", b / 1000000000.0}')"
        echo "BAM detected; computed with samtools view (stats.sh does not support BAM)."
        echo "Read count: $reads"
        echo "Total bases: $total_bases"
        echo "Total sequence (Gbp): $total_gbp"
      else
        singularity exec -B "$WORKDIR:$WORKDIR,$RUNS_DIR:$RUNS_DIR" "$CONTAINER" \
          stats.sh in="$run_file"
        fq_counts="$(count_fastq_reads_bases "$run_file")"
        reads="${fq_counts%%$'\t'*}"
        total_bases="${fq_counts##*$'\t'}"
        total_gbp="$(awk -v b="$total_bases" 'BEGIN {printf "%.6f", b / 1000000000.0}')"
        echo "Computed from FASTQ sequences:"
        echo "Read count: $reads"
        echo "Total bases: $total_bases"
        echo "Total sequence (Gbp): $total_gbp"
      fi

      echo
    } >> "$LOG_FILE" 2>&1 || status="error"

    if [[ "$status" != "ok" ]]; then
      reads=0
      total_bases=0
      total_gbp="0.000000"
    fi

    echo -e "$timestamp\t$sample\t$platform\t$run_file\t$rel_path\t$file_type\t$reads\t$total_bases\t$total_gbp\t$status" >> "$CLEAN_TABLE_FILE"

    processed_files=$((processed_files + 1))
  done
done

if [[ $processed_files -eq 0 ]]; then
  echo "No sequencing run files processed under: $RUNS_DIR" >&2
  exit 1
fi

echo "Appended stats for $processed_files sequencing files across $processed_samples samples to $LOG_FILE"
echo "Wrote clean per-run summary table to $CLEAN_TABLE_FILE"

#!/bin/bash
#SBATCH --job-name=filter_subsample_contigs
#SBATCH --output=logs/filter_subsample_contigs_%j.out
#SBATCH --error=logs/filter_subsample_contigs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=math-alderaan

# Usage:
#   sbatch filter_subsampling_contigs_min_length.sh [SUBSAMPLING_ROOT]
# Example:
#   sbatch filter_subsampling_contigs_min_length.sh /storage/biology/projects/miller-lowry/beitner/data/subsampling

set -euo pipefail

echo "Starting filter_subsampling_contigs_min_length at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/assembly_stats"
SUBSAMPLING_ROOT="${1:-/storage/biology/projects/miller-lowry/beitner/data/subsampling}"
CONTAINER_PATH=""

container_candidates=(
  "${QC_TOOLS_CONTAINER:-}"
  "$WORKDIR/containers/qc_tools_miniconda.sif"
  "$WORKDIR/../containers/qc_tools_miniconda.sif"
  "$WORKDIR/../Lowry-assemblies/containers/qc_tools_miniconda.sif"
  "/storage/biology/projects/miller-lowry/beitner/Lowry-assemblies/containers/qc_tools_miniconda.sif"
)

for candidate in "${container_candidates[@]}"; do
  if [[ -n "$candidate" && -f "$candidate" ]]; then
    CONTAINER_PATH="$candidate"
    break
  fi
done

if [[ -z "$CONTAINER_PATH" ]]; then
  echo "Container not found. Set QC_TOOLS_CONTAINER or place qc_tools_miniconda.sif in a standard containers path." >&2
  exit 1
fi

if [[ ! -d "$SUBSAMPLING_ROOT" ]]; then
  echo "Subsampling root not found: $SUBSAMPLING_ROOT" >&2
  exit 1
fi

mapfile -t assembly_files < <(
  find "$SUBSAMPLING_ROOT" -type f -path '*/S1_subsample_*/assembly.*/*assembly.fasta' | sort -u
)

if [[ ${#assembly_files[@]} -eq 0 ]]; then
  echo "No subsampling assembly FASTA files found under: $SUBSAMPLING_ROOT" >&2
  exit 1
fi

processed_files=0

for assembly_file in "${assembly_files[@]}"; do
  assembly_dir="$(dirname "$assembly_file")"
  output_filtered="${assembly_dir}/contigs.ge1000.fa"

  if [[ -s "$output_filtered" ]]; then
    echo "Skipping $(basename "$(dirname "$assembly_dir")")/$(basename "$assembly_dir"): output already exists"
    continue
  fi

  echo ""
  echo "=== Filtering $assembly_file at $(date) ==="
  echo "Output: $output_filtered"

  singularity exec -B "$SUBSAMPLING_ROOT:$SUBSAMPLING_ROOT" -B "$WORKDIR:$WORKDIR" "$CONTAINER_PATH" reformat.sh \
    in="$assembly_file" \
    out="$output_filtered" \
    minlength=1000

  processed_files=$((processed_files + 1))
done

echo "Finished filtering $processed_files subsampling FASTA files to minlength=1000"
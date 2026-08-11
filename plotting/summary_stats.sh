#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSEMBLIES_DIR="${1:-$SCRIPT_DIR/../data/assemblies}"
CONTAINER=""
LOG_FILE="summary_stats_log.txt"
WORKDIR="$(pwd)"

container_candidates=(
  "${QC_TOOLS_CONTAINER:-}"
  "$SCRIPT_DIR/containers/qc_tools_miniconda.sif"
  "$SCRIPT_DIR/../containers/qc_tools_miniconda.sif"
  "$SCRIPT_DIR/../Lowry-assemblies/containers/qc_tools_miniconda.sif"
)

for candidate in "${container_candidates[@]}"; do
  if [[ -n "$candidate" && -f "$candidate" ]]; then
    CONTAINER="$candidate"
    break
  fi
done

pick_final_assembly() {
  local assembly_dir="$1"
  local assembler="$2"
  local candidates=()

  case "$assembler" in
    flye)
      candidates=("$assembly_dir/flye.assembly.fasta" "$assembly_dir/assembly.fasta")
      ;;
    megahit)
      candidates=("$assembly_dir/megahit.final.contigs.fa" "$assembly_dir/final.contigs.fa")
      ;;
    metamdbg)
      candidates=("$assembly_dir/metamdbg.contigs.fasta" "$assembly_dir/contigs.fasta")
      ;;
    metaspades|metaspades_hybrid|spades)
      candidates=("$assembly_dir/contigs.fasta")
      ;;
    idbaud)
      candidates=("$assembly_dir/assembly.fasta")
      ;;
    *)
      candidates=(
        "$assembly_dir/assembly.fasta"
        "$assembly_dir/contigs.fasta"
        "$assembly_dir/final.contigs.fa"
      )
      ;;
  esac

  for candidate in "${candidates[@]}"; do
    if [[ -f "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done

  return 1
}

if [[ ! -f "$CONTAINER" ]]; then
  echo "Container not found. Set QC_TOOLS_CONTAINER or place qc_tools_miniconda.sif in a standard containers path." >&2
  exit 1
fi

if [[ ! -d "$ASSEMBLIES_DIR" ]]; then
  echo "Assemblies directory not found: $ASSEMBLIES_DIR" >&2
  exit 1
fi

allowed_assemblers=(
  flye
  idbaud
  megahit
  metaconnet
  metamdbg
  metaspades
  metaspades_hybrid
  opera-ms
)

shopt -s nullglob
assemblies=()
samples=()
assemblers=()
declare -A seen_sample_assembler=()

infer_sample_from_flat_file() {
  local path="$1"
  local base sample_part

  base="$(basename "$path")"
  sample_part="${base%%.assembly.*}"

  # Normalize s1/S1 -> S1
  if [[ "$sample_part" =~ ^[sS]([0-9]+)$ ]]; then
    echo "S${BASH_REMATCH[1]}"
    return 0
  fi

  return 1
}

# Only scan assembler-named top-level folders (ignore S# directories entirely).
for assembler in "${allowed_assemblers[@]}"; do
  assembler_root="$ASSEMBLIES_DIR/$assembler"
  [[ -d "$assembler_root" ]] || continue

  for assembly in \
    "$assembler_root"/s*.assembly.fasta \
    "$assembler_root"/s*.assembly.fa \
    "$assembler_root"/S*.assembly.fasta \
    "$assembler_root"/S*.assembly.fa; do
    [[ -f "$assembly" ]] || continue

    if ! sample="$(infer_sample_from_flat_file "$assembly")"; then
      echo "Skipping $assembly (unable to infer sample name from filename)" >&2
      continue
    fi

    key="${sample}|${assembler}"
    if [[ -n "${seen_sample_assembler[$key]+x}" ]]; then
      continue
    fi

    assemblies+=("$assembly")
    samples+=("$sample")
    assemblers+=("$assembler")
    seen_sample_assembler["$key"]=1
  done
done
shopt -u nullglob

if [[ ${#assemblies[@]} -eq 0 ]]; then
  echo "No canonical assembly FASTA files found under: $ASSEMBLIES_DIR" >&2
  exit 1
fi

# Recreate the log from scratch on each run.
: > "$LOG_FILE"

for idx in "${!assemblies[@]}"; do
  assembly="${assemblies[$idx]}"
  sample="${samples[$idx]}"
  assembler="${assemblers[$idx]}"
  rel_path="${assembly#"$ASSEMBLIES_DIR"/}"

  {
    echo "============================================================"
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Sample: $sample"
    echo "Assembler: $assembler"
    echo "Assembly: $assembly"
    echo "Relative path: $rel_path"
    echo "============================================================"

    singularity exec -B "$WORKDIR:$WORKDIR,$ASSEMBLIES_DIR:$ASSEMBLIES_DIR" "$CONTAINER" \
      stats.sh in="$assembly"

    echo
  } >> "$LOG_FILE" 2>&1

done

echo "Appended stats for ${#assemblies[@]} assemblies to $LOG_FILE"
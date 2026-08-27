#!/usr/bin/env bash
#SBATCH --job-name=metawrap_binner_set
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/metawrap_binner_set_%j.out
#SBATCH --error=logs/metawrap_binner_set_%j.err
#SBATCH --chdir=/storage/biology/projects/miller-lowry/beitner/binning-classification-visualization

set -euo pipefail

# Use absolute paths to work correctly with sbatch
PROJECT_ROOT="/storage/biology/projects/miller-lowry/beitner"
BINNING_DIR="${PROJECT_ROOT}/binning-classification-visualization/binning"
CONFIG_DIR="${PROJECT_ROOT}/binning-classification-visualization/configs"

if [[ -f "$CONFIG_DIR/configs_master.conf" ]]; then
    # shellcheck source=/dev/null
    source "$CONFIG_DIR/configs_master.conf"
fi

BASE_DIR="${PROJECT_BASE_DIR:-$PROJECT_ROOT}"
WORK_DIR="${PROJECT_WORK_DIR:-${PROJECT_ROOT}/binning-classification-visualization}"

ASSEMBLY_NAME=""
SAMPLE=""
SET_NAME=""
ASSEMBLY_FASTA=""
OUT_ROOT=""
METAWRAP_CONTAINER="${WORK_DIR}/containers/metawrap.sif"
CPUS="${SLURM_CPUS_PER_TASK:-32}"
MEM_GB=360

SHORT_READS=()
LONG_READS=()

usage() {
    cat <<EOF
Usage:
  bash scripts/metawrap_binner_single_set.sh \
    --assembly-name NAME --sample NAME --set-name NAME \
    --assembly-fasta PATH --out-root PATH \
    [--short-read PATH ...] [--long-read PATH ...] [--container PATH] [--cpus N] [--mem-gb N]

Notes:
  - All reads are staged as single-end FASTQ for metawrap binning --single-end.
  - Provide at least one --short-read or --long-read.
EOF
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --assembly-name)
            ASSEMBLY_NAME="$2"
            shift 2
            ;;
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --set-name)
            SET_NAME="$2"
            shift 2
            ;;
        --assembly-fasta)
            ASSEMBLY_FASTA="$2"
            shift 2
            ;;
        --out-root)
            OUT_ROOT="$2"
            shift 2
            ;;
        --short-read)
            SHORT_READS+=("$2")
            shift 2
            ;;
        --long-read)
            LONG_READS+=("$2")
            shift 2
            ;;
        --container)
            METAWRAP_CONTAINER="$2"
            shift 2
            ;;
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --mem-gb)
            MEM_GB="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$ASSEMBLY_NAME" || -z "$SAMPLE" || -z "$SET_NAME" || -z "$ASSEMBLY_FASTA" || -z "$OUT_ROOT" ]]; then
    echo "ERROR: missing required arguments" >&2
    usage >&2
    exit 1
fi

if [[ ! -f "$ASSEMBLY_FASTA" ]]; then
    echo "ERROR: assembly fasta not found: $ASSEMBLY_FASTA" >&2
    exit 1
fi
if [[ ! -f "$METAWRAP_CONTAINER" ]]; then
    echo "ERROR: metawrap container not found: $METAWRAP_CONTAINER" >&2
    exit 1
fi
if [[ "${#SHORT_READS[@]}" -eq 0 && "${#LONG_READS[@]}" -eq 0 ]]; then
    echo "ERROR: no reads provided" >&2
    exit 1
fi

for read in "${SHORT_READS[@]}" "${LONG_READS[@]}"; do
    if [[ -n "$read" && ! -f "$read" ]]; then
        echo "ERROR: read file not found: $read" >&2
        exit 1
    fi
done

mkdir -p "$OUT_ROOT" "$WORK_DIR/slurm_logs"
RUN_DIR="${OUT_ROOT}/${SET_NAME}"
STAGE_ROOT="${OUT_ROOT}/staging"
LOG_DIR="${OUT_ROOT}/logs"
mkdir -p "$RUN_DIR" "$STAGE_ROOT" "$LOG_DIR"

BIN_FASTA_DIR="${RUN_DIR}/bin_fasta"

STAGE_DIR="$(mktemp -d -p "$STAGE_ROOT" "${SET_NAME}_metawrap_XXXXXX")"
trap 'rm -rf "$STAGE_DIR"' EXIT

stage_read() {
    local src="$1"
    local prefix="$2"
    local base
    local dest
    base="$(basename "$src")"
    if [[ "$base" == *.fastq.gz ]]; then
        dest="${STAGE_DIR}/${prefix}_${base%.gz}"
        gzip -dc "$src" > "$dest"
    elif [[ "$base" == *.fq.gz ]]; then
        dest="${STAGE_DIR}/${prefix}_${base%.gz}"
        gzip -dc "$src" > "$dest"
    elif [[ "$base" == *.fastq || "$base" == *.fq ]]; then
        dest="${STAGE_DIR}/${prefix}_${base}"
        cp "$src" "$dest"
    else
        dest="${STAGE_DIR}/${prefix}_${base}.fastq"
        cp "$src" "$dest"
    fi
    printf '%s\n' "$dest"
}

STAGED_READS=()
for read in "${SHORT_READS[@]}"; do
    STAGED_READS+=("$(stage_read "$read" "short")")
done
for read in "${LONG_READS[@]}"; do
    STAGED_READS+=("$(stage_read "$read" "long")")
done

if [[ "${#STAGED_READS[@]}" -eq 0 ]]; then
    echo "ERROR: no staged reads created" >&2
    exit 1
fi

rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"

echo "[$(date)] Starting metaWRAP binning for $SET_NAME"

singularity exec "$METAWRAP_CONTAINER" metawrap binning \
    -a "$ASSEMBLY_FASTA" \
    -o "$RUN_DIR" \
    -t "$CPUS" \
    -m "$MEM_GB" \
    --single-end \
    --metabat2 --maxbin2 --concoct \
    "${STAGED_READS[@]}"

rm -rf "$BIN_FASTA_DIR"
mkdir -p "$BIN_FASTA_DIR"

copy_bins_from_dir() {
    local src_dir="$1"
    local prefix="$2"
    [[ -d "$src_dir" ]] || return 0

    local f
    local n=0
    for f in "$src_dir"/*.fa "$src_dir"/*.fna "$src_dir"/*.fasta; do
        [[ -f "$f" ]] || continue
        n=$((n + 1))
        cp "$f" "${BIN_FASTA_DIR}/${prefix}_${n}.fa"
    done
}

copy_bins_from_dir "${RUN_DIR}/metabat2_bins" "metabat2"
copy_bins_from_dir "${RUN_DIR}/maxbin2_bins" "maxbin2"
copy_bins_from_dir "${RUN_DIR}/concoct_bins" "concoct"

BIN_FASTA_COUNT=$(find "$BIN_FASTA_DIR" -maxdepth 1 -type f -name '*.fa' | wc -l)
if [[ "$BIN_FASTA_COUNT" -eq 0 ]]; then
    echo "ERROR: no FASTA bins were produced for standardized output in ${BIN_FASTA_DIR}" >&2
    exit 1
fi

SUMMARY_TSV="${OUT_ROOT}/metawrap_binner_summary.tsv"
if [[ ! -f "$SUMMARY_TSV" ]]; then
    printf 'assembly\tsample\tset_name\tassembly_fasta\tshort_reads\tlong_reads\trun_dir\tbin_fasta_dir\tbin_fasta_bins\tstatus\n' > "$SUMMARY_TSV"
fi

TMP_SUMMARY="${SUMMARY_TSV}.tmp.$$"
awk -F'\t' -v s="$SET_NAME" 'NR==1 || $3 != s' "$SUMMARY_TSV" > "$TMP_SUMMARY"
mv "$TMP_SUMMARY" "$SUMMARY_TSV"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$ASSEMBLY_NAME" "$SAMPLE" "$SET_NAME" "$ASSEMBLY_FASTA" \
    "${#SHORT_READS[@]}" "${#LONG_READS[@]}" "$RUN_DIR" "$BIN_FASTA_DIR" "$BIN_FASTA_COUNT" "completed" >> "$SUMMARY_TSV"

echo "[$(date)] metaWRAP complete for $SET_NAME"
echo "[$(date)] Output: $RUN_DIR"
echo "[$(date)] Standardized bin FASTA dir: $BIN_FASTA_DIR (n=${BIN_FASTA_COUNT})"

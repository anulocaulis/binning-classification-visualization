#!/usr/bin/env bash
#SBATCH --job-name=vamb_binner_set
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/vamb_binner_set_%j.out
#SBATCH --error=logs/vamb_binner_set_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/../configs"
if [[ -f "$CONFIG_DIR/configs_master.conf" ]]; then
    # shellcheck source=/dev/null
    source "$CONFIG_DIR/configs_master.conf"
fi

BASE_DIR="${PROJECT_BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
WORK_DIR="${PROJECT_WORK_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"

ASSEMBLY_NAME=""
SAMPLE=""
SET_NAME=""
ASSEMBLY_FASTA=""
OUT_ROOT=""
QC_CONTAINER="${WORK_DIR}/containers/qc_tools_miniconda.sif"
VAMB_CONTAINER="${WORK_DIR}/containers/vamb.sif"
CPUS="${SLURM_CPUS_PER_TASK:-32}"
MIN_CONTIG_LEN=2000
EXTRACT_MIN_BP=100000

SHORT_READS=()
LONG_READS=()

usage() {
    cat <<EOF
Usage:
  bash scripts/vamb_binner_single_set.sh \
    --assembly-name NAME --sample NAME --set-name NAME \
    --assembly-fasta PATH --out-root PATH \
        [--short-read PATH ...] [--long-read PATH ...] [--qc-container PATH] [--vamb-container PATH] [--extract-min-bp N]
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
        --qc-container)
            QC_CONTAINER="$2"
            shift 2
            ;;
        --vamb-container)
            VAMB_CONTAINER="$2"
            shift 2
            ;;
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --min-contig-len)
            MIN_CONTIG_LEN="$2"
            shift 2
            ;;
        --extract-min-bp)
            EXTRACT_MIN_BP="$2"
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
if [[ ! -f "$QC_CONTAINER" ]]; then
    echo "ERROR: qc container not found: $QC_CONTAINER" >&2
    exit 1
fi
if [[ ! -f "$VAMB_CONTAINER" ]]; then
    echo "ERROR: vamb container not found: $VAMB_CONTAINER" >&2
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

ROOT="${OUT_ROOT}/${SET_NAME}"
MAP_DIR="${ROOT}/mapping"
OUT_DIR="${ROOT}/vamb_bins"
ABUNDANCE_TSV="${ROOT}/combined.abundance.tsv"
BIN_FASTA_DIR="${ROOT}/bin_fasta"
mkdir -p "$MAP_DIR" "$WORK_DIR/slurm_logs"

BAMS=()
LABELS=()

if [[ "${#SHORT_READS[@]}" -gt 0 ]]; then
    singularity exec "$QC_CONTAINER" bwa index "$ASSEMBLY_FASTA"
    idx=1
    for read in "${SHORT_READS[@]}"; do
        bam="${MAP_DIR}/short_${idx}.sorted.bam"
        singularity exec "$QC_CONTAINER" bash -o pipefail -c \
            "bwa mem -t ${CPUS} '$ASSEMBLY_FASTA' '$read' | samtools sort -@ 8 -o '$bam'"
        singularity exec "$QC_CONTAINER" samtools index "$bam"
        BAMS+=("$bam")
        LABELS+=("short_${idx}")
        idx=$((idx + 1))
    done
fi

if [[ "${#LONG_READS[@]}" -gt 0 ]]; then
    idx=1
    for read in "${LONG_READS[@]}"; do
        bam="${MAP_DIR}/long_${idx}.sorted.bam"
        singularity exec "$QC_CONTAINER" bash -o pipefail -c \
            "minimap2 -ax map-ont -t ${CPUS} '$ASSEMBLY_FASTA' '$read' | samtools sort -@ 8 -o '$bam'"
        singularity exec "$QC_CONTAINER" samtools index "$bam"
        BAMS+=("$bam")
        LABELS+=("long_${idx}")
        idx=$((idx + 1))
    done
fi

if [[ "${#BAMS[@]}" -eq 0 ]]; then
    echo "ERROR: no BAM files produced" >&2
    exit 1
fi

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
for i in "${!BAMS[@]}"; do
    singularity exec "$QC_CONTAINER" samtools idxstats "${BAMS[$i]}" > "${TMPDIR}/${LABELS[$i]}.idxstats.tsv"
done

LABELS_JOINED="$(IFS=,; echo "${LABELS[*]}")"
python3 - "$ASSEMBLY_FASTA" "$ABUNDANCE_TSV" "$TMPDIR" "$LABELS_JOINED" <<'PY'
import csv
import os
import sys

fasta_path, out_tsv, idx_dir, labels_joined = sys.argv[1:5]
labels = [x for x in labels_joined.split(',') if x]

contigs = []
with open(fasta_path, 'r', encoding='utf-8', errors='replace') as handle:
    for line in handle:
        if line.startswith('>'):
            contigs.append(line[1:].split()[0])

sample_maps = []
for label in labels:
    path = os.path.join(idx_dir, f'{label}.idxstats.tsv')
    abundance = {}
    with open(path, newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if len(row) < 4 or row[0] == '*':
                continue
            contig = row[0]
            try:
                length = float(row[1])
                mapped = float(row[2])
            except ValueError:
                continue
            abundance[contig] = '0' if length <= 0 else f'{mapped/length:.12f}'
    sample_maps.append(abundance)

with open(out_tsv, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['contigname', *labels])
    for contig in contigs:
        writer.writerow([contig, *[m.get(contig, '0') for m in sample_maps]])
PY

rm -rf "$OUT_DIR"
singularity exec "$VAMB_CONTAINER" vamb bin default \
    --outdir "$OUT_DIR" \
    --fasta "$ASSEMBLY_FASTA" \
    --abundance_tsv "$ABUNDANCE_TSV" \
    -m "$MIN_CONTIG_LEN" \
    -p "$CPUS"

if [[ -s "${OUT_DIR}/vae_clusters_unsplit.tsv" && ! -s "${OUT_DIR}/clusters.tsv" ]]; then
    cp "${OUT_DIR}/vae_clusters_unsplit.tsv" "${OUT_DIR}/clusters.tsv"
fi

CLUSTERS_TSV="${OUT_DIR}/clusters.tsv"
if [[ ! -s "$CLUSTERS_TSV" && -s "${OUT_DIR}/vae_clusters_unsplit.tsv" ]]; then
    CLUSTERS_TSV="${OUT_DIR}/vae_clusters_unsplit.tsv"
fi
if [[ ! -s "$CLUSTERS_TSV" ]]; then
    echo "ERROR: no VAMB clusters file found in ${OUT_DIR}" >&2
    exit 1
fi

extract_vamb_bins() {
    local min_bp="$1"
    python3 - "$ASSEMBLY_FASTA" "$CLUSTERS_TSV" "$BIN_FASTA_DIR" "$min_bp" <<'PY'
import os
import sys
from collections import defaultdict

assembly_fasta, clusters_tsv, outdir, min_bp = sys.argv[1:5]
min_bp = int(min_bp)

if not os.path.exists(assembly_fasta):
    raise SystemExit(f"Missing assembly FASTA: {assembly_fasta}")
if not os.path.exists(clusters_tsv):
    raise SystemExit(f"Missing clusters TSV: {clusters_tsv}")

contig_to_cluster = {}
with open(clusters_tsv, "r", encoding="utf-8", errors="replace") as fh:
    header = None
    for line in fh:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue

        if header is None and ("contigname" in parts or "clustername" in parts):
            header = [p.strip().lower() for p in parts]
            continue

        if header and "contigname" in header and "clustername" in header:
            cidx = header.index("contigname")
            zidx = header.index("clustername")
            if cidx >= len(parts) or zidx >= len(parts):
                continue
            contig = parts[cidx]
            cluster = parts[zidx]
        else:
            contig, cluster = parts[0].strip(), parts[1].strip()

        if contig and cluster:
            contig_to_cluster[contig] = cluster

clusters = defaultdict(list)
for contig, cluster in contig_to_cluster.items():
    clusters[cluster].append(contig)

seqs = {}
current = None
chunks = []
with open(assembly_fasta, "r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith(">"):
            if current is not None:
                seqs[current] = "".join(chunks)
            current = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
if current is not None:
    seqs[current] = "".join(chunks)

os.makedirs(outdir, exist_ok=True)
for existing in os.listdir(outdir):
    if existing.endswith('.fa'):
        os.remove(os.path.join(outdir, existing))

written = 0
for cluster in sorted(clusters):
    total_len = 0
    records = []
    for contig in clusters[cluster]:
        seq = seqs.get(contig)
        if seq is None:
            continue
        total_len += len(seq)
        records.append((contig, seq))

    if total_len < min_bp or not records:
        continue

    written += 1
    out_path = os.path.join(outdir, f"Bin_{written}.fa")
    with open(out_path, "w", encoding="utf-8") as out:
        for contig, seq in records:
            out.write(f">{contig}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

print(written)
PY
}

mkdir -p "$BIN_FASTA_DIR"
BIN_FASTA_COUNT="$(extract_vamb_bins "$EXTRACT_MIN_BP")"
if [[ "$BIN_FASTA_COUNT" -eq 0 ]]; then
    for retry_bp in 50000 10000 1000; do
        BIN_FASTA_COUNT="$(extract_vamb_bins "$retry_bp")"
        if [[ "$BIN_FASTA_COUNT" -gt 0 ]]; then
            EXTRACT_MIN_BP="$retry_bp"
            break
        fi
    done
fi

if [[ "$BIN_FASTA_COUNT" -eq 0 ]]; then
    echo "ERROR: VAMB extraction produced zero FASTA bins in ${BIN_FASTA_DIR}" >&2
    exit 1
fi

SUMMARY_TSV="${OUT_ROOT}/vamb_binner_summary.tsv"
if [[ ! -f "$SUMMARY_TSV" ]]; then
    printf 'assembly\tsample\tset_name\tassembly_fasta\tshort_reads\tlong_reads\trun_dir\tbin_fasta_dir\tbin_fasta_bins\textract_min_bp\tstatus\n' > "$SUMMARY_TSV"
fi

TMP_SUMMARY="${SUMMARY_TSV}.tmp.$$"
awk -F'\t' -v s="$SET_NAME" 'NR==1 || $3 != s' "$SUMMARY_TSV" > "$TMP_SUMMARY"
mv "$TMP_SUMMARY" "$SUMMARY_TSV"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$ASSEMBLY_NAME" "$SAMPLE" "$SET_NAME" "$ASSEMBLY_FASTA" \
    "${#SHORT_READS[@]}" "${#LONG_READS[@]}" "$ROOT" "$BIN_FASTA_DIR" "$BIN_FASTA_COUNT" "$EXTRACT_MIN_BP" "completed" >> "$SUMMARY_TSV"

echo "[$(date)] VAMB complete for $SET_NAME"
echo "[$(date)] Output: $ROOT"
echo "[$(date)] Standardized bin FASTA dir: $BIN_FASTA_DIR (n=${BIN_FASTA_COUNT}, min_bp=${EXTRACT_MIN_BP})"

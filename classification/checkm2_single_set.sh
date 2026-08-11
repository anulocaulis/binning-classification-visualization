#!/usr/bin/env bash
#SBATCH --job-name=checkm2_single_set
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/checkm2_single_set_%j.out
#SBATCH --error=logs/checkm2_single_set_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/../configs"
if [[ -f "$CONFIG_DIR/configs_master.conf" ]]; then
    # shellcheck source=/dev/null
    source "$CONFIG_DIR/configs_master.conf"
fi

BASE_DIR="${PROJECT_BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
WORK_DIR="${PROJECT_WORK_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
DATA_DIR="${PROJECT_DATA_DIR:-${BASE_DIR}/data}"

ASSEMBLY=""
SAMPLE=""
SET_NAME=""
SET_SUFFIX="metawrap_bins_refined"
BIN_DIR=""

OUT_ROOT="${DATA_DIR}/binning_outputs/checkm2_single_runs"
CONTAINER="${WORK_DIR}/containers/checkm2-fastani-drep.sif"
CHECKM2_DB="${BASE_DIR}/databases/CheckM2_database/uniref100.KO.1.dmnd"

HQ_COMPLETENESS=90
HQ_CONTAMINATION=5
MQ_COMPLETENESS=50
MQ_CONTAMINATION=10

usage() {
    cat <<EOF
Usage:
  bash scripts/checkm2_single_set.sh --assembly <ASSEMBLY> --sample <SAMPLE> --bin-dir <BIN_DIR> [options]

Required:
  --assembly NAME          Assembly label (e.g., flye, idbaud)
  --sample NAME            Sample label (e.g., S1)
  --bin-dir PATH           Directory containing bin FASTA files

Optional:
  --set-name NAME          Explicit set name (default: <assembly>_<sample_lower>_<set_suffix>)
  --set-suffix NAME        Suffix used for default set-name (default: metawrap_bins_refined)
  --out-root PATH          Output root (default: data/binning_outputs/checkm2_single_runs)
  --container PATH         CheckM2 container (default: containers/checkm2-fastani-drep.sif)
  --checkm2-db PATH        CheckM2 .dmnd database path
  --hq-completeness FLOAT  HQ completeness cutoff (default: 90)
  --hq-contamination FLOAT HQ contamination cutoff (default: 5)
  --mq-completeness FLOAT  MQ completeness cutoff (default: 50)
  --mq-contamination FLOAT MQ contamination cutoff (default: 10)
  -h, --help               Show help
EOF
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --assembly)
            ASSEMBLY="$2"
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
        --set-suffix)
            SET_SUFFIX="$2"
            shift 2
            ;;
        --bin-dir)
            BIN_DIR="$2"
            shift 2
            ;;
        --out-root)
            OUT_ROOT="$2"
            shift 2
            ;;
        --container)
            CONTAINER="$2"
            shift 2
            ;;
        --checkm2-db)
            CHECKM2_DB="$2"
            shift 2
            ;;
        --hq-completeness)
            HQ_COMPLETENESS="$2"
            shift 2
            ;;
        --hq-contamination)
            HQ_CONTAMINATION="$2"
            shift 2
            ;;
        --mq-completeness)
            MQ_COMPLETENESS="$2"
            shift 2
            ;;
        --mq-contamination)
            MQ_CONTAMINATION="$2"
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

if [[ -z "$ASSEMBLY" || -z "$SAMPLE" || -z "$BIN_DIR" ]]; then
    echo "ERROR: --assembly, --sample, and --bin-dir are required" >&2
    usage >&2
    exit 1
fi

sample_lower="$(printf '%s' "$SAMPLE" | tr '[:upper:]' '[:lower:]')"
if [[ -z "$SET_NAME" ]]; then
    SET_NAME="${ASSEMBLY}_${sample_lower}_${SET_SUFFIX}"
fi

if [[ ! -d "$BIN_DIR" ]]; then
    echo "ERROR: bin directory not found: $BIN_DIR" >&2
    exit 1
fi
if [[ ! -f "$CONTAINER" ]]; then
    echo "ERROR: container not found: $CONTAINER" >&2
    exit 1
fi
if [[ ! -f "$CHECKM2_DB" ]]; then
    echo "ERROR: CheckM2 database not found: $CHECKM2_DB" >&2
    exit 1
fi

mkdir -p "$OUT_ROOT" "$OUT_ROOT/runs" "$OUT_ROOT/summary" "$WORK_DIR/slurm_logs"

infer_ext() {
    local dir="$1"
    if find "$dir" -maxdepth 1 -type f -name '*.fa' | grep -q .; then
        echo "fa"
    elif find "$dir" -maxdepth 1 -type f -name '*.fna' | grep -q .; then
        echo "fna"
    elif find "$dir" -maxdepth 1 -type f -name '*.fasta' | grep -q .; then
        echo "fasta"
    else
        echo ""
    fi
}

BIN_EXT="$(infer_ext "$BIN_DIR")"
if [[ -z "$BIN_EXT" ]]; then
    echo "ERROR: no .fa/.fna/.fasta files found in: $BIN_DIR" >&2
    exit 1
fi

RUN_OUT="${OUT_ROOT}/runs/${SET_NAME}"
rm -rf "$RUN_OUT"
mkdir -p "$RUN_OUT"

echo "[$(date)] Running CheckM2 for ${SET_NAME}"
echo "[$(date)] BIN_DIR=$BIN_DIR"
echo "[$(date)] RUN_OUT=$RUN_OUT"

INPUT_DIR="$BIN_DIR" \
OUT_DIR="$RUN_OUT" \
THREADS="${SLURM_CPUS_PER_TASK:-32}" \
DB_PATH="$CHECKM2_DB" \
BIN_EXT="$BIN_EXT" \
/usr/local/bin/singularity exec "$CONTAINER" bash -lc '
set -euo pipefail

if command -v checkm2 >/dev/null 2>&1; then
    CHECKM2_CMD="checkm2"
else
    CHECKM2_CMD=""
    for c in /usr/local/miniconda3/envs/checkm2/bin/checkm2 /opt/conda/envs/checkm2/bin/checkm2 /opt/conda/bin/checkm2; do
        if [[ -x "$c" ]]; then
            CHECKM2_CMD="$c"
            break
        fi
    done
fi

if [[ -z "$CHECKM2_CMD" ]]; then
    echo "ERROR: checkm2 not found in container" >&2
    exit 1
fi

"$CHECKM2_CMD" predict \
    --input "$INPUT_DIR" \
    --output-directory "$OUT_DIR" \
    --threads "$THREADS" \
    --general \
    --extension "$BIN_EXT" \
    --database_path "$DB_PATH"
'

QUALITY_REPORT="${RUN_OUT}/quality_report.tsv"
if [[ ! -f "$QUALITY_REPORT" ]]; then
    echo "ERROR: quality_report.tsv missing after CheckM2 run: $QUALITY_REPORT" >&2
    exit 1
fi

COUNTS="$(python3 - "$QUALITY_REPORT" "$HQ_COMPLETENESS" "$HQ_CONTAMINATION" "$MQ_COMPLETENESS" "$MQ_CONTAMINATION" <<'PY'
import csv
import sys

report, hq_c, hq_x, mq_c, mq_x = sys.argv[1:6]
hq_c = float(hq_c)
hq_x = float(hq_x)
mq_c = float(mq_c)
mq_x = float(mq_x)

def pick(fields, names):
    for n in names:
        if n in fields:
            return n
    return None

total = hq = mq = 0
with open(report, newline="", encoding="utf-8", errors="replace") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []
    comp_col = pick(fields, ["Completeness", "Completeness_General", "Completeness_Specific"])
    cont_col = pick(fields, ["Contamination"])
    for row in reader:
        total += 1
        try:
            comp = float(row.get(comp_col, 0) or 0) if comp_col else 0.0
        except Exception:
            comp = 0.0
        try:
            cont = float(row.get(cont_col, 0) or 0) if cont_col else 0.0
        except Exception:
            cont = 0.0

        if comp >= hq_c and cont <= hq_x:
            hq += 1
        elif comp >= mq_c and cont <= mq_x:
            mq += 1

lq = total - hq - mq
print(f"{total}\t{hq}\t{mq}\t{lq}")
PY
)"

TOTAL_BINS="$(printf '%s' "$COUNTS" | cut -f1)"
HQ_BINS="$(printf '%s' "$COUNTS" | cut -f2)"
MQ_BINS="$(printf '%s' "$COUNTS" | cut -f3)"
LQ_BINS="$(printf '%s' "$COUNTS" | cut -f4)"

SUMMARY_TSV="${OUT_ROOT}/summary/checkm2_single_summary.tsv"
if [[ ! -f "$SUMMARY_TSV" ]]; then
    printf "assembly\tsample\tset_name\tbin_dir\tquality_report\ttotal_bins\thq_bins\tmq_bins\tlq_bins\thq_completeness\thq_contamination\tmq_completeness\tmq_contamination\tstatus\tnotes\n" > "$SUMMARY_TSV"
fi

TMP_SUMMARY="${SUMMARY_TSV}.tmp.$$"
awk -F'\t' -v set_name="$SET_NAME" 'NR==1 || $3 != set_name' "$SUMMARY_TSV" > "$TMP_SUMMARY"
mv "$TMP_SUMMARY" "$SUMMARY_TSV"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$ASSEMBLY" "$SAMPLE" "$SET_NAME" "$BIN_DIR" "$QUALITY_REPORT" \
    "$TOTAL_BINS" "$HQ_BINS" "$MQ_BINS" "$LQ_BINS" \
    "$HQ_COMPLETENESS" "$HQ_CONTAMINATION" "$MQ_COMPLETENESS" "$MQ_CONTAMINATION" \
    "completed" "single-set-checkm2" >> "$SUMMARY_TSV"

echo "[$(date)] CheckM2 complete: ${SET_NAME}"
echo "[$(date)] Quality report: $QUALITY_REPORT"
echo "[$(date)] Summary TSV: $SUMMARY_TSV"

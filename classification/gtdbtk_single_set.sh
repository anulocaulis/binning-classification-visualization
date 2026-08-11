#!/usr/bin/env bash
#SBATCH --job-name=gtdbtk_single_set
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=4-00:00:00
#SBATCH --output=/storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper/slurm_logs/gtdbtk_single_set_%j.out
#SBATCH --error=/storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper/slurm_logs/gtdbtk_single_set_%j.err

set -euo pipefail

BASE_DIR="/storage/biology/projects/miller-lowry/beitner"
WORK_DIR="${BASE_DIR}/binning-classification-wrapper"
DATA_DIR="${BASE_DIR}/data"

ASSEMBLY=""
SAMPLE=""
SET_NAME=""
SET_SUFFIX="metawrap_bins_refined"

BIN_DIR=""
CHECKM2_QUALITY_REPORT=""

OUT_ROOT="${DATA_DIR}/binning_outputs/gtdbtk_single_runs"
CONTAINER="${WORK_DIR}/containers/classification.sif"
GTDBTK_DB="${BASE_DIR}/databases/gtdbtk_db"

HQ_COMPLETENESS=90
HQ_CONTAMINATION=5
MQ_COMPLETENESS=50
MQ_CONTAMINATION=10
INCLUDE_MQ=1

usage() {
    cat <<EOF
Usage:
  bash scripts/gtdbtk_single_set.sh --assembly <ASSEMBLY> --sample <SAMPLE> --bin-dir <BIN_DIR> --checkm2-quality-report <QUALITY_REPORT> [options]

Required:
  --assembly NAME                 Assembly label (e.g., flye)
  --sample NAME                   Sample label (e.g., S1)
  --bin-dir PATH                  Directory containing source bin FASTA files
  --checkm2-quality-report PATH   CheckM2 quality_report.tsv for this set

Optional:
  --set-name NAME                 Explicit set name (default: <assembly>_<sample_lower>_<set_suffix>)
  --set-suffix NAME               Suffix used for default set-name (default: metawrap_bins_refined)
  --out-root PATH                 Output root (default: data/binning_outputs/gtdbtk_single_runs)
  --container PATH                GTDB-Tk container (default: containers/classification.sif)
  --gtdbtk-db PATH                GTDB-Tk database directory path
  --hq-completeness FLOAT         HQ completeness cutoff (default: 90)
  --hq-contamination FLOAT        HQ contamination cutoff (default: 5)
  --mq-completeness FLOAT         MQ completeness cutoff (default: 50)
  --mq-contamination FLOAT        MQ contamination cutoff (default: 10)
  --hq-only                       Select HQ bins only (default is HQ+MQ)
  -h, --help                      Show help
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
        --checkm2-quality-report)
            CHECKM2_QUALITY_REPORT="$2"
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
        --gtdbtk-db)
            GTDBTK_DB="$2"
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
        --hq-only)
            INCLUDE_MQ=0
            shift
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

if [[ -z "$ASSEMBLY" || -z "$SAMPLE" || -z "$BIN_DIR" || -z "$CHECKM2_QUALITY_REPORT" ]]; then
    echo "ERROR: --assembly, --sample, --bin-dir, and --checkm2-quality-report are required" >&2
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
if [[ ! -f "$CHECKM2_QUALITY_REPORT" ]]; then
    echo "ERROR: CheckM2 quality report not found: $CHECKM2_QUALITY_REPORT" >&2
    exit 1
fi
if [[ ! -f "$CONTAINER" ]]; then
    echo "ERROR: container not found: $CONTAINER" >&2
    exit 1
fi
if [[ ! -d "$GTDBTK_DB" ]]; then
    echo "ERROR: GTDB-Tk DB directory not found: $GTDBTK_DB" >&2
    exit 1
fi

mkdir -p "$OUT_ROOT" "$OUT_ROOT/runs" "$OUT_ROOT/staging" "$OUT_ROOT/summary" "$WORK_DIR/slurm_logs"

if ! /usr/local/bin/singularity exec "$CONTAINER" bash -lc "command -v gtdbtk >/dev/null 2>&1"; then
    echo "ERROR: gtdbtk not found in container: $CONTAINER" >&2
    exit 1
fi

STAGING_DIR="${OUT_ROOT}/staging/${SET_NAME}"
RUN_DIR="${OUT_ROOT}/runs/${SET_NAME}"
rm -rf "$STAGING_DIR" "$RUN_DIR"
mkdir -p "$STAGING_DIR" "$RUN_DIR"

SELECTION_COUNTS="$(python3 - "$CHECKM2_QUALITY_REPORT" "$BIN_DIR" "$STAGING_DIR" "$HQ_COMPLETENESS" "$HQ_CONTAMINATION" "$MQ_COMPLETENESS" "$MQ_CONTAMINATION" "$INCLUDE_MQ" <<'PY'
import csv
import glob
import os
import sys

report_path, bin_dir, staging_dir, hq_c, hq_x, mq_c, mq_x, include_mq = sys.argv[1:9]
hq_c = float(hq_c)
hq_x = float(hq_x)
mq_c = float(mq_c)
mq_x = float(mq_x)
include_mq = int(include_mq)

def normalize_name(value: str) -> str:
    return os.path.splitext((value or "").strip())[0]

def normalize_key(value: str) -> str:
    return "".join(ch.lower() for ch in normalize_name(value) if ch.isalnum())

def pick_column(fieldnames, candidates):
    if not fieldnames:
        return None
    norm_map = {normalize_key(name): name for name in fieldnames}
    candidate_norms = [normalize_key(c) for c in candidates]
    for cand in candidate_norms:
        if cand in norm_map:
            return norm_map[cand]
    for name in fieldnames:
        norm_name = normalize_key(name)
        for cand in candidate_norms:
            if cand in norm_name or norm_name in cand:
                return name
    return None

def find_source_fasta(bin_name):
    base = os.path.join(bin_dir, bin_name)
    candidates = [base, base + ".fa", base + ".fna", base + ".fasta"]
    for c in candidates:
        if os.path.isfile(c):
            return c
    matches = sorted(glob.glob(base + ".*"))
    for m in matches:
        if os.path.isfile(m):
            return m
    return None

def quality_tier(c, x):
    if c >= hq_c and x <= hq_x:
        return "HQ"
    if include_mq and c >= mq_c and x <= mq_x:
        return "MQ"
    return "LQ"

selected = 0
hq = 0
mq = 0
total = 0

selected_rows = []
with open(report_path, newline="", encoding="utf-8", errors="replace") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []
    bin_col = pick_column(fields, ["Name", "Genome", "Genome Name", "Genome_Name", "Bin", "Bin ID", "bin_id", "bin_name", "GenomeID"])
    comp_col = pick_column(fields, ["Completeness", "Completeness_General", "Completeness_Specific"])
    cont_col = pick_column(fields, ["Contamination"])

    if bin_col is None:
        raise SystemExit("ERROR: could not identify bin name column in CheckM2 report")

    for row in reader:
        total += 1
        name = (row.get(bin_col) or "").strip()
        if not name:
            continue
        try:
            comp = float(row.get(comp_col) or 0) if comp_col else 0.0
        except Exception:
            comp = 0.0
        try:
            cont = float(row.get(cont_col) or 0) if cont_col else 0.0
        except Exception:
            cont = 0.0

        tier = quality_tier(comp, cont)
        if tier == "LQ":
            continue

        source_fasta = find_source_fasta(name)
        if source_fasta is None:
            raise SystemExit(f"ERROR: unable to locate source FASTA for selected bin {name} in {bin_dir}")

        staged_name = normalize_name(os.path.basename(source_fasta)) + ".fa"
        staged_path = os.path.join(staging_dir, staged_name)
        if os.path.lexists(staged_path):
            os.remove(staged_path)
        os.symlink(os.path.abspath(source_fasta), staged_path)

        selected += 1
        if tier == "HQ":
            hq += 1
        elif tier == "MQ":
            mq += 1

        selected_rows.append((staged_name, os.path.abspath(source_fasta), name, f"{comp:.6f}", f"{cont:.6f}", tier))

selected_manifest = os.path.join(staging_dir, "selected_bins.tsv")
with open(selected_manifest, "w", encoding="utf-8") as out:
    out.write("staged_bin\tsource_fasta\tcheckm2_bin_name\tcheckm2_completeness\tcheckm2_contamination\tcheckm2_quality_tier\n")
    for row in sorted(selected_rows):
        out.write("\t".join(row) + "\n")

print(f"{selected}\t{hq}\t{mq}\t{total}")
PY
)"

SELECTED_BINS="$(printf '%s' "$SELECTION_COUNTS" | cut -f1)"
SELECTED_HQ="$(printf '%s' "$SELECTION_COUNTS" | cut -f2)"
SELECTED_MQ="$(printf '%s' "$SELECTION_COUNTS" | cut -f3)"
TOTAL_BINS="$(printf '%s' "$SELECTION_COUNTS" | cut -f4)"

if [[ "$SELECTED_BINS" -eq 0 ]]; then
    echo "ERROR: no bins selected from CheckM2 report for ${SET_NAME}" >&2
    exit 1
fi

echo "[$(date)] Running GTDB-Tk for ${SET_NAME} (selected=${SELECTED_BINS}, hq=${SELECTED_HQ}, mq=${SELECTED_MQ})"
/usr/local/bin/singularity exec "$CONTAINER" bash -lc "set -euo pipefail; export GTDBTK_DATA_PATH='$GTDBTK_DB'; gtdbtk classify_wf --genome_dir '$STAGING_DIR' --out_dir '$RUN_DIR' --cpus '${SLURM_CPUS_PER_TASK:-32}' --extension .fa --skip_ani_screen"

GTDB_SUMMARY="${RUN_DIR}/gtdbtk.bac120.summary.tsv"
if [[ ! -f "$GTDB_SUMMARY" ]]; then
    GTDB_SUMMARY="$(find "$RUN_DIR" -type f -name '*.summary.tsv' | sort | head -n 1 || true)"
fi
if [[ -z "$GTDB_SUMMARY" || ! -f "$GTDB_SUMMARY" ]]; then
    echo "ERROR: no GTDB summary TSV found in: $RUN_DIR" >&2
    exit 1
fi

MERGED_SET_TSV="${OUT_ROOT}/summary/${SET_NAME}.merged_bins.tsv"
python3 - "$CHECKM2_QUALITY_REPORT" "$GTDB_SUMMARY" "$MERGED_SET_TSV" "$ASSEMBLY" "$SAMPLE" "$SET_NAME" "$HQ_COMPLETENESS" "$HQ_CONTAMINATION" "$MQ_COMPLETENESS" "$MQ_CONTAMINATION" <<'PY'
import csv
import os
import sys

checkm2_report, gtdb_summary, merged_tsv, assembly, sample, set_name, hq_c, hq_x, mq_c, mq_x = sys.argv[1:11]
hq_c = float(hq_c)
hq_x = float(hq_x)
mq_c = float(mq_c)
mq_x = float(mq_x)

def normalize_name(value: str) -> str:
    return os.path.splitext((value or "").strip())[0]

def normalize_key(value: str) -> str:
    return "".join(ch.lower() for ch in normalize_name(value) if ch.isalnum())

def pick_column(fieldnames, candidates):
    if not fieldnames:
        return None
    norm_map = {normalize_key(name): name for name in fieldnames}
    candidate_norms = [normalize_key(c) for c in candidates]
    for cand in candidate_norms:
        if cand in norm_map:
            return norm_map[cand]
    for name in fieldnames:
        norm_name = normalize_key(name)
        for cand in candidate_norms:
            if cand in norm_name or norm_name in cand:
                return name
    return None

def quality_tier(comp, cont):
    if comp >= hq_c and cont <= hq_x:
        return "HQ"
    if comp >= mq_c and cont <= mq_x:
        return "MQ"
    return "LQ"

checkm2_map = {}
with open(checkm2_report, newline="", encoding="utf-8", errors="replace") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []
    bin_col = pick_column(fields, ["Name", "Genome", "Genome Name", "Genome_Name", "Bin", "Bin ID", "bin_id", "bin_name", "GenomeID"])
    comp_col = pick_column(fields, ["Completeness", "Completeness_General", "Completeness_Specific"])
    cont_col = pick_column(fields, ["Contamination"])
    if bin_col:
        for row in reader:
            name = (row.get(bin_col) or "").strip()
            if not name:
                continue
            try:
                comp = float(row.get(comp_col) or 0) if comp_col else 0.0
            except Exception:
                comp = 0.0
            try:
                cont = float(row.get(cont_col) or 0) if cont_col else 0.0
            except Exception:
                cont = 0.0
            checkm2_map[normalize_key(name)] = {
                "checkm2_bin_name": name,
                "checkm2_completeness": f"{comp:.6f}",
                "checkm2_contamination": f"{cont:.6f}",
                "checkm2_quality_tier": quality_tier(comp, cont),
            }

rows = []
with open(gtdb_summary, newline="", encoding="utf-8", errors="replace") as gh:
    gtdb_reader = csv.DictReader(gh, delimiter="\t")
    gtdb_fields = gtdb_reader.fieldnames or []
    for grow in gtdb_reader:
        bin_id = (grow.get("user_genome") or grow.get("genome") or "").strip()
        c2 = checkm2_map.get(normalize_key(bin_id), {})
        out = {
            "assembly": assembly,
            "sample": sample,
            "set_name": set_name,
            "bin_id": bin_id,
            "checkm2_bin_name": c2.get("checkm2_bin_name", ""),
            "checkm2_completeness": c2.get("checkm2_completeness", ""),
            "checkm2_contamination": c2.get("checkm2_contamination", ""),
            "checkm2_quality_tier": c2.get("checkm2_quality_tier", ""),
        }
        for field in gtdb_fields:
            out[field] = grow.get(field, "")
        rows.append(out)

meta_fields = [
    "assembly",
    "sample",
    "set_name",
    "bin_id",
    "checkm2_bin_name",
    "checkm2_completeness",
    "checkm2_contamination",
    "checkm2_quality_tier",
]
header = meta_fields + gtdb_fields

with open(merged_tsv, "w", newline="", encoding="utf-8") as out:
    writer = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(len(rows))
PY

SUMMARY_TSV="${OUT_ROOT}/summary/gtdbtk_single_summary.tsv"
if [[ ! -f "$SUMMARY_TSV" ]]; then
    printf "assembly\tsample\tset_name\tcheckm2_quality_report\tbin_dir\tselected_bins\tselected_hq\tselected_mq\ttotal_bins\tgtdbtk_run_dir\tgtdbtk_summary\tmerged_bins_tsv\thq_completeness\thq_contamination\tmq_completeness\tmq_contamination\tstatus\tnotes\n" > "$SUMMARY_TSV"
fi

TMP_SUMMARY="${SUMMARY_TSV}.tmp.$$"
awk -F'\t' -v set_name="$SET_NAME" 'NR==1 || $3 != set_name' "$SUMMARY_TSV" > "$TMP_SUMMARY"
mv "$TMP_SUMMARY" "$SUMMARY_TSV"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$ASSEMBLY" "$SAMPLE" "$SET_NAME" "$CHECKM2_QUALITY_REPORT" "$BIN_DIR" \
    "$SELECTED_BINS" "$SELECTED_HQ" "$SELECTED_MQ" "$TOTAL_BINS" \
    "$RUN_DIR" "$GTDB_SUMMARY" "$MERGED_SET_TSV" \
    "$HQ_COMPLETENESS" "$HQ_CONTAMINATION" "$MQ_COMPLETENESS" "$MQ_CONTAMINATION" \
    "completed" "single-set-gtdbtk" >> "$SUMMARY_TSV"

echo "[$(date)] GTDB-Tk complete: ${SET_NAME}"
echo "[$(date)] GTDB summary: $GTDB_SUMMARY"
echo "[$(date)] Merged bins TSV: $MERGED_SET_TSV"
echo "[$(date)] Summary TSV: $SUMMARY_TSV"

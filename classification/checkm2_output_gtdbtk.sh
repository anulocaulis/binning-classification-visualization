#!/usr/bin/env bash
#SBATCH --job-name=checkm2_output_gtdbtk
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/checkm2_output_gtdbtk_%j.out
#SBATCH --error=logs/checkm2_output_gtdbtk_%j.err

set -euo pipefail

# Classify CheckM2 HQ/MQ bins with GTDB-Tk.
# Accepts one or more CheckM2 summary TSVs or CheckM2 output roots.
#
# Usage:
#   sbatch scripts/checkm2_output_gtdbtk.sh
#   sbatch scripts/checkm2_output_gtdbtk.sh <CHECKM2_INPUT_1> [CHECKM2_INPUT_2 ...]
#   sbatch scripts/checkm2_output_gtdbtk.sh <CHECKM2_INPUT_1> [CHECKM2_INPUT_2 ...] --out-root <OUT_ROOT> --gtdbtk-db <GTDBTK_DB>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/../configs"
if [[ -f "$CONFIG_DIR/configs_master.conf" ]]; then
    # shellcheck source=/dev/null
    source "$CONFIG_DIR/configs_master.conf"
fi

BASE_DIR="${PROJECT_BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
WORK_DIR="${PROJECT_WORK_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
DATA_DIR="${PROJECT_DATA_DIR:-${BASE_DIR}/data}"

DEFAULT_INPUTS=(
    "${DATA_DIR}/binning_outputs/checkm2_refined_bins"
    "${DATA_DIR}/binning_outputs/checkm2_vamb_bins_full"
)
CHECKM2_INPUTS=()
OUT_ROOT="${DATA_DIR}/binning_outputs/gtdbtk_checkm2_all"
GTDBTK_DB="${BASE_DIR}/databases/gtdbtk_db"

while [ "$#" -gt 0 ]; do
    case "$1" in
        --out-root)
            OUT_ROOT="$2"
            shift 2
            ;;
        --gtdbtk-db)
            GTDBTK_DB="$2"
            shift 2
            ;;
        -h|--help)
            cat <<EOF
Usage:
  sbatch scripts/checkm2_output_gtdbtk.sh [CHECKM2_INPUT ...] [--out-root OUT_ROOT] [--gtdbtk-db GTDBTK_DB]

Examples:
  sbatch scripts/checkm2_output_gtdbtk.sh
  sbatch scripts/checkm2_output_gtdbtk.sh \
        /path/to/data/binning_outputs/checkm2_refined_bins \
        /path/to/data/binning_outputs/checkm2_vamb_bins
EOF
            exit 0
            ;;
        --*)
            echo "ERROR: Unknown option: $1" >&2
            exit 1
            ;;
        *)
            CHECKM2_INPUTS+=("$1")
            shift
            ;;
    esac
done

if [ "${#CHECKM2_INPUTS[@]}" -eq 0 ]; then
    CHECKM2_INPUTS=("${DEFAULT_INPUTS[@]}")
fi

CONTAINER="${WORK_DIR}/containers/classification.sif"
GTDBTK_BIN="gtdbtk"

mkdir -p "${OUT_ROOT}" "${OUT_ROOT}/summary" "${OUT_ROOT}/runs" "${OUT_ROOT}/staging" "${WORK_DIR}/slurm_logs"

resolve_summary() {
    local input="$1"
    if [ -f "$input" ]; then
        echo "$input"
        return 0
    fi
    if [ -d "$input" ]; then
        for candidate in \
            "$input/summary/checkm2_refined_bins_summary.tsv" \
            "$input/summary/metawrap_fastmode_checkm2_summary.tsv"; do
            if [ -f "$candidate" ]; then
                echo "$candidate"
                return 0
            fi
        done

        if find "$input" -mindepth 1 -maxdepth 2 -type f -name 'quality_report.tsv' | grep -q .; then
            echo "$input"
            return 0
        fi
    fi
    return 1
}

build_manifest_from_quality_reports() {
    local input_root="$1"
    local manifest_path="$2"

    python3 - "$input_root" "$manifest_path" "$DATA_DIR" <<'PY'
import csv
import os
import sys
from pathlib import Path

input_root, manifest_path, data_dir = sys.argv[1:4]

def split_set_name(set_name: str):
    stem = set_name
    for suffix in ("_metawrap_bins_refined", "_metawrap_bins", "_vamb"):
        if stem.endswith(suffix):
            stem = stem[:-len(suffix)]
            break
    parts = stem.rsplit("_", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return stem, ""

def candidate_bin_dirs(root_name: str, set_name: str):
    candidates = []

    def add_if_dir(path):
        if os.path.isdir(path):
            candidates.append(path)

    root_lc = root_name.lower()
    bo = os.path.join(data_dir, "binning_outputs")

    # Normalize possible naming variants for VAMB-style sets.
    vamb_names = {set_name}
    if set_name.endswith("_metawrap_bins_refined"):
        vamb_names.add(set_name[: -len("_metawrap_bins_refined")] + "_vamb")
    if set_name.endswith("_metawrap_bins"):
        vamb_names.add(set_name[: -len("_metawrap_bins")] + "_vamb")

    # Refined-bin candidates (non-gold-standard metawrap-style bins).
    add_if_dir(os.path.join(bo, "refined_bins", set_name, "fast_mode_bins"))
    add_if_dir(os.path.join(bo, "refined_bins", set_name))

    # VAMB candidates from long-mapped outputs.
    for name in sorted(vamb_names):
        add_if_dir(os.path.join(bo, "long_mapped_binning", name, "vamb_bins_extracted"))
        add_if_dir(os.path.join(bo, "long_mapped_binning", name, "vamb_bins"))
        add_if_dir(os.path.join(bo, "long_mapped_binning", name))

    # VAMB candidates from pre-existing CheckM2 prepared bins.
    for name in sorted(vamb_names):
        add_if_dir(os.path.join(bo, "checkm2_existing_bins", "prepared_vamb_bins", name))

    # Gold-standard VAMB outputs can be under run_*/vamb_bins/<set_name>/vamb_bins_extracted|vamb_bins.
    gsb_root = os.path.join(bo, "gold_standard_binner_runs")
    if os.path.isdir(gsb_root):
        for run_dir in sorted(os.listdir(gsb_root), reverse=True):
            run_path = os.path.join(gsb_root, run_dir)
            if not os.path.isdir(run_path):
                continue
            for name in sorted(vamb_names):
                add_if_dir(os.path.join(run_path, "vamb_bins", name, "vamb_bins_extracted"))
                add_if_dir(os.path.join(run_path, "vamb_bins", name, "vamb_bins"))

    # Prefer root-specific intent when possible.
    if "vamb" in root_lc:
        ordered = [p for p in candidates if "vamb" in p.lower()] + [p for p in candidates if "vamb" not in p.lower()]
    else:
        ordered = [p for p in candidates if "refined_bins" in p] + [p for p in candidates if "refined_bins" not in p]

    # Deduplicate while preserving order.
    seen = set()
    unique = []
    for p in ordered:
        if p in seen:
            continue
        seen.add(p)
        unique.append(p)
    return unique

rows = []
root_name = os.path.basename(os.path.normpath(input_root))
for report in sorted(Path(input_root).glob("*/quality_report.tsv")):
    set_name = report.parent.name
    assembly, sample = split_set_name(set_name)
    bin_dir = ""
    for candidate in candidate_bin_dirs(root_name, set_name):
        if os.path.isdir(candidate):
            bin_dir = candidate
            break

    if not bin_dir:
        continue

    total_bins = 0
    hq_bins = 0
    mq_bins = 0
    lq_bins = 0
    with open(report, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        comp_col = "Completeness_General" if "Completeness_General" in fields else "Completeness"
        cont_col = "Contamination"
        for row in reader:
            total_bins += 1
            try:
                completeness = float(row.get(comp_col) or 0)
            except Exception:
                completeness = 0.0
            try:
                contamination = float(row.get(cont_col) or 0)
            except Exception:
                contamination = 0.0

            if completeness >= 90 and contamination <= 5:
                hq_bins += 1
            elif completeness >= 50 and contamination <= 10:
                mq_bins += 1
            else:
                lq_bins += 1

    rows.append([
        assembly,
        sample,
        set_name,
        bin_dir,
        str(total_bins),
        str(hq_bins),
        str(mq_bins),
        str(lq_bins),
        str(report),
        "completed",
        "generated-from-quality-report",
    ])

with open(manifest_path, "w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["assembly", "sample", "set_name", "bin_dir", "total_bins", "hq_bins", "mq_bins", "lq_bins", "quality_report", "status", "notes"])
    writer.writerows(rows)
PY
}

SUMMARY_INPUTS=()
GENERATED_MANIFESTS=()
for input in "${CHECKM2_INPUTS[@]}"; do
    if resolved=$(resolve_summary "$input"); then
        if [ -d "$resolved" ]; then
            manifest_name="$(basename "$resolved")_gtdbtk_manifest.tsv"
            manifest_path="${OUT_ROOT}/summary/${manifest_name}"
            build_manifest_from_quality_reports "$resolved" "$manifest_path"
            GENERATED_MANIFESTS+=("$manifest_path")
            SUMMARY_INPUTS+=("$manifest_path")
        else
            SUMMARY_INPUTS+=("$resolved")
        fi
    else
        echo "ERROR: CheckM2 summary TSV not found from input: $input" >&2
        exit 1
    fi
done

if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: container not found: $CONTAINER" >&2
    exit 1
fi

if [ ! -d "$GTDBTK_DB" ]; then
    echo "ERROR: GTDB-Tk database directory not found: $GTDBTK_DB" >&2
    exit 1
fi

if ! /usr/local/bin/singularity exec "$CONTAINER" bash -lc "command -v $GTDBTK_BIN >/dev/null 2>&1"; then
    echo "ERROR: gtdbtk not found inside container: $CONTAINER" >&2
    exit 1
fi

echo "[$(date)] Starting GTDB-Tk classification of CheckM2 HQ/MQ bins"
printf '[%s] SUMMARY_INPUT=%s\n' "$(date)" "${SUMMARY_INPUTS[@]}"
echo "[$(date)] OUT_ROOT=$OUT_ROOT"
echo "[$(date)] GTDBTK_DB=$GTDBTK_DB"
echo "[$(date)] CONTAINER=$CONTAINER"

SUMMARY_TSV="${OUT_ROOT}/summary/checkm2_gtdbtk_summary.tsv"
MERGED_BINS_TSV="${OUT_ROOT}/summary/checkm2_gtdbtk_merged_bins.tsv"
TMP_SUMMARY="${SUMMARY_TSV}.tmp.$$"
printf "source_summary\tassembly\tsample\tset_name\tcheckm2_quality_report\ttotal_bins\thq_bins\tmq_bins\tselected_bins\tselection_dir\tgtdbtk_run_dir\tgtdbtk_summary\n" > "$TMP_SUMMARY"

SELECTED_ANY=0

for CHECKM2_SUMMARY in "${SUMMARY_INPUTS[@]}"; do
    while IFS=$'\t' read -r assembly sample set_name bin_dir total_bins hq_bins mq_bins lq_bins quality_report status notes; do
        if [ "$assembly" = "assembly" ] || [ -z "$set_name" ]; then
            continue
        fi

        if [ "$status" != "completed" ] && [ "$notes" != "generated-from-quality-report" ]; then
            echo "[$(date)] [skip] ${set_name}: status=${status}"
            continue
        fi

        if [ ! -d "$bin_dir" ]; then
            echo "[$(date)] [skip] ${set_name}: bin_dir missing: $bin_dir"
            continue
        fi

        if [ -z "$quality_report" ] || [ ! -f "$quality_report" ]; then
            echo "[$(date)] [skip] ${set_name}: quality_report missing: $quality_report"
            continue
        fi

        selection_dir="${OUT_ROOT}/staging/${set_name}"
        run_dir="${OUT_ROOT}/runs/${set_name}"
        rm -rf "$selection_dir" "$run_dir"
        mkdir -p "$selection_dir" "$run_dir"

        IFS=$'\t' read -r selected_bins selected_hq selected_mq selected_total < <(python3 - "$quality_report" "$bin_dir" "$selection_dir" <<'PY'
import csv
import glob
import os
import sys

report_path, bin_dir, selection_dir = sys.argv[1:4]

def normalize(value: str) -> str:
    return ''.join(ch.lower() for ch in value if ch.isalnum())

def pick_column(fieldnames, candidates):
    if not fieldnames:
        return None
    norm_map = {normalize(name): name for name in fieldnames}
    candidate_norms = [normalize(c) for c in candidates]
    for cand in candidate_norms:
        if cand in norm_map:
            return norm_map[cand]
    for name in fieldnames:
        norm_name = normalize(name)
        for cand in candidate_norms:
            if cand in norm_name or norm_name in cand:
                return name
    return None

def find_source_fasta(bin_name):
    base = os.path.join(bin_dir, bin_name)
    candidates = [base, base + '.fa', base + '.fasta', base + '.fna']
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    matches = sorted(glob.glob(base + '.*'))
    return matches[0] if matches else None

with open(report_path, 'r', encoding='utf-8', errors='replace') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    fieldnames = reader.fieldnames or []
    bin_col = pick_column(fieldnames, ['Name', 'Genome', 'Genome Name', 'Genome_Name', 'Bin', 'Bin ID', 'bin_id', 'bin_name', 'GenomeID'])
    comp_col = pick_column(fieldnames, ['Completeness', 'Completeness_General', 'Completeness_Specific'])
    cont_col = pick_column(fieldnames, ['Contamination'])

    if comp_col is None or cont_col is None:
        raise SystemExit(f'ERROR: could not locate completeness/contamination columns in {report_path}')

    total = 0
    hq = 0
    mq = 0
    selected = 0

    for row in reader:
        total += 1
        try:
            completeness = float(row.get(comp_col, 0) or 0)
        except Exception:
            completeness = 0.0
        try:
            contamination = float(row.get(cont_col, 0) or 0)
        except Exception:
            contamination = 0.0

        is_hq = completeness >= 90 and contamination <= 5
        is_mq = completeness >= 50 and completeness < 90 and contamination <= 10
        if not (is_hq or is_mq):
            continue

        bin_name = (row.get(bin_col, '') or '').strip()
        if not bin_name:
            continue

        source_fasta = find_source_fasta(bin_name)
        if source_fasta is None:
            raise SystemExit(f'ERROR: could not find FASTA for selected bin {bin_name} in {bin_dir}')

        staged_name = os.path.splitext(os.path.basename(source_fasta))[0] + '.fa'
        staged_path = os.path.join(selection_dir, staged_name)
        if os.path.lexists(staged_path):
            os.remove(staged_path)
        os.symlink(os.path.abspath(source_fasta), staged_path)

        selected += 1
        if is_hq:
            hq += 1
        elif is_mq:
            mq += 1

    selected_list_path = os.path.join(selection_dir, 'selected_bins.tsv')
    with open(selected_list_path, 'w', encoding='utf-8') as out:
        out.write('staged_bin\tsource_fasta\n')
        for staged in sorted(os.listdir(selection_dir)):
            if staged.endswith('.fa'):
                source = os.path.realpath(os.path.join(selection_dir, staged))
                out.write(f'{staged}\t{source}\n')

    print(f'{selected}\t{hq}\t{mq}\t{total}')
PY
)

        if [ "$selected_bins" -eq 0 ]; then
            echo "[$(date)] [skip] ${set_name}: no HQ/MQ bins selected from CheckM2"

            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$CHECKM2_SUMMARY" "$assembly" "$sample" "$set_name" "$quality_report" \
                "$total_bins" "$selected_hq" "$selected_mq" "$selected_bins" \
                "$selection_dir" "" "skipped-no-hqmq-bins" >> "$TMP_SUMMARY"
            continue
        fi

        SELECTED_ANY=1
        set_outdir="${OUT_ROOT}/runs/${set_name}"
        rm -rf "$set_outdir"
        mkdir -p "$set_outdir"

        echo "[$(date)] Running GTDB-Tk for ${set_name} (selected=${selected_bins}, HQ=${selected_hq}, MQ=${selected_mq}, total_checkm2=${selected_total})"
        /usr/local/bin/singularity exec "$CONTAINER" bash -lc "set -euo pipefail; export GTDBTK_DATA_PATH='$GTDBTK_DB'; gtdbtk classify_wf --genome_dir '$selection_dir' --out_dir '$set_outdir' --cpus '${SLURM_CPUS_PER_TASK:-32}' --extension .fa --skip_ani_screen"

        gtdb_summary="${set_outdir}/gtdbtk.bac120.summary.tsv"
        if [ ! -f "$gtdb_summary" ]; then
            fallback_summary="$(find "$set_outdir" -maxdepth 1 -type f -name '*.summary.tsv' | sort | head -n 1 || true)"
            if [ -n "$fallback_summary" ]; then
                gtdb_summary="$fallback_summary"
            fi
        fi

        if [ ! -f "$gtdb_summary" ]; then
            echo "[$(date)] ERROR: GTDB-Tk did not produce a summary TSV for ${set_name}" >&2
            exit 1
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$CHECKM2_SUMMARY" "$assembly" "$sample" "$set_name" "$quality_report" \
            "$total_bins" "$selected_hq" "$selected_mq" "$selected_bins" \
            "$selection_dir" "$set_outdir" "$gtdb_summary" >> "$TMP_SUMMARY"

        echo "[$(date)] Completed ${set_name}: GTDB-Tk summary=$gtdb_summary"
    done < "$CHECKM2_SUMMARY"
done

if [ "$SELECTED_ANY" -eq 0 ]; then
    echo "[$(date)] ERROR: No HQ/MQ bins selected from CheckM2 summaries" >&2
    rm -f "$TMP_SUMMARY"
    exit 1
fi

mv "$TMP_SUMMARY" "$SUMMARY_TSV"

python3 - "$SUMMARY_TSV" "$MERGED_BINS_TSV" <<'PY'
import csv
import os
import sys

summary_tsv, merged_tsv = sys.argv[1:3]

meta_fields = [
    "source_summary",
    "assembly",
    "sample",
    "set_name",
    "checkm2_quality_report",
    "total_bins",
    "hq_bins",
    "mq_bins",
    "selected_bins",
    "gtdbtk_run_dir",
    "gtdbtk_summary",
    "bin_id",
]

merged_rows = []
gtdb_fields = []
seen_fields = set()

with open(summary_tsv, newline="", encoding="utf-8", errors="replace") as sh:
    summary_reader = csv.DictReader(sh, delimiter="\t")
    for summary_row in summary_reader:
        gtdb_summary = (summary_row.get("gtdbtk_summary") or "").strip()
        if not gtdb_summary or gtdb_summary == "skipped-no-hqmq-bins" or not os.path.isfile(gtdb_summary):
            continue

        with open(gtdb_summary, newline="", encoding="utf-8", errors="replace") as gh:
            gtdb_reader = csv.DictReader(gh, delimiter="\t")
            fields = gtdb_reader.fieldnames or []
            for field in fields:
                if field not in seen_fields:
                    seen_fields.add(field)
                    gtdb_fields.append(field)

            for gtdb_row in gtdb_reader:
                out_row = {k: summary_row.get(k, "") for k in meta_fields if k != "bin_id"}
                out_row["bin_id"] = gtdb_row.get("user_genome", "") or gtdb_row.get("genome", "") or ""
                for field in fields:
                    out_row[field] = gtdb_row.get(field, "")
                merged_rows.append(out_row)

with open(merged_tsv, "w", newline="", encoding="utf-8") as out:
    header = meta_fields + gtdb_fields
    writer = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    for row in merged_rows:
        writer.writerow(row)
PY

echo ""
echo "[$(date)] GTDB-Tk classification finished"
echo "[$(date)] Summary TSV: $SUMMARY_TSV"
echo "[$(date)] Merged bins TSV: $MERGED_BINS_TSV"

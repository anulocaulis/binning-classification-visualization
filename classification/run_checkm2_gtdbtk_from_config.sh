#!/usr/bin/env bash
#SBATCH --job-name=checkm2_gtdbtk_wrapper
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/checkm2_gtdbtk_wrapper_%j.out
#SBATCH --error=logs/checkm2_gtdbtk_wrapper_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/../configs"
DEFAULT_CONFIG="${CONFIG_DIR}/checkm2_gtdbtk_wrapper.conf"
CONFIG_PATH="$DEFAULT_CONFIG"

SAMPLE_NAME_OVERRIDE=""
SAMPLE_NUMBER_OVERRIDE=""
ASSEMBLY_NAME_OVERRIDE=""
ASSEMBLY_NUMBER_OVERRIDE=""

usage() {
    cat <<EOF
Usage:
    bash classification/run_checkm2_gtdbtk_from_config.sh [options]

Options:
    --config PATH              Config file (default: configs/checkm2_gtdbtk_wrapper.conf)
  --sample NAME              Run only one sample by name (e.g., S1)
  --sample-number N          Run only one sample by 1-based index from config SAMPLES array
  --assembly NAME            Run only one assembly by name (e.g., flye)
  --assembly-number N        Run only one assembly by 1-based index from config ASSEMBLIES array
  -h, --help                 Show help

Notes:
  - Without explicit sample/assembly selectors, all configured combinations are run.
  - Final combined TSV path is controlled by COMBINED_TSV in the config file.
EOF
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --sample)
            SAMPLE_NAME_OVERRIDE="$2"
            shift 2
            ;;
        --sample-number)
            SAMPLE_NUMBER_OVERRIDE="$2"
            shift 2
            ;;
        --assembly)
            ASSEMBLY_NAME_OVERRIDE="$2"
            shift 2
            ;;
        --assembly-number)
            ASSEMBLY_NUMBER_OVERRIDE="$2"
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

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "ERROR: config file not found: $CONFIG_PATH" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG_PATH"

required_vars=(
    DATA_DIR
    CHECKM2_SCRIPT
    GTDBTK_SCRIPT
    CHECKM2_OUT_ROOT
    GTDBTK_OUT_ROOT
    COMBINED_TSV
    BIN_DIR_TEMPLATE
    CHECKM2_CONTAINER
    CHECKM2_DB
    GTDBTK_CONTAINER
    GTDBTK_DB
    SET_SUFFIX
)

for var_name in "${required_vars[@]}"; do
    if [[ -z "${!var_name:-}" ]]; then
        echo "ERROR: missing required config variable: $var_name" >&2
        exit 1
    fi
done

if [[ -z "${SAMPLES+x}" || "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: SAMPLES array is empty in config" >&2
    exit 1
fi
if [[ -z "${ASSEMBLIES+x}" || "${#ASSEMBLIES[@]}" -eq 0 ]]; then
    echo "ERROR: ASSEMBLIES array is empty in config" >&2
    exit 1
fi

if [[ ! -f "$CHECKM2_SCRIPT" ]]; then
    echo "ERROR: CHECKM2_SCRIPT not found: $CHECKM2_SCRIPT" >&2
    exit 1
fi
if [[ ! -f "$GTDBTK_SCRIPT" ]]; then
    echo "ERROR: GTDBTK_SCRIPT not found: $GTDBTK_SCRIPT" >&2
    exit 1
fi

resolve_index() {
    local idx="$1"
    local max="$2"
    if [[ ! "$idx" =~ ^[0-9]+$ ]]; then
        echo "ERROR"
        return
    fi
    if (( idx < 1 || idx > max )); then
        echo "ERROR"
        return
    fi
    echo $((idx - 1))
}

TARGET_SAMPLES=()
if [[ -n "$SAMPLE_NAME_OVERRIDE" ]]; then
    TARGET_SAMPLES=("$SAMPLE_NAME_OVERRIDE")
elif [[ -n "$SAMPLE_NUMBER_OVERRIDE" ]]; then
    idx="$(resolve_index "$SAMPLE_NUMBER_OVERRIDE" "${#SAMPLES[@]}")"
    if [[ "$idx" == "ERROR" ]]; then
        echo "ERROR: --sample-number out of range for configured SAMPLES" >&2
        exit 1
    fi
    TARGET_SAMPLES=("${SAMPLES[$idx]}")
else
    TARGET_SAMPLES=("${SAMPLES[@]}")
fi

TARGET_ASSEMBLIES=()
if [[ -n "$ASSEMBLY_NAME_OVERRIDE" ]]; then
    TARGET_ASSEMBLIES=("$ASSEMBLY_NAME_OVERRIDE")
elif [[ -n "$ASSEMBLY_NUMBER_OVERRIDE" ]]; then
    idx="$(resolve_index "$ASSEMBLY_NUMBER_OVERRIDE" "${#ASSEMBLIES[@]}")"
    if [[ "$idx" == "ERROR" ]]; then
        echo "ERROR: --assembly-number out of range for configured ASSEMBLIES" >&2
        exit 1
    fi
    TARGET_ASSEMBLIES=("${ASSEMBLIES[$idx]}")
else
    TARGET_ASSEMBLIES=("${ASSEMBLIES[@]}")
fi

mkdir -p "$CHECKM2_OUT_ROOT" "$GTDBTK_OUT_ROOT" "$(dirname "$COMBINED_TSV")"

expand_bin_template() {
    local tmpl="$1"
    local assembly="$2"
    local sample="$3"
    local sample_lower="$4"
    local set_suffix="$5"

    local out="$tmpl"
    out="${out//\{data_dir\}/$DATA_DIR}"
    out="${out//\{assembly\}/$assembly}"
    out="${out//\{sample\}/$sample}"
    out="${out//\{sample_lower\}/$sample_lower}"
    out="${out//\{set_suffix\}/$set_suffix}"
    printf '%s\n' "$out"
}

MERGED_SET_TABLES=()
RUN_LOG_TSV="${GTDBTK_OUT_ROOT}/summary/wrapper_run_log.tsv"
printf "assembly\tsample\tset_name\tbin_dir\tcheckm2_quality_report\tmerged_set_tsv\tstatus\tnotes\n" > "$RUN_LOG_TSV"

for sample in "${TARGET_SAMPLES[@]}"; do
    sample_lower="$(printf '%s' "$sample" | tr '[:upper:]' '[:lower:]')"
    for assembly in "${TARGET_ASSEMBLIES[@]}"; do
        set_name="${assembly}_${sample_lower}_${SET_SUFFIX}"
        bin_dir="$(expand_bin_template "$BIN_DIR_TEMPLATE" "$assembly" "$sample" "$sample_lower" "$SET_SUFFIX")"

        echo "[$(date)] Processing ${set_name}"

        if [[ ! -d "$bin_dir" ]]; then
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$assembly" "$sample" "$set_name" "$bin_dir" "" "" "skipped" "bin_dir_missing" >> "$RUN_LOG_TSV"
            echo "[$(date)] Skipping ${set_name}: bin_dir missing (${bin_dir})"
            continue
        fi

        bash "$CHECKM2_SCRIPT" \
            --assembly "$assembly" \
            --sample "$sample" \
            --set-name "$set_name" \
            --set-suffix "$SET_SUFFIX" \
            --bin-dir "$bin_dir" \
            --out-root "$CHECKM2_OUT_ROOT" \
            --container "$CHECKM2_CONTAINER" \
            --checkm2-db "$CHECKM2_DB" \
            --hq-completeness "$HQ_COMPLETENESS" \
            --hq-contamination "$HQ_CONTAMINATION" \
            --mq-completeness "$MQ_COMPLETENESS" \
            --mq-contamination "$MQ_CONTAMINATION"

        quality_report="${CHECKM2_OUT_ROOT}/runs/${set_name}/quality_report.tsv"
        if [[ ! -f "$quality_report" ]]; then
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$assembly" "$sample" "$set_name" "$bin_dir" "" "" "failed" "missing_checkm2_quality_report" >> "$RUN_LOG_TSV"
            echo "[$(date)] ERROR: missing quality report for ${set_name}" >&2
            continue
        fi

        gtdb_args=(
            --assembly "$assembly"
            --sample "$sample"
            --set-name "$set_name"
            --set-suffix "$SET_SUFFIX"
            --bin-dir "$bin_dir"
            --checkm2-quality-report "$quality_report"
            --out-root "$GTDBTK_OUT_ROOT"
            --container "$GTDBTK_CONTAINER"
            --gtdbtk-db "$GTDBTK_DB"
            --hq-completeness "$HQ_COMPLETENESS"
            --hq-contamination "$HQ_CONTAMINATION"
            --mq-completeness "$MQ_COMPLETENESS"
            --mq-contamination "$MQ_CONTAMINATION"
        )
        if [[ "${GTDBTK_HQ_ONLY:-0}" == "1" ]]; then
            gtdb_args+=(--hq-only)
        fi

        bash "$GTDBTK_SCRIPT" "${gtdb_args[@]}"

        merged_set_tsv="${GTDBTK_OUT_ROOT}/summary/${set_name}.merged_bins.tsv"
        if [[ -f "$merged_set_tsv" ]]; then
            MERGED_SET_TABLES+=("$merged_set_tsv")
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$assembly" "$sample" "$set_name" "$bin_dir" "$quality_report" "$merged_set_tsv" "completed" "ok" >> "$RUN_LOG_TSV"
        else
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$assembly" "$sample" "$set_name" "$bin_dir" "$quality_report" "" "failed" "missing_merged_set_tsv" >> "$RUN_LOG_TSV"
            echo "[$(date)] ERROR: missing merged set TSV for ${set_name}" >&2
        fi
    done
done

if [[ "${#MERGED_SET_TABLES[@]}" -eq 0 ]]; then
    echo "ERROR: no merged set TSVs were produced" >&2
    exit 1
fi

python3 - "$COMBINED_TSV" "${MERGED_SET_TABLES[@]}" <<'PY'
import csv
import sys

out_path = sys.argv[1]
in_paths = sys.argv[2:]

header = []
seen = set()
rows = []

for path in in_paths:
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        for field in fields:
            if field not in seen:
                seen.add(field)
                header.append(field)
        for row in reader:
            rows.append(row)

with open(out_path, "w", newline="", encoding="utf-8") as out:
    writer = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(f"WROTE {out_path} rows={len(rows)}")
PY

echo "[$(date)] Wrapper complete"
echo "[$(date)] Combined TSV: $COMBINED_TSV"
echo "[$(date)] Run log: $RUN_LOG_TSV"

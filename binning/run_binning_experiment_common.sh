#!/usr/bin/env bash
#SBATCH --job-name=binning_wrapper
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/binning_wrapper_%j.out
#SBATCH --error=logs/binning_wrapper_%j.err

set -euo pipefail

CONFIG_PATH=""
BINNER_OVERRIDE=""
SAMPLE_OVERRIDE=""
ASSEMBLY_OVERRIDE=""
SAMPLE_NUMBER_OVERRIDE=""
ASSEMBLY_NUMBER_OVERRIDE=""

usage() {
    cat <<EOF
Usage:
  bash scripts/run_binning_experiment_common.sh --config PATH [options]

Options:
  --binner NAME         metawrap or vamb (override config BINNER)
  --sample NAME         single sample override
  --sample-number N     1-based sample index from config TARGET_SAMPLES
  --assembly NAME       single assembly override
  --assembly-number N   1-based assembly index from config ASSEMBLIES
EOF
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --binner)
            BINNER_OVERRIDE="$2"
            shift 2
            ;;
        --sample)
            SAMPLE_OVERRIDE="$2"
            shift 2
            ;;
        --sample-number)
            SAMPLE_NUMBER_OVERRIDE="$2"
            shift 2
            ;;
        --assembly)
            ASSEMBLY_OVERRIDE="$2"
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

if [[ -z "$CONFIG_PATH" || ! -f "$CONFIG_PATH" ]]; then
    echo "ERROR: --config is required and must exist" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG_PATH"

required_vars=(
    EXPERIMENT_NAME
    BINNER
    ASSEMBLIES
    TARGET_SAMPLES
    ASSEMBLY_TEMPLATE
    SHORT_READ_TEMPLATE
    LONG_READ_TEMPLATE
    OUTPUT_ROOT
    SET_SUFFIX
    METAWRAP_BINNER_SCRIPT
    VAMB_BINNER_SCRIPT
)

for var in "${required_vars[@]}"; do
    if [[ "$var" == "ASSEMBLIES" || "$var" == "TARGET_SAMPLES" ]]; then
        continue
    fi
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: required config variable missing: $var" >&2
        exit 1
    fi
done

if [[ -z "${ASSEMBLIES+x}" || "${#ASSEMBLIES[@]}" -eq 0 ]]; then
    echo "ERROR: ASSEMBLIES array is empty" >&2
    exit 1
fi
if [[ -z "${TARGET_SAMPLES+x}" || "${#TARGET_SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: TARGET_SAMPLES array is empty" >&2
    exit 1
fi

if [[ -n "$BINNER_OVERRIDE" ]]; then
    BINNER="$BINNER_OVERRIDE"
fi
if [[ "$BINNER" != "metawrap" && "$BINNER" != "vamb" ]]; then
    echo "ERROR: BINNER must be metawrap or vamb" >&2
    exit 1
fi

CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-32}}"
MEM_GB="${MEM_GB:-360}"
MIN_CONTIG_LEN="${MIN_CONTIG_LEN:-2000}"
COBINNING_MODE="${COBINNING_MODE:-0}"

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

TARGET_SAMPLES_RUN=()
if [[ -n "$SAMPLE_OVERRIDE" ]]; then
    TARGET_SAMPLES_RUN=("$SAMPLE_OVERRIDE")
elif [[ -n "$SAMPLE_NUMBER_OVERRIDE" ]]; then
    sidx="$(resolve_index "$SAMPLE_NUMBER_OVERRIDE" "${#TARGET_SAMPLES[@]}")"
    if [[ "$sidx" == "ERROR" ]]; then
        echo "ERROR: sample-number out of range" >&2
        exit 1
    fi
    TARGET_SAMPLES_RUN=("${TARGET_SAMPLES[$sidx]}")
else
    TARGET_SAMPLES_RUN=("${TARGET_SAMPLES[@]}")
fi

ASSEMBLIES_RUN=()
if [[ -n "$ASSEMBLY_OVERRIDE" ]]; then
    ASSEMBLIES_RUN=("$ASSEMBLY_OVERRIDE")
elif [[ -n "$ASSEMBLY_NUMBER_OVERRIDE" ]]; then
    aidx="$(resolve_index "$ASSEMBLY_NUMBER_OVERRIDE" "${#ASSEMBLIES[@]}")"
    if [[ "$aidx" == "ERROR" ]]; then
        echo "ERROR: assembly-number out of range" >&2
        exit 1
    fi
    ASSEMBLIES_RUN=("${ASSEMBLIES[$aidx]}")
else
    ASSEMBLIES_RUN=("${ASSEMBLIES[@]}")
fi

expand_template() {
    local tmpl="$1"
    local assembly="$2"
    local sample="$3"
    local sample_lower="$4"
    local set_suffix="$5"

    local out="$tmpl"
    out="${out//\{assembly\}/$assembly}"
    out="${out//\{sample\}/$sample}"
    out="${out//\{sample_lower\}/$sample_lower}"
    out="${out//\{set_suffix\}/$set_suffix}"
    printf '%s\n' "$out"
}

if [[ -n "${READ_SAMPLES+x}" && "${#READ_SAMPLES[@]}" -gt 0 ]]; then
    READ_SAMPLES_DEFAULT=("${READ_SAMPLES[@]}")
else
    READ_SAMPLES_DEFAULT=()
fi

mkdir -p "$OUTPUT_ROOT"
RUN_LOG="${OUTPUT_ROOT}/${EXPERIMENT_NAME}_${BINNER}_run_log.tsv"
printf 'experiment\tbinner\tassembly\ttarget_sample\tset_name\tassembly_fasta\tshort_reads\tlong_reads\tbin_fasta_dir\tbin_fasta_bins\tstatus\tnotes\n' > "$RUN_LOG"

for sample in "${TARGET_SAMPLES_RUN[@]}"; do
    sample_lower="$(printf '%s' "$sample" | tr '[:upper:]' '[:lower:]')"
    for assembly in "${ASSEMBLIES_RUN[@]}"; do
        set_name="${assembly}_${sample_lower}_${SET_SUFFIX}"
        assembly_fasta="$(expand_template "$ASSEMBLY_TEMPLATE" "$assembly" "$sample" "$sample_lower" "$SET_SUFFIX")"

        if [[ ! -f "$assembly_fasta" ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t%s\t%s\tskipped\tassembly_missing\n' \
                "$EXPERIMENT_NAME" "$BINNER" "$assembly" "$sample" "$set_name" "$assembly_fasta" "" "" >> "$RUN_LOG"
            continue
        fi

        read_samples=()
        if [[ "$COBINNING_MODE" == "1" ]]; then
            if [[ "${#READ_SAMPLES_DEFAULT[@]}" -eq 0 ]]; then
                echo "ERROR: COBINNING_MODE=1 requires READ_SAMPLES array in config" >&2
                exit 1
            fi
            read_samples=("${READ_SAMPLES_DEFAULT[@]}")
        else
            read_samples=("$sample")
        fi

        short_reads=()
        long_reads=()
        for read_sample in "${read_samples[@]}"; do
            read_sample_lower="$(printf '%s' "$read_sample" | tr '[:upper:]' '[:lower:]')"

            sr_path="$(expand_template "$SHORT_READ_TEMPLATE" "$assembly" "$read_sample" "$read_sample_lower" "$SET_SUFFIX")"
            if [[ -f "$sr_path" ]]; then
                short_reads+=("$sr_path")
            fi

            lr_path="$(expand_template "$LONG_READ_TEMPLATE" "$assembly" "$read_sample" "$read_sample_lower" "$SET_SUFFIX")"
            if [[ -f "$lr_path" ]]; then
                long_reads+=("$lr_path")
            fi
        done

        if [[ "${#short_reads[@]}" -eq 0 && "${#long_reads[@]}" -eq 0 ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t%s\t%s\tskipped\tno_reads_found\n' \
                "$EXPERIMENT_NAME" "$BINNER" "$assembly" "$sample" "$set_name" "$assembly_fasta" "" "" >> "$RUN_LOG"
            continue
        fi

        if [[ "$BINNER" == "metawrap" ]]; then
            cmd=(
                bash "$METAWRAP_BINNER_SCRIPT"
                --assembly-name "$assembly"
                --sample "$sample"
                --set-name "$set_name"
                --assembly-fasta "$assembly_fasta"
                --out-root "$OUTPUT_ROOT/metawrap"
                --container "$METAWRAP_CONTAINER"
                --cpus "$CPUS"
                --mem-gb "$MEM_GB"
            )
            for r in "${short_reads[@]}"; do cmd+=(--short-read "$r"); done
            for r in "${long_reads[@]}"; do cmd+=(--long-read "$r"); done
        else
            cmd=(
                bash "$VAMB_BINNER_SCRIPT"
                --assembly-name "$assembly"
                --sample "$sample"
                --set-name "$set_name"
                --assembly-fasta "$assembly_fasta"
                --out-root "$OUTPUT_ROOT/vamb"
                --qc-container "$QC_CONTAINER"
                --vamb-container "$VAMB_CONTAINER"
                --cpus "$CPUS"
                --min-contig-len "$MIN_CONTIG_LEN"
            )
            for r in "${short_reads[@]}"; do cmd+=(--short-read "$r"); done
            for r in "${long_reads[@]}"; do cmd+=(--long-read "$r"); done
        fi

        if "${cmd[@]}"; then
            if [[ "$BINNER" == "metawrap" ]]; then
                bin_fasta_dir="${OUTPUT_ROOT}/metawrap/${set_name}/bin_fasta"
            else
                bin_fasta_dir="${OUTPUT_ROOT}/vamb/${set_name}/bin_fasta"
            fi

            if [[ -d "$bin_fasta_dir" ]]; then
                bin_fasta_bins="$(find "$bin_fasta_dir" -maxdepth 1 \( -type f -o -type l \) -name '*.fa' | wc -l)"
            else
                bin_fasta_bins=0
            fi

            if [[ "$bin_fasta_bins" -gt 0 ]]; then
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tcompleted\tok\n' \
                    "$EXPERIMENT_NAME" "$BINNER" "$assembly" "$sample" "$set_name" "$assembly_fasta" \
                    "${#short_reads[@]}" "${#long_reads[@]}" "$bin_fasta_dir" "$bin_fasta_bins" >> "$RUN_LOG"
            else
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tfailed\tmissing_bin_fasta_dir\n' \
                    "$EXPERIMENT_NAME" "$BINNER" "$assembly" "$sample" "$set_name" "$assembly_fasta" \
                    "${#short_reads[@]}" "${#long_reads[@]}" "$bin_fasta_dir" "$bin_fasta_bins" >> "$RUN_LOG"
            fi
        else
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tfailed\tbinner_failed\n' \
                "$EXPERIMENT_NAME" "$BINNER" "$assembly" "$sample" "$set_name" "$assembly_fasta" \
                "${#short_reads[@]}" "${#long_reads[@]}" "" "0" >> "$RUN_LOG"
        fi
    done
done

echo "[$(date)] Wrapper complete"
echo "[$(date)] Run log: $RUN_LOG"

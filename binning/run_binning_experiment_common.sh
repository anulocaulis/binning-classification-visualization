#!/usr/bin/env bash
#SBATCH --job-name=binning_wrapper
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
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
    bash binning/run_binning_experiment_common.sh --config PATH [options]

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

resolve_assembly_template() {
    local assembly="$1"
    local key override_var

    key="$(printf '%s' "$assembly" | tr '[:lower:]-.' '[:upper:]__')"
    override_var="ASSEMBLY_TEMPLATE_${key}"

    if [[ -n "${!override_var:-}" ]]; then
        printf '%s\n' "${!override_var}"
    else
        printf '%s\n' "$ASSEMBLY_TEMPLATE"
    fi
}

resolve_read_sample_token() {
    local token="$1"
    local target_sample="$2"

    if [[ "$token" == "__TARGET_SAMPLE__" ]]; then
        printf '%s\n' "$target_sample"
    else
        printf '%s\n' "$token"
    fi
}

READ_SAMPLES_SHORT_DEFAULT=()
READ_SAMPLES_LONG_DEFAULT=()

# Backward compatibility: READ_SAMPLES applies to both read types unless
# per-type arrays are provided.
if [[ -n "${READ_SAMPLES+x}" && "${#READ_SAMPLES[@]}" -gt 0 ]]; then
    READ_SAMPLES_SHORT_DEFAULT=("${READ_SAMPLES[@]}")
    READ_SAMPLES_LONG_DEFAULT=("${READ_SAMPLES[@]}")
fi
if [[ -n "${READ_SAMPLES_SHORT+x}" && "${#READ_SAMPLES_SHORT[@]}" -gt 0 ]]; then
    READ_SAMPLES_SHORT_DEFAULT=("${READ_SAMPLES_SHORT[@]}")
fi
if [[ -n "${READ_SAMPLES_LONG+x}" && "${#READ_SAMPLES_LONG[@]}" -gt 0 ]]; then
    READ_SAMPLES_LONG_DEFAULT=("${READ_SAMPLES_LONG[@]}")
fi

mkdir -p "$OUTPUT_ROOT"
RUN_LOG="${OUTPUT_ROOT}/${EXPERIMENT_NAME}_${BINNER}_run_log.tsv"
printf 'experiment\tbinner\tassembly\ttarget_sample\tset_name\tassembly_fasta\tshort_reads\tlong_reads\tbin_fasta_dir\tbin_fasta_bins\tstatus\tnotes\n' > "$RUN_LOG"

for sample in "${TARGET_SAMPLES_RUN[@]}"; do
    sample_lower="$(printf '%s' "$sample" | tr '[:upper:]' '[:lower:]')"
    for assembly in "${ASSEMBLIES_RUN[@]}"; do
        set_name="${assembly}_${sample_lower}_${SET_SUFFIX}"
        assembly_template_use="$(resolve_assembly_template "$assembly")"
        assembly_fasta="$(expand_template "$assembly_template_use" "$assembly" "$sample" "$sample_lower" "$SET_SUFFIX")"

        if [[ ! -f "$assembly_fasta" ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t%s\t%s\tskipped\tassembly_missing\n' \
                "$EXPERIMENT_NAME" "$BINNER" "$assembly" "$sample" "$set_name" "$assembly_fasta" "" "" >> "$RUN_LOG"
            continue
        fi

        read_samples_short=()
        read_samples_long=()
        if [[ "$COBINNING_MODE" == "1" ]]; then
            if [[ "${#READ_SAMPLES_SHORT_DEFAULT[@]}" -eq 0 && "${#READ_SAMPLES_LONG_DEFAULT[@]}" -eq 0 ]]; then
                echo "ERROR: COBINNING_MODE=1 requires read sample arrays in config" >&2
                exit 1
            fi
            read_samples_short=("${READ_SAMPLES_SHORT_DEFAULT[@]}")
            read_samples_long=("${READ_SAMPLES_LONG_DEFAULT[@]}")
        else
            read_samples_short=("$sample")
            read_samples_long=("$sample")
        fi

        for i in "${!read_samples_short[@]}"; do
            read_samples_short[$i]="$(resolve_read_sample_token "${read_samples_short[$i]}" "$sample")"
        done
        for i in "${!read_samples_long[@]}"; do
            read_samples_long[$i]="$(resolve_read_sample_token "${read_samples_long[$i]}" "$sample")"
        done

        short_reads=()
        long_reads=()
        for read_sample in "${read_samples_short[@]}"; do
            read_sample_lower="$(printf '%s' "$read_sample" | tr '[:upper:]' '[:lower:]')"

            sr_path="$(expand_template "$SHORT_READ_TEMPLATE" "$assembly" "$read_sample" "$read_sample_lower" "$SET_SUFFIX")"
            if [[ -f "$sr_path" ]]; then
                short_reads+=("$sr_path")
            fi
        done

        for read_sample in "${read_samples_long[@]}"; do
            read_sample_lower="$(printf '%s' "$read_sample" | tr '[:upper:]' '[:lower:]')"

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
                sbatch
                "--job-name=${EXPERIMENT_NAME}_${set_name}"
                "--cpus-per-task=$CPUS"
                "--mem=${MEM_GB}G"
                "--time=4-00:00:00"
                "$METAWRAP_BINNER_SCRIPT"
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
                sbatch
                "--job-name=${EXPERIMENT_NAME}_${set_name}"
                "--cpus-per-task=$CPUS"
                "--mem=${MEM_GB}G"
                "--time=4-00:00:00"
                "$VAMB_BINNER_SCRIPT"
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

        job_output=$("${cmd[@]}" 2>&1)
        job_id=$(echo "$job_output" | grep -oP 'Submitted batch job \K[0-9]+' || echo "")

        if [[ -n "$job_id" ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tsubmitted\tjob_%s\n' \
                "$EXPERIMENT_NAME" "$BINNER" "$assembly" "$sample" "$set_name" "$assembly_fasta" \
                "${#short_reads[@]}" "${#long_reads[@]}" "" "" "$job_id" >> "$RUN_LOG"
        else
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tfailed\tsubmission_failed\n' \
                "$EXPERIMENT_NAME" "$BINNER" "$assembly" "$sample" "$set_name" "$assembly_fasta" \
                "${#short_reads[@]}" "${#long_reads[@]}" "" "0" >> "$RUN_LOG"
            echo "ERROR: sbatch submission failed for $set_name" >&2
            echo "Output: $job_output" >&2
        fi
    done
done

echo "[$(date)] Wrapper complete"
echo "[$(date)] Run log: $RUN_LOG"

#!/bin/bash
#SBATCH --job-name=filter_contigs
#SBATCH --output=logs/filter_contigs_%j.out
#SBATCH --error=logs/filter_contigs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=math-alderaan

# Usage:
#   sbatch filter_contigs_min_length.sh <SAMPLE> <ASSEMBLER>
# Example:
#   sbatch filter_contigs_min_length.sh S1 metaspades
# If no args are passed, process all existing assemblies with contigs files.

set -euo pipefail

echo "Starting filter_contigs_min_length at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper"
CONTAINER_PATH="${WORKDIR}/containers/qc_tools_miniconda.sif"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

cd "${WORKDIR}"

resolve_input_contigs() {
    local sample="$1"
    local assembler="$2"
    local base="${WORKDIR}/assemblies/${sample}/assembly.${assembler}"

    if [[ "${assembler}" == "metamdbg" ]]; then
        if [[ -s "${base}/metamdbg.contigs.fasta" ]]; then
            echo "${base}/metamdbg.contigs.fasta"
            return 0
        fi
        if [[ -s "${base}/contigs.fasta" ]]; then
            echo "${base}/contigs.fasta"
            return 0
        fi
        return 1
    fi

    if [[ "${assembler}" == "megahit" ]]; then
        if [[ -s "${base}/megahit.final.contigs.fa" ]]; then
            echo "${base}/megahit.final.contigs.fa"
            return 0
        fi
        if [[ -s "${base}/final.contigs.fa" ]]; then
            echo "${base}/final.contigs.fa"
            return 0
        fi
        return 1
    fi

    if [[ "${assembler}" == "idbaud" ]]; then
        if [[ -s "${base}/scaffold.fa" ]]; then
            echo "${base}/scaffold.fa"
            return 0
        fi
        return 1
    fi

    if [[ "${assembler}" == "flye" ]]; then
        if [[ -s "${base}/assembly.fasta" ]]; then
            echo "${base}/assembly.fasta"
            return 0
        fi
        if [[ -s "${base}/flye.assembly.fasta" ]]; then
            echo "${base}/flye.assembly.fasta"
            return 0
        fi
        return 1
    fi

    [[ -s "${base}/contigs.fasta" ]] && echo "${base}/contigs.fasta"
}

run_filter() {
    local sample="$1"
    local assembler="$2"
    local input_contigs

    input_contigs="$(resolve_input_contigs "${sample}" "${assembler}" || true)"
    if [[ -z "${input_contigs}" || ! -s "${input_contigs}" ]]; then
        echo "Skipping ${sample}/${assembler}: no contigs input found" >&2
        return 0
    fi

    local output_filtered="${WORKDIR}/assemblies/${sample}/assembly.${assembler}/contigs.ge1000.fa"

    if [[ -s "${output_filtered}" ]]; then
        echo "Skipping ${sample}/${assembler}: output already exists" >&2
        return 0
    fi

    echo ""
    echo "=== Filtering ${sample}/${assembler} at $(date) ==="
    echo "Input:  ${input_contigs}"
    echo "Output: ${output_filtered}"

    singularity exec -B "$PWD" "${CONTAINER_PATH}" reformat.sh \
        in="${input_contigs}" \
        out="${output_filtered}" \
        minlength=1000
}

if [[ $# -eq 2 ]]; then
    run_filter "$1" "$2"
elif [[ $# -eq 0 ]]; then
    mapfile -t targets < <(
        find "${WORKDIR}/assemblies" -maxdepth 3 -type d -name "assembly.*" \
            | sed -E 's#.*/assemblies/([^/]+)/assembly\.([^/]+)#\1 \2#' \
            | sort -u
    )

    if [[ ${#targets[@]} -eq 0 ]]; then
        echo "No assembly.* directories found under ${WORKDIR}/assemblies" >&2
        exit 1
    fi

    for pair in "${targets[@]}"; do
        sample="${pair%% *}"
        assembler="${pair##* }"
        run_filter "${sample}" "${assembler}"
    done
else
    echo "Usage: sbatch filter_contigs_min_length.sh [SAMPLE ASSEMBLER]" >&2
    exit 1
fi
echo "Finished at $(date)"

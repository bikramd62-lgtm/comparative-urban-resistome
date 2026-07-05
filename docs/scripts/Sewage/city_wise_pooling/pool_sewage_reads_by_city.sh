#!/usr/bin/env bash

# 01_pool_sewage_reads_by_city.sh
#
# Purpose:
# Pool human-depleted paired-end sewage reads by city.
#
# Workflow:
# 1. Search for human-depleted paired FASTQ files for each city/run.
# 2. Concatenate all R1 files for a city into one pooled R1 file.
# 3. Concatenate all R2 files for a city into one pooled R2 file.
#
# Expected input structure:
#   <input_dir>/<City>/<Run_accession>/*_paired_1.fastq.gz
#   <input_dir>/<City>/<Run_accession>/*_paired_2.fastq.gz
#
# Example:
#   data/processed_reads/sewage/human_depleted/Berlin/ERRxxxxxxx/ERRxxxxxxx_paired_1.fastq.gz
#   data/processed_reads/sewage/human_depleted/Berlin/ERRxxxxxxx/ERRxxxxxxx_paired_2.fastq.gz
#
# Output structure:
#   <output_dir>/<City>/<City>_pooled_1.fastq.gz
#   <output_dir>/<City>/<City>_pooled_2.fastq.gz
#
# Example run:
#   bash scripts/05_pooling/01_pool_sewage_reads_by_city.sh \
#     --input-dir data/processed_reads/sewage/human_depleted \
#     --output-dir data/processed_reads/sewage/pooled \
#     --logdir logs/05_pooling_sewage

set -euo pipefail

INPUT_DIR="data/processed_reads/sewage/human_depleted"
OUTPUT_DIR="data/processed_reads/sewage/pooled"
LOGDIR="logs/05_pooling_sewage"

usage() {
    echo ""
    echo "Usage:"
    echo "  bash scripts/05_pooling/01_pool_sewage_reads_by_city.sh \\"
    echo "    --input-dir <human_depleted_input_dir> \\"
    echo "    --output-dir <pooled_output_dir> \\"
    echo "    --logdir <log_output_dir>"
    echo ""
    echo "Optional:"
    echo "  --input-dir    Input directory containing human-depleted reads"
    echo "                 Default: data/processed_reads/sewage/human_depleted"
    echo "  --output-dir   Output directory for city-wise pooled reads"
    echo "                 Default: data/processed_reads/sewage/pooled"
    echo "  --logdir       Log directory"
    echo "                 Default: logs/05_pooling_sewage"
    echo ""
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --logdir)
            LOGDIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if ! command -v cat >/dev/null 2>&1; then
    echo "ERROR: cat is not available in PATH."
    exit 1
fi

if ! command -v find >/dev/null 2>&1; then
    echo "ERROR: find is not available in PATH."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOGDIR"

VERSION_FILE="${LOGDIR}/software_versions_pooling.txt"
MANIFEST="${LOGDIR}/sewage_city_pooling_manifest.tsv"
FAILED="${LOGDIR}/failed_pooling.tsv"

{
    echo "City-wise pooling script run date: $(date)"
    echo ""
    echo "bash:"
    bash --version | head -n 1
    echo ""
    echo "cat:"
    cat --version 2>/dev/null | head -n 1 || echo "cat version not available"
    echo ""
    echo "find:"
    find --version 2>/dev/null | head -n 1 || echo "find version not available"
} > "$VERSION_FILE"

echo -e "city\tn_runs_pooled\tn_r1_files\tn_r2_files\tpooled_r1\tpooled_r2\tstatus" > "$MANIFEST"
echo -e "city\trun_accession\treason" > "$FAILED"

sanitize_name() {
    echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

find_human_depleted_pair() {
    local run_dir="$1"

    local r1=""
    local r2=""

    r1=$(find "$run_dir" -maxdepth 1 -type f \( \
        -name "*_paired_1.fastq.gz" -o \
        -name "*_paired_1.fq.gz" -o \
        -name "*_clean_1.fastq.gz" -o \
        -name "*_clean_1.fq.gz" -o \
        -name "*_decontaminated_1.fastq.gz" -o \
        -name "*_decontaminated_1.fq.gz" \
    \) | sort | head -n 1)

    r2=$(find "$run_dir" -maxdepth 1 -type f \( \
        -name "*_paired_2.fastq.gz" -o \
        -name "*_paired_2.fq.gz" -o \
        -name "*_clean_2.fastq.gz" -o \
        -name "*_clean_2.fq.gz" -o \
        -name "*_decontaminated_2.fastq.gz" -o \
        -name "*_decontaminated_2.fq.gz" \
    \) | sort | head -n 1)

    echo -e "${r1}\t${r2}"
}

echo ""
echo "Starting city-wise pooling of sewage reads"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $LOGDIR"
echo ""

while IFS= read -r CITY_DIR; do

    CITY=$(basename "$CITY_DIR")
    CITY_SAFE=$(sanitize_name "$CITY")

    echo "------------------------------------------------------------"
    echo "City: $CITY"
    echo "Input: $CITY_DIR"

    CITY_OUTPUT_DIR="${OUTPUT_DIR}/${CITY_SAFE}"
    mkdir -p "$CITY_OUTPUT_DIR"

    R1_LIST="${LOGDIR}/${CITY_SAFE}_R1_files.txt"
    R2_LIST="${LOGDIR}/${CITY_SAFE}_R2_files.txt"

    : > "$R1_LIST"
    : > "$R2_LIST"

    COMPLETE_RUNS=0

    while IFS= read -r RUN_DIR; do

        RUN_ACCESSION=$(basename "$RUN_DIR")

        PAIR_INFO=$(find_human_depleted_pair "$RUN_DIR")
        R1_FILE=$(echo "$PAIR_INFO" | cut -f1)
        R2_FILE=$(echo "$PAIR_INFO" | cut -f2)

        if [[ -z "$R1_FILE" || -z "$R2_FILE" ]]; then
            echo "WARNING: Missing paired human-depleted FASTQ for ${CITY}/${RUN_ACCESSION}"
            echo -e "${CITY}\t${RUN_ACCESSION}\thuman_depleted_pair_not_found" >> "$FAILED"
            continue
        fi

        echo "$R1_FILE" >> "$R1_LIST"
        echo "$R2_FILE" >> "$R2_LIST"

        COMPLETE_RUNS=$((COMPLETE_RUNS + 1))

    done < <(find "$CITY_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

    N_R1=$(wc -l < "$R1_LIST" | tr -d ' ')
    N_R2=$(wc -l < "$R2_LIST" | tr -d ' ')

    if [[ "$N_R1" -eq 0 || "$N_R2" -eq 0 ]]; then
        echo "WARNING: No complete read pairs found for city: $CITY"
        echo -e "${CITY}\tNA\tno_complete_pairs_for_city" >> "$FAILED"
        echo -e "${CITY}\t0\t${N_R1}\t${N_R2}\tNA\tNA\tfailed_no_pairs" >> "$MANIFEST"
        continue
    fi

    if [[ "$N_R1" -ne "$N_R2" ]]; then
        echo "ERROR: Unequal number of R1 and R2 files for city: $CITY"
        echo -e "${CITY}\tNA\tunequal_r1_r2_file_counts" >> "$FAILED"
        echo -e "${CITY}\t${COMPLETE_RUNS}\t${N_R1}\t${N_R2}\tNA\tNA\tfailed_unequal_pairs" >> "$MANIFEST"
        continue
    fi

    POOLED_R1="${CITY_OUTPUT_DIR}/${CITY_SAFE}_pooled_1.fastq.gz"
    POOLED_R2="${CITY_OUTPUT_DIR}/${CITY_SAFE}_pooled_2.fastq.gz"

    echo "Number of complete runs to pool: $COMPLETE_RUNS"
    echo "Pooling R1 files into: $POOLED_R1"
    cat $(cat "$R1_LIST") > "$POOLED_R1"

    echo "Pooling R2 files into: $POOLED_R2"
    cat $(cat "$R2_LIST") > "$POOLED_R2"

    echo -e "${CITY}\t${COMPLETE_RUNS}\t${N_R1}\t${N_R2}\t${POOLED_R1}\t${POOLED_R2}\tcompleted" >> "$MANIFEST"

done < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

echo ""
echo "City-wise pooling complete."
echo "Pooled reads: $OUTPUT_DIR"
echo "Manifest: $MANIFEST"
echo "Failed pooling file: $FAILED"
echo "Software versions: $VERSION_FILE"
echo ""

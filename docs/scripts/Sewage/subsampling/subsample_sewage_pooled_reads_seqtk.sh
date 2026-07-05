#!/usr/bin/env bash

# 01_subsample_sewage_pooled_reads_seqtk.sh
#
# Purpose:
# Subsample city-wise pooled sewage paired-end reads to 5 million read pairs
# per city using seqtk.
#
# Workflow:
# 1. Read city-wise pooled R1 and R2 FASTQ files.
# 2. Subsample each mate file to the same number of reads using the same seed.
# 3. Compress the subsampled FASTQ output files.
#
# Expected input structure:
#   <input_dir>/<City>/<City>_pooled_1.fastq.gz
#   <input_dir>/<City>/<City>_pooled_2.fastq.gz
#
# Output structure:
#   <output_dir>/<City>/<City>_subsampled_5M_1.fastq.gz
#   <output_dir>/<City>/<City>_subsampled_5M_2.fastq.gz
#
# Example run:
#   bash scripts/06_subsampling/01_subsample_sewage_pooled_reads_seqtk.sh \
#     --input-dir data/processed_reads/sewage/pooled \
#     --output-dir data/processed_reads/sewage/subsampled_5M \
#     --logdir logs/06_subsampling_sewage \
#     --reads 5000000 \
#     --seed 42

set -euo pipefail

INPUT_DIR="data/processed_reads/sewage/pooled"
OUTPUT_DIR="data/processed_reads/sewage/subsampled_5M"
LOGDIR="logs/06_subsampling_sewage"
N_READS=5000000
SEED=42
CHECK_COUNTS=1

usage() {
    echo ""
    echo "Usage:"
    echo "  bash scripts/06_subsampling/01_subsample_sewage_pooled_reads_seqtk.sh \\"
    echo "    --input-dir <pooled_reads_dir> \\"
    echo "    --output-dir <subsampled_output_dir> \\"
    echo "    --logdir <log_output_dir> \\"
    echo "    --reads 5000000 \\"
    echo "    --seed 42"
    echo ""
    echo "Optional:"
    echo "  --input-dir          Input directory containing city-wise pooled reads"
    echo "                       Default: data/processed_reads/sewage/pooled"
    echo "  --output-dir         Output directory for subsampled reads"
    echo "                       Default: data/processed_reads/sewage/subsampled_5M"
    echo "  --logdir             Log directory"
    echo "                       Default: logs/06_subsampling_sewage"
    echo "  --reads              Number of reads to sample from each mate file"
    echo "                       Default: 5000000"
    echo "  --seed               Random seed for seqtk"
    echo "                       Default: 42"
    echo "  --skip-count-check   Skip read-count checking before and after subsampling"
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
        --reads)
            N_READS="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --skip-count-check)
            CHECK_COUNTS=0
            shift 1
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

if ! command -v seqtk >/dev/null 2>&1; then
    echo "ERROR: seqtk is not installed or not available in PATH."
    exit 1
fi

if ! command -v gzip >/dev/null 2>&1; then
    echo "ERROR: gzip is not installed or not available in PATH."
    exit 1
fi

if [[ "$CHECK_COUNTS" -eq 1 ]]; then
    if ! command -v zcat >/dev/null 2>&1; then
        echo "ERROR: zcat is required for read-count checking."
        echo "Run again with --skip-count-check if you want to skip count checks."
        exit 1
    fi
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOGDIR"

VERSION_FILE="${LOGDIR}/software_versions_subsampling.txt"
MANIFEST="${LOGDIR}/sewage_subsampling_manifest.tsv"
FAILED="${LOGDIR}/failed_subsampling.tsv"

{
    echo "Subsampling script run date: $(date)"
    echo ""
    echo "seqtk:"
    seqtk 2>&1 | head -n 1 || true
    echo ""
    echo "gzip:"
    gzip --version | head -n 1
    echo ""
    echo "bash:"
    bash --version | head -n 1
    echo ""
    echo "Subsampling settings:"
    echo "reads_per_mate=${N_READS}"
    echo "seed=${SEED}"
} > "$VERSION_FILE"

echo -e "city\tpooled_r1\tpooled_r2\tinput_r1_reads\tinput_r2_reads\toutput_r1_reads\toutput_r2_reads\tsubsampled_r1\tsubsampled_r2\treads_requested\tseed\tstatus" > "$MANIFEST"
echo -e "city\treason" > "$FAILED"

sanitize_name() {
    echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

count_reads_fastq_gz() {
    local fastq_file="$1"
    zcat -f "$fastq_file" | awk 'END {print NR / 4}'
}

find_pooled_pair() {
    local city_dir="$1"
    local city_name="$2"

    local r1=""
    local r2=""

    r1=$(find "$city_dir" -maxdepth 1 -type f \( \
        -name "${city_name}_pooled_1.fastq.gz" -o \
        -name "${city_name}_pooled_1.fq.gz" -o \
        -name "*_pooled_1.fastq.gz" -o \
        -name "*_pooled_1.fq.gz" -o \
        -name "*_1.fastq.gz" -o \
        -name "*_1.fq.gz" \
    \) | sort | head -n 1)

    r2=$(find "$city_dir" -maxdepth 1 -type f \( \
        -name "${city_name}_pooled_2.fastq.gz" -o \
        -name "${city_name}_pooled_2.fq.gz" -o \
        -name "*_pooled_2.fastq.gz" -o \
        -name "*_pooled_2.fq.gz" -o \
        -name "*_2.fastq.gz" -o \
        -name "*_2.fq.gz" \
    \) | sort | head -n 1)

    echo -e "${r1}\t${r2}"
}

echo ""
echo "Starting sewage read subsampling"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $LOGDIR"
echo "Reads per mate: $N_READS"
echo "Seed: $SEED"
echo ""

while IFS= read -r CITY_DIR; do

    CITY=$(basename "$CITY_DIR")
    CITY_SAFE=$(sanitize_name "$CITY")

    echo "------------------------------------------------------------"
    echo "City: $CITY"
    echo "Input: $CITY_DIR"

    PAIR_INFO=$(find_pooled_pair "$CITY_DIR" "$CITY")
    POOLED_R1=$(echo "$PAIR_INFO" | cut -f1)
    POOLED_R2=$(echo "$PAIR_INFO" | cut -f2)

    if [[ -z "$POOLED_R1" || -z "$POOLED_R2" ]]; then
        echo "WARNING: Could not find pooled paired FASTQ files for city: $CITY"
        echo -e "${CITY}\tpooled_pair_not_found" >> "$FAILED"
        continue
    fi

    CITY_OUTPUT_DIR="${OUTPUT_DIR}/${CITY_SAFE}"
    mkdir -p "$CITY_OUTPUT_DIR"

    OUT_R1="${CITY_OUTPUT_DIR}/${CITY_SAFE}_subsampled_5M_1.fastq.gz"
    OUT_R2="${CITY_OUTPUT_DIR}/${CITY_SAFE}_subsampled_5M_2.fastq.gz"

    INPUT_R1_READS="not_checked"
    INPUT_R2_READS="not_checked"
    OUTPUT_R1_READS="not_checked"
    OUTPUT_R2_READS="not_checked"

    echo "Pooled R1: $POOLED_R1"
    echo "Pooled R2: $POOLED_R2"

    if [[ "$CHECK_COUNTS" -eq 1 ]]; then
        echo "Counting input reads..."
        INPUT_R1_READS=$(count_reads_fastq_gz "$POOLED_R1")
        INPUT_R2_READS=$(count_reads_fastq_gz "$POOLED_R2")

        echo "Input R1 reads: $INPUT_R1_READS"
        echo "Input R2 reads: $INPUT_R2_READS"

        if [[ "$INPUT_R1_READS" != "$INPUT_R2_READS" ]]; then
            echo "WARNING: Unequal input R1/R2 read counts for city: $CITY"
            echo -e "${CITY}\tunequal_input_r1_r2_counts" >> "$FAILED"
            echo -e "${CITY}\t${POOLED_R1}\t${POOLED_R2}\t${INPUT_R1_READS}\t${INPUT_R2_READS}\tNA\tNA\tNA\tNA\t${N_READS}\t${SEED}\tfailed_unequal_input_counts" >> "$MANIFEST"
            continue
        fi

        if [[ "$INPUT_R1_READS" -lt "$N_READS" ]]; then
            echo "WARNING: City has fewer reads than requested subsampling depth: $CITY"
            echo -e "${CITY}\tinput_reads_less_than_requested_depth" >> "$FAILED"
            echo -e "${CITY}\t${POOLED_R1}\t${POOLED_R2}\t${INPUT_R1_READS}\t${INPUT_R2_READS}\tNA\tNA\tNA\tNA\t${N_READS}\t${SEED}\tfailed_insufficient_reads" >> "$MANIFEST"
            continue
        fi
    fi

    echo "Subsampling R1 to ${N_READS} reads using seed ${SEED}..."
    seqtk sample -s"${SEED}" "$POOLED_R1" "$N_READS" | gzip -c > "$OUT_R1"

    echo "Subsampling R2 to ${N_READS} reads using seed ${SEED}..."
    seqtk sample -s"${SEED}" "$POOLED_R2" "$N_READS" | gzip -c > "$OUT_R2"

    if [[ "$CHECK_COUNTS" -eq 1 ]]; then
        echo "Counting output reads..."
        OUTPUT_R1_READS=$(count_reads_fastq_gz "$OUT_R1")
        OUTPUT_R2_READS=$(count_reads_fastq_gz "$OUT_R2")

        echo "Output R1 reads: $OUTPUT_R1_READS"
        echo "Output R2 reads: $OUTPUT_R2_READS"

        if [[ "$OUTPUT_R1_READS" -ne "$N_READS" || "$OUTPUT_R2_READS" -ne "$N_READS" ]]; then
            echo "WARNING: Output read count does not match requested depth for city: $CITY"
            echo -e "${CITY}\toutput_read_count_not_equal_requested_depth" >> "$FAILED"
            echo -e "${CITY}\t${POOLED_R1}\t${POOLED_R2}\t${INPUT_R1_READS}\t${INPUT_R2_READS}\t${OUTPUT_R1_READS}\t${OUTPUT_R2_READS}\t${OUT_R1}\t${OUT_R2}\t${N_READS}\t${SEED}\tfailed_output_count_check" >> "$MANIFEST"
            continue
        fi
    fi

    echo -e "${CITY}\t${POOLED_R1}\t${POOLED_R2}\t${INPUT_R1_READS}\t${INPUT_R2_READS}\t${OUTPUT_R1_READS}\t${OUTPUT_R2_READS}\t${OUT_R1}\t${OUT_R2}\t${N_READS}\t${SEED}\tcompleted" >> "$MANIFEST"

done < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

echo ""
echo "Sewage read subsampling complete."
echo "Subsampled reads: $OUTPUT_DIR"
echo "Manifest: $MANIFEST"
echo "Failed subsampling file: $FAILED"
echo "Software versions: $VERSION_FILE"
echo ""

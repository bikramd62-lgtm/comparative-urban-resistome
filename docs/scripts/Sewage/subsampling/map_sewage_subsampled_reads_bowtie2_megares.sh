#!/usr/bin/env bash

# 01_map_sewage_subsampled_reads_bowtie2_megares.sh
#
# Purpose:
# Map 5M subsampled sewage paired-end reads against MEGARes v3.00
# using Bowtie2 with the --very-sensitive setting.
#
# Workflow:
# 1. Read city-wise subsampled paired-end FASTQ files.
# 2. Map each city against the MEGARes v3.00 Bowtie2 index.
# 3. Save alignments as SAM files.
# 4. Save Bowtie2 mapping logs and a run manifest.
#
# Expected input structure:
#   <input_dir>/<City>/<City>_subsampled_5M_1.fastq.gz
#   <input_dir>/<City>/<City>_subsampled_5M_2.fastq.gz
#
# Output structure:
#   <output_dir>/<City>/<City>_vs_MEGARes_v3.sam
#
# Example run:
#   bash scripts/07_amr_mapping/01_map_sewage_subsampled_reads_bowtie2_megares.sh \
#     --input-dir data/processed_reads/sewage/subsampled_5M \
#     --megares-index path/to/megares_v3 \
#     --output-dir results/amr_mapping/sewage/megares_v3/sam \
#     --logdir logs/07_amr_mapping_sewage \
#     --threads 4

set -euo pipefail

INPUT_DIR="data/processed_reads/sewage/subsampled_5M"
MEGARES_INDEX=""
OUTPUT_DIR="results/amr_mapping/sewage/megares_v3/sam"
LOGDIR="logs/07_amr_mapping_sewage"
THREADS=4

usage() {
    echo ""
    echo "Usage:"
    echo "  bash scripts/07_amr_mapping/01_map_sewage_subsampled_reads_bowtie2_megares.sh \\"
    echo "    --input-dir <subsampled_reads_dir> \\"
    echo "    --megares-index <bowtie2_index_prefix> \\"
    echo "    --output-dir <sam_output_dir> \\"
    echo "    --logdir <log_output_dir> \\"
    echo "    --threads <number_of_threads>"
    echo ""
    echo "Required:"
    echo "  --megares-index   Bowtie2 index prefix for MEGARes v3.00"
    echo "                    Example: database/megares_v3/index/megares_v3"
    echo ""
    echo "Optional:"
    echo "  --input-dir       Input directory containing 5M subsampled sewage reads"
    echo "                    Default: data/processed_reads/sewage/subsampled_5M"
    echo "  --output-dir      Output directory for SAM files"
    echo "                    Default: results/amr_mapping/sewage/megares_v3/sam"
    echo "  --logdir          Directory for Bowtie2 logs and manifest files"
    echo "                    Default: logs/07_amr_mapping_sewage"
    echo "  --threads         Number of Bowtie2 threads"
    echo "                    Default: 4"
    echo ""
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --megares-index)
            MEGARES_INDEX="$2"
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
        --threads)
            THREADS="$2"
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

if [[ -z "$MEGARES_INDEX" ]]; then
    echo "ERROR: --megares-index is required."
    usage
    exit 1
fi

if ! command -v bowtie2 >/dev/null 2>&1; then
    echo "ERROR: Bowtie2 is not installed or not available in PATH."
    exit 1
fi

if ! command -v bowtie2-inspect >/dev/null 2>&1; then
    echo "ERROR: bowtie2-inspect is not installed or not available in PATH."
    exit 1
fi

if ! bowtie2-inspect -s "$MEGARES_INDEX" >/dev/null 2>&1; then
    echo "ERROR: Bowtie2 index could not be inspected:"
    echo "$MEGARES_INDEX"
    echo ""
    echo "Make sure --megares-index points to the index prefix, not only the folder."
    echo "Example:"
    echo "  --megares-index database/megares_v3/index/megares_v3"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOGDIR"

VERSION_FILE="${LOGDIR}/software_versions_bowtie2_mapping.txt"
MANIFEST="${LOGDIR}/sewage_megares_mapping_manifest.tsv"
FAILED="${LOGDIR}/failed_mapping.tsv"

{
    echo "Bowtie2 MEGARes mapping script run date: $(date)"
    echo ""
    echo "Bowtie2:"
    bowtie2 --version | head -n 1
    echo ""
    echo "Bowtie2 index prefix:"
    echo "$MEGARES_INDEX"
    echo ""
    echo "Mapping settings:"
    echo "mode=--very-sensitive"
    echo "threads=${THREADS}"
} > "$VERSION_FILE"

echo -e "city\tsubsampled_r1\tsubsampled_r2\tsam_file\tbowtie2_log\tthreads\tsetting\tstatus" > "$MANIFEST"
echo -e "city\treason" > "$FAILED"

sanitize_name() {
    echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

find_subsampled_pair() {
    local city_dir="$1"
    local city_name="$2"

    local r1=""
    local r2=""

    r1=$(find "$city_dir" -maxdepth 1 -type f \( \
        -name "${city_name}_subsampled_5M_1.fastq.gz" -o \
        -name "${city_name}_subsampled_5M_1.fq.gz" -o \
        -name "*_subsampled_5M_1.fastq.gz" -o \
        -name "*_subsampled_5M_1.fq.gz" -o \
        -name "*_1.fastq.gz" -o \
        -name "*_1.fq.gz" \
    \) | sort | head -n 1)

    r2=$(find "$city_dir" -maxdepth 1 -type f \( \
        -name "${city_name}_subsampled_5M_2.fastq.gz" -o \
        -name "${city_name}_subsampled_5M_2.fq.gz" -o \
        -name "*_subsampled_5M_2.fastq.gz" -o \
        -name "*_subsampled_5M_2.fq.gz" -o \
        -name "*_2.fastq.gz" -o \
        -name "*_2.fq.gz" \
    \) | sort | head -n 1)

    echo -e "${r1}\t${r2}"
}

echo ""
echo "Starting Bowtie2 mapping of sewage reads against MEGARes v3.00"
echo "Input directory: $INPUT_DIR"
echo "MEGARes Bowtie2 index prefix: $MEGARES_INDEX"
echo "SAM output directory: $OUTPUT_DIR"
echo "Log directory: $LOGDIR"
echo "Threads: $THREADS"
echo "Bowtie2 setting: --very-sensitive"
echo ""

while IFS= read -r CITY_DIR; do

    CITY=$(basename "$CITY_DIR")
    CITY_SAFE=$(sanitize_name "$CITY")

    echo "------------------------------------------------------------"
    echo "City: $CITY"
    echo "Input: $CITY_DIR"

    PAIR_INFO=$(find_subsampled_pair "$CITY_DIR" "$CITY")
    SUB_R1=$(echo "$PAIR_INFO" | cut -f1)
    SUB_R2=$(echo "$PAIR_INFO" | cut -f2)

    if [[ -z "$SUB_R1" || -z "$SUB_R2" ]]; then
        echo "WARNING: Could not find subsampled paired FASTQ files for city: $CITY"
        echo -e "${CITY}\tsubsampled_pair_not_found" >> "$FAILED"
        continue
    fi

    CITY_OUTPUT_DIR="${OUTPUT_DIR}/${CITY_SAFE}"
    mkdir -p "$CITY_OUTPUT_DIR"

    SAM_FILE="${CITY_OUTPUT_DIR}/${CITY_SAFE}_vs_MEGARes_v3.sam"
    BOWTIE2_LOG="${LOGDIR}/${CITY_SAFE}_bowtie2_megares_v3.log"

    echo "Subsampled R1: $SUB_R1"
    echo "Subsampled R2: $SUB_R2"
    echo "SAM output: $SAM_FILE"
    echo "Bowtie2 log: $BOWTIE2_LOG"

    bowtie2 \
        --very-sensitive \
        -x "$MEGARES_INDEX" \
        -1 "$SUB_R1" \
        -2 "$SUB_R2" \
        -S "$SAM_FILE" \
        -p "$THREADS" \
        2> "$BOWTIE2_LOG"

    if [[ ! -s "$SAM_FILE" ]]; then
        echo "WARNING: SAM file was not created or is empty for city: $CITY"
        echo -e "${CITY}\tsam_file_missing_or_empty" >> "$FAILED"
        echo -e "${CITY}\t${SUB_R1}\t${SUB_R2}\t${SAM_FILE}\t${BOWTIE2_LOG}\t${THREADS}\t--very-sensitive\tfailed_empty_sam" >> "$MANIFEST"
        continue
    fi

    echo -e "${CITY}\t${SUB_R1}\t${SUB_R2}\t${SAM_FILE}\t${BOWTIE2_LOG}\t${THREADS}\t--very-sensitive\tcompleted" >> "$MANIFEST"

done < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

echo ""
echo "Bowtie2 MEGARes mapping complete."
echo "SAM files: $OUTPUT_DIR"
echo "Mapping manifest: $MANIFEST"
echo "Failed mapping file: $FAILED"
echo "Software versions: $VERSION_FILE"
echo ""

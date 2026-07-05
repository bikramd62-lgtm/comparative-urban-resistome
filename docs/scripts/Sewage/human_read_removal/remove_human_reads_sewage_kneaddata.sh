#!/usr/bin/env bash

# 01_remove_human_reads_sewage_kneaddata.sh
#
# Purpose:
# Remove human-derived reads from trimmed sewage paired-end FASTQ files
# using KneadData with Bowtie2.
#
# Workflow:
# 1. Read trimmed paired-end FASTQ files.
# 2. Screen reads against a human Bowtie2 database using KneadData.
# 3. Retain non-human paired reads.
# 4. Compress final paired non-human FASTQ files.
#
# Expected input structure:
#   <trimmed_dir>/<City>/<Run_accession>/<Run_accession>_trimmed_1.fastq.gz
#   <trimmed_dir>/<City>/<Run_accession>/<Run_accession>_trimmed_2.fastq.gz
#
# Example run:
#   bash scripts/04_human_removal/01_remove_human_reads_sewage_kneaddata.sh \
#     --trimmed-dir data/processed_reads/sewage/trimmed \
#     --human-db /path/to/Homo_sapiens_hg39_T2T_Bowtie2_v0.1 \
#     --output-dir data/processed_reads/sewage/human_depleted \
#     --logdir logs/04_human_removal_sewage \
#     --threads 4

set -euo pipefail

TRIMMED_DIR="data/processed_reads/sewage/trimmed"
HUMAN_DB=""
OUTPUT_DIR="data/processed_reads/sewage/human_depleted"
LOGDIR="logs/04_human_removal_sewage"
THREADS=4

usage() {
    echo ""
    echo "Usage:"
    echo "  bash scripts/04_human_removal/01_remove_human_reads_sewage_kneaddata.sh \\"
    echo "    --trimmed-dir <trimmed_fastq_dir> \\"
    echo "    --human-db <kneaddata_human_bowtie2_db_dir> \\"
    echo "    --output-dir <human_depleted_output_dir> \\"
    echo "    --logdir <log_output_dir> \\"
    echo "    --threads <number_of_threads>"
    echo ""
    echo "Required:"
    echo "  --human-db     Path to KneadData/Bowtie2 human reference database"
    echo ""
    echo "Optional:"
    echo "  --trimmed-dir  Input directory containing fastp-trimmed sewage reads"
    echo "                 Default: data/processed_reads/sewage/trimmed"
    echo "  --output-dir   Output directory for human-depleted reads"
    echo "                 Default: data/processed_reads/sewage/human_depleted"
    echo "  --logdir       Log directory"
    echo "                 Default: logs/04_human_removal_sewage"
    echo "  --threads      Number of threads"
    echo "                 Default: 4"
    echo ""
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --trimmed-dir)
            TRIMMED_DIR="$2"
            shift 2
            ;;
        --human-db)
            HUMAN_DB="$2"
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

if [[ ! -d "$TRIMMED_DIR" ]]; then
    echo "ERROR: Trimmed read directory not found: $TRIMMED_DIR"
    exit 1
fi

if [[ -z "$HUMAN_DB" ]]; then
    echo "ERROR: --human-db is required."
    usage
    exit 1
fi

if [[ ! -d "$HUMAN_DB" ]]; then
    echo "ERROR: Human database directory not found: $HUMAN_DB"
    exit 1
fi

if ! command -v kneaddata >/dev/null 2>&1; then
    echo "ERROR: KneadData is not installed or not available in PATH."
    exit 1
fi

if ! command -v bowtie2 >/dev/null 2>&1; then
    echo "ERROR: Bowtie2 is not installed or not available in PATH."
    exit 1
fi

if ! command -v gzip >/dev/null 2>&1; then
    echo "ERROR: gzip is not installed or not available in PATH."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOGDIR"

VERSION_FILE="${LOGDIR}/software_versions_human_removal.txt"
MANIFEST="${LOGDIR}/sewage_human_removal_manifest.tsv"
FAILED="${LOGDIR}/failed_human_removal.tsv"

{
    echo "Human read removal script run date: $(date)"
    echo ""
    echo "KneadData:"
    kneaddata --version 2>&1 || true
    echo ""
    echo "Bowtie2:"
    bowtie2 --version | head -n 1
    echo ""
    echo "gzip:"
    gzip --version | head -n 1
    echo ""
    echo "Human database:"
    echo "$HUMAN_DB"
} > "$VERSION_FILE"

echo -e "city\trun_accession\ttrimmed_r1\ttrimmed_r2\thuman_depleted_r1\thuman_depleted_r2\tkneaddata_output_dir\tstatus" > "$MANIFEST"
echo -e "city\trun_accession\treason" > "$FAILED"

sanitize_name() {
    echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

find_trimmed_pair() {
    local run_dir="$1"
    local run_accession="$2"

    local r1=""
    local r2=""

    r1=$(find "$run_dir" -maxdepth 1 -type f \( \
        -name "${run_accession}_trimmed_1.fastq.gz" -o \
        -name "${run_accession}_trimmed_1.fq.gz" -o \
        -name "*_trimmed_1.fastq.gz" -o \
        -name "*_trimmed_1.fq.gz" -o \
        -name "*_1.fastq.gz" -o \
        -name "*_1.fq.gz" -o \
        -name "*_R1.fastq.gz" -o \
        -name "*_R1.fq.gz" \
    \) | sort | head -n 1)

    r2=$(find "$run_dir" -maxdepth 1 -type f \( \
        -name "${run_accession}_trimmed_2.fastq.gz" -o \
        -name "${run_accession}_trimmed_2.fq.gz" -o \
        -name "*_trimmed_2.fastq.gz" -o \
        -name "*_trimmed_2.fq.gz" -o \
        -name "*_2.fastq.gz" -o \
        -name "*_2.fq.gz" -o \
        -name "*_R2.fastq.gz" -o \
        -name "*_R2.fq.gz" \
    \) | sort | head -n 1)

    echo -e "${r1}\t${r2}"
}

echo ""
echo "Starting human read removal for sewage reads"
echo "Trimmed input directory: $TRIMMED_DIR"
echo "Human database: $HUMAN_DB"
echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $LOGDIR"
echo "Threads: $THREADS"
echo ""

while IFS= read -r RUN_DIR; do

    RUN_ACCESSION=$(basename "$RUN_DIR")
    CITY=$(basename "$(dirname "$RUN_DIR")")
    CITY_SAFE=$(sanitize_name "$CITY")

    echo "------------------------------------------------------------"
    echo "City: $CITY"
    echo "Run:  $RUN_ACCESSION"
    echo "Input: $RUN_DIR"

    PAIR_INFO=$(find_trimmed_pair "$RUN_DIR" "$RUN_ACCESSION")
    TRIMMED_R1=$(echo "$PAIR_INFO" | cut -f1)
    TRIMMED_R2=$(echo "$PAIR_INFO" | cut -f2)

    if [[ -z "$TRIMMED_R1" || -z "$TRIMMED_R2" ]]; then
        echo "WARNING: Could not find trimmed paired FASTQ files for $RUN_ACCESSION"
        echo -e "${CITY}\t${RUN_ACCESSION}\ttrimmed_pair_not_found" >> "$FAILED"
        continue
    fi

    RUN_OUTPUT_DIR="${OUTPUT_DIR}/${CITY_SAFE}/${RUN_ACCESSION}"
    mkdir -p "$RUN_OUTPUT_DIR"

    echo "Trimmed R1: $TRIMMED_R1"
    echo "Trimmed R2: $TRIMMED_R2"
    echo "KneadData output: $RUN_OUTPUT_DIR"

    kneaddata \
        --input1 "$TRIMMED_R1" \
        --input2 "$TRIMMED_R2" \
        --reference-db "$HUMAN_DB" \
        --output "$RUN_OUTPUT_DIR" \
        --output-prefix "$RUN_ACCESSION" \
        --threads "$THREADS" \
        --bypass-trim \
        --bypass-trf \
        --decontaminate-pairs strict \
        --remove-intermediate-output

    PAIRED_R1=$(find "$RUN_OUTPUT_DIR" -maxdepth 1 -type f -name "*_paired_1.fastq" | sort | head -n 1)
    PAIRED_R2=$(find "$RUN_OUTPUT_DIR" -maxdepth 1 -type f -name "*_paired_2.fastq" | sort | head -n 1)

    if [[ -z "$PAIRED_R1" || -z "$PAIRED_R2" ]]; then
        echo "WARNING: KneadData paired output files not found for $RUN_ACCESSION"
        echo -e "${CITY}\t${RUN_ACCESSION}\tkneaddata_paired_output_not_found" >> "$FAILED"
        continue
    fi

    echo "Compressing human-depleted paired reads..."

    gzip -f "$PAIRED_R1"
    gzip -f "$PAIRED_R2"

    HUMAN_DEPLETED_R1="${PAIRED_R1}.gz"
    HUMAN_DEPLETED_R2="${PAIRED_R2}.gz"

    echo -e "${CITY}\t${RUN_ACCESSION}\t${TRIMMED_R1}\t${TRIMMED_R2}\t${HUMAN_DEPLETED_R1}\t${HUMAN_DEPLETED_R2}\t${RUN_OUTPUT_DIR}\tcompleted" >> "$MANIFEST"

done < <(find "$TRIMMED_DIR" -mindepth 2 -maxdepth 2 -type d | sort)

echo ""
echo "Human read removal step complete."
echo "Human-depleted reads: $OUTPUT_DIR"
echo "Manifest: $MANIFEST"
echo "Failed runs: $FAILED"
echo "Software versions: $VERSION_FILE"
echo ""

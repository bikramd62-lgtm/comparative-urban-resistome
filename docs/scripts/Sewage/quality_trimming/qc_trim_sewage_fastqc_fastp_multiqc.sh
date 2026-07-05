#!/usr/bin/env bash

# qc_trim_sewage_fastqc_fastp_multiqc.sh
#
# Purpose:
# Perform read quality assessment and trimming for sewage paired-end FASTQ files.
#
# Workflow:
# 1. Run FastQC on raw downloaded reads.
# 2. Trim reads using fastp.
# 3. Run FastQC on trimmed reads.
# 4. Summarize raw and trimmed QC reports using MultiQC.
#
# Expected input structure:
#   <raw_dir>/<City>/<Run_accession>/*_1.fastq.gz
#   <raw_dir>/<City>/<Run_accession>/*_2.fastq.gz
#
# Example:
#   data/raw_reads/sewage/Berlin/ERRxxxxxxx/ERRxxxxxxx_1.fastq.gz
#   data/raw_reads/sewage/Berlin/ERRxxxxxxx/ERRxxxxxxx_2.fastq.gz
#
# Example run:
#   bash scripts/03_quality_trimming/01_qc_trim_sewage_fastqc_fastp_multiqc.sh \
#     --raw-dir data/raw_reads/sewage \
#     --trimmed-dir data/processed_reads/sewage/trimmed \
#     --qc-dir results/qc/sewage \
#     --logdir logs/03_quality_trimming_sewage \
#     --threads 4

set -euo pipefail

RAW_DIR="data/raw_reads/sewage"
TRIMMED_DIR="data/processed_reads/sewage/trimmed"
QC_DIR="results/qc/sewage"
LOGDIR="logs/03_quality_trimming_sewage"
THREADS=4
FASTP_QUAL=20
MIN_LENGTH=50

usage() {
    echo ""
    echo "Usage:"
    echo "  bash scripts/03_quality_trimming/01_qc_trim_sewage_fastqc_fastp_multiqc.sh \\"
    echo "    --raw-dir <raw_fastq_dir> \\"
    echo "    --trimmed-dir <trimmed_output_dir> \\"
    echo "    --qc-dir <qc_output_dir> \\"
    echo "    --logdir <log_output_dir> \\"
    echo "    --threads <number_of_threads>"
    echo ""
    echo "Optional:"
    echo "  --raw-dir       Input raw sewage FASTQ directory"
    echo "                  Default: data/raw_reads/sewage"
    echo "  --trimmed-dir   Output directory for fastp-trimmed reads"
    echo "                  Default: data/processed_reads/sewage/trimmed"
    echo "  --qc-dir        Output directory for FastQC and MultiQC reports"
    echo "                  Default: results/qc/sewage"
    echo "  --logdir        Log directory"
    echo "                  Default: logs/03_quality_trimming_sewage"
    echo "  --threads       Number of threads"
    echo "                  Default: 4"
    echo "  --fastp-qual    fastp qualified quality threshold"
    echo "                  Default: 20"
    echo "  --min-length    Minimum read length retained by fastp"
    echo "                  Default: 50"
    echo ""
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --raw-dir)
            RAW_DIR="$2"
            shift 2
            ;;
        --trimmed-dir)
            TRIMMED_DIR="$2"
            shift 2
            ;;
        --qc-dir)
            QC_DIR="$2"
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
        --fastp-qual)
            FASTP_QUAL="$2"
            shift 2
            ;;
        --min-length)
            MIN_LENGTH="$2"
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

if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: Raw read directory not found: $RAW_DIR"
    exit 1
fi

if ! command -v fastqc >/dev/null 2>&1; then
    echo "ERROR: FastQC is not installed or not available in PATH."
    exit 1
fi

if ! command -v multiqc >/dev/null 2>&1; then
    echo "ERROR: MultiQC is not installed or not available in PATH."
    exit 1
fi

if ! command -v fastp >/dev/null 2>&1; then
    echo "ERROR: fastp is not installed or not available in PATH."
    exit 1
fi

mkdir -p "$TRIMMED_DIR"
mkdir -p "$QC_DIR/raw_fastqc"
mkdir -p "$QC_DIR/trimmed_fastqc"
mkdir -p "$QC_DIR/multiqc_raw"
mkdir -p "$QC_DIR/multiqc_trimmed"
mkdir -p "$QC_DIR/fastp_reports"
mkdir -p "$LOGDIR"

VERSION_FILE="${LOGDIR}/software_versions_qc_trimming.txt"
MANIFEST="${LOGDIR}/sewage_qc_trimming_manifest.tsv"
FAILED="${LOGDIR}/failed_qc_trimming.tsv"

{
    echo "QC and trimming script run date: $(date)"
    echo ""
    echo "FastQC:"
    fastqc --version
    echo ""
    echo "MultiQC:"
    multiqc --version
    echo ""
    echo "fastp:"
    fastp --version 2>&1 || true
} > "$VERSION_FILE"

echo -e "city\trun_accession\traw_r1\traw_r2\ttrimmed_r1\ttrimmed_r2\tfastp_html\tfastp_json\tstatus" > "$MANIFEST"
echo -e "city\trun_accession\treason" > "$FAILED"

sanitize_name() {
    echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

find_read_pair() {
    local run_dir="$1"
    local run_accession="$2"

    local r1=""
    local r2=""

    r1=$(find "$run_dir" -maxdepth 1 -type f \( \
        -name "${run_accession}_1.fastq.gz" -o \
        -name "${run_accession}_1.fq.gz" -o \
        -name "*_1.fastq.gz" -o \
        -name "*_1.fq.gz" -o \
        -name "*_R1.fastq.gz" -o \
        -name "*_R1.fq.gz" -o \
        -name "*_R1_*.fastq.gz" -o \
        -name "*_R1_*.fq.gz" \
    \) | sort | head -n 1)

    r2=$(find "$run_dir" -maxdepth 1 -type f \( \
        -name "${run_accession}_2.fastq.gz" -o \
        -name "${run_accession}_2.fq.gz" -o \
        -name "*_2.fastq.gz" -o \
        -name "*_2.fq.gz" -o \
        -name "*_R2.fastq.gz" -o \
        -name "*_R2.fq.gz" -o \
        -name "*_R2_*.fastq.gz" -o \
        -name "*_R2_*.fq.gz" \
    \) | sort | head -n 1)

    echo -e "${r1}\t${r2}"
}

echo ""
echo "Starting sewage read QC and trimming"
echo "Raw input directory: $RAW_DIR"
echo "Trimmed output directory: $TRIMMED_DIR"
echo "QC output directory: $QC_DIR"
echo "Log directory: $LOGDIR"
echo "Threads: $THREADS"
echo "fastp quality threshold: Q${FASTP_QUAL}"
echo "fastp minimum read length: ${MIN_LENGTH} bp"
echo ""

while IFS= read -r RUN_DIR; do

    RUN_ACCESSION=$(basename "$RUN_DIR")
    CITY=$(basename "$(dirname "$RUN_DIR")")
    CITY_SAFE=$(sanitize_name "$CITY")

    echo "------------------------------------------------------------"
    echo "City: $CITY"
    echo "Run:  $RUN_ACCESSION"
    echo "Input: $RUN_DIR"

    PAIR_INFO=$(find_read_pair "$RUN_DIR" "$RUN_ACCESSION")
    RAW_R1=$(echo "$PAIR_INFO" | cut -f1)
    RAW_R2=$(echo "$PAIR_INFO" | cut -f2)

    if [[ -z "$RAW_R1" || -z "$RAW_R2" ]]; then
        echo "WARNING: Could not find paired FASTQ files for $RUN_ACCESSION"
        echo -e "${CITY}\t${RUN_ACCESSION}\tpaired_fastq_not_found" >> "$FAILED"
        continue
    fi

    RUN_TRIMMED_DIR="${TRIMMED_DIR}/${CITY_SAFE}/${RUN_ACCESSION}"
    RUN_FASTP_DIR="${QC_DIR}/fastp_reports/${CITY_SAFE}/${RUN_ACCESSION}"

    mkdir -p "$RUN_TRIMMED_DIR"
    mkdir -p "$RUN_FASTP_DIR"

    TRIMMED_R1="${RUN_TRIMMED_DIR}/${RUN_ACCESSION}_trimmed_1.fastq.gz"
    TRIMMED_R2="${RUN_TRIMMED_DIR}/${RUN_ACCESSION}_trimmed_2.fastq.gz"
    FASTP_HTML="${RUN_FASTP_DIR}/${RUN_ACCESSION}_fastp.html"
    FASTP_JSON="${RUN_FASTP_DIR}/${RUN_ACCESSION}_fastp.json"

    echo "Raw R1: $RAW_R1"
    echo "Raw R2: $RAW_R2"

    echo "Running FastQC on raw reads..."
    fastqc \
        --threads "$THREADS" \
        --outdir "$QC_DIR/raw_fastqc" \
        "$RAW_R1" "$RAW_R2"

    echo "Running fastp trimming..."
    fastp \
        --in1 "$RAW_R1" \
        --in2 "$RAW_R2" \
        --out1 "$TRIMMED_R1" \
        --out2 "$TRIMMED_R2" \
        --detect_adapter_for_pe \
        --qualified_quality_phred "$FASTP_QUAL" \
        --length_required "$MIN_LENGTH" \
        --thread "$THREADS" \
        --html "$FASTP_HTML" \
        --json "$FASTP_JSON"

    echo "Running FastQC on trimmed reads..."
    fastqc \
        --threads "$THREADS" \
        --outdir "$QC_DIR/trimmed_fastqc" \
        "$TRIMMED_R1" "$TRIMMED_R2"

    echo -e "${CITY}\t${RUN_ACCESSION}\t${RAW_R1}\t${RAW_R2}\t${TRIMMED_R1}\t${TRIMMED_R2}\t${FASTP_HTML}\t${FASTP_JSON}\tcompleted" >> "$MANIFEST"

done < <(find "$RAW_DIR" -mindepth 2 -maxdepth 2 -type d | sort)

echo ""
echo "Generating MultiQC report for raw FastQC results..."
multiqc \
    "$QC_DIR/raw_fastqc" \
    --outdir "$QC_DIR/multiqc_raw" \
    --filename "multiqc_raw_sewage.html" \
    --force

echo ""
echo "Generating MultiQC report for trimmed FastQC and fastp results..."
multiqc \
    "$QC_DIR/trimmed_fastqc" \
    "$QC_DIR/fastp_reports" \
    --outdir "$QC_DIR/multiqc_trimmed" \
    --filename "multiqc_trimmed_sewage.html" \
    --force

echo ""
echo "QC and trimming step complete."
echo "Trimmed reads: $TRIMMED_DIR"
echo "Raw FastQC reports: $QC_DIR/raw_fastqc"
echo "Trimmed FastQC reports: $QC_DIR/trimmed_fastqc"
echo "fastp reports: $QC_DIR/fastp_reports"
echo "Raw MultiQC report: $QC_DIR/multiqc_raw/multiqc_raw_sewage.html"
echo "Trimmed MultiQC report: $QC_DIR/multiqc_trimmed/multiqc_trimmed_sewage.html"
echo "Manifest: $MANIFEST"
echo "Failed runs: $FAILED"
echo "Software versions: $VERSION_FILE"
echo ""

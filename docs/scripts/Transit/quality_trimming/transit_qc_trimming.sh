#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Script: 01_transit_qc_trimming.sh
#
# Purpose:
#   Perform quality control and trimming for public transit
#   metagenomic paired-end reads.
#
# Workflow:
#   1. Run FastQC on raw reads
#   2. Summarize raw-read QC with MultiQC
#   3. Trim paired-end reads using fastp
#   4. Run FastQC on trimmed reads
#   5. Summarize trimmed-read QC with MultiQC
#
# Input:
#   City-wise raw paired-end reads under:
#   bioinformatics/raw_reads/Transit/<City>/
#
# Output:
#   Trimmed reads:
#   bioinformatics/trimmed_reads/Transit/<City>/
#
#   QC reports:
#   bioinformatics/qc/Transit/pre_trim/
#   bioinformatics/qc/Transit/post_trim/
#
# Notes:
#   This script assumes paired-end reads with filenames ending in:
#   *_1.fastq.gz / *_2.fastq.gz
#   or
#   *_R1.fastq.gz / *_R2.fastq.gz
# ============================================================

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/Thesis_AMR_Project}"

RAW_DIR="${RAW_DIR:-$PROJECT_ROOT/bioinformatics/raw_reads/Transit}"
TRIM_DIR="${TRIM_DIR:-$PROJECT_ROOT/bioinformatics/trimmed_reads/Transit}"

QC_DIR="${QC_DIR:-$PROJECT_ROOT/bioinformatics/qc/Transit}"

PRE_FASTQC_DIR="$QC_DIR/pre_trim/fastqc"
PRE_MULTIQC_DIR="$QC_DIR/pre_trim/multiqc"

POST_FASTQC_DIR="$QC_DIR/post_trim/fastqc"
POST_MULTIQC_DIR="$QC_DIR/post_trim/multiqc"

FASTP_REPORT_DIR="$QC_DIR/post_trim/fastp_reports"
LOG_DIR="$QC_DIR/post_trim/logs"

THREADS="${THREADS:-8}"

# fastp trimming parameters
QUAL="${QUAL:-20}"
MIN_LEN="${MIN_LEN:-50}"

# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p "$PRE_FASTQC_DIR" "$PRE_MULTIQC_DIR"
mkdir -p "$POST_FASTQC_DIR" "$POST_MULTIQC_DIR"
mkdir -p "$FASTP_REPORT_DIR" "$LOG_DIR"
mkdir -p "$TRIM_DIR"

# ============================================================
# CHECK REQUIRED TOOLS
# ============================================================

for tool in fastqc multiqc fastp; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool is not available in PATH."
        echo "Please activate the correct conda/micromamba environment."
        exit 1
    fi
done

# ============================================================
# CHECK INPUT DIRECTORY
# ============================================================

if [ ! -d "$RAW_DIR" ]; then
    echo "ERROR: Raw read directory not found:"
    echo "$RAW_DIR"
    exit 1
fi

echo "============================================================"
echo "Transit QC and trimming pipeline"
echo "============================================================"
echo "Project root       : $PROJECT_ROOT"
echo "Raw reads          : $RAW_DIR"
echo "Trimmed reads      : $TRIM_DIR"
echo "QC directory       : $QC_DIR"
echo "Threads            : $THREADS"
echo "fastp Q threshold  : $QUAL"
echo "fastp min length   : $MIN_LEN"
echo "============================================================"
echo

# ============================================================
# STEP 1: PRE-TRIM FASTQC
# ============================================================

echo "STEP 1: Running pre-trim FastQC..."

mapfile -t RAW_FASTQ_FILES < <(
    find "$RAW_DIR" -type f \
    \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
    | sort
)

if [ "${#RAW_FASTQ_FILES[@]}" -eq 0 ]; then
    echo "ERROR: No raw FASTQ files found under:"
    echo "$RAW_DIR"
    exit 1
fi

echo "Raw FASTQ files found: ${#RAW_FASTQ_FILES[@]}"

fastqc \
    -t "$THREADS" \
    -o "$PRE_FASTQC_DIR" \
    "${RAW_FASTQ_FILES[@]}" \
    2>&1 | tee "$QC_DIR/pre_trim_fastqc.log"

echo "Pre-trim FastQC completed."
echo

# ============================================================
# STEP 2: PRE-TRIM MULTIQC
# ============================================================

echo "STEP 2: Running pre-trim MultiQC..."

multiqc "$PRE_FASTQC_DIR" \
    -o "$PRE_MULTIQC_DIR" \
    -n "transit_pretrim_multiqc_report.html" \
    2>&1 | tee "$QC_DIR/pre_trim_multiqc.log"

echo "Pre-trim MultiQC completed."
echo

# ============================================================
# STEP 3: FASTP TRIMMING
# ============================================================

echo "STEP 3: Running fastp trimming..."

mapfile -t CITY_DIRS < <(
    find "$RAW_DIR" -mindepth 1 -maxdepth 1 -type d | sort
)

if [ "${#CITY_DIRS[@]}" -eq 0 ]; then
    echo "ERROR: No city folders found under:"
    echo "$RAW_DIR"
    exit 1
fi

echo "City folders found: ${#CITY_DIRS[@]}"
echo

for CITY_DIR in "${CITY_DIRS[@]}"; do
    CITY="$(basename "$CITY_DIR")"

    echo "------------------------------------------------------------"
    echo "Processing city: $CITY"
    echo "------------------------------------------------------------"

    mkdir -p "$TRIM_DIR/$CITY"
    mkdir -p "$FASTP_REPORT_DIR/$CITY"

    mapfile -t R1_FILES < <(
        find "$CITY_DIR" -maxdepth 1 -type f \
        \( -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" -o -name "*_R1.fq.gz" -o -name "*_1.fq.gz" \) \
        | sort
    )

    echo "R1 files found for $CITY: ${#R1_FILES[@]}"

    if [ "${#R1_FILES[@]}" -eq 0 ]; then
        echo "WARNING: No R1 files found for $CITY. Skipping city."
        echo
        continue
    fi

    for R1 in "${R1_FILES[@]}"; do
        FNAME="$(basename "$R1")"
        SAMPLE=""
        R2=""

        if [[ "$FNAME" == *_R1.fastq.gz ]]; then
            SAMPLE="${FNAME%_R1.fastq.gz}"
            R2="$CITY_DIR/${SAMPLE}_R2.fastq.gz"
        elif [[ "$FNAME" == *_1.fastq.gz ]]; then
            SAMPLE="${FNAME%_1.fastq.gz}"
            R2="$CITY_DIR/${SAMPLE}_2.fastq.gz"
        elif [[ "$FNAME" == *_R1.fq.gz ]]; then
            SAMPLE="${FNAME%_R1.fq.gz}"
            R2="$CITY_DIR/${SAMPLE}_R2.fq.gz"
        elif [[ "$FNAME" == *_1.fq.gz ]]; then
            SAMPLE="${FNAME%_1.fq.gz}"
            R2="$CITY_DIR/${SAMPLE}_2.fq.gz"
        else
            echo "WARNING: Could not parse R1 filename: $FNAME"
            continue
        fi

        if [ ! -f "$R2" ]; then
            echo "WARNING: Matching R2 file not found for sample: $SAMPLE"
            echo "Expected R2: $R2"
            continue
        fi

        OUT_R1="$TRIM_DIR/$CITY/${SAMPLE}_trimmed_R1.fastq.gz"
        OUT_R2="$TRIM_DIR/$CITY/${SAMPLE}_trimmed_R2.fastq.gz"

        HTML_REPORT="$FASTP_REPORT_DIR/$CITY/${SAMPLE}_fastp.html"
        JSON_REPORT="$FASTP_REPORT_DIR/$CITY/${SAMPLE}_fastp.json"
        LOG_FILE="$LOG_DIR/${CITY}_${SAMPLE}_fastp.log"

        echo "Trimming sample: $SAMPLE"

        fastp \
            -i "$R1" \
            -I "$R2" \
            -o "$OUT_R1" \
            -O "$OUT_R2" \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality "$QUAL" \
            --length_required "$MIN_LEN" \
            --thread "$THREADS" \
            --html "$HTML_REPORT" \
            --json "$JSON_REPORT" \
            > "$LOG_FILE" 2>&1

        echo "Finished trimming: $SAMPLE"
        echo
    done
done

echo "fastp trimming completed."
echo

# ============================================================
# STEP 4: POST-TRIM FASTQC
# ============================================================

echo "STEP 4: Running post-trim FastQC..."

mapfile -t TRIMMED_FASTQ_FILES < <(
    find "$TRIM_DIR" -type f \
    \( -name "*_trimmed_R1.fastq.gz" -o -name "*_trimmed_R2.fastq.gz" \) \
    | sort
)

if [ "${#TRIMMED_FASTQ_FILES[@]}" -eq 0 ]; then
    echo "ERROR: No trimmed FASTQ files found under:"
    echo "$TRIM_DIR"
    exit 1
fi

echo "Trimmed FASTQ files found: ${#TRIMMED_FASTQ_FILES[@]}"

fastqc \
    -t "$THREADS" \
    -o "$POST_FASTQC_DIR" \
    "${TRIMMED_FASTQ_FILES[@]}" \
    2>&1 | tee "$QC_DIR/post_trim_fastqc.log"

echo "Post-trim FastQC completed."
echo

# ============================================================
# STEP 5: POST-TRIM MULTIQC
# ============================================================

echo "STEP 5: Running post-trim MultiQC..."

multiqc "$POST_FASTQC_DIR" \
    -o "$POST_MULTIQC_DIR" \
    -n "transit_posttrim_multiqc_report.html" \
    2>&1 | tee "$QC_DIR/post_trim_multiqc.log"

echo "Post-trim MultiQC completed."
echo

# ============================================================
# FINISH
# ============================================================

echo "============================================================"
echo "Transit QC and trimming pipeline completed successfully."
echo "============================================================"
echo "Trimmed reads:"
echo "$TRIM_DIR"
echo
echo "Pre-trim MultiQC report:"
echo "$PRE_MULTIQC_DIR/transit_pretrim_multiqc_report.html"
echo
echo "Post-trim MultiQC report:"
echo "$POST_MULTIQC_DIR/transit_posttrim_multiqc_report.html"
echo
echo "fastp reports:"
echo "$FASTP_REPORT_DIR"
echo "============================================================"

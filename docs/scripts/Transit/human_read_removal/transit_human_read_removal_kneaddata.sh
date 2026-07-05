#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Script: 02_transit_human_read_removal_kneaddata.sh
#
# Purpose:
#   Remove human-derived reads from trimmed public transit
#   metagenomic paired-end reads using KneadData.
#
# Workflow:
#   1. Take fastp-trimmed paired-end reads as input
#   2. Run KneadData against the human Bowtie2 database
#   3. Bypass KneadData trimming because fastp was already used
#   4. Bypass TRF to avoid unnecessary runtime
#   5. Keep strictly decontaminated paired reads
#   6. Compress final clean paired reads with gzip
#
# Input:
#   bioinformatics/trimmed_reads/Transit/<City>/
#
# Expected input filenames:
#   *_trimmed_R1.fastq.gz / *_trimmed_R2.fastq.gz
#
# Output:
#   bioinformatics/decontaminated_reads/Transit/<City>/<Sample>/
#
# Final clean paired reads:
#   *_paired_1.fastq.gz
#   *_paired_2.fastq.gz
#
# Notes:
#   The KneadData database is not included in GitHub.
#   Users should download/build the database separately and update DB_DIR.
# ============================================================

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/Thesis_AMR_Project}"

TRIM_DIR="${TRIM_DIR:-$PROJECT_ROOT/bioinformatics/trimmed_reads/Transit}"
OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/bioinformatics/decontaminated_reads/Transit}"
LOG_DIR="${LOG_DIR:-$PROJECT_ROOT/bioinformatics/decontamination_logs/Transit}"

# KneadData-compatible human Bowtie2 database directory
DB_DIR="${DB_DIR:-$PROJECT_ROOT/bioinformatics/db/kneaddata}"

THREADS="${THREADS:-4}"

# If true, samples with existing gzipped paired outputs are skipped
SKIP_EXISTING="${SKIP_EXISTING:-true}"

# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p "$OUT_DIR" "$LOG_DIR"

# ============================================================
# CHECK REQUIRED TOOLS
# ============================================================

if ! command -v kneaddata >/dev/null 2>&1; then
    echo "ERROR: kneaddata is not available in PATH."
    echo "Please activate the correct conda/micromamba environment."
    exit 1
fi

if ! command -v gzip >/dev/null 2>&1; then
    echo "ERROR: gzip is not available in PATH."
    exit 1
fi

# ============================================================
# CHECK INPUTS
# ============================================================

if [ ! -d "$TRIM_DIR" ]; then
    echo "ERROR: Trimmed read directory not found:"
    echo "$TRIM_DIR"
    exit 1
fi

if [ ! -d "$DB_DIR" ]; then
    echo "ERROR: KneadData database directory not found:"
    echo "$DB_DIR"
    echo
    echo "Update DB_DIR or run with:"
    echo "DB_DIR=/path/to/kneaddata_db bash scripts/01_preprocessing/02_transit_human_read_removal_kneaddata.sh"
    exit 1
fi

echo "============================================================"
echo "Transit human read removal with KneadData"
echo "============================================================"
echo "Project root        : $PROJECT_ROOT"
echo "Trimmed reads       : $TRIM_DIR"
echo "Decontaminated reads: $OUT_DIR"
echo "Logs                : $LOG_DIR"
echo "KneadData database  : $DB_DIR"
echo "Threads             : $THREADS"
echo "Skip existing       : $SKIP_EXISTING"
echo "============================================================"
echo

# ============================================================
# FIND CITY DIRECTORIES
# ============================================================

mapfile -t CITY_DIRS < <(
    find "$TRIM_DIR" -mindepth 1 -maxdepth 1 -type d | sort
)

if [ "${#CITY_DIRS[@]}" -eq 0 ]; then
    echo "ERROR: No city folders found under:"
    echo "$TRIM_DIR"
    exit 1
fi

echo "City folders found: ${#CITY_DIRS[@]}"
echo

# ============================================================
# PROCESS EACH CITY
# ============================================================

for CITY_DIR in "${CITY_DIRS[@]}"; do
    CITY="$(basename "$CITY_DIR")"

    echo "------------------------------------------------------------"
    echo "Processing city: $CITY"
    echo "------------------------------------------------------------"

    mkdir -p "$OUT_DIR/$CITY"

    mapfile -t R1_FILES < <(
        find "$CITY_DIR" -maxdepth 1 -type f \
        \( -name "*_trimmed_R1.fastq.gz" -o -name "*_trimmed_1.fastq.gz" \) \
        | sort
    )

    echo "Trimmed R1 files found for $CITY: ${#R1_FILES[@]}"

    if [ "${#R1_FILES[@]}" -eq 0 ]; then
        echo "WARNING: No trimmed R1 files found for $CITY. Skipping city."
        echo
        continue
    fi

    for R1 in "${R1_FILES[@]}"; do
        FNAME="$(basename "$R1")"
        SAMPLE=""
        R2=""

        if [[ "$FNAME" == *_trimmed_R1.fastq.gz ]]; then
            SAMPLE="${FNAME%_trimmed_R1.fastq.gz}"
            R2="$CITY_DIR/${SAMPLE}_trimmed_R2.fastq.gz"
        elif [[ "$FNAME" == *_trimmed_1.fastq.gz ]]; then
            SAMPLE="${FNAME%_trimmed_1.fastq.gz}"
            R2="$CITY_DIR/${SAMPLE}_trimmed_2.fastq.gz"
        else
            echo "WARNING: Could not parse R1 filename: $FNAME"
            continue
        fi

        if [ ! -f "$R2" ]; then
            echo "WARNING: Matching R2 file not found for sample: $SAMPLE"
            echo "Expected R2: $R2"
            continue
        fi

        SAMPLE_OUT_DIR="$OUT_DIR/$CITY/$SAMPLE"
        LOG_FILE="$LOG_DIR/${CITY}_${SAMPLE}_kneaddata.log"

        mkdir -p "$SAMPLE_OUT_DIR"

        EXPECTED_CLEAN_R1="$SAMPLE_OUT_DIR/${SAMPLE}_paired_1.fastq.gz"
        EXPECTED_CLEAN_R2="$SAMPLE_OUT_DIR/${SAMPLE}_paired_2.fastq.gz"

        if [ "$SKIP_EXISTING" = "true" ] && [ -f "$EXPECTED_CLEAN_R1" ] && [ -f "$EXPECTED_CLEAN_R2" ]; then
            echo "Skipping existing sample: $SAMPLE"
            echo
            continue
        fi

        echo "Running KneadData human read removal for sample: $SAMPLE"

        kneaddata \
            -i1 "$R1" \
            -i2 "$R2" \
            -db "$DB_DIR" \
            -o "$SAMPLE_OUT_DIR" \
            --output-prefix "$SAMPLE" \
            --bypass-trim \
            --bypass-trf \
            --decontaminate-pairs strict \
            -t "$THREADS" \
            --remove-intermediate-output \
            --log "$LOG_FILE"

        echo "Compressing final clean paired reads for sample: $SAMPLE"

        if [ -f "$SAMPLE_OUT_DIR/${SAMPLE}_paired_1.fastq" ]; then
            gzip -f "$SAMPLE_OUT_DIR/${SAMPLE}_paired_1.fastq"
        fi

        if [ -f "$SAMPLE_OUT_DIR/${SAMPLE}_paired_2.fastq" ]; then
            gzip -f "$SAMPLE_OUT_DIR/${SAMPLE}_paired_2.fastq"
        fi

        if [ ! -f "$EXPECTED_CLEAN_R1" ] || [ ! -f "$EXPECTED_CLEAN_R2" ]; then
            echo "WARNING: Expected clean paired outputs were not found after KneadData:"
            echo "$EXPECTED_CLEAN_R1"
            echo "$EXPECTED_CLEAN_R2"
            echo "Check log file:"
            echo "$LOG_FILE"
        else
            echo "Finished human read removal for sample: $SAMPLE"
        fi

        echo
    done
done

echo "============================================================"
echo "Transit human read removal completed."
echo "============================================================"
echo "Clean non-human paired reads:"
echo "$OUT_DIR"
echo
echo "KneadData logs:"
echo "$LOG_DIR"
echo "============================================================"

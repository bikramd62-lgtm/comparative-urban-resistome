#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Script: 03_transit_citywise_pooling.sh
#
# Purpose:
#   Pool decontaminated public transit metagenomic reads
#   city-wise after human read removal.
#
# Workflow:
#   1. Search each city folder for KneadData clean paired reads
#   2. Concatenate all R1 reads for each city
#   3. Concatenate all R2 reads for each city
#   4. Count pooled read pairs per city
#   5. Report whether each city has at least 5 million pairs
#
# Input:
#   bioinformatics/decontaminated_reads/Transit/<City>/<Sample>/
#
# Expected input files:
#   *_paired_1.fastq.gz
#   *_paired_2.fastq.gz
#
# Output:
#   bioinformatics/pooled_reads/Transit/<City>/
#
# Pooled files:
#   <City>_pooled_R1.fastq.gz
#   <City>_pooled_R2.fastq.gz
#
# Summary:
#   bioinformatics/pooled_reads/Transit/summary/
#   transit_city_pooled_read_summary.tsv
# ============================================================

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/Thesis_AMR_Project}"

INPUT_DIR="${INPUT_DIR:-$PROJECT_ROOT/bioinformatics/decontaminated_reads/Transit}"
POOL_DIR="${POOL_DIR:-$PROJECT_ROOT/bioinformatics/pooled_reads/Transit}"

SUMMARY_DIR="$POOL_DIR/summary"
SUMMARY_FILE="$SUMMARY_DIR/transit_city_pooled_read_summary.tsv"

MIN_PAIRS="${MIN_PAIRS:-5000000}"

# If true, existing pooled files are overwritten
OVERWRITE="${OVERWRITE:-true}"

# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p "$POOL_DIR" "$SUMMARY_DIR"

# ============================================================
# CHECK INPUT DIRECTORY
# ============================================================

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found:"
    echo "$INPUT_DIR"
    exit 1
fi

# ============================================================
# CHECK REQUIRED TOOLS
# ============================================================

if ! command -v zcat >/dev/null 2>&1; then
    echo "ERROR: zcat is not available in PATH."
    exit 1
fi

# ============================================================
# START SUMMARY FILE
# ============================================================

echo -e "City\tNum_Samples\tPooled_R1\tPooled_R2\tR1_Reads\tR2_Reads\tRead_Pairs\tSupports_5M_Pairs" > "$SUMMARY_FILE"

echo "============================================================"
echo "Transit city-wise pooling"
echo "============================================================"
echo "Project root          : $PROJECT_ROOT"
echo "Decontaminated reads  : $INPUT_DIR"
echo "Pooled reads          : $POOL_DIR"
echo "Summary file          : $SUMMARY_FILE"
echo "Minimum required pairs: $MIN_PAIRS"
echo "Overwrite existing    : $OVERWRITE"
echo "============================================================"
echo

# ============================================================
# FIND CITY DIRECTORIES
# ============================================================

mapfile -t CITY_DIRS < <(
    find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort
)

if [ "${#CITY_DIRS[@]}" -eq 0 ]; then
    echo "ERROR: No city folders found under:"
    echo "$INPUT_DIR"
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

    OUT_CITY_DIR="$POOL_DIR/$CITY"
    mkdir -p "$OUT_CITY_DIR"

    POOLED_R1="$OUT_CITY_DIR/${CITY}_pooled_R1.fastq.gz"
    POOLED_R2="$OUT_CITY_DIR/${CITY}_pooled_R2.fastq.gz"

    if [ "$OVERWRITE" != "true" ] && [ -f "$POOLED_R1" ] && [ -f "$POOLED_R2" ]; then
        echo "Skipping $CITY because pooled files already exist."
        echo
        continue
    fi

    mapfile -t R1_FILES < <(
        find "$CITY_DIR" -type f -name "*_paired_1.fastq.gz" | sort
    )

    if [ "${#R1_FILES[@]}" -eq 0 ]; then
        echo "WARNING: No *_paired_1.fastq.gz files found for $CITY. Skipping."
        echo
        continue
    fi

    R2_FILES=()

    for R1 in "${R1_FILES[@]}"; do
        R2="${R1%_paired_1.fastq.gz}_paired_2.fastq.gz"

        if [ ! -f "$R2" ]; then
            echo "WARNING: Missing R2 pair for:"
            echo "$R1"
            echo "Expected:"
            echo "$R2"
            echo "Skipping this sample."
            continue
        fi

        R2_FILES+=("$R2")
    done

    NUM_SAMPLES="${#R2_FILES[@]}"

    if [ "$NUM_SAMPLES" -eq 0 ]; then
        echo "WARNING: No complete paired samples found for $CITY. Skipping."
        echo
        continue
    fi

    echo "Complete paired samples found: $NUM_SAMPLES"

    echo "Pooling R1 reads..."
    cat "${R1_FILES[@]}" > "$POOLED_R1"

    echo "Pooling R2 reads..."
    cat "${R2_FILES[@]}" > "$POOLED_R2"

    echo "Counting pooled reads..."

    R1_READS=$(zcat "$POOLED_R1" | awk 'END {print NR/4}')
    R2_READS=$(zcat "$POOLED_R2" | awk 'END {print NR/4}')

    if [ "$R1_READS" -le "$R2_READS" ]; then
        READ_PAIRS="$R1_READS"
    else
        READ_PAIRS="$R2_READS"
    fi

    if [ "$R1_READS" -eq "$R2_READS" ] && [ "$READ_PAIRS" -ge "$MIN_PAIRS" ]; then
        SUPPORTS_5M="YES"
    else
        SUPPORTS_5M="NO"
    fi

    echo "R1 reads       : $R1_READS"
    echo "R2 reads       : $R2_READS"
    echo "Read pairs     : $READ_PAIRS"
    echo "Supports 5M    : $SUPPORTS_5M"
    echo

    echo -e "${CITY}\t${NUM_SAMPLES}\t${POOLED_R1}\t${POOLED_R2}\t${R1_READS}\t${R2_READS}\t${READ_PAIRS}\t${SUPPORTS_5M}" >> "$SUMMARY_FILE"
done

echo "============================================================"
echo "Transit city-wise pooling completed."
echo "============================================================"
echo "Pooled reads:"
echo "$POOL_DIR"
echo
echo "Summary table:"
echo "$SUMMARY_FILE"
echo "============================================================"

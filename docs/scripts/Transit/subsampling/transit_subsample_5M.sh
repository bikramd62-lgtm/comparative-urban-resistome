#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Script: 04_transit_subsample_5M.sh
#
# Purpose:
#   Subsample city-wise pooled public transit metagenomic reads
#   to 5 million paired reads per city.
#
# Workflow:
#   1. Read pooled city-wise paired FASTQ files
#   2. Confirm that each city has at least 5 million read pairs
#   3. Subsample R1 and R2 using seqtk with the same random seed
#   4. Write standardized 5M paired-end FASTQ files
#   5. Count output reads and generate a summary table
#
# Input:
#   bioinformatics/pooled_reads/Transit/<City>/
#
# Expected input files:
#   <City>_pooled_R1.fastq.gz
#   <City>_pooled_R2.fastq.gz
#
# Output:
#   bioinformatics/subsampled_reads/Transit_5M/<City>/
#
# Output files:
#   <City>_5M_R1.fastq.gz
#   <City>_5M_R2.fastq.gz
#
# Important:
#   5M means 5 million paired reads:
#   5,000,000 records in R1 and 5,000,000 records in R2.
# ============================================================

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/Thesis_AMR_Project}"

INPUT_DIR="${INPUT_DIR:-$PROJECT_ROOT/bioinformatics/pooled_reads/Transit}"
OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/bioinformatics/subsampled_reads/Transit_5M}"

SUMMARY_DIR="$OUT_DIR/summary"
SUMMARY_FILE="$SUMMARY_DIR/transit_5M_subsampling_summary.tsv"

TARGET_PAIRS="${TARGET_PAIRS:-5000000}"
SEED="${SEED:-100}"

# If true, existing subsampled files are overwritten
OVERWRITE="${OVERWRITE:-true}"

# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p "$OUT_DIR" "$SUMMARY_DIR"

# ============================================================
# CHECK REQUIRED TOOLS
# ============================================================

for tool in seqtk gzip zcat; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool is not available in PATH."
        echo "Please activate the correct conda/micromamba environment."
        exit 1
    fi
done

# ============================================================
# CHECK INPUT DIRECTORY
# ============================================================

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found:"
    echo "$INPUT_DIR"
    exit 1
fi

# ============================================================
# START SUMMARY FILE
# ============================================================

echo -e "City\tInput_R1\tInput_R2\tInput_R1_Reads\tInput_R2_Reads\tInput_Read_Pairs\tOutput_R1\tOutput_R2\tOutput_R1_Reads\tOutput_R2_Reads\tTarget_Pairs\tSeed\tStatus" > "$SUMMARY_FILE"

echo "============================================================"
echo "Transit 5M paired-read subsampling"
echo "============================================================"
echo "Project root     : $PROJECT_ROOT"
echo "Input pooled reads: $INPUT_DIR"
echo "Output directory : $OUT_DIR"
echo "Target pairs     : $TARGET_PAIRS"
echo "Random seed      : $SEED"
echo "Overwrite existing: $OVERWRITE"
echo "Summary file     : $SUMMARY_FILE"
echo "============================================================"
echo

# ============================================================
# FIND CITY DIRECTORIES
# ============================================================

mapfile -t CITY_DIRS < <(
    find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d ! -name "summary" | sort
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

    INPUT_R1="$CITY_DIR/${CITY}_pooled_R1.fastq.gz"
    INPUT_R2="$CITY_DIR/${CITY}_pooled_R2.fastq.gz"

    if [ ! -f "$INPUT_R1" ] || [ ! -f "$INPUT_R2" ]; then
        echo "WARNING: Missing pooled paired files for $CITY. Skipping."
        echo "Expected:"
        echo "$INPUT_R1"
        echo "$INPUT_R2"
        echo

        echo -e "${CITY}\t${INPUT_R1}\t${INPUT_R2}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t${TARGET_PAIRS}\t${SEED}\tMISSING_INPUT" >> "$SUMMARY_FILE"
        continue
    fi

    OUT_CITY_DIR="$OUT_DIR/$CITY"
    mkdir -p "$OUT_CITY_DIR"

    OUTPUT_R1="$OUT_CITY_DIR/${CITY}_5M_R1.fastq.gz"
    OUTPUT_R2="$OUT_CITY_DIR/${CITY}_5M_R2.fastq.gz"

    if [ "$OVERWRITE" != "true" ] && [ -f "$OUTPUT_R1" ] && [ -f "$OUTPUT_R2" ]; then
        echo "Skipping $CITY because output files already exist."
        echo

        echo -e "${CITY}\t${INPUT_R1}\t${INPUT_R2}\tNA\tNA\tNA\t${OUTPUT_R1}\t${OUTPUT_R2}\tNA\tNA\t${TARGET_PAIRS}\t${SEED}\tSKIPPED_EXISTING" >> "$SUMMARY_FILE"
        continue
    fi

    echo "Counting input reads..."

    INPUT_R1_READS=$(zcat "$INPUT_R1" | awk 'END {print NR/4}')
    INPUT_R2_READS=$(zcat "$INPUT_R2" | awk 'END {print NR/4}')

    if [ "$INPUT_R1_READS" -le "$INPUT_R2_READS" ]; then
        INPUT_PAIRS="$INPUT_R1_READS"
    else
        INPUT_PAIRS="$INPUT_R2_READS"
    fi

    echo "Input R1 reads : $INPUT_R1_READS"
    echo "Input R2 reads : $INPUT_R2_READS"
    echo "Input pairs    : $INPUT_PAIRS"

    if [ "$INPUT_R1_READS" -ne "$INPUT_R2_READS" ]; then
        echo "WARNING: R1 and R2 read counts are unequal for $CITY. Skipping."
        echo

        echo -e "${CITY}\t${INPUT_R1}\t${INPUT_R2}\t${INPUT_R1_READS}\t${INPUT_R2_READS}\t${INPUT_PAIRS}\t${OUTPUT_R1}\t${OUTPUT_R2}\tNA\tNA\t${TARGET_PAIRS}\t${SEED}\tFAILED_UNEQUAL_INPUT_PAIRS" >> "$SUMMARY_FILE"
        continue
    fi

    if [ "$INPUT_PAIRS" -lt "$TARGET_PAIRS" ]; then
        echo "WARNING: $CITY has fewer than $TARGET_PAIRS paired reads. Skipping."
        echo

        echo -e "${CITY}\t${INPUT_R1}\t${INPUT_R2}\t${INPUT_R1_READS}\t${INPUT_R2_READS}\t${INPUT_PAIRS}\t${OUTPUT_R1}\t${OUTPUT_R2}\tNA\tNA\t${TARGET_PAIRS}\t${SEED}\tFAILED_INSUFFICIENT_READS" >> "$SUMMARY_FILE"
        continue
    fi

    echo "Subsampling R1 to $TARGET_PAIRS reads..."
    seqtk sample -s"$SEED" "$INPUT_R1" "$TARGET_PAIRS" | gzip > "$OUTPUT_R1"

    echo "Subsampling R2 to $TARGET_PAIRS reads..."
    seqtk sample -s"$SEED" "$INPUT_R2" "$TARGET_PAIRS" | gzip > "$OUTPUT_R2"

    echo "Counting output reads..."

    OUTPUT_R1_READS=$(zcat "$OUTPUT_R1" | awk 'END {print NR/4}')
    OUTPUT_R2_READS=$(zcat "$OUTPUT_R2" | awk 'END {print NR/4}')

    echo "Output R1 reads: $OUTPUT_R1_READS"
    echo "Output R2 reads: $OUTPUT_R2_READS"

    if [ "$OUTPUT_R1_READS" -eq "$TARGET_PAIRS" ] && [ "$OUTPUT_R2_READS" -eq "$TARGET_PAIRS" ]; then
        STATUS="SUCCESS"
        echo "Subsampling completed successfully for $CITY."
    else
        STATUS="FAILED_OUTPUT_COUNT_CHECK"
        echo "WARNING: Output read count check failed for $CITY."
    fi

    echo

    echo -e "${CITY}\t${INPUT_R1}\t${INPUT_R2}\t${INPUT_R1_READS}\t${INPUT_R2_READS}\t${INPUT_PAIRS}\t${OUTPUT_R1}\t${OUTPUT_R2}\t${OUTPUT_R1_READS}\t${OUTPUT_R2_READS}\t${TARGET_PAIRS}\t${SEED}\t${STATUS}" >> "$SUMMARY_FILE"
done

echo "============================================================"
echo "Transit 5M paired-read subsampling completed."
echo "============================================================"
echo "Subsampled reads:"
echo "$OUT_DIR"
echo
echo "Summary table:"
echo "$SUMMARY_FILE"
echo "============================================================"

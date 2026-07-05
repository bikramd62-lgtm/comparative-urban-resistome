#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Script: 05_transit_bowtie2_mapping_megares_v3.sh
#
# Purpose:
#   Map subsampled public transit metagenomic reads against
#   MEGARes v3.00 using Bowtie2 very-sensitive settings.
#
# Workflow:
#   1. Read city-wise 5M paired-end transit reads
#   2. Map each city against the MEGARes v3.00 Bowtie2 index
#   3. Save one SAM alignment file per city
#   4. Save Bowtie2 mapping logs
#   5. Generate a mapping summary table
#
# Input:
#   bioinformatics/subsampled_reads/Transit_5M/<City>/
#
# Expected input files:
#   <City>_5M_R1.fastq.gz
#   <City>_5M_R2.fastq.gz
#
# Output:
#   bioinformatics/amr_results/Transit_5M/sam/
#
# Output SAM files:
#   <City>_vs_MEGARes.sam
#
# MEGARes Bowtie2 index:
#   bioinformatics/amr_db/megares/index/megares_v3.00_bt2
#
# Notes:
#   The MEGARes database/index files are not included in GitHub.
#   Users should download/build the database separately and update
#   MEGARES_INDEX if needed.
# ============================================================

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/Thesis_AMR_Project}"

INPUT_DIR="${INPUT_DIR:-$PROJECT_ROOT/bioinformatics/subsampled_reads/Transit_5M}"

OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/bioinformatics/amr_results/Transit_5M}"
SAM_DIR="$OUT_DIR/sam"
LOG_DIR="$OUT_DIR/logs"
SUMMARY_DIR="$OUT_DIR/summary"
SUMMARY_FILE="$SUMMARY_DIR/transit_bowtie2_megares_mapping_summary.tsv"

# Bowtie2 index prefix, not the .bt2 file itself
MEGARES_INDEX="${MEGARES_INDEX:-$PROJECT_ROOT/bioinformatics/amr_db/megares/index/megares_v3.00_bt2}"

THREADS="${THREADS:-4}"

# If true, existing SAM files are overwritten
OVERWRITE="${OVERWRITE:-true}"

# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p "$SAM_DIR" "$LOG_DIR" "$SUMMARY_DIR"

# ============================================================
# CHECK REQUIRED TOOLS
# ============================================================

if ! command -v bowtie2 >/dev/null 2>&1; then
    echo "ERROR: bowtie2 is not available in PATH."
    echo "Please activate the correct conda/micromamba environment."
    exit 1
fi

# ============================================================
# CHECK INPUTS
# ============================================================

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found:"
    echo "$INPUT_DIR"
    exit 1
fi

if [ ! -f "${MEGARES_INDEX}.1.bt2" ] && [ ! -f "${MEGARES_INDEX}.1.bt2l" ]; then
    echo "ERROR: MEGARes Bowtie2 index not found with prefix:"
    echo "$MEGARES_INDEX"
    echo
    echo "Expected files like:"
    echo "${MEGARES_INDEX}.1.bt2"
    echo "${MEGARES_INDEX}.2.bt2"
    echo "${MEGARES_INDEX}.rev.1.bt2"
    echo
    echo "Update MEGARES_INDEX if your index prefix is different."
    exit 1
fi

# ============================================================
# START SUMMARY FILE
# ============================================================

echo -e "City\tInput_R1\tInput_R2\tSAM_File\tBowtie2_Log\tThreads\tStatus" > "$SUMMARY_FILE"

echo "============================================================"
echo "Transit Bowtie2 mapping against MEGARes v3.00"
echo "============================================================"
echo "Project root    : $PROJECT_ROOT"
echo "Input reads     : $INPUT_DIR"
echo "MEGARes index   : $MEGARES_INDEX"
echo "SAM output dir  : $SAM_DIR"
echo "Log dir         : $LOG_DIR"
echo "Threads         : $THREADS"
echo "Overwrite       : $OVERWRITE"
echo "Summary file    : $SUMMARY_FILE"
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

    INPUT_R1="$CITY_DIR/${CITY}_5M_R1.fastq.gz"
    INPUT_R2="$CITY_DIR/${CITY}_5M_R2.fastq.gz"

    SAM_FILE="$SAM_DIR/${CITY}_vs_MEGARes.sam"
    BOWTIE2_LOG="$LOG_DIR/${CITY}_bowtie2_megares.log"

    if [ ! -f "$INPUT_R1" ] || [ ! -f "$INPUT_R2" ]; then
        echo "WARNING: Missing subsampled paired reads for $CITY. Skipping."
        echo "Expected:"
        echo "$INPUT_R1"
        echo "$INPUT_R2"
        echo

        echo -e "${CITY}\t${INPUT_R1}\t${INPUT_R2}\t${SAM_FILE}\t${BOWTIE2_LOG}\t${THREADS}\tMISSING_INPUT" >> "$SUMMARY_FILE"
        continue
    fi

    if [ "$OVERWRITE" != "true" ] && [ -f "$SAM_FILE" ]; then
        echo "Skipping $CITY because SAM file already exists."
        echo

        echo -e "${CITY}\t${INPUT_R1}\t${INPUT_R2}\t${SAM_FILE}\t${BOWTIE2_LOG}\t${THREADS}\tSKIPPED_EXISTING" >> "$SUMMARY_FILE"
        continue
    fi

    echo "Input R1 : $INPUT_R1"
    echo "Input R2 : $INPUT_R2"
    echo "SAM file : $SAM_FILE"
    echo "Log file : $BOWTIE2_LOG"
    echo

    echo "Running Bowtie2 very-sensitive mapping..."

    bowtie2 \
        --very-sensitive \
        -x "$MEGARES_INDEX" \
        -1 "$INPUT_R1" \
        -2 "$INPUT_R2" \
        -S "$SAM_FILE" \
        -p "$THREADS" \
        2> "$BOWTIE2_LOG"

    if [ $? -eq 0 ]; then
        STATUS="SUCCESS"
        echo "Finished mapping: $CITY"
    else
        STATUS="FAILED"
        echo "WARNING: Bowtie2 failed for $CITY. Check log:"
        echo "$BOWTIE2_LOG"
    fi

    echo -e "${CITY}\t${INPUT_R1}\t${INPUT_R2}\t${SAM_FILE}\t${BOWTIE2_LOG}\t${THREADS}\t${STATUS}" >> "$SUMMARY_FILE"
    echo
done

echo "============================================================"
echo "Transit Bowtie2 mapping completed."
echo "============================================================"
echo "SAM files:"
echo "$SAM_DIR"
echo
echo "Bowtie2 logs:"
echo "$LOG_DIR"
echo
echo "Summary table:"
echo "$SUMMARY_FILE"
echo "============================================================"

#!/usr/bin/env bash

# 01_download_sewage_reads_wget.sh
#
# Purpose:
# Download sewage paired-end FASTQ files from ENA using run accession IDs.
#
# Input:
# A TSV file containing at least these columns:
#   city
#   run_accession
#
# Example input:
#   R_analysis/input_data/metadata/derived_runlists_munk2022/sewage_city_run_map_matched_capped5_seed42.tsv
#
# Output structure:
#   <outdir>/<City>/<Run_accession>/*.fastq.gz


set -euo pipefail

RUN_MAP=""
OUTDIR="data/raw_reads/sewage"
LOGDIR="logs/02_download_sewage_reads"
CHECK_MD5=1

usage() {
    echo ""
    echo "Usage:"
    echo "  bash scripts/02_download/01_download_sewage_reads_wget.sh \\"
    echo "    --run-map <city_run_map.tsv> \\"
    echo "    --outdir <download_output_dir> \\"
    echo "    --logdir <log_output_dir>"
    echo ""
    echo "Required:"
    echo "  --run-map   TSV file with columns: city and run_accession"
    echo ""
    echo "Optional:"
    echo "  --outdir    Output directory for downloaded FASTQ files"
    echo "              Default: data/raw_reads/sewage"
    echo "  --logdir    Directory for logs and ENA file reports"
    echo "              Default: logs/02_download_sewage_reads"
    echo "  --no-md5    Skip MD5 checksum verification"
    echo ""
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-map)
            RUN_MAP="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --logdir)
            LOGDIR="$2"
            shift 2
            ;;
        --no-md5)
            CHECK_MD5=0
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

if [[ -z "$RUN_MAP" ]]; then
    echo "ERROR: --run-map is required."
    usage
    exit 1
fi

if [[ ! -f "$RUN_MAP" ]]; then
    echo "ERROR: run map file not found: $RUN_MAP"
    exit 1
fi

if ! command -v wget >/dev/null 2>&1; then
    echo "ERROR: wget is not installed or not available in PATH."
    exit 1
fi

if ! command -v awk >/dev/null 2>&1; then
    echo "ERROR: awk is not installed or not available in PATH."
    exit 1
fi

if [[ "$CHECK_MD5" -eq 1 ]]; then
    if ! command -v md5sum >/dev/null 2>&1; then
        echo "ERROR: md5sum is required for checksum verification."
        echo "Run again with --no-md5 to skip checksum checks."
        exit 1
    fi
fi

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

VERSION_FILE="${LOGDIR}/software_versions_download.txt"
PAIR_FILE="${LOGDIR}/run_city_pairs.tsv"
MANIFEST="${LOGDIR}/download_manifest.tsv"
FAILED="${LOGDIR}/failed_downloads.tsv"

{
    echo "Download script run date: $(date)"
    echo ""
    echo "wget:"
    wget --version | head -n 1
    echo ""
    echo "awk:"
    awk --version 2>/dev/null | head -n 1 || echo "awk version not available"
    echo ""
    if [[ "$CHECK_MD5" -eq 1 ]]; then
        echo "md5sum:"
        md5sum --version 2>/dev/null | head -n 1 || echo "md5sum version not available"
    fi
} > "$VERSION_FILE"

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
NR == 1 {
    city_col = 0
    run_col = 0

    for (i = 1; i <= NF; i++) {
        if ($i == "city") {
            city_col = i
        }
        if ($i == "run_accession") {
            run_col = i
        }
    }

    if (city_col == 0 || run_col == 0) {
        print "ERROR: Input TSV must contain columns named city and run_accession." > "/dev/stderr"
        exit 1
    }

    next
}
NR > 1 {
    if ($city_col != "" && $run_col != "") {
        print $city_col, $run_col
    }
}
' "$RUN_MAP" > "$PAIR_FILE"

if [[ ! -s "$PAIR_FILE" ]]; then
    echo "ERROR: No city-run pairs were extracted from: $RUN_MAP"
    exit 1
fi

echo -e "city\trun_accession\tfile_name\tfile_path\tdownload_url\texpected_md5\tmd5_status" > "$MANIFEST"
echo -e "city\trun_accession\treason" > "$FAILED"

sanitize_city() {
    echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

echo ""
echo "Starting sewage FASTQ download"
echo "Run map: $RUN_MAP"
echo "Output directory: $OUTDIR"
echo "Log directory: $LOGDIR"
echo ""

while IFS=$'\t' read -r CITY RUN; do

    CITY_SAFE=$(sanitize_city "$CITY")
    RUN_DIR="${OUTDIR}/${CITY_SAFE}/${RUN}"
    mkdir -p "$RUN_DIR"

    echo "------------------------------------------------------------"
    echo "City: $CITY"
    echo "Run:  $RUN"
    echo "Output: $RUN_DIR"

    ENA_REPORT="${LOGDIR}/ena_filereport_${RUN}.tsv"

    ENA_QUERY_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"

    wget -qO "$ENA_REPORT" "$ENA_QUERY_URL"

    FASTQ_FTP=$(awk -F '\t' 'NR == 2 {print $2}' "$ENA_REPORT")
    FASTQ_MD5=$(awk -F '\t' 'NR == 2 {print $3}' "$ENA_REPORT")

    if [[ -z "$FASTQ_FTP" ]]; then
        echo "WARNING: No FASTQ FTP links found for $RUN"
        echo -e "${CITY}\t${RUN}\tno_fastq_ftp_links_found" >> "$FAILED"
        continue
    fi

    IFS=';' read -r -a URLS <<< "$FASTQ_FTP"
    IFS=';' read -r -a MD5S <<< "$FASTQ_MD5"

    if [[ "${#URLS[@]}" -lt 2 ]]; then
        echo "WARNING: Fewer than two FASTQ files found for $RUN. Check whether the run is paired-end."
    fi

    for IDX in "${!URLS[@]}"; do

        RAW_URL="${URLS[$IDX]}"

        if [[ "$RAW_URL" =~ ^ftp:// || "$RAW_URL" =~ ^http:// || "$RAW_URL" =~ ^https:// ]]; then
            DOWNLOAD_URL="$RAW_URL"
        else
            DOWNLOAD_URL="https://${RAW_URL}"
        fi

        FILE_NAME=$(basename "$RAW_URL")
        FILE_PATH="${RUN_DIR}/${FILE_NAME}"

        EXPECTED_MD5=""
        if [[ "${#MD5S[@]}" -gt "$IDX" ]]; then
            EXPECTED_MD5="${MD5S[$IDX]}"
        fi

        echo "Downloading: $FILE_NAME"

        wget \
            -c \
            --tries=5 \
            --timeout=60 \
            --waitretry=10 \
            --retry-connrefused \
            -P "$RUN_DIR" \
            "$DOWNLOAD_URL"

        MD5_STATUS="not_checked"

        if [[ "$CHECK_MD5" -eq 1 && -n "$EXPECTED_MD5" ]]; then
            echo "Checking MD5 for: $FILE_NAME"

            if (
                cd "$RUN_DIR"
                echo "${EXPECTED_MD5}  ${FILE_NAME}" | md5sum -c -
            ); then
                MD5_STATUS="passed"
            else
                MD5_STATUS="failed"
                echo "WARNING: MD5 check failed for $FILE_PATH"
                echo -e "${CITY}\t${RUN}\tmd5_failed_${FILE_NAME}" >> "$FAILED"
            fi
        fi

        echo -e "${CITY}\t${RUN}\t${FILE_NAME}\t${FILE_PATH}\t${DOWNLOAD_URL}\t${EXPECTED_MD5}\t${MD5_STATUS}" >> "$MANIFEST"

    done

done < "$PAIR_FILE"

echo ""
echo "Download step complete."
echo "Downloaded files are under: $OUTDIR"
echo "Download manifest: $MANIFEST"
echo "Failed downloads file: $FAILED"
echo "Software versions: $VERSION_FILE"
echo ""

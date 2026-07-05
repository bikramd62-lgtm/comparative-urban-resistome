#!/usr/bin/env python3
"""
01_parse_metadata_identify_matched_cities.py

Purpose
-------
Parse sewage and transit metadata to identify the matched city panel used in the
thesis. The script:

1. Reads sewage metadata files from a metadata directory.
2. Extracts city names and ENA/SRA run accessions.
3. Reads transit/MetaSUB/GeoSeeq metadata.
4. Standardizes city names across both environments.
5. Identifies cities present in both sewage and transit datasets.
6. Adds Bogota sewage run ERR1713346 manually if not already present.
7. Applies a maximum of five sewage runs per city using reproducible random
   sampling with NumPy default_rng(seed=42).
8. Writes derived metadata tables for downstream read download and processing.

This script does not download reads. It only creates the matched-city and run
selection tables.

Author: Bikram Dutta
Project: Comparative metagenomic profiling of urban resistomes in sewage and transit
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import unicodedata
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Project-specific constants
# ---------------------------------------------------------------------

EXPECTED_FINAL_CITIES = [
    "Barcelona",
    "Berlin",
    "Bogota",
    "Hanoi",
    "Hong_Kong",
    "Ilorin",
    "Kuala_Lumpur",
    "Lisbon",
    "Oslo",
    "Porto",
    "Rio_de_Janeiro",
    "Santiago",
    "Singapore",
    "Sofia",
    "Taipei",
    "Vienna",
]

# City harmonization aliases.
# Add more aliases here if metadata contains alternative spellings.
CITY_ALIASES = {
    "barcelona": "Barcelona",
    "berlin": "Berlin",
    "bogota": "Bogota",
    "bogotá": "Bogota",
    "hanoi": "Hanoi",
    "ha noi": "Hanoi",
    "hong kong": "Hong_Kong",
    "hong_kong": "Hong_Kong",
    "ilorin": "Ilorin",
    "kuala lumpur": "Kuala_Lumpur",
    "kuala_lumpur": "Kuala_Lumpur",
    "lisbon": "Lisbon",
    "lisboa": "Lisbon",
    "oslo": "Oslo",
    "porto": "Porto",
    "rio de janeiro": "Rio_de_Janeiro",
    "rio_de_janeiro": "Rio_de_Janeiro",
    "santiago": "Santiago",
    "santiago de chile": "Santiago",
    "singapore": "Singapore",
    "sofia": "Sofia",
    "taipei": "Taipei",
    "vienna": "Vienna",
    "wien": "Vienna",
}

RUN_ACCESSION_RE = re.compile(r"\b(?:ERR|SRR|DRR)\d+\b", flags=re.IGNORECASE)


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def normalize_text(value: object) -> str:
    """Normalize text for matching."""
    if pd.isna(value):
        return ""

    text = str(value).strip()
    text = unicodedata.normalize("NFKD", text)
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = text.lower()
    text = text.replace("_", " ")
    text = text.replace("-", " ")
    text = re.sub(r"[^a-z0-9\s]", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def harmonize_city(value: object) -> Optional[str]:
    """
    Convert city names from different metadata sources into one canonical format.
    Returns None if no expected city is detected.
    """
    text = normalize_text(value)

    if not text:
        return None

    # Exact alias match first
    if text in CITY_ALIASES:
        return CITY_ALIASES[text]

    # Then allow city names embedded in longer metadata strings,
    # e.g. "Rio de Janeiro, Brazil" or "Barcelona_Spain".
    for alias in sorted(CITY_ALIASES, key=len, reverse=True):
        alias_norm = normalize_text(alias)
        pattern = r"\b" + re.escape(alias_norm) + r"\b"
        if re.search(pattern, text):
            return CITY_ALIASES[alias]

    return None


def extract_run_accession(value: object) -> Optional[str]:
    """Extract the first ERR/SRR/DRR accession from a value."""
    if pd.isna(value):
        return None

    match = RUN_ACCESSION_RE.search(str(value))
    if match:
        return match.group(0).upper()

    return None


def extract_run_accession_from_row(row: pd.Series) -> Optional[str]:
    """Extract run accession by scanning all fields in one metadata row."""
    joined = " ".join(row.dropna().astype(str).tolist())
    return extract_run_accession(joined)


def read_table_file(path: Path) -> list[tuple[pd.DataFrame, str]]:
    """
    Read CSV, TSV, TXT, XLS, or XLSX metadata files.

    Returns a list of (dataframe, sheet_or_source_name).
    """
    tables: list[tuple[pd.DataFrame, str]] = []

    suffix = path.suffix.lower()

    try:
        if suffix in {".xlsx", ".xls"}:
            excel = pd.ExcelFile(path)
            for sheet in excel.sheet_names:
                df = pd.read_excel(path, sheet_name=sheet, dtype=str)
                tables.append((df, sheet))

        elif suffix in {".tsv"}:
            df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
            tables.append((df, path.name))

        elif suffix in {".csv"}:
            df = pd.read_csv(path, dtype=str, low_memory=False)
            tables.append((df, path.name))

        elif suffix in {".txt"}:
            # Try tab-delimited first, then allow auto separator detection.
            try:
                df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
            except Exception:
                df = pd.read_csv(path, sep=None, engine="python", dtype=str)
            tables.append((df, path.name))

    except Exception as exc:
        print(f"[WARNING] Could not read {path}: {exc}", file=sys.stderr)

    cleaned_tables: list[tuple[pd.DataFrame, str]] = []

    for df, source_name in tables:
        if df.empty:
            continue

        df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
        df.columns = [str(c).strip() for c in df.columns]

        if not df.empty:
            cleaned_tables.append((df, source_name))

    return cleaned_tables


def iter_metadata_tables(input_path: Path) -> Iterable[tuple[pd.DataFrame, Path, str]]:
    """Yield metadata tables from a file or all readable files in a directory."""
    supported = {".csv", ".tsv", ".txt", ".xlsx", ".xls"}

    if input_path.is_file():
        files = [input_path]
    elif input_path.is_dir():
        files = sorted(
            p for p in input_path.rglob("*")
            if p.is_file() and p.suffix.lower() in supported
        )
    else:
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    for file_path in files:
        for df, source_name in read_table_file(file_path):
            yield df, file_path, source_name


def detect_city_column(df: pd.DataFrame) -> Optional[str]:
    """
    Detect the column most likely to contain city names.

    The function scores columns by:
    - how many values match expected city aliases;
    - whether the column name contains words like city, location, site.
    """
    best_column = None
    best_score = 0

    city_name_keywords = [
        "city",
        "sampling_city",
        "sample_city",
        "location",
        "site",
        "urban",
        "metasub_name",
        "geo_loc_name",
    ]

    for col in df.columns:
        series = df[col].dropna().astype(str)

        if series.empty:
            continue

        values = series.head(10000).tolist()
        matched_count = sum(harmonize_city(v) is not None for v in values)

        col_norm = normalize_text(col)
        name_bonus = 0
        if any(keyword in col_norm for keyword in city_name_keywords):
            name_bonus = 5

        score = matched_count * 10 + name_bonus

        if score > best_score:
            best_score = score
            best_column = col

    if best_score == 0:
        return None

    return best_column


def detect_sample_column(df: pd.DataFrame) -> Optional[str]:
    """Detect a likely sample identifier column in transit metadata."""
    candidates = [
        "sample",
        "sample_id",
        "sampleid",
        "metasub_id",
        "geoseeq",
        "run",
        "accession",
        "library",
    ]

    for col in df.columns:
        col_norm = normalize_text(col)
        if any(candidate in col_norm for candidate in candidates):
            return col

    return None


# ---------------------------------------------------------------------
# Sewage metadata parsing
# ---------------------------------------------------------------------

def parse_sewage_metadata(
    sewage_metadata_dir: Path,
    bogota_run: str = "ERR1713346",
    add_bogota: bool = True,
) -> pd.DataFrame:
    """
    Parse sewage metadata and extract city-to-run mapping.

    This function scans all metadata files for rows containing both:
    - a recognized city name;
    - an ERR/SRR/DRR run accession.
    """
    records: list[dict[str, str]] = []

    for df, file_path, source_name in iter_metadata_tables(sewage_metadata_dir):
        city_col = detect_city_column(df)

        if city_col is None:
            continue

        for _, row in df.iterrows():
            city = harmonize_city(row.get(city_col))
            run_accession = extract_run_accession_from_row(row)

            if city is None or run_accession is None:
                continue

            records.append(
                {
                    "city": city,
                    "run_accession": run_accession,
                    "source_file": str(file_path),
                    "source_sheet_or_table": source_name,
                    "selection_note": "parsed_from_metadata",
                }
            )

    sewage = pd.DataFrame(records)

    if sewage.empty:
        raise RuntimeError(
            "No sewage city-run records were detected. "
            "Check metadata files and city/run column names."
        )

    sewage = sewage.drop_duplicates(subset=["city", "run_accession"]).copy()

    if add_bogota:
        bogota_run = bogota_run.upper()

        already_present = (
            (sewage["city"] == "Bogota")
            & (sewage["run_accession"] == bogota_run)
        ).any()

        if not already_present:
            bogota_record = pd.DataFrame(
                [
                    {
                        "city": "Bogota",
                        "run_accession": bogota_run,
                        "source_file": "manual_entry",
                        "source_sheet_or_table": "Hendriksen_2019",
                        "selection_note": "manual_bogota_hendriksen2019",
                    }
                ]
            )
            sewage = pd.concat([sewage, bogota_record], ignore_index=True)

    sewage = sewage.sort_values(["city", "run_accession"]).reset_index(drop=True)
    return sewage


# ---------------------------------------------------------------------
# Transit metadata parsing
# ---------------------------------------------------------------------

def parse_transit_metadata(transit_metadata_path: Path) -> pd.DataFrame:
    """
    Parse transit metadata and extract city-level sample information.

    The script only needs city names to identify the matched panel, but it also
    retains a likely sample identifier column when available.
    """
    records: list[dict[str, str]] = []

    for df, file_path, source_name in iter_metadata_tables(transit_metadata_path):
        city_col = detect_city_column(df)

        if city_col is None:
            continue

        sample_col = detect_sample_column(df)

        for idx, row in df.iterrows():
            city = harmonize_city(row.get(city_col))

            if city is None:
                continue

            sample_id = ""
            if sample_col is not None and not pd.isna(row.get(sample_col)):
                sample_id = str(row.get(sample_col)).strip()
            else:
                sample_id = f"{file_path.stem}_row_{idx}"

            records.append(
                {
                    "city": city,
                    "sample_id": sample_id,
                    "source_file": str(file_path),
                    "source_sheet_or_table": source_name,
                    "city_column_used": city_col,
                    "sample_column_used": sample_col if sample_col else "",
                }
            )

    transit = pd.DataFrame(records)

    if transit.empty:
        raise RuntimeError(
            "No transit city records were detected. "
            "Check transit metadata path and city column names."
        )

    transit = transit.drop_duplicates(subset=["city", "sample_id"]).copy()
    transit = transit.sort_values(["city", "sample_id"]).reset_index(drop=True)
    return transit


# ---------------------------------------------------------------------
# Matching and sewage run capping
# ---------------------------------------------------------------------

def order_cities(cities: Iterable[str]) -> list[str]:
    """Order cities according to EXPECTED_FINAL_CITIES, then alphabetically."""
    cities = list(cities)
    expected_order = {city: i for i, city in enumerate(EXPECTED_FINAL_CITIES)}

    return sorted(
        cities,
        key=lambda c: (expected_order.get(c, 999), c)
    )


def cap_sewage_runs_per_city(
    sewage_map: pd.DataFrame,
    max_runs_per_city: int = 5,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Retain all runs for cities with <= max_runs_per_city.
    For cities with more runs, randomly select max_runs_per_city runs using
    NumPy default_rng(seed).
    """
    rng = np.random.default_rng(seed)
    selected_frames: list[pd.DataFrame] = []

    for city in order_cities(sewage_map["city"].unique()):
        city_df = sewage_map[sewage_map["city"] == city].copy()
        city_df = city_df.sort_values("run_accession").reset_index(drop=True)

        n_available = len(city_df)

        if n_available <= max_runs_per_city:
            city_df["selected_for_final"] = True
            city_df["selection_reason"] = (
                f"all_available_runs_retained_n_{n_available}"
            )
            selected_frames.append(city_df)
            continue

        selected_idx = rng.choice(
            city_df.index.to_numpy(),
            size=max_runs_per_city,
            replace=False,
        )

        city_df["selected_for_final"] = city_df.index.isin(selected_idx)
        city_df["selection_reason"] = np.where(
            city_df["selected_for_final"],
            f"randomly_selected_max_{max_runs_per_city}_seed_{seed}",
            f"excluded_by_random_cap_max_{max_runs_per_city}_seed_{seed}",
        )

        selected_frames.append(city_df)

    capped = pd.concat(selected_frames, ignore_index=True)
    capped = capped.sort_values(["city", "run_accession"]).reset_index(drop=True)
    return capped


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Parse sewage and transit metadata, identify matched cities, "
            "and generate capped sewage run lists."
        )
    )

    parser.add_argument(
        "--sewage-metadata-dir",
        required=True,
        type=Path,
        help="Directory containing sewage metadata files.",
    )

    parser.add_argument(
        "--transit-metadata",
        required=True,
        type=Path,
        help="Transit metadata file or directory.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for derived metadata tables.",
    )

    parser.add_argument(
        "--max-sewage-runs-per-city",
        type=int,
        default=5,
        help="Maximum number of sewage runs retained per city. Default: 5.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for capped sewage run selection. Default: 42.",
    )

    parser.add_argument(
        "--bogota-run",
        default="ERR1713346",
        help="Bogota sewage run accession from Hendriksen et al. Default: ERR1713346.",
    )

    parser.add_argument(
        "--no-add-bogota",
        action="store_true",
        help="Disable manual addition of Bogota run.",
    )

    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # 1. Parse metadata
    sewage_all = parse_sewage_metadata(
        sewage_metadata_dir=args.sewage_metadata_dir,
        bogota_run=args.bogota_run,
        add_bogota=not args.no_add_bogota,
    )

    transit_all = parse_transit_metadata(args.transit_metadata)

    # 2. Identify matched cities
    sewage_cities = set(sewage_all["city"].unique())
    transit_cities = set(transit_all["city"].unique())

    matched_cities = order_cities(sewage_cities.intersection(transit_cities))

    matched_df = pd.DataFrame({"city": matched_cities})

    expected_set = set(EXPECTED_FINAL_CITIES)
    matched_set = set(matched_cities)

    missing_expected = order_cities(expected_set - matched_set)
    extra_matched = order_cities(matched_set - expected_set)

    # 3. Keep sewage records only for matched cities
    sewage_matched_uncapped = sewage_all[
        sewage_all["city"].isin(matched_cities)
    ].copy()

    # 4. Cap sewage runs to max 5 per city
    sewage_capped_all_status = cap_sewage_runs_per_city(
        sewage_matched_uncapped,
        max_runs_per_city=args.max_sewage_runs_per_city,
        seed=args.seed,
    )

    sewage_capped_selected = sewage_capped_all_status[
        sewage_capped_all_status["selected_for_final"]
    ].copy()

    # 5. Summaries
    sewage_counts = (
        sewage_capped_selected
        .groupby("city", as_index=False)
        .agg(n_selected_sewage_runs=("run_accession", "nunique"))
    )

    transit_counts = (
        transit_all[transit_all["city"].isin(matched_cities)]
        .groupby("city", as_index=False)
        .agg(n_transit_records=("sample_id", "nunique"))
    )

    matched_summary = (
        matched_df
        .merge(sewage_counts, on="city", how="left")
        .merge(transit_counts, on="city", how="left")
        .sort_values("city")
        .reset_index(drop=True)
    )

    # 6. Write outputs
    sewage_all.to_csv(
        args.outdir / "sewage_city_run_map_all_parsed.tsv",
        sep="\t",
        index=False,
    )

    transit_all.to_csv(
        args.outdir / "transit_city_sample_map_all_parsed.tsv",
        sep="\t",
        index=False,
    )

    matched_df.to_csv(
        args.outdir / "matched_cities.tsv",
        sep="\t",
        index=False,
    )

    matched_summary.to_csv(
        args.outdir / "matched_cities_summary.tsv",
        sep="\t",
        index=False,
    )

    sewage_matched_uncapped.to_csv(
        args.outdir / "sewage_city_run_map_matched_uncapped.tsv",
        sep="\t",
        index=False,
    )

    sewage_capped_all_status.to_csv(
        args.outdir / (
            f"sewage_city_run_map_matched_capped"
            f"{args.max_sewage_runs_per_city}_seed{args.seed}_all_status.tsv"
        ),
        sep="\t",
        index=False,
    )

    sewage_capped_selected.to_csv(
        args.outdir / (
            f"sewage_city_run_map_matched_capped"
            f"{args.max_sewage_runs_per_city}_seed{args.seed}.tsv"
        ),
        sep="\t",
        index=False,
    )

    run_list_path = args.outdir / (
        f"sewage_run_accessions_matched_capped"
        f"{args.max_sewage_runs_per_city}_seed{args.seed}.txt"
    )

    with open(run_list_path, "w", encoding="utf-8") as handle:
        for run in sewage_capped_selected["run_accession"].drop_duplicates():
            handle.write(f"{run}\n")

    summary = {
        "n_sewage_cities_detected": len(sewage_cities),
        "n_transit_cities_detected": len(transit_cities),
        "n_matched_cities": len(matched_cities),
        "matched_cities": matched_cities,
        "expected_final_cities": EXPECTED_FINAL_CITIES,
        "missing_expected_cities": missing_expected,
        "extra_matched_cities_not_in_expected_list": extra_matched,
        "max_sewage_runs_per_city": args.max_sewage_runs_per_city,
        "random_seed": args.seed,
        "bogota_run": args.bogota_run,
        "bogota_added_manually": not args.no_add_bogota,
    }

    with open(args.outdir / "metadata_matching_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    # 7. Console report
    print("\nMetadata parsing complete.")
    print(f"Sewage cities detected: {len(sewage_cities)}")
    print(f"Transit cities detected: {len(transit_cities)}")
    print(f"Matched cities detected: {len(matched_cities)}")
    print("\nMatched cities:")
    for city in matched_cities:
        print(f"  - {city}")

    if missing_expected:
        print("\n[WARNING] Expected final cities missing from matched set:")
        for city in missing_expected:
            print(f"  - {city}")

    if extra_matched:
        print("\n[WARNING] Extra matched cities not in expected final list:")
        for city in extra_matched:
            print(f"  - {city}")

    print("\nSelected sewage run counts:")
    print(matched_summary.to_string(index=False))

    print(f"\nOutputs written to: {args.outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env Rscript

# ============================================================
# Script: 07_transit_prevalence_filter_prev10.R
#
# Purpose:
#   Generate qualitative presence/absence matrices from
#   ResistomeAnalyzer ARG profiling outputs and apply a 10%
#   prevalence filter for public transit samples.
#
# Workflow:
#   1. Read per-city ResistomeAnalyzer outputs:
#      - *_gene.tsv
#      - *_group.tsv
#      - *_mechanism.tsv
#      - *_class.tsv
#   2. Combine files into long-format tables
#   3. Generate wide count matrices
#   4. Convert counts to binary presence/absence:
#        count > 0  -> 1
#        count == 0 -> 0
#   5. Calculate city-level prevalence for each feature
#   6. Apply 10% prevalence filter:
#        retain features present in at least 10% of cities
#
# Important:
#   The 80% gene-fraction threshold was already applied upstream
#   during ResistomeAnalyzer processing using:
#        -t 80
#
#   No additional count cutoff is applied here.
#   Binary conversion is based only on count > 0.
#
# Input:
#   bioinformatics/amr_results/Transit_5M/tables/<City>/
#
# Expected input files:
#   <City>_gene.tsv
#   <City>_group.tsv
#   <City>_mechanism.tsv
#   <City>_class.tsv
#
# Output:
#   bioinformatics/downstream_analysis/Transit/amr_geneFrac80_prev10/
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# USER SETTINGS
# ============================================================

project_root <- Sys.getenv(
  "PROJECT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "Thesis_AMR_Project")
)

input_dir <- Sys.getenv(
  "INPUT_DIR",
  unset = file.path(project_root, "bioinformatics/amr_results/Transit_5M/tables")
)

out_root <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(project_root, "bioinformatics/downstream_analysis/Transit/amr_geneFrac80_prev10")
)

environment_name <- "Transit"

levels_to_process <- c("gene", "group", "mechanism", "class")

prevalence_threshold <- 0.10

# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

dirs <- list(
  long_tables      = file.path(out_root, "long_tables"),
  matrices_full    = file.path(out_root, "matrices_full_unfiltered"),
  matrices_prev10  = file.path(out_root, "matrices_prev10_filtered"),
  frequency        = file.path(out_root, "frequency"),
  richness         = file.path(out_root, "richness"),
  summary          = file.path(out_root, "summary")
)

walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

read_ra_table <- function(file_path) {
  readr::read_tsv(file_path, show_col_types = FALSE, progress = FALSE)
}

detect_feature_column <- function(df, level_name) {
  candidates <- c(
    str_to_title(level_name),
    level_name,
    "Gene", "Group", "Mechanism", "Class",
    "gene", "group", "mechanism", "class",
    "feature", "Feature"
  )

  hit <- intersect(candidates, names(df))

  if (length(hit) > 0) {
    return(hit[1])
  }

  character_columns <- names(df)[vapply(df, is.character, logical(1))]

  if (length(character_columns) > 0) {
    return(character_columns[1])
  }

  stop("Could not detect feature column.")
}

detect_count_column <- function(df) {
  candidates <- c(
    "Hits", "hits",
    "Count", "count",
    "Counts", "counts",
    "read_count", "Read_Count"
  )

  hit <- intersect(candidates, names(df))

  if (length(hit) > 0) {
    return(hit[1])
  }

  numeric_columns <- names(df)[vapply(df, is.numeric, logical(1))]

  if (length(numeric_columns) > 0) {
    return(numeric_columns[1])
  }

  stop("Could not detect count column.")
}

standardize_ra_table <- function(file_path, city, level_name) {
  df <- read_ra_table(file_path)

  feature_col <- detect_feature_column(df, level_name)
  count_col <- detect_count_column(df)

  df %>%
    transmute(
      city = city,
      environment = environment_name,
      level = level_name,
      feature = as.character(.data[[feature_col]]),
      count = as.numeric(.data[[count_col]])
    ) %>%
    filter(!is.na(feature), feature != "") %>%
    mutate(count = replace_na(count, 0)) %>%
    group_by(city, environment, level, feature) %>%
    summarise(count = sum(count), .groups = "drop")
}

make_count_matrix <- function(long_df) {
  long_df %>%
    select(feature, city, count) %>%
    pivot_wider(
      names_from = city,
      values_from = count,
      values_fill = 0
    ) %>%
    arrange(feature)
}

make_binary_matrix <- function(count_matrix) {
  count_matrix %>%
    mutate(across(-feature, ~ ifelse(.x > 0, 1L, 0L)))
}

make_frequency_table <- function(binary_matrix) {
  city_columns <- setdiff(names(binary_matrix), "feature")
  n_cities <- length(city_columns)

  binary_matrix %>%
    mutate(
      n_cities_present = rowSums(across(all_of(city_columns))),
      proportion_cities_present = n_cities_present / n_cities
    ) %>%
    arrange(desc(n_cities_present), feature)
}

apply_prev10_filter <- function(count_matrix, binary_matrix, frequency_table, min_cities_required) {
  retained_features <- frequency_table %>%
    filter(n_cities_present >= min_cities_required) %>%
    pull(feature)

  list(
    count_matrix_prev10 = count_matrix %>%
      filter(feature %in% retained_features),

    binary_matrix_prev10 = binary_matrix %>%
      filter(feature %in% retained_features),

    frequency_prev10 = frequency_table %>%
      filter(feature %in% retained_features)
  )
}

make_richness_table <- function(binary_matrix, level_name) {
  city_columns <- setdiff(names(binary_matrix), "feature")

  tibble(
    city = city_columns,
    environment = environment_name,
    level = level_name,
    richness = map_int(city_columns, ~ sum(binary_matrix[[.x]] > 0))
  ) %>%
    arrange(city)
}

# ============================================================
# CHECK INPUT
# ============================================================

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

city_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)

if (length(city_dirs) == 0) {
  stop("No city folders found in: ", input_dir)
}

city_names <- basename(city_dirs)
n_cities <- length(city_names)

min_cities_required <- ceiling(prevalence_threshold * n_cities)

message("============================================================")
message("Transit 10% prevalence filtering")
message("============================================================")
message("Input directory             : ", input_dir)
message("Output directory            : ", out_root)
message("Environment                 : ", environment_name)
message("Number of cities detected   : ", n_cities)
message("Prevalence threshold        : ", prevalence_threshold)
message("Minimum cities required     : ", min_cities_required)
message("Binary rule                 : count > 0")
message("Extra count cutoff          : none")
message("Gene-fraction threshold     : already applied upstream with ResistomeAnalyzer -t 80")
message("============================================================")

# ============================================================
# MAIN PROCESSING
# ============================================================

all_richness <- list()
summary_rows <- list()

for (level_name in levels_to_process) {

  message("------------------------------------------------------------")
  message("Processing level: ", level_name)
  message("------------------------------------------------------------")

  level_tables <- list()

  for (city_dir in city_dirs) {

    city <- basename(city_dir)
    file_path <- file.path(city_dir, paste0(city, "_", level_name, ".tsv"))

    if (!file.exists(file_path)) {
      warning("Missing file: ", file_path)
      next
    }

    level_tables[[city]] <- standardize_ra_table(
      file_path = file_path,
      city = city,
      level_name = level_name
    )
  }

  long_df <- bind_rows(level_tables)

  if (nrow(long_df) == 0) {
    warning("No data found for level: ", level_name)
    next
  }

  # ----------------------------------------------------------
  # Long table
  # ----------------------------------------------------------

  long_out <- file.path(
    dirs$long_tables,
    paste0("transit_", level_name, "_long.tsv")
  )

  write_tsv(long_df, long_out)

  # ----------------------------------------------------------
  # Full count and binary matrices
  # ----------------------------------------------------------

  count_matrix <- make_count_matrix(long_df)
  binary_matrix <- make_binary_matrix(count_matrix)
  frequency_table <- make_frequency_table(binary_matrix)

  count_full_out <- file.path(
    dirs$matrices_full,
    paste0("transit_", level_name, "_count_matrix_full.tsv")
  )

  binary_full_out <- file.path(
    dirs$matrices_full,
    paste0("transit_", level_name, "_binary_matrix_full.tsv")
  )

  freq_full_out <- file.path(
    dirs$frequency,
    paste0("transit_", level_name, "_frequency_full.tsv")
  )

  write_tsv(count_matrix, count_full_out)
  write_tsv(binary_matrix, binary_full_out)
  write_tsv(frequency_table, freq_full_out)

  # ----------------------------------------------------------
  # 10% prevalence filtering
  # ----------------------------------------------------------

  filtered <- apply_prev10_filter(
    count_matrix = count_matrix,
    binary_matrix = binary_matrix,
    frequency_table = frequency_table,
    min_cities_required = min_cities_required
  )

  count_matrix_prev10 <- filtered$count_matrix_prev10
  binary_matrix_prev10 <- filtered$binary_matrix_prev10
  frequency_prev10 <- filtered$frequency_prev10

  count_prev10_out <- file.path(
    dirs$matrices_prev10,
    paste0("transit_", level_name, "_count_matrix_prev10.tsv")
  )

  binary_prev10_out <- file.path(
    dirs$matrices_prev10,
    paste0("transit_", level_name, "_binary_matrix_prev10.tsv")
  )

  freq_prev10_out <- file.path(
    dirs$frequency,
    paste0("transit_", level_name, "_frequency_prev10.tsv")
  )

  write_tsv(count_matrix_prev10, count_prev10_out)
  write_tsv(binary_matrix_prev10, binary_prev10_out)
  write_tsv(frequency_prev10, freq_prev10_out)

  # ----------------------------------------------------------
  # Richness after prevalence filtering
  # ----------------------------------------------------------

  richness_prev10 <- make_richness_table(
    binary_matrix = binary_matrix_prev10,
    level_name = level_name
  )

  richness_out <- file.path(
    dirs$richness,
    paste0("transit_", level_name, "_richness_prev10.tsv")
  )

  write_tsv(richness_prev10, richness_out)

  all_richness[[level_name]] <- richness_prev10

  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------

  summary_rows[[level_name]] <- tibble(
    environment = environment_name,
    level = level_name,
    n_cities = n_cities,
    prevalence_threshold = prevalence_threshold,
    min_cities_required = min_cities_required,
    n_features_full = nrow(binary_matrix),
    n_features_retained_prev10 = nrow(binary_matrix_prev10),
    n_features_removed_prev10 = nrow(binary_matrix) - nrow(binary_matrix_prev10),
    binary_rule = "count > 0",
    gene_fraction_threshold = "80% applied upstream with ResistomeAnalyzer -t 80"
  )

  message("Features before filtering : ", nrow(binary_matrix))
  message("Features after prev10     : ", nrow(binary_matrix_prev10))
}

# ============================================================
# COMBINED RICHNESS SUMMARY
# ============================================================

if (length(all_richness) > 0) {

  combined_richness <- bind_rows(all_richness)

  combined_richness_out <- file.path(
    dirs$richness,
    "transit_all_levels_richness_prev10.tsv"
  )

  write_tsv(combined_richness, combined_richness_out)

  city_richness_summary <- combined_richness %>%
    select(city, level, richness) %>%
    pivot_wider(
      names_from = level,
      values_from = richness,
      values_fill = 0
    ) %>%
    arrange(city)

  city_richness_summary_out <- file.path(
    dirs$richness,
    "transit_city_richness_summary_prev10.tsv"
  )

  write_tsv(city_richness_summary, city_richness_summary_out)
}

# ============================================================
# WRITE FILTERING SUMMARY
# ============================================================

filtering_summary <- bind_rows(summary_rows)

summary_out <- file.path(
  dirs$summary,
  "transit_prev10_filtering_summary.tsv"
)

write_tsv(filtering_summary, summary_out)

message("============================================================")
message("Transit 10% prevalence filtering completed.")
message("============================================================")
message("Output directory:")
message(out_root)
message("")
message("Filtering summary:")
message(summary_out)
message("============================================================")

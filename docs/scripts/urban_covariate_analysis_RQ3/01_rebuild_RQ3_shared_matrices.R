# ============================================================
# 01_rebuild_RQ3_shared_matrices.R
#
# Purpose:
# Rebuild the corrected RQ3 city-level ARG-group overlap tables.
#
# Inputs:
#   results/RQ3_output_tables/sewage_group_binary_matrix_filtered10pct.tsv
#   results/RQ3_output_tables/transit_group_binary_matrix_prev10.tsv
#   results/RQ3_output_tables/RQ3_master_table_with_coordinates_and_NASA_climate.csv
#
# Outputs:
#   results/RQ3_output_tables/RQ3_city_environment_combined_matrix.csv
#   results/RQ3_output_tables/RQ3_shared_ARG_group_matrix.csv
#   results/RQ3_output_tables/RQ3_overlap_metrics_corrected.csv
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_dir  <- "results/RQ3_output_tables"
output_dir <- "results/RQ3_output_tables"

sewage_file <- file.path(input_dir, "sewage_group_binary_matrix_filtered10pct.tsv")
transit_file <- file.path(input_dir, "transit_group_binary_matrix_prev10.tsv")
master_file <- file.path(input_dir, "RQ3_master_table_with_coordinates_and_NASA_climate.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(sewage_file, transit_file, master_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required input files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

standardize_city <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[[:space:]-]+", "_", x)
  x
}

detect_city_column <- function(df) {
  candidates <- c(
    "city", "City", "CITY",
    "city_name", "City_name", "city_clean", "City_clean",
    "sample", "Sample", "sample_id", "Sample_ID", "SampleID"
  )

  hit <- candidates[candidates %in% names(df)]

  if (length(hit) > 0) {
    return(hit[1])
  }

  names(df)[1]
}

to_binary <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.numeric(x) || is.integer(x)) {
    return(as.integer(x > 0))
  }

  x_chr <- trimws(as.character(x))
  x_lower <- tolower(x_chr)

  out <- rep(NA_real_, length(x_chr))

  out[x_lower %in% c("1", "true", "present", "yes")] <- 1
  out[x_lower %in% c("0", "false", "absent", "no")] <- 0

  suppressWarnings({
    numeric_values <- as.numeric(x_chr)
  })

  out[is.na(out) & !is.na(numeric_values)] <- numeric_values[is.na(out) & !is.na(numeric_values)]

  if (any(is.na(out))) {
    warning("Some non-numeric/non-binary values were converted to 0.")
    out[is.na(out)] <- 0
  }

  as.integer(out > 0)
}

make_unique_feature_names <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- paste0("unnamed_feature_", seq_len(sum(is.na(x) | x == "")))
  make.unique(x)
}

read_binary_matrix <- function(file, expected_cities, label) {
  message("Reading ", label, " matrix: ", file)

  df <- read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )

  if (ncol(df) < 2) {
    stop(label, " matrix must contain at least one city column and one feature column.")
  }

  first_col <- names(df)[1]

  row_city_hits <- sum(standardize_city(df[[first_col]]) %in% expected_cities)
  col_city_hits <- sum(standardize_city(names(df)[-1]) %in% expected_cities)

  # Case 1: rows are cities and columns are ARG groups
  if (row_city_hits >= col_city_hits) {
    city_col <- detect_city_column(df)

    city_values <- standardize_city(df[[city_col]])
    feature_cols <- setdiff(names(df), city_col)

    feature_df <- df[, feature_cols, drop = FALSE]
    feature_df[] <- lapply(feature_df, to_binary)

    mat <- as.matrix(feature_df)
    storage.mode(mat) <- "integer"

    rownames(mat) <- city_values
    colnames(mat) <- make_unique_feature_names(feature_cols)

  } else {
    # Case 2: rows are ARG groups and columns are cities; transpose
    feature_names <- make_unique_feature_names(df[[first_col]])
    value_df <- df[, -1, drop = FALSE]

    city_values <- standardize_city(names(value_df))

    value_df[] <- lapply(value_df, to_binary)

    mat <- t(as.matrix(value_df))
    storage.mode(mat) <- "integer"

    rownames(mat) <- city_values
    colnames(mat) <- feature_names
  }

  if (any(duplicated(rownames(mat)))) {
    duplicated_cities <- unique(rownames(mat)[duplicated(rownames(mat))])
    stop(
      "Duplicated city names detected in ", label, " matrix:\n",
      paste(duplicated_cities, collapse = ", ")
    )
  }

  message(label, " matrix dimensions: ", nrow(mat), " cities x ", ncol(mat), " ARG groups")

  mat
}

align_matrix_features <- function(mat, all_features) {
  aligned <- matrix(
    0L,
    nrow = nrow(mat),
    ncol = length(all_features),
    dimnames = list(rownames(mat), all_features)
  )

  shared_features <- intersect(colnames(mat), all_features)

  aligned[, shared_features] <- mat[, shared_features, drop = FALSE]

  aligned
}

write_matrix_csv <- function(mat, file) {
  out <- data.frame(
    city = rownames(mat),
    mat,
    check.names = FALSE
  )

  write_csv(out, file)
}

# ------------------------------------------------------------
# 3. Read master covariate table
# ------------------------------------------------------------

master <- read_csv(
  master_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

master_city_col <- detect_city_column(master)

master[[master_city_col]] <- standardize_city(master[[master_city_col]])

if (master_city_col != "city") {
  master <- master %>%
    rename(city = all_of(master_city_col))
}

if (any(duplicated(master$city))) {
  duplicated_cities <- unique(master$city[duplicated(master$city)])
  stop(
    "Duplicated city names detected in master table:\n",
    paste(duplicated_cities, collapse = ", ")
  )
}

expected_cities <- master$city

message("Master table cities: ", length(expected_cities))

# ------------------------------------------------------------
# 4. Read sewage and transit prevalence-filtered matrices
# ------------------------------------------------------------

sewage_mat <- read_binary_matrix(
  file = sewage_file,
  expected_cities = expected_cities,
  label = "sewage"
)

transit_mat <- read_binary_matrix(
  file = transit_file,
  expected_cities = expected_cities,
  label = "transit"
)

# ------------------------------------------------------------
# 5. Match cities and align ARG-group columns
# ------------------------------------------------------------

common_cities <- expected_cities[
  expected_cities %in% rownames(sewage_mat) &
    expected_cities %in% rownames(transit_mat)
]

if (length(common_cities) == 0) {
  stop("No matching cities were found between sewage, transit and master table.")
}

if (length(common_cities) < length(expected_cities)) {
  warning(
    "Some cities from the master table were not present in both matrices:\n",
    paste(setdiff(expected_cities, common_cities), collapse = ", ")
  )
}

message("Matched cities used for RQ3: ", length(common_cities))

sewage_mat  <- sewage_mat[common_cities, , drop = FALSE]
transit_mat <- transit_mat[common_cities, , drop = FALSE]

all_features <- sort(union(colnames(sewage_mat), colnames(transit_mat)))

sewage_aligned  <- align_matrix_features(sewage_mat, all_features)
transit_aligned <- align_matrix_features(transit_mat, all_features)

# ------------------------------------------------------------
# 6. Build combined city-environment matrix
# ------------------------------------------------------------

combined_matrix <- bind_rows(
  data.frame(
    city = rownames(sewage_aligned),
    environment = "sewage",
    sewage_aligned,
    check.names = FALSE
  ),
  data.frame(
    city = rownames(transit_aligned),
    environment = "transit",
    transit_aligned,
    check.names = FALSE
  )
)

write_csv(
  combined_matrix,
  file.path(output_dir, "RQ3_city_environment_combined_matrix.csv")
)

# ------------------------------------------------------------
# 7. Build city-level shared ARG-group matrix
# ------------------------------------------------------------

shared_mat <- (sewage_aligned == 1L) & (transit_aligned == 1L)
shared_mat <- matrix(
  as.integer(shared_mat),
  nrow = nrow(shared_mat),
  ncol = ncol(shared_mat),
  dimnames = dimnames(shared_mat)
)

# Remove ARG groups not shared in any city
shared_mat <- shared_mat[, colSums(shared_mat) > 0, drop = FALSE]

write_matrix_csv(
  shared_mat,
  file.path(output_dir, "RQ3_shared_ARG_group_matrix.csv")
)

# ------------------------------------------------------------
# 8. Calculate corrected overlap metrics
# ------------------------------------------------------------

union_mat <- (sewage_aligned == 1L) | (transit_aligned == 1L)

sewage_ARG_group_count <- rowSums(sewage_aligned)
transit_ARG_group_count <- rowSums(transit_aligned)
shared_ARG_group_count <- rowSums(shared_mat)
union_ARG_group_count <- rowSums(union_mat)

jaccard_similarity <- ifelse(
  union_ARG_group_count > 0,
  shared_ARG_group_count / union_ARG_group_count,
  NA_real_
)

overlap_metrics <- tibble(
  city = common_cities,
  sewage_ARG_group_count = as.integer(sewage_ARG_group_count),
  transit_ARG_group_count = as.integer(transit_ARG_group_count),
  shared_ARG_group_count = as.integer(shared_ARG_group_count),
  union_ARG_group_count = as.integer(union_ARG_group_count),
  jaccard_similarity_sewage_transit_ARG_group = jaccard_similarity,
  jaccard_distance_sewage_transit_ARG_group = 1 - jaccard_similarity
)

write_csv(
  overlap_metrics,
  file.path(output_dir, "RQ3_overlap_metrics_corrected.csv")
)

# ------------------------------------------------------------
# 9. Create final corrected RQ3 master table
# ------------------------------------------------------------

metric_columns_to_replace <- c(
  "sewage_ARG_group_count",
  "transit_ARG_group_count",
  "shared_ARG_group_count",
  "union_ARG_group_count",
  "jaccard_similarity_sewage_transit_ARG_group",
  "jaccard_distance_sewage_transit_ARG_group"
)

master_clean <- master %>%
  select(-any_of(metric_columns_to_replace))

master_corrected <- master_clean %>%
  inner_join(overlap_metrics, by = "city")

write_csv(
  master_corrected,
  file.path(output_dir, "RQ3_master_table_FINAL_corrected.csv")
)

# ------------------------------------------------------------
# 10. Summary printed to console
# ------------------------------------------------------------

message("\nRQ3 matrix reconstruction completed successfully.")
message("Number of matched cities: ", nrow(master_corrected))
message("Total ARG groups in aligned sewage/transit union: ", length(all_features))
message("ARG groups shared in at least one city: ", ncol(shared_mat))

message("\nCorrected overlap metric summary:")
print(
  overlap_metrics %>%
    summarise(
      min_jaccard = min(jaccard_similarity_sewage_transit_ARG_group, na.rm = TRUE),
      median_jaccard = median(jaccard_similarity_sewage_transit_ARG_group, na.rm = TRUE),
      max_jaccard = max(jaccard_similarity_sewage_transit_ARG_group, na.rm = TRUE),
      min_shared_count = min(shared_ARG_group_count, na.rm = TRUE),
      median_shared_count = median(shared_ARG_group_count, na.rm = TRUE),
      max_shared_count = max(shared_ARG_group_count, na.rm = TRUE)
    )
)

message("\nFiles written to: ", output_dir)
message(" - RQ3_city_environment_combined_matrix.csv")
message(" - RQ3_shared_ARG_group_matrix.csv")
message(" - RQ3_overlap_metrics_corrected.csv")
message(" - RQ3_master_table_FINAL_corrected.csv")

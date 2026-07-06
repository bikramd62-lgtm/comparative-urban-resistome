# ============================================================
# 01_build_rq1_binary_matrices.R
#
# Purpose:
# Build corrected RQ1 combined binary matrices for comparative
# sewage-transit resistome analysis.
#
# Key correction:
# Prevalence filtering is used only to define eligible features.
# The final combined matrix is rebuilt from the unfiltered
# sewage and transit matrices so that genuine detections in the
# opposite environment are preserved.
#
# Expected input orientation:
# Rows    = ARG features
# Columns = cities
#
# Example input files:
# data/interim/rq1/unfiltered/sewage_group_binary_matrix.tsv
# data/interim/rq1/unfiltered/transit_group_binary_matrix.tsv
# data/interim/rq1/unfiltered/sewage_class_binary_matrix.tsv
# data/interim/rq1/unfiltered/transit_class_binary_matrix.tsv
#
# Main outputs:
# results/rq1/matrices/rq1_group_combined_binary_matrix.tsv
# results/rq1/matrices/rq1_class_combined_binary_matrix.tsv
# ============================================================

# ------------------------------------------------------------
# 1. User-adjustable settings
# ------------------------------------------------------------

prevalence_threshold <- 0.10
expected_number_of_cities <- 16

min_city_prevalence <- ceiling(
  prevalence_threshold * expected_number_of_cities
)

input_dir <- file.path(
  "data",
  "interim",
  "rq1",
  "unfiltered"
)

output_dir <- file.path(
  "results",
  "rq1",
  "matrices"
)

# ARG Group input files
sewage_group_file <- file.path(
  input_dir,
  "sewage_group_binary_matrix.tsv"
)

transit_group_file <- file.path(
  input_dir,
  "transit_group_binary_matrix.tsv"
)

# Resistance Class input files
sewage_class_file <- file.path(
  input_dir,
  "sewage_class_binary_matrix.tsv"
)

transit_class_file <- file.path(
  input_dir,
  "transit_class_binary_matrix.tsv"
)

message(
  "\nPrevalence threshold: ",
  prevalence_threshold * 100,
  "%"
)

message(
  "Expected number of matched cities per environment: ",
  expected_number_of_cities
)

message(
  "Minimum city prevalence used for filtering: ",
  min_city_prevalence,
  " cities"
)

# ------------------------------------------------------------
# 2. Create output directory
# ------------------------------------------------------------

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# 3. Check required input files
# ------------------------------------------------------------

required_input_files <- c(
  sewage_group_file,
  transit_group_file,
  sewage_class_file,
  transit_class_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "The following required input files are missing:\n",
      paste(missing_input_files, collapse = "\n")
    )
  )
}

# ------------------------------------------------------------
# 4. Helper function: read unfiltered binary matrix
# ------------------------------------------------------------

read_binary_matrix <- function(file_path, matrix_name) {

  if (!file.exists(file_path)) {
    stop(
      paste0(
        matrix_name,
        " file not found:\n",
        file_path
      )
    )
  }

  message(
    "\nReading ",
    matrix_name,
    ":\n",
    file_path
  )

  raw_df <- read.delim(
    file = file_path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(raw_df) < 2) {
    stop(
      paste0(
        matrix_name,
        " must contain one feature column and at least one city column."
      )
    )
  }

  feature_ids <- trimws(
    as.character(raw_df[[1]])
  )

  if (any(is.na(feature_ids)) || any(feature_ids == "")) {
    stop(
      paste0(
        matrix_name,
        " contains missing or empty feature names in the first column."
      )
    )
  }

  value_df <- raw_df[
    ,
    -1,
    drop = FALSE
  ]

  value_df[] <- lapply(
    value_df,
    function(x) {
      suppressWarnings(
        as.numeric(x)
      )
    }
  )

  value_matrix <- as.matrix(
    value_df
  )

  rownames(value_matrix) <- feature_ids

  if (any(is.na(value_matrix))) {
    stop(
      paste0(
        matrix_name,
        " contains missing or non-numeric values."
      )
    )
  }

  # Convert any positive value to presence.
  # This allows the function to accept either binary matrices
  # or count-like matrices.
  value_matrix <- ifelse(
    value_matrix > 0,
    1,
    0
  )

  storage.mode(value_matrix) <- "numeric"

  # If duplicated feature rows exist, collapse them by maximum presence.
  if (any(duplicated(rownames(value_matrix)))) {

    message(
      matrix_name,
      ": duplicated feature names detected; collapsing duplicates using maximum presence."
    )

    collapsed_df <- data.frame(
      Feature = rownames(value_matrix),
      value_matrix,
      check.names = FALSE
    )

    collapsed_df <- aggregate(
      . ~ Feature,
      data = collapsed_df,
      FUN = max
    )

    feature_ids_collapsed <- collapsed_df$Feature

    value_matrix <- as.matrix(
      collapsed_df[
        ,
        -1,
        drop = FALSE
      ]
    )

    rownames(value_matrix) <- feature_ids_collapsed

    storage.mode(value_matrix) <- "numeric"
  }

  colnames(value_matrix) <- trimws(
    colnames(value_matrix)
  )

  if (any(colnames(value_matrix) == "")) {
    stop(
      paste0(
        matrix_name,
        " contains empty city names."
      )
    )
  }

  message(
    matrix_name,
    ": ",
    nrow(value_matrix),
    " features × ",
    ncol(value_matrix),
    " cities"
  )

  return(value_matrix)
}

# ------------------------------------------------------------
# 5. Helper function: standardise city names
# ------------------------------------------------------------

standardise_city_names <- function(city_names) {

  city_names <- trimws(
    city_names
  )

  # Remove environment suffixes if they accidentally exist.
  city_names <- sub(
    "_Sewage$",
    "",
    city_names
  )

  city_names <- sub(
    "_Transit$",
    "",
    city_names
  )

  return(city_names)
}

# ------------------------------------------------------------
# 6. Helper function: align matched cities
# ------------------------------------------------------------

align_matched_cities <- function(
  sewage_matrix,
  transit_matrix,
  analysis_label
) {

  colnames(sewage_matrix) <- standardise_city_names(
    colnames(sewage_matrix)
  )

  colnames(transit_matrix) <- standardise_city_names(
    colnames(transit_matrix)
  )

  sewage_cities <- colnames(
    sewage_matrix
  )

  transit_cities <- colnames(
    transit_matrix
  )

  matched_cities <- intersect(
    sewage_cities,
    transit_cities
  )

  matched_cities <- sort(
    matched_cities
  )

  if (length(matched_cities) == 0) {
    stop(
      paste0(
        analysis_label,
        ": no matched cities were found between sewage and transit matrices."
      )
    )
  }

  if (length(matched_cities) != expected_number_of_cities) {
    warning(
      paste0(
        analysis_label,
        ": expected ",
        expected_number_of_cities,
        " matched cities, but found ",
        length(matched_cities),
        ".\nMatched cities:\n",
        paste(matched_cities, collapse = ", ")
      )
    )
  }

  sewage_only_cities <- setdiff(
    sewage_cities,
    matched_cities
  )

  transit_only_cities <- setdiff(
    transit_cities,
    matched_cities
  )

  if (length(sewage_only_cities) > 0) {
    warning(
      paste0(
        analysis_label,
        ": sewage-only cities excluded: ",
        paste(sewage_only_cities, collapse = ", ")
      )
    )
  }

  if (length(transit_only_cities) > 0) {
    warning(
      paste0(
        analysis_label,
        ": transit-only cities excluded: ",
        paste(transit_only_cities, collapse = ", ")
      )
    )
  }

  sewage_matrix <- sewage_matrix[
    ,
    matched_cities,
    drop = FALSE
  ]

  transit_matrix <- transit_matrix[
    ,
    matched_cities,
    drop = FALSE
  ]

  return(
    list(
      sewage = sewage_matrix,
      transit = transit_matrix,
      cities = matched_cities
    )
  )
}

# ------------------------------------------------------------
# 7. Helper function: build corrected combined matrix
# ------------------------------------------------------------

build_corrected_rq1_matrix <- function(
  sewage_file,
  transit_file,
  level_name,
  output_prefix
) {

  message(
    "\n============================================================"
  )

  message(
    "Building corrected RQ1 matrix: ",
    level_name
  )

  message(
    "============================================================"
  )

  sewage_matrix <- read_binary_matrix(
    file_path = sewage_file,
    matrix_name = paste0(level_name, " sewage matrix")
  )

  transit_matrix <- read_binary_matrix(
    file_path = transit_file,
    matrix_name = paste0(level_name, " transit matrix")
  )

  aligned <- align_matched_cities(
    sewage_matrix = sewage_matrix,
    transit_matrix = transit_matrix,
    analysis_label = level_name
  )

  sewage_matrix <- aligned$sewage
  transit_matrix <- aligned$transit
  matched_cities <- aligned$cities

  # ----------------------------------------------------------
  # 7.1 Determine feature prevalence in each environment
  # ----------------------------------------------------------

  all_features <- union(
    rownames(sewage_matrix),
    rownames(transit_matrix)
  )

  all_features <- sort(
    all_features
  )

  sewage_full <- matrix(
    0,
    nrow = length(all_features),
    ncol = length(matched_cities),
    dimnames = list(
      all_features,
      matched_cities
    )
  )

  transit_full <- matrix(
    0,
    nrow = length(all_features),
    ncol = length(matched_cities),
    dimnames = list(
      all_features,
      matched_cities
    )
  )

  sewage_full[
    rownames(sewage_matrix),
    matched_cities
  ] <- sewage_matrix[
    rownames(sewage_matrix),
    matched_cities,
    drop = FALSE
  ]

  transit_full[
    rownames(transit_matrix),
    matched_cities
  ] <- transit_matrix[
    rownames(transit_matrix),
    matched_cities,
    drop = FALSE
  ]

  sewage_prevalence <- rowSums(
    sewage_full
  )

  transit_prevalence <- rowSums(
    transit_full
  )

  eligible_in_sewage <- sewage_prevalence >= min_city_prevalence
  eligible_in_transit <- transit_prevalence >= min_city_prevalence

  eligible_union <- eligible_in_sewage | eligible_in_transit

  retained_features <- names(
    eligible_union
  )[
    eligible_union
  ]

  retained_features <- sort(
    retained_features
  )

  message(
    level_name,
    ": total unfiltered feature union = ",
    length(all_features)
  )

  message(
    level_name,
    ": features passing sewage prevalence threshold = ",
    sum(eligible_in_sewage)
  )

  message(
    level_name,
    ": features passing transit prevalence threshold = ",
    sum(eligible_in_transit)
  )

  message(
    level_name,
    ": retained union after prevalence filtering = ",
    length(retained_features)
  )

  # ----------------------------------------------------------
  # 7.2 Rebuild final matrix from unfiltered data
  #
  # This is the main correction:
  # Do not combine already-filtered matrices.
  # Return to the unfiltered sewage_full and transit_full matrices.
  # ----------------------------------------------------------

  sewage_retained <- sewage_full[
    retained_features,
    ,
    drop = FALSE
  ]

  transit_retained <- transit_full[
    retained_features,
    ,
    drop = FALSE
  ]

  colnames(sewage_retained) <- paste0(
    matched_cities,
    "_Sewage"
  )

  colnames(transit_retained) <- paste0(
    matched_cities,
    "_Transit"
  )

  combined_features_by_profile <- cbind(
    sewage_retained,
    transit_retained
  )

  # Reorder profiles by city:
  # City_Sewage, City_Transit
  ordered_profiles <- as.vector(
    rbind(
      paste0(matched_cities, "_Sewage"),
      paste0(matched_cities, "_Transit")
    )
  )

  combined_features_by_profile <- combined_features_by_profile[
    retained_features,
    ordered_profiles,
    drop = FALSE
  ]

  combined_profiles_by_feature <- t(
    combined_features_by_profile
  )

  # ----------------------------------------------------------
  # 7.3 Feature classification after corrected filtering
  # ----------------------------------------------------------

  corrected_sewage_prevalence <- rowSums(
    sewage_retained
  )

  corrected_transit_prevalence <- rowSums(
    transit_retained
  )

  classification <- ifelse(
    corrected_sewage_prevalence > 0 &
      corrected_transit_prevalence > 0,
    "Shared",
    ifelse(
      corrected_sewage_prevalence > 0 &
        corrected_transit_prevalence == 0,
      "Sewage_only",
      ifelse(
        corrected_sewage_prevalence == 0 &
          corrected_transit_prevalence > 0,
        "Transit_only",
        "Absent"
      )
    )
  )

  prevalence_summary <- data.frame(
    Feature = retained_features,
    Sewage_prevalence_unfiltered = sewage_prevalence[retained_features],
    Transit_prevalence_unfiltered = transit_prevalence[retained_features],
    Eligible_in_sewage = eligible_in_sewage[retained_features],
    Eligible_in_transit = eligible_in_transit[retained_features],
    Eligible_in_union = eligible_union[retained_features],
    Final_sewage_prevalence = corrected_sewage_prevalence,
    Final_transit_prevalence = corrected_transit_prevalence,
    Classification = classification,
    stringsAsFactors = FALSE
  )

  overlap_summary <- data.frame(
    Level = level_name,
    Prevalence_threshold = prevalence_threshold,
    Min_city_prevalence = min_city_prevalence,
    Matched_cities = length(matched_cities),
    Total_unfiltered_feature_union = length(all_features),
    Retained_feature_union = length(retained_features),
    Sewage_detected_features = sum(corrected_sewage_prevalence > 0),
    Transit_detected_features = sum(corrected_transit_prevalence > 0),
    Shared_features = sum(classification == "Shared"),
    Sewage_only_features = sum(classification == "Sewage_only"),
    Transit_only_features = sum(classification == "Transit_only"),
    stringsAsFactors = FALSE
  )

  message(
    "\n",
    level_name,
    " corrected overlap summary:"
  )

  print(
    overlap_summary
  )

  # ----------------------------------------------------------
  # 7.4 Metadata table
  # ----------------------------------------------------------

  metadata <- data.frame(
    Profile = rownames(combined_profiles_by_feature),
    City = sub(
      "_(Sewage|Transit)$",
      "",
      rownames(combined_profiles_by_feature)
    ),
    Environment = sub(
      "^.*_(Sewage|Transit)$",
      "\\1",
      rownames(combined_profiles_by_feature)
    ),
    stringsAsFactors = FALSE
  )

  metadata$Environment <- factor(
    metadata$Environment,
    levels = c(
      "Sewage",
      "Transit"
    )
  )

  city_environment_table <- table(
    metadata$City,
    metadata$Environment
  )

  if (!all(city_environment_table == 1)) {
    print(
      city_environment_table
    )

    stop(
      paste0(
        level_name,
        ": final metadata does not contain exactly one sewage and one transit profile per city."
      )
    )
  }

  # ----------------------------------------------------------
  # 7.5 Write outputs
  # ----------------------------------------------------------

  features_by_profile_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_combined_binary_matrix_features_by_profile.tsv"
    )
  )

  profiles_by_feature_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_combined_binary_matrix.tsv"
    )
  )

  metadata_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_metadata.tsv"
    )
  )

  prevalence_summary_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_feature_prevalence_summary.tsv"
    )
  )

  overlap_summary_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_overlap_summary.tsv"
    )
  )

  shared_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_shared_features.tsv"
    )
  )

  sewage_only_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_sewage_only_features.tsv"
    )
  )

  transit_only_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_transit_only_features.tsv"
    )
  )

  features_by_profile_df <- data.frame(
    Feature = rownames(combined_features_by_profile),
    combined_features_by_profile,
    check.names = FALSE
  )

  profiles_by_feature_df <- data.frame(
    Profile = rownames(combined_profiles_by_feature),
    combined_profiles_by_feature,
    check.names = FALSE
  )

  write.table(
    features_by_profile_df,
    file = features_by_profile_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    profiles_by_feature_df,
    file = profiles_by_feature_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    metadata,
    file = metadata_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    prevalence_summary,
    file = prevalence_summary_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    overlap_summary,
    file = overlap_summary_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    prevalence_summary[
      prevalence_summary$Classification == "Shared",
      ,
      drop = FALSE
    ],
    file = shared_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    prevalence_summary[
      prevalence_summary$Classification == "Sewage_only",
      ,
      drop = FALSE
    ],
    file = sewage_only_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    prevalence_summary[
      prevalence_summary$Classification == "Transit_only",
      ,
      drop = FALSE
    ],
    file = transit_only_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message(
    "\nFiles written for ",
    level_name,
    ":"
  )

  message(
    features_by_profile_output
  )

  message(
    profiles_by_feature_output
  )

  message(
    metadata_output
  )

  message(
    prevalence_summary_output
  )

  message(
    overlap_summary_output
  )

  return(
    list(
      combined_features_by_profile = combined_features_by_profile,
      combined_profiles_by_feature = combined_profiles_by_feature,
      metadata = metadata,
      prevalence_summary = prevalence_summary,
      overlap_summary = overlap_summary
    )
  )
}

# ------------------------------------------------------------
# 8. Build ARG Group matrix
# ------------------------------------------------------------

group_results <- build_corrected_rq1_matrix(
  sewage_file = sewage_group_file,
  transit_file = transit_group_file,
  level_name = "ARG Group",
  output_prefix = "rq1_group"
)

# ------------------------------------------------------------
# 9. Build Resistance Class matrix
# ------------------------------------------------------------

class_results <- build_corrected_rq1_matrix(
  sewage_file = sewage_class_file,
  transit_file = transit_class_file,
  level_name = "Resistance Class",
  output_prefix = "rq1_class"
)

# ------------------------------------------------------------
# 10. Save combined run summary
# ------------------------------------------------------------

run_summary <- rbind(
  group_results$overlap_summary,
  class_results$overlap_summary
)

run_summary_output <- file.path(
  output_dir,
  "rq1_binary_matrix_build_summary.tsv"
)

write.table(
  run_summary,
  file = run_summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "\n============================================================"
)

message(
  "RQ1 binary matrix construction completed."
)

message(
  "Summary written to:"
)

message(
  run_summary_output
)

message(
  "============================================================"
)

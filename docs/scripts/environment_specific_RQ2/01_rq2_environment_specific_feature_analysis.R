# ============================================================
# 01_rq2_environment_specific_feature_analysis.R
#
# Purpose:
# Identify strict environment-specific and shared resistome
# features for RQ2.
#
# This script is descriptive. It does not perform formal
# statistical testing.
#
# Strict classification:
# - Shared: detected in at least one sewage city and at least
#   one transit city
# - Sewage-only: detected in sewage and absent from all transit
#   profiles
# - Transit-only: detected in transit and absent from all sewage
#   profiles
#
# Input:
# results/rq1/matrices/rq1_group_combined_binary_matrix.tsv
# results/rq1/matrices/rq1_class_combined_binary_matrix.tsv
#
# Expected input matrix orientation:
# Rows    = city-environment profiles
# Columns = ARG features
#
# Example profile names:
# Berlin_Sewage
# Berlin_Transit
#
# Output:
# results/rq2/environment_specific/
# ============================================================

# ------------------------------------------------------------
# 1. User-adjustable settings
# ------------------------------------------------------------

input_dir <- file.path(
  "results",
  "rq1",
  "matrices"
)

output_dir <- file.path(
  "results",
  "rq2",
  "environment_specific"
)

expected_number_of_cities <- 16

group_matrix_file <- file.path(
  input_dir,
  "rq1_group_combined_binary_matrix.tsv"
)

class_matrix_file <- file.path(
  input_dir,
  "rq1_class_combined_binary_matrix.tsv"
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
  group_matrix_file,
  class_matrix_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "The following required input files are missing:\n",
      paste(missing_input_files, collapse = "\n"),
      "\nRun 01_build_rq1_binary_matrices.R before running this RQ2 script."
    )
  )
}

# ------------------------------------------------------------
# 4. Helper function: read combined binary matrix
# ------------------------------------------------------------

read_combined_binary_matrix <- function(file_path, level_name) {

  message(
    "\nReading ",
    level_name,
    " combined binary matrix:\n",
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
        level_name,
        ": matrix must contain one profile column and at least one feature column."
      )
    )
  }

  profile_ids <- trimws(
    as.character(raw_df[[1]])
  )

  if (any(is.na(profile_ids)) || any(profile_ids == "")) {
    stop(
      paste0(
        level_name,
        ": missing or empty profile names detected in the first column."
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

  rownames(value_matrix) <- profile_ids

  if (any(is.na(value_matrix))) {
    stop(
      paste0(
        level_name,
        ": matrix contains missing or non-numeric values."
      )
    )
  }

  value_matrix <- ifelse(
    value_matrix > 0,
    1,
    0
  )

  storage.mode(value_matrix) <- "numeric"

  if (any(duplicated(rownames(value_matrix)))) {
    stop(
      paste0(
        level_name,
        ": duplicated profile names detected."
      )
    )
  }

  if (any(duplicated(colnames(value_matrix)))) {
    stop(
      paste0(
        level_name,
        ": duplicated feature names detected."
      )
    )
  }

  message(
    level_name,
    ": ",
    nrow(value_matrix),
    " profiles × ",
    ncol(value_matrix),
    " features"
  )

  return(value_matrix)
}

# ------------------------------------------------------------
# 5. Helper function: build metadata from profile names
# ------------------------------------------------------------

build_profile_metadata <- function(binary_matrix, level_name) {

  metadata <- data.frame(
    Profile = rownames(binary_matrix),
    stringsAsFactors = FALSE
  )

  metadata$Environment <- sub(
    "^.*_(Sewage|Transit)$",
    "\\1",
    metadata$Profile
  )

  metadata$City <- sub(
    "_(Sewage|Transit)$",
    "",
    metadata$Profile
  )

  if (!all(metadata$Environment %in% c("Sewage", "Transit"))) {
    stop(
      paste0(
        level_name,
        ": environment could not be extracted from some profile names:\n",
        paste(metadata$Profile, collapse = "\n")
      )
    )
  }

  metadata$Environment <- factor(
    metadata$Environment,
    levels = c(
      "Sewage",
      "Transit"
    )
  )

  rownames(metadata) <- metadata$Profile

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
        ": each city must contain exactly one sewage and one transit profile."
      )
    )
  }

  if (length(unique(metadata$City)) != expected_number_of_cities) {
    warning(
      paste0(
        level_name,
        ": expected ",
        expected_number_of_cities,
        " matched cities, but found ",
        length(unique(metadata$City)),
        "."
      )
    )
  }

  metadata <- metadata[
    order(
      metadata$City,
      metadata$Environment
    ),
    ,
    drop = FALSE
  ]

  return(metadata)
}

# ------------------------------------------------------------
# 6. Helper function: make short feature labels
# ------------------------------------------------------------

make_short_feature_label <- function(feature_name) {

  label <- feature_name

  # If labels contain comma-separated MEGARes hierarchy,
  # use the final element as the display label.
  label <- sub(
    "^.*,",
    "",
    label
  )

  label <- trimws(
    label
  )

  return(label)
}

# ------------------------------------------------------------
# 7. Helper function: write feature list
# ------------------------------------------------------------

write_feature_list <- function(feature_df, file_path) {

  write.table(
    feature_df,
    file = file_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 8. Main function: strict environment-specific analysis
# ------------------------------------------------------------

run_environment_specific_analysis <- function(
  matrix_file,
  level_name,
  output_prefix
) {

  message(
    "\n============================================================"
  )

  message(
    "Running RQ2 strict environment-specific analysis: ",
    level_name
  )

  message(
    "============================================================"
  )

  binary_matrix <- read_combined_binary_matrix(
    file_path = matrix_file,
    level_name = level_name
  )

  metadata <- build_profile_metadata(
    binary_matrix = binary_matrix,
    level_name = level_name
  )

  binary_matrix <- binary_matrix[
    metadata$Profile,
    ,
    drop = FALSE
  ]

  sewage_profiles <- metadata$Profile[
    metadata$Environment == "Sewage"
  ]

  transit_profiles <- metadata$Profile[
    metadata$Environment == "Transit"
  ]

  sewage_matrix <- binary_matrix[
    sewage_profiles,
    ,
    drop = FALSE
  ]

  transit_matrix <- binary_matrix[
    transit_profiles,
    ,
    drop = FALSE
  ]

  number_of_cities <- length(
    unique(
      metadata$City
    )
  )

  # ----------------------------------------------------------
  # 8.1 Calculate detection counts and strict classification
  # ----------------------------------------------------------

  sewage_city_count <- colSums(
    sewage_matrix
  )

  transit_city_count <- colSums(
    transit_matrix
  )

  present_in_sewage <- sewage_city_count > 0
  present_in_transit <- transit_city_count > 0

  strict_classification <- ifelse(
    present_in_sewage & present_in_transit,
    "Shared",
    ifelse(
      present_in_sewage & !present_in_transit,
      "Sewage_only",
      ifelse(
        !present_in_sewage & present_in_transit,
        "Transit_only",
        "Absent"
      )
    )
  )

  feature_classification <- data.frame(
    Feature = colnames(binary_matrix),
    Feature_label = vapply(
      colnames(binary_matrix),
      make_short_feature_label,
      FUN.VALUE = character(1)
    ),
    Level = level_name,
    Sewage_city_count = as.numeric(
      sewage_city_count
    ),
    Transit_city_count = as.numeric(
      transit_city_count
    ),
    Sewage_prevalence = as.numeric(
      sewage_city_count / number_of_cities
    ),
    Transit_prevalence = as.numeric(
      transit_city_count / number_of_cities
    ),
    Present_in_sewage = as.logical(
      present_in_sewage
    ),
    Present_in_transit = as.logical(
      present_in_transit
    ),
    Strict_environment_classification = strict_classification,
    stringsAsFactors = FALSE
  )

  feature_classification <- feature_classification[
    order(
      feature_classification$Strict_environment_classification,
      -feature_classification$Sewage_city_count,
      -feature_classification$Transit_city_count,
      feature_classification$Feature
    ),
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 8.2 Extract feature groups
  # ----------------------------------------------------------

  shared_features <- feature_classification[
    feature_classification$Strict_environment_classification == "Shared",
    ,
    drop = FALSE
  ]

  sewage_only_features <- feature_classification[
    feature_classification$Strict_environment_classification == "Sewage_only",
    ,
    drop = FALSE
  ]

  transit_only_features <- feature_classification[
    feature_classification$Strict_environment_classification == "Transit_only",
    ,
    drop = FALSE
  ]

  absent_features <- feature_classification[
    feature_classification$Strict_environment_classification == "Absent",
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 8.3 Environment-level summary
  # ----------------------------------------------------------

  summary_df <- data.frame(
    Level = level_name,
    Matched_cities = number_of_cities,
    Total_features = ncol(binary_matrix),
    Sewage_detected_features = sum(
      present_in_sewage
    ),
    Transit_detected_features = sum(
      present_in_transit
    ),
    Shared_features = nrow(
      shared_features
    ),
    Strict_sewage_only_features = nrow(
      sewage_only_features
    ),
    Strict_transit_only_features = nrow(
      transit_only_features
    ),
    Absent_features = nrow(
      absent_features
    ),
    Fraction_shared_of_all_features = nrow(shared_features) /
      ncol(binary_matrix),
    Fraction_sewage_only_of_all_features = nrow(sewage_only_features) /
      ncol(binary_matrix),
    Fraction_transit_only_of_all_features = nrow(transit_only_features) /
      ncol(binary_matrix),
    stringsAsFactors = FALSE
  )

  message(
    "\nStrict environment-specific summary for ",
    level_name,
    ":"
  )

  print(
    summary_df
  )

  # ----------------------------------------------------------
  # 8.4 City-level summary
  # ----------------------------------------------------------

  city_level_summary <- data.frame(
    City = sort(
      unique(
        metadata$City
      )
    ),
    stringsAsFactors = FALSE
  )

  sewage_by_city <- sewage_matrix
  transit_by_city <- transit_matrix

  rownames(sewage_by_city) <- metadata[
    rownames(sewage_by_city),
    "City"
  ]

  rownames(transit_by_city) <- metadata[
    rownames(transit_by_city),
    "City"
  ]

  sewage_by_city <- sewage_by_city[
    city_level_summary$City,
    ,
    drop = FALSE
  ]

  transit_by_city <- transit_by_city[
    city_level_summary$City,
    ,
    drop = FALSE
  ]

  city_level_summary$Sewage_feature_richness <- rowSums(
    sewage_by_city
  )

  city_level_summary$Transit_feature_richness <- rowSums(
    transit_by_city
  )

  city_level_summary$Shared_within_city_features <- rowSums(
    sewage_by_city * transit_by_city
  )

  city_level_summary$Sewage_only_within_city_features <- rowSums(
    sewage_by_city == 1 & transit_by_city == 0
  )

  city_level_summary$Transit_only_within_city_features <- rowSums(
    sewage_by_city == 0 & transit_by_city == 1
  )

  city_level_summary$Within_city_jaccard_similarity <- with(
    city_level_summary,
    Shared_within_city_features /
      (
        Sewage_feature_richness +
          Transit_feature_richness -
          Shared_within_city_features
      )
  )

  city_level_summary$Level <- level_name

  city_level_summary <- city_level_summary[
    ,
    c(
      "Level",
      "City",
      "Sewage_feature_richness",
      "Transit_feature_richness",
      "Shared_within_city_features",
      "Sewage_only_within_city_features",
      "Transit_only_within_city_features",
      "Within_city_jaccard_similarity"
    ),
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 8.5 Write outputs
  # ----------------------------------------------------------

  classification_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_environment_specific_feature_classification.tsv"
    )
  )

  summary_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_environment_specific_summary.tsv"
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
      "_strict_sewage_only_features.tsv"
    )
  )

  transit_only_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_strict_transit_only_features.tsv"
    )
  )

  city_summary_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_city_level_environment_specific_summary.tsv"
    )
  )

  write.table(
    feature_classification,
    file = classification_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    summary_df,
    file = summary_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write_feature_list(
    feature_df = shared_features,
    file_path = shared_output
  )

  write_feature_list(
    feature_df = sewage_only_features,
    file_path = sewage_only_output
  )

  write_feature_list(
    feature_df = transit_only_features,
    file_path = transit_only_output
  )

  write.table(
    city_level_summary,
    file = city_summary_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message(
    "\nFiles written for ",
    level_name,
    ":"
  )

  message(classification_output)
  message(summary_output)
  message(shared_output)
  message(sewage_only_output)
  message(transit_only_output)
  message(city_summary_output)

  return(
    list(
      feature_classification = feature_classification,
      summary = summary_df,
      city_level_summary = city_level_summary
    )
  )
}

# ------------------------------------------------------------
# 9. Run ARG Group analysis
# ------------------------------------------------------------

group_results <- run_environment_specific_analysis(
  matrix_file = group_matrix_file,
  level_name = "ARG Group",
  output_prefix = "rq2_group"
)

# ------------------------------------------------------------
# 10. Run Resistance Class analysis
# ------------------------------------------------------------

class_results <- run_environment_specific_analysis(
  matrix_file = class_matrix_file,
  level_name = "Resistance Class",
  output_prefix = "rq2_class"
)

# ------------------------------------------------------------
# 11. Save combined summaries
# ------------------------------------------------------------

combined_summary <- rbind(
  group_results$summary,
  class_results$summary
)

combined_city_summary <- rbind(
  group_results$city_level_summary,
  class_results$city_level_summary
)

combined_summary_output <- file.path(
  output_dir,
  "rq2_combined_environment_specific_summary.tsv"
)

combined_city_summary_output <- file.path(
  output_dir,
  "rq2_combined_city_level_environment_specific_summary.tsv"
)

write.table(
  combined_summary,
  file = combined_summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  combined_city_summary,
  file = combined_city_summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "\n============================================================"
)

message(
  "RQ2 strict environment-specific feature analysis completed."
)

message(
  "Combined outputs:"
)

message(combined_summary_output)
message(combined_city_summary_output)

message(
  "============================================================"
)

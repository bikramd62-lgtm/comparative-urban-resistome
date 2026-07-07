# ============================================================
# 02_rq2_prevalence_difference_analysis.R
#
# Purpose:
# Calculate descriptive prevalence differences between sewage
# and transit resistome features for RQ2.
#
# This script is descriptive. It does not perform formal
# statistical testing
#
# Prevalence definition:
# Prevalence is calculated as the proportion of matched cities
# in which a feature is detected in each environment.
#
# Prevalence difference:
# Transit prevalence - Sewage prevalence
#
# Interpretation:
# - Positive value: higher city-level prevalence in transit
# - Negative value: higher city-level prevalence in sewage
# - Zero value: equal city-level prevalence
#
# Input:
# results/rq2/environment_specific/rq2_group_environment_specific_feature_classification.tsv
# results/rq2/environment_specific/rq2_class_environment_specific_feature_classification.tsv
#
# Output:
# results/rq2/prevalence_difference/
# ============================================================

# ------------------------------------------------------------
# 1. User-adjustable settings
# ------------------------------------------------------------

input_dir <- file.path(
  "results",
  "rq2",
  "environment_specific"
)

output_dir <- file.path(
  "results",
  "rq2",
  "prevalence_difference"
)

expected_number_of_cities <- 16

top_n_features <- 30

group_input_file <- file.path(
  input_dir,
  "rq2_group_environment_specific_feature_classification.tsv"
)

class_input_file <- file.path(
  input_dir,
  "rq2_class_environment_specific_feature_classification.tsv"
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
  group_input_file,
  class_input_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "The following required input files are missing:\n",
      paste(missing_input_files, collapse = "\n"),
      "\nRun 01_rq2_environment_specific_feature_analysis.R before running this script."
    )
  )
}

# ------------------------------------------------------------
# 4. Helper function: read RQ2 feature-classification table
# ------------------------------------------------------------

read_feature_classification <- function(file_path, level_name) {

  message(
    "\nReading ",
    level_name,
    " feature-classification table:\n",
    file_path
  )

  feature_table <- read.delim(
    file = file_path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  required_columns <- c(
    "Feature",
    "Feature_label",
    "Level",
    "Sewage_city_count",
    "Transit_city_count",
    "Sewage_prevalence",
    "Transit_prevalence",
    "Present_in_sewage",
    "Present_in_transit",
    "Strict_environment_classification"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(feature_table)
  )

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        level_name,
        ": the following required columns are missing:\n",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  numeric_columns <- c(
    "Sewage_city_count",
    "Transit_city_count",
    "Sewage_prevalence",
    "Transit_prevalence"
  )

  for (column_name in numeric_columns) {
    feature_table[[column_name]] <- suppressWarnings(
      as.numeric(
        feature_table[[column_name]]
      )
    )

    if (any(is.na(feature_table[[column_name]]))) {
      stop(
        paste0(
          level_name,
          ": column contains missing or non-numeric values: ",
          column_name
        )
      )
    }
  }

  if (any(feature_table$Sewage_prevalence < 0) ||
      any(feature_table$Sewage_prevalence > 1) ||
      any(feature_table$Transit_prevalence < 0) ||
      any(feature_table$Transit_prevalence > 1)) {
    stop(
      paste0(
        level_name,
        ": prevalence values must be between 0 and 1."
      )
    )
  }

  message(
    level_name,
    ": ",
    nrow(feature_table),
    " features"
  )

  return(feature_table)
}

# ------------------------------------------------------------
# 5. Helper function: classify prevalence-difference direction
# ------------------------------------------------------------

classify_prevalence_direction <- function(prevalence_difference) {

  direction <- ifelse(
    prevalence_difference > 0,
    "Higher_in_transit",
    ifelse(
      prevalence_difference < 0,
      "Higher_in_sewage",
      "Equal_prevalence"
    )
  )

  return(direction)
}

# ------------------------------------------------------------
# 6. Helper function: label prevalence direction for figures
# ------------------------------------------------------------

make_direction_label <- function(prevalence_direction) {

  label <- ifelse(
    prevalence_direction == "Higher_in_transit",
    "Higher in transit",
    ifelse(
      prevalence_direction == "Higher_in_sewage",
      "Higher in sewage",
      "Equal prevalence"
    )
  )

  return(label)
}

# ------------------------------------------------------------
# 7. Main function: prevalence-difference analysis
# ------------------------------------------------------------

run_prevalence_difference_analysis <- function(
  input_file,
  level_name,
  output_prefix
) {

  message(
    "\n============================================================"
  )

  message(
    "Running RQ2 descriptive prevalence-difference analysis: ",
    level_name
  )

  message(
    "============================================================"
  )

  feature_table <- read_feature_classification(
    file_path = input_file,
    level_name = level_name
  )

  # ----------------------------------------------------------
  # 7.1 Calculate prevalence differences
  # ----------------------------------------------------------

  feature_table$Prevalence_difference <- feature_table$Transit_prevalence -
    feature_table$Sewage_prevalence

  feature_table$Absolute_prevalence_difference <- abs(
    feature_table$Prevalence_difference
  )

  feature_table$Prevalence_direction <- classify_prevalence_direction(
    feature_table$Prevalence_difference
  )

  feature_table$Prevalence_direction_label <- make_direction_label(
    feature_table$Prevalence_direction
  )

  feature_table$Maximum_environment_prevalence <- pmax(
    feature_table$Sewage_prevalence,
    feature_table$Transit_prevalence
  )

  feature_table$Minimum_environment_prevalence <- pmin(
    feature_table$Sewage_prevalence,
    feature_table$Transit_prevalence
  )

  feature_table$Both_environment_prevalence_sum <- feature_table$Sewage_prevalence +
    feature_table$Transit_prevalence

  # ----------------------------------------------------------
  # 7.2 Add readable percentage columns
  # ----------------------------------------------------------

  feature_table$Sewage_prevalence_percent <- feature_table$Sewage_prevalence * 100

  feature_table$Transit_prevalence_percent <- feature_table$Transit_prevalence * 100

  feature_table$Prevalence_difference_percent_points <- feature_table$Prevalence_difference * 100

  feature_table$Absolute_prevalence_difference_percent_points <-
    feature_table$Absolute_prevalence_difference * 100

  # ----------------------------------------------------------
  # 7.3 Reorder output columns
  # ----------------------------------------------------------

  output_columns <- c(
    "Feature",
    "Feature_label",
    "Level",
    "Sewage_city_count",
    "Transit_city_count",
    "Sewage_prevalence",
    "Transit_prevalence",
    "Sewage_prevalence_percent",
    "Transit_prevalence_percent",
    "Prevalence_difference",
    "Prevalence_difference_percent_points",
    "Absolute_prevalence_difference",
    "Absolute_prevalence_difference_percent_points",
    "Prevalence_direction",
    "Prevalence_direction_label",
    "Maximum_environment_prevalence",
    "Minimum_environment_prevalence",
    "Both_environment_prevalence_sum",
    "Present_in_sewage",
    "Present_in_transit",
    "Strict_environment_classification"
  )

  feature_table <- feature_table[
    ,
    output_columns,
    drop = FALSE
  ]

  feature_table <- feature_table[
    order(
      -feature_table$Absolute_prevalence_difference,
      feature_table$Prevalence_direction,
      -feature_table$Maximum_environment_prevalence,
      feature_table$Feature_label,
      feature_table$Feature
    ),
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # 7.4 Create direction-specific tables
  # ----------------------------------------------------------

  higher_in_sewage <- feature_table[
    feature_table$Prevalence_direction == "Higher_in_sewage",
    ,
    drop = FALSE
  ]

  higher_in_transit <- feature_table[
    feature_table$Prevalence_direction == "Higher_in_transit",
    ,
    drop = FALSE
  ]

  equal_prevalence <- feature_table[
    feature_table$Prevalence_direction == "Equal_prevalence",
    ,
    drop = FALSE
  ]

  top_prevalence_difference <- feature_table[
    order(
      -feature_table$Absolute_prevalence_difference,
      -feature_table$Maximum_environment_prevalence,
      feature_table$Feature_label,
      feature_table$Feature
    ),
    ,
    drop = FALSE
  ]

  if (nrow(top_prevalence_difference) > top_n_features) {
    top_prevalence_difference <- top_prevalence_difference[
      seq_len(top_n_features),
      ,
      drop = FALSE
    ]
  }

  # ----------------------------------------------------------
  # 7.5 Summary table
  # ----------------------------------------------------------

  summary_df <- data.frame(
    Level = level_name,
    Total_features = nrow(feature_table),
    Higher_in_sewage_features = nrow(higher_in_sewage),
    Higher_in_transit_features = nrow(higher_in_transit),
    Equal_prevalence_features = nrow(equal_prevalence),
    Strict_shared_features = sum(
      feature_table$Strict_environment_classification == "Shared"
    ),
    Strict_sewage_only_features = sum(
      feature_table$Strict_environment_classification == "Sewage_only"
    ),
    Strict_transit_only_features = sum(
      feature_table$Strict_environment_classification == "Transit_only"
    ),
    Maximum_absolute_prevalence_difference =
      max(feature_table$Absolute_prevalence_difference),
    Median_absolute_prevalence_difference =
      median(feature_table$Absolute_prevalence_difference),
    Mean_absolute_prevalence_difference =
      mean(feature_table$Absolute_prevalence_difference),
    stringsAsFactors = FALSE
  )

  message(
    "\nDescriptive prevalence-difference summary for ",
    level_name,
    ":"
  )

  print(
    summary_df
  )

  # ----------------------------------------------------------
  # 7.6 Write outputs
  # ----------------------------------------------------------

  prevalence_comparison_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_prevalence_comparison.tsv"
    )
  )

  summary_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_prevalence_difference_summary.tsv"
    )
  )

  higher_in_sewage_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_higher_prevalence_in_sewage.tsv"
    )
  )

  higher_in_transit_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_higher_prevalence_in_transit.tsv"
    )
  )

  equal_prevalence_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_equal_prevalence_features.tsv"
    )
  )

  top_prevalence_difference_output <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_top",
      top_n_features,
      "_absolute_prevalence_difference_features.tsv"
    )
  )

  write.table(
    feature_table,
    file = prevalence_comparison_output,
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

  write.table(
    higher_in_sewage,
    file = higher_in_sewage_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    higher_in_transit,
    file = higher_in_transit_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    equal_prevalence,
    file = equal_prevalence_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    top_prevalence_difference,
    file = top_prevalence_difference_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message(
    "\nFiles written for ",
    level_name,
    ":"
  )

  message(prevalence_comparison_output)
  message(summary_output)
  message(higher_in_sewage_output)
  message(higher_in_transit_output)
  message(equal_prevalence_output)
  message(top_prevalence_difference_output)

  return(
    list(
      prevalence_comparison = feature_table,
      summary = summary_df,
      higher_in_sewage = higher_in_sewage,
      higher_in_transit = higher_in_transit,
      equal_prevalence = equal_prevalence,
      top_prevalence_difference = top_prevalence_difference
    )
  )
}

# ------------------------------------------------------------
# 8. Run ARG Group analysis
# ------------------------------------------------------------

group_results <- run_prevalence_difference_analysis(
  input_file = group_input_file,
  level_name = "ARG Group",
  output_prefix = "rq2_group"
)

# ------------------------------------------------------------
# 9. Run Resistance Class analysis
# ------------------------------------------------------------

class_results <- run_prevalence_difference_analysis(
  input_file = class_input_file,
  level_name = "Resistance Class",
  output_prefix = "rq2_class"
)

# ------------------------------------------------------------
# 10. Save combined summaries
# ------------------------------------------------------------

combined_summary <- rbind(
  group_results$summary,
  class_results$summary
)

combined_prevalence_comparison <- rbind(
  group_results$prevalence_comparison,
  class_results$prevalence_comparison
)

combined_summary_output <- file.path(
  output_dir,
  "rq2_combined_prevalence_difference_summary.tsv"
)

combined_prevalence_comparison_output <- file.path(
  output_dir,
  "rq2_combined_prevalence_comparison.tsv"
)

write.table(
  combined_summary,
  file = combined_summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  combined_prevalence_comparison,
  file = combined_prevalence_comparison_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "\n============================================================"
)

message(
  "RQ2 descriptive prevalence-difference analysis completed."
)

message(
  "Combined outputs:"
)

message(combined_summary_output)
message(combined_prevalence_comparison_output)

message(
  "============================================================"
)

# ============================================================
# 03_pcoa1_richness_count_sample_number_correlations.R
#
# Purpose:
#   Test whether ARG-group PCoA axes are associated with
#   ARG-group richness, total ARG count, and sample number.
#
# Analyses:
#   1. Combined sewage-transit profiles:
#      PCoA1 and PCoA2 versus:
#        - ARG-group richness
#        - total ARG count
#        - number of contributing samples
#
#   2. Sewage-only profiles:
#      PCoA1 and PCoA2 versus:
#        - ARG-group richness
#        - total ARG count
#        - number of contributing samples
#
#   3. Transit-only profiles:
#      PCoA1 and PCoA2 versus:
#        - ARG-group richness
#        - total ARG count
#
# Multiple testing:
#   All valid correlations are treated as one testing family.
#   Benjamini-Hochberg correction is applied across all tests.
#
# Output:
#   - Spearman correlation table
#   - BH-adjusted correlation table
#   - Input validation table
#   - Scatterplots for all valid correlations
#
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following required R package(s) are missing:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nInstall them before running this script."
    )
  )
}

library(ggplot2)


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------


project_root <- Sys.getenv(
  "PROJECT_ROOT",
  unset = getwd()
)

input_dir <- file.path(
  project_root,
  "data",
  "processed",
  "sensitivity_analysis",
  "pcoa1_correlations"
)

output_dir <- file.path(
  project_root,
  "results",
  "sensitivity_analysis",
  "pcoa1_correlations"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

figures_dir <- file.path(
  output_dir,
  "figures"
)

diagnostics_dir <- file.path(
  output_dir,
  "diagnostics"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_dir)) {
  stop(
    paste0(
      "Input directory not found:\n",
      input_dir,
      "\n\nPlace the PCoA1 correlation input files in this folder."
    )
  )
}


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

combined_file <- file.path(
  input_dir,
  "Combined_ARG_group_PCoA_coordinates_metrics.csv"
)

sewage_only_file <- file.path(
  input_dir,
  "Sewage_only_ARG_group_PCoA_coordinates_metrics.csv"
)

sample_number_file <- file.path(
  input_dir,
  "sample_number_metadata_used.csv"
)

all_coordinates_file <- file.path(
  input_dir,
  "ALL_ARG_group_PCoA_coordinates_metrics.csv"
)

transit_only_file <- file.path(
  input_dir,
  "Transit_only_ARG_group_PCoA_coordinates_metrics.csv"
)

required_input_files <- c(
  combined_file,
  sewage_only_file,
  sample_number_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "The following required input file(s) were not found:\n",
      paste(missing_input_files, collapse = "\n")
    )
  )
}


# ------------------------------------------------------------
# 4. General settings
# ------------------------------------------------------------

significance_threshold <- 0.05

make_scatterplots <- TRUE

point_size <- 3

plot_width <- 6.5
plot_height <- 5.5
plot_dpi <- 300


# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

clean_column_name <- function(x) {

  tolower(
    gsub(
      pattern = "[^a-zA-Z0-9]",
      replacement = "",
      x = x
    )
  )
}


find_column <- function(
    data,
    candidate_names,
    required = TRUE,
    data_label = "data"
) {

  cleaned_names <- clean_column_name(
    colnames(data)
  )

  cleaned_candidates <- clean_column_name(
    candidate_names
  )

  matches <- which(
    cleaned_names %in% cleaned_candidates
  )

  if (length(matches) == 0) {

    if (required) {
      stop(
        paste0(
          data_label,
          ": could not find any of the following columns:\n",
          paste(candidate_names, collapse = ", "),
          "\n\nAvailable columns:\n",
          paste(colnames(data), collapse = ", ")
        )
      )
    } else {
      return(NULL)
    }
  }

  if (length(matches) > 1) {
    warning(
      paste0(
        data_label,
        ": multiple candidate columns were found for ",
        paste(candidate_names, collapse = "/"),
        ". Using: ",
        colnames(data)[matches[1]]
      )
    )
  }

  colnames(data)[matches[1]]
}


read_csv_base <- function(file_path) {

  read.csv(
    file_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


standardise_coordinate_table <- function(
    data,
    data_label,
    expected_environment = NULL
) {

  sample_col <- find_column(
    data = data,
    candidate_names = c(
      "SampleID",
      "Sample_ID",
      "sample_id",
      "profile_id",
      "ProfileID",
      "Profile",
      "ID"
    ),
    required = TRUE,
    data_label = data_label
  )

  city_col <- find_column(
    data = data,
    candidate_names = c(
      "City",
      "city"
    ),
    required = FALSE,
    data_label = data_label
  )

  environment_col <- find_column(
    data = data,
    candidate_names = c(
      "Environment",
      "environment"
    ),
    required = FALSE,
    data_label = data_label
  )

  pcoa1_col <- find_column(
    data = data,
    candidate_names = c(
      "PCoA1",
      "PC1",
      "Axis1",
      "Axis_1",
      "MDS1",
      "Dim1"
    ),
    required = TRUE,
    data_label = data_label
  )

  pcoa2_col <- find_column(
    data = data,
    candidate_names = c(
      "PCoA2",
      "PC2",
      "Axis2",
      "Axis_2",
      "MDS2",
      "Dim2"
    ),
    required = TRUE,
    data_label = data_label
  )

  richness_col <- find_column(
    data = data,
    candidate_names = c(
      "ARG_richness",
      "ARG_group_richness",
      "ARG_Group_richness",
      "ARGGroupRichness",
      "ARG_groups_richness",
      "ARG_group_count",
      "Richness",
      "richness"
    ),
    required = FALSE,
    data_label = data_label
  )

  total_col <- find_column(
    data = data,
    candidate_names = c(
      "Total_ARG_count",
      "total_ARG_count",
      "TotalARGcount",
      "ARG_count",
      "ARG_counts",
      "Total_count",
      "total_count",
      "Total_ARG_aligned_count",
      "Total_aligned_ARG_count"
    ),
    required = FALSE,
    data_label = data_label
  )

  sample_number_col <- find_column(
    data = data,
    candidate_names = c(
      "Number_of_samples",
      "number_of_samples",
      "NumberOfSamples",
      "Sample_number",
      "sample_number",
      "n_samples",
      "N_samples",
      "Sample_count",
      "contributing_samples",
      "Sewage_sample_number"
    ),
    required = FALSE,
    data_label = data_label
  )

  output <- data.frame(
    SampleID = trimws(as.character(data[[sample_col]])),
    PCoA1 = as.numeric(as.character(data[[pcoa1_col]])),
    PCoA2 = as.numeric(as.character(data[[pcoa2_col]])),
    stringsAsFactors = FALSE
  )

  if (!is.null(city_col)) {
    output$City <- trimws(as.character(data[[city_col]]))
  } else {
    output$City <- NA_character_
  }

  if (!is.null(environment_col)) {

    environment_values <- trimws(
      as.character(data[[environment_col]])
    )

    environment_values[
      tolower(environment_values) == "sewage"
    ] <- "Sewage"

    environment_values[
      tolower(environment_values) == "transit"
    ] <- "Transit"

    output$Environment <- environment_values

  } else if (!is.null(expected_environment)) {

    output$Environment <- expected_environment

  } else {

    output$Environment <- NA_character_
  }

  if (!is.null(richness_col)) {
    output$ARG_richness <- as.numeric(
      as.character(data[[richness_col]])
    )
  } else {
    output$ARG_richness <- NA_real_
  }

  if (!is.null(total_col)) {
    output$Total_ARG_count <- as.numeric(
      as.character(data[[total_col]])
    )
  } else {
    output$Total_ARG_count <- NA_real_
  }

  if (!is.null(sample_number_col)) {
    output$Number_of_samples <- as.numeric(
      as.character(data[[sample_number_col]])
    )
  } else {
    output$Number_of_samples <- NA_real_
  }

  if (anyNA(output$PCoA1) || anyNA(output$PCoA2)) {
    stop(
      paste0(
        data_label,
        ": PCoA1 or PCoA2 contains missing/non-numeric values."
      )
    )
  }

  output
}


merge_missing_metrics <- function(
    coordinate_data,
    metadata_data,
    data_label
) {

  metric_columns <- c(
    "ARG_richness",
    "Total_ARG_count",
    "Number_of_samples"
  )

  metadata_available <- metric_columns[
    metric_columns %in% colnames(metadata_data)
  ]

  if (length(metadata_available) == 0) {
    return(coordinate_data)
  }

  # Prefer merging by SampleID when available.

  if (
    "SampleID" %in% colnames(coordinate_data) &&
      "SampleID" %in% colnames(metadata_data)
  ) {

    metadata_unique <- metadata_data[
      !duplicated(metadata_data$SampleID),
      ,
      drop = FALSE
    ]

    matched_index <- match(
      coordinate_data$SampleID,
      metadata_unique$SampleID
    )

  } else if (
    all(c("City", "Environment") %in% colnames(coordinate_data)) &&
      all(c("City", "Environment") %in% colnames(metadata_data))
  ) {

    coordinate_key <- paste(
      coordinate_data$City,
      coordinate_data$Environment,
      sep = "___"
    )

    metadata_key <- paste(
      metadata_data$City,
      metadata_data$Environment,
      sep = "___"
    )

    metadata_unique <- metadata_data[
      !duplicated(metadata_key),
      ,
      drop = FALSE
    ]

    metadata_key_unique <- metadata_key[
      !duplicated(metadata_key)
    ]

    matched_index <- match(
      coordinate_key,
      metadata_key_unique
    )

  } else {

    warning(
      paste0(
        data_label,
        ": could not merge external metadata because no suitable ",
        "SampleID or City+Environment key was available."
      )
    )

    return(coordinate_data)
  }

  for (metric in metric_columns) {

    if (
      metric %in% colnames(metadata_unique) &&
        (
          !metric %in% colnames(coordinate_data) ||
            all(is.na(coordinate_data[[metric]]))
        )
    ) {

      coordinate_data[[metric]] <- metadata_unique[
        matched_index,
        metric
      ]
    }
  }

  coordinate_data
}


validate_analysis_data <- function(
    data,
    analysis_label
) {

  required_columns <- c(
    "SampleID",
    "PCoA1",
    "PCoA2",
    "ARG_richness",
    "Total_ARG_count"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        analysis_label,
        ": missing required columns after standardisation:\n",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  if (anyDuplicated(data$SampleID) > 0) {
    warning(
      paste0(
        analysis_label,
        ": duplicate SampleID values were found."
      )
    )
  }

  validation <- data.frame(
    analysis = analysis_label,
    n_profiles = nrow(data),
    n_missing_PCoA1 = sum(is.na(data$PCoA1)),
    n_missing_PCoA2 = sum(is.na(data$PCoA2)),
    n_missing_ARG_richness = sum(is.na(data$ARG_richness)),
    n_missing_Total_ARG_count = sum(is.na(data$Total_ARG_count)),
    n_missing_Number_of_samples = sum(is.na(data$Number_of_samples)),
    ARG_richness_unique_values = length(unique(data$ARG_richness[!is.na(data$ARG_richness)])),
    Total_ARG_count_unique_values = length(unique(data$Total_ARG_count[!is.na(data$Total_ARG_count)])),
    Number_of_samples_unique_values = length(unique(data$Number_of_samples[!is.na(data$Number_of_samples)])),
    stringsAsFactors = FALSE
  )

  validation
}


run_spearman_test <- function(
    data,
    analysis_label,
    axis_column,
    descriptor_column
) {

  complete_rows <- complete.cases(
    data[
      ,
      c(
        axis_column,
        descriptor_column
      )
    ]
  )

  test_data <- data[
    complete_rows,
    ,
    drop = FALSE
  ]

  n_valid <- nrow(test_data)

  if (n_valid < 3) {

    return(
      data.frame(
        analysis = analysis_label,
        axis = axis_column,
        descriptor = descriptor_column,
        n = n_valid,
        spearman_rho = NA_real_,
        p_value = NA_real_,
        status = "not_tested_too_few_observations",
        stringsAsFactors = FALSE
      )
    )
  }

  if (
    length(unique(test_data[[axis_column]])) < 2 ||
      length(unique(test_data[[descriptor_column]])) < 2
  ) {

    return(
      data.frame(
        analysis = analysis_label,
        axis = axis_column,
        descriptor = descriptor_column,
        n = n_valid,
        spearman_rho = NA_real_,
        p_value = NA_real_,
        status = "not_tested_no_variation",
        stringsAsFactors = FALSE
      )
    )
  }

  test_result <- suppressWarnings(
    stats::cor.test(
      x = test_data[[axis_column]],
      y = test_data[[descriptor_column]],
      method = "spearman",
      exact = FALSE,
      alternative = "two.sided"
    )
  )

  data.frame(
    analysis = analysis_label,
    axis = axis_column,
    descriptor = descriptor_column,
    n = n_valid,
    spearman_rho = unname(test_result$estimate),
    p_value = test_result$p.value,
    status = "tested",
    stringsAsFactors = FALSE
  )
}


descriptor_label <- function(descriptor) {

  labels <- c(
    ARG_richness = "ARG-group richness",
    Total_ARG_count = "Total ARG count",
    Number_of_samples = "Number of samples"
  )

  unname(labels[descriptor])
}


axis_label <- function(axis) {

  labels <- c(
    PCoA1 = "PCoA1",
    PCoA2 = "PCoA2"
  )

  unname(labels[axis])
}


safe_filename <- function(x) {

  gsub(
    pattern = "[^A-Za-z0-9]+",
    replacement = "_",
    x = x
  )
}


make_correlation_plot <- function(
    data,
    result_row,
    output_file
) {

  axis_column <- result_row$axis
  descriptor_column <- result_row$descriptor

  complete_rows <- complete.cases(
    data[
      ,
      c(
        axis_column,
        descriptor_column
      )
    ]
  )

  plot_data <- data[
    complete_rows,
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) < 3) {
    return(invisible(NULL))
  }

  subtitle_text <- paste0(
    "Spearman rho = ",
    round(result_row$spearman_rho, 3),
    ", raw p = ",
    signif(result_row$p_value, 3),
    ", BH-adjusted p = ",
    signif(result_row$p_adjusted_BH, 3)
  )

  if (
    "Environment" %in% colnames(plot_data) &&
      length(unique(plot_data$Environment[!is.na(plot_data$Environment)])) > 1
  ) {

    p <- ggplot(
      plot_data,
      aes(
        x = .data[[descriptor_column]],
        y = .data[[axis_column]],
        shape = Environment
      )
    ) +
      geom_point(size = point_size, alpha = 0.85) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.6,
        linetype = "dashed"
      )

  } else {

    p <- ggplot(
      plot_data,
      aes(
        x = .data[[descriptor_column]],
        y = .data[[axis_column]]
      )
    ) +
      geom_point(size = point_size, alpha = 0.85) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        linewidth = 0.6,
        linetype = "dashed"
      )
  }

  p <- p +
    labs(
      title = paste0(
        gsub("_", " ", result_row$analysis),
        ": ",
        axis_label(axis_column),
        " versus ",
        descriptor_label(descriptor_column)
      ),
      subtitle = subtitle_text,
      x = descriptor_label(descriptor_column),
      y = axis_label(axis_column)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.title = element_blank()
    )

  ggsave(
    filename = output_file,
    plot = p,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )

  invisible(NULL)
}


# ------------------------------------------------------------
# 6. Read and standardise input tables
# ------------------------------------------------------------

combined_raw <- read_csv_base(
  combined_file
)

sewage_only_raw <- read_csv_base(
  sewage_only_file
)

sample_number_raw <- read_csv_base(
  sample_number_file
)

combined_data <- standardise_coordinate_table(
  data = combined_raw,
  data_label = "Combined_ARG_group_PCoA_coordinates_metrics.csv",
  expected_environment = NULL
)

sewage_only_data <- standardise_coordinate_table(
  data = sewage_only_raw,
  data_label = "Sewage_only_ARG_group_PCoA_coordinates_metrics.csv",
  expected_environment = "Sewage"
)

sample_number_data <- standardise_coordinate_table(
  data = sample_number_raw,
  data_label = "sample_number_metadata_used.csv",
  expected_environment = NULL
)

combined_data <- merge_missing_metrics(
  coordinate_data = combined_data,
  metadata_data = sample_number_data,
  data_label = "Combined analysis"
)

sewage_only_data <- merge_missing_metrics(
  coordinate_data = sewage_only_data,
  metadata_data = sample_number_data,
  data_label = "Sewage-only analysis"
)


# ------------------------------------------------------------
# 7. Prepare transit-only data
# ------------------------------------------------------------

if (file.exists(transit_only_file)) {

  transit_only_raw <- read_csv_base(
    transit_only_file
  )

  transit_only_data <- standardise_coordinate_table(
    data = transit_only_raw,
    data_label = "Transit_only_ARG_group_PCoA_coordinates_metrics.csv",
    expected_environment = "Transit"
  )

  transit_only_data <- merge_missing_metrics(
    coordinate_data = transit_only_data,
    metadata_data = sample_number_data,
    data_label = "Transit-only analysis"
  )

} else {

  if (!"Environment" %in% colnames(combined_data)) {
    stop(
      paste0(
        "No transit-only file was found and the combined coordinate table ",
        "does not contain an Environment column."
      )
    )
  }

  transit_only_data <- combined_data[
    combined_data$Environment == "Transit",
    ,
    drop = FALSE
  ]

  if (nrow(transit_only_data) == 0) {
    stop(
      "No transit profiles could be extracted from the combined coordinate table."
    )
  }

  warning(
    paste0(
      "Transit_only_ARG_group_PCoA_coordinates_metrics.csv was not found. ",
      "Transit-only correlations will be calculated using transit profiles ",
      "extracted from the combined PCoA coordinate table."
    )
  )
}


# ------------------------------------------------------------
# 8. Optional cross-check with ALL coordinate table
# ------------------------------------------------------------

if (file.exists(all_coordinates_file)) {

  all_coordinates_raw <- read_csv_base(
    all_coordinates_file
  )

  all_coordinates_data <- standardise_coordinate_table(
    data = all_coordinates_raw,
    data_label = "ALL_ARG_group_PCoA_coordinates_metrics.csv",
    expected_environment = NULL
  )

  write.csv(
    data.frame(
      file = "ALL_ARG_group_PCoA_coordinates_metrics.csv",
      n_profiles = nrow(all_coordinates_data),
      n_sewage = sum(all_coordinates_data$Environment == "Sewage", na.rm = TRUE),
      n_transit = sum(all_coordinates_data$Environment == "Transit", na.rm = TRUE),
      stringsAsFactors = FALSE
    ),
    file = file.path(
      diagnostics_dir,
      "ALL_coordinate_table_crosscheck.csv"
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 9. Validate analysis datasets
# ------------------------------------------------------------

validation_table <- do.call(
  rbind,
  list(
    validate_analysis_data(
      data = combined_data,
      analysis_label = "Combined"
    ),
    validate_analysis_data(
      data = sewage_only_data,
      analysis_label = "Sewage_only"
    ),
    validate_analysis_data(
      data = transit_only_data,
      analysis_label = "Transit_only"
    )
  )
)

write.csv(
  validation_table,
  file = file.path(
    diagnostics_dir,
    "pcoa_correlation_input_validation.csv"
  ),
  row.names = FALSE
)

write.table(
  validation_table,
  file = file.path(
    diagnostics_dir,
    "pcoa_correlation_input_validation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 10. Define planned tests
# ------------------------------------------------------------

planned_tests <- rbind(

  expand.grid(
    analysis = "Combined",
    axis = c("PCoA1", "PCoA2"),
    descriptor = c(
      "ARG_richness",
      "Total_ARG_count",
      "Number_of_samples"
    ),
    stringsAsFactors = FALSE
  ),

  expand.grid(
    analysis = "Sewage_only",
    axis = c("PCoA1", "PCoA2"),
    descriptor = c(
      "ARG_richness",
      "Total_ARG_count",
      "Number_of_samples"
    ),
    stringsAsFactors = FALSE
  ),

  expand.grid(
    analysis = "Transit_only",
    axis = c("PCoA1", "PCoA2"),
    descriptor = c(
      "ARG_richness",
      "Total_ARG_count"
    ),
    stringsAsFactors = FALSE
  )
)


analysis_data_list <- list(
  Combined = combined_data,
  Sewage_only = sewage_only_data,
  Transit_only = transit_only_data
)


# ------------------------------------------------------------
# 11. Run Spearman correlations
# ------------------------------------------------------------

correlation_results_list <- vector(
  mode = "list",
  length = nrow(planned_tests)
)

for (i in seq_len(nrow(planned_tests))) {

  current_analysis <- planned_tests$analysis[i]
  current_axis <- planned_tests$axis[i]
  current_descriptor <- planned_tests$descriptor[i]

  current_data <- analysis_data_list[[current_analysis]]

  correlation_results_list[[i]] <- run_spearman_test(
    data = current_data,
    analysis_label = current_analysis,
    axis_column = current_axis,
    descriptor_column = current_descriptor
  )
}

correlation_results <- do.call(
  rbind,
  correlation_results_list
)

tested_rows <- correlation_results$status == "tested"

correlation_results$p_adjusted_BH <- NA_real_

correlation_results$p_adjusted_BH[tested_rows] <- p.adjust(
  correlation_results$p_value[tested_rows],
  method = "BH"
)

correlation_results$significant_raw_p <- ifelse(
  tested_rows,
  correlation_results$p_value < significance_threshold,
  NA
)

correlation_results$significant_BH <- ifelse(
  tested_rows,
  correlation_results$p_adjusted_BH < significance_threshold,
  NA
)

correlation_results$descriptor_label <- vapply(
  correlation_results$descriptor,
  descriptor_label,
  character(1)
)


# ------------------------------------------------------------
# 12. Save correlation results
# ------------------------------------------------------------

write.csv(
  correlation_results,
  file = file.path(
    tables_dir,
    "pcoa_axis_spearman_correlations_BH_adjusted.csv"
  ),
  row.names = FALSE
)

write.table(
  correlation_results,
  file = file.path(
    tables_dir,
    "pcoa_axis_spearman_correlations_BH_adjusted.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

tested_results <- correlation_results[
  correlation_results$status == "tested",
  ,
  drop = FALSE
]

write.csv(
  tested_results,
  file = file.path(
    tables_dir,
    "pcoa_axis_spearman_correlations_tested_only.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Save thesis-focused PCoA1-only table
# ------------------------------------------------------------

pcoa1_results <- correlation_results[
  correlation_results$axis == "PCoA1",
  ,
  drop = FALSE
]

write.csv(
  pcoa1_results,
  file = file.path(
    tables_dir,
    "pcoa1_spearman_correlations_BH_adjusted.csv"
  ),
  row.names = FALSE
)

write.table(
  pcoa1_results,
  file = file.path(
    tables_dir,
    "pcoa1_spearman_correlations_BH_adjusted.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 14. Generate scatterplots
# ------------------------------------------------------------

if (make_scatterplots) {

  for (i in seq_len(nrow(correlation_results))) {

    result_row <- correlation_results[i, , drop = FALSE]

    if (result_row$status != "tested") {
      next
    }

    current_data <- analysis_data_list[
      [result_row$analysis]
    ]

    output_file <- file.path(
      figures_dir,
      paste0(
        safe_filename(result_row$analysis),
        "_",
        result_row$axis,
        "_vs_",
        safe_filename(result_row$descriptor),
        ".png"
      )
    )

    make_correlation_plot(
      data = current_data,
      result_row = result_row,
      output_file = output_file
    )
  }
}


# ------------------------------------------------------------
# 15. Save concise thesis report
# ------------------------------------------------------------

number_of_valid_tests <- sum(
  correlation_results$status == "tested"
)

number_significant_BH <- sum(
  correlation_results$significant_BH,
  na.rm = TRUE
)

strongest_correlations <- tested_results[
  order(
    abs(tested_results$spearman_rho),
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

strongest_correlations <- head(
  strongest_correlations,
  10
)

report_lines <- c(
  "PCoA AXIS CORRELATION ANALYSIS",
  "",
  paste0("Number of planned correlations: ", nrow(planned_tests)),
  paste0("Number of valid tested correlations: ", number_of_valid_tests),
  paste0(
    "Multiple-testing correction: Benjamini-Hochberg across all valid tests"
  ),
  paste0(
    "Number of BH-significant correlations at FDR < ",
    significance_threshold,
    ": ",
    number_significant_BH
  ),
  "",
  "Strongest correlations by absolute Spearman rho:",
  capture.output(
    print(
      strongest_correlations[
        ,
        c(
          "analysis",
          "axis",
          "descriptor",
          "n",
          "spearman_rho",
          "p_value",
          "p_adjusted_BH"
        )
      ],
      row.names = FALSE
    )
  )
)

writeLines(
  report_lines,
  con = file.path(
    tables_dir,
    "pcoa_axis_correlation_thesis_report.txt"
  )
)


# ------------------------------------------------------------
# 16. Save session information
# ------------------------------------------------------------

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("R version: ", R.version.string),
  "Correlation method: two-sided Spearman rank correlation",
  "P-value adjustment: Benjamini-Hochberg across all valid tests",
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_text,
  con = file.path(
    diagnostics_dir,
    "pcoa_axis_correlation_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 17. Print final summary
# ------------------------------------------------------------

message("")
message("====================================================")
message("PCoA AXIS CORRELATION ANALYSIS COMPLETED")
message("====================================================")

message("")
message("Input validation:")
print(validation_table, row.names = FALSE)

message("")
message("Correlation results:")
print(
  correlation_results[
    ,
    c(
      "analysis",
      "axis",
      "descriptor",
      "n",
      "spearman_rho",
      "p_value",
      "p_adjusted_BH",
      "status"
    )
  ],
  row.names = FALSE
)

message("")
message("Valid tested correlations: ", number_of_valid_tests)
message("BH-significant correlations: ", number_significant_BH)

message("")
message("Outputs saved in:")
message(output_dir)

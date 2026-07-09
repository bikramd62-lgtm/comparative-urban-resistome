# ============================================================
# 01_reduced_depth_city_restricted_permanova.R
#
# Purpose:
#   Reanalyse reduced-depth sewage-versus-transit ARG-group
#   sensitivity tests using city-restricted PERMANOVA.
#
# Comparisons:
#   1. Sewage 1M versus Transit 5M
#   2. Sewage 2M versus Transit 5M
#   3. Sewage 3M versus Transit 5M
#
# Input:
#   Binary ARG-group presence/absence matrices and matching
#   metadata tables for each depth comparison.
#
# Methods:
#   - Binary Jaccard dissimilarity
#   - PCoA using stats::cmdscale()
#   - City-restricted PERMANOVA using vegan::adonis2()
#   - PERMDISP using vegan::betadisper() and vegan::permutest()
#   - Negative-eigenvalue diagnostics for each PCoA
#
# Important:
#   PERMANOVA permutations are restricted within city using:
#     strata = metadata$city
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c("vegan")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
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

library(vegan)


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
  "reduced_depth"
)

output_dir <- file.path(
  project_root,
  "results",
  "sensitivity_analysis",
  "reduced_depth"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

distance_dir <- file.path(
  output_dir,
  "distance_matrices"
)

ordination_dir <- file.path(
  output_dir,
  "ordination"
)

diagnostics_dir <- file.path(
  output_dir,
  "diagnostics"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(distance_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ordination_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_dir)) {
  stop(
    paste0(
      "Input directory not found:\n",
      input_dir,
      "\n\nPlace the reduced-depth input files in this folder."
    )
  )
}


# ------------------------------------------------------------
# 3. Analysis settings
# ------------------------------------------------------------

n_permutations <- 9999
analysis_seed <- 123


# ------------------------------------------------------------
# 4. Input file specification
# ------------------------------------------------------------

comparisons <- list(

  Sewage_1M_vs_Transit_5M = list(
    sewage_depth = "1M",
    transit_depth = "5M",
    binary_matrix = "combined_arg_group_binary_matrix_sewage1M_transit5M_t80_prev10.tsv",
    metadata = "metadata_sewage1M_transit5M_arg_group.tsv"
  ),

  Sewage_2M_vs_Transit_5M = list(
    sewage_depth = "2M",
    transit_depth = "5M",
    binary_matrix = "combined_arg_group_binary_matrix_sewage2M_transit5M_t80_prev10.tsv",
    metadata = "metadata_sewage2M_transit5M_arg_group.tsv"
  ),

  Sewage_3M_vs_Transit_5M = list(
    sewage_depth = "3M",
    transit_depth = "5M",
    binary_matrix = "combined_arg_group_binary_matrix_sewage3M_transit5M_t80_prev10.tsv",
    metadata = "metadata_sewage3M_transit5M_arg_group.tsv"
  )
)


# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

standardise_metadata_names <- function(metadata) {

  original_names <- colnames(metadata)
  lower_names <- tolower(original_names)

  profile_col <- which(lower_names %in% c("profile_id", "sampleid", "sample_id"))
  city_col <- which(lower_names == "city")
  environment_col <- which(lower_names == "environment")

  if (length(profile_col) != 1) {
    stop("Metadata must contain one profile ID column: profile_id, SampleID, or sample_id.")
  }

  if (length(city_col) != 1) {
    stop("Metadata must contain one city column.")
  }

  if (length(environment_col) != 1) {
    stop("Metadata must contain one environment column.")
  }

  colnames(metadata)[profile_col] <- "profile_id"
  colnames(metadata)[city_col] <- "city"
  colnames(metadata)[environment_col] <- "environment"

  metadata
}


read_binary_matrix <- function(file_path, comparison_name) {

  binary_table <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  lower_names <- tolower(colnames(binary_table))
  profile_col <- which(lower_names %in% c("profile_id", "sampleid", "sample_id"))

  if (length(profile_col) != 1) {
    stop(
      paste0(
        comparison_name,
        ": binary matrix must contain one profile ID column named profile_id, SampleID, or sample_id."
      )
    )
  }

  colnames(binary_table)[profile_col] <- "profile_id"

  if (anyDuplicated(binary_table$profile_id) > 0) {
    stop(
      paste0(
        comparison_name,
        ": duplicate profile IDs found in binary matrix."
      )
    )
  }

  profile_ids <- trimws(as.character(binary_table$profile_id))
  binary_table$profile_id <- NULL

  binary_table[] <- lapply(
    binary_table,
    function(x) as.numeric(as.character(x))
  )

  binary_matrix <- as.matrix(binary_table)
  rownames(binary_matrix) <- profile_ids
  storage.mode(binary_matrix) <- "numeric"

  if (anyNA(binary_matrix)) {
    stop(
      paste0(
        comparison_name,
        ": binary matrix contains missing, non-numeric, or invalid values."
      )
    )
  }

  observed_values <- sort(unique(as.vector(binary_matrix)))

  if (!all(observed_values %in% c(0, 1))) {
    stop(
      paste0(
        comparison_name,
        ": binary matrix contains values other than 0 and 1.\n",
        "Observed values: ",
        paste(observed_values, collapse = ", ")
      )
    )
  }

  empty_profiles <- rownames(binary_matrix)[rowSums(binary_matrix) == 0]

  if (length(empty_profiles) > 0) {
    stop(
      paste0(
        comparison_name,
        ": the following profiles contain no ARG groups:\n",
        paste(empty_profiles, collapse = ", ")
      )
    )
  }

  binary_matrix
}


read_and_align_metadata <- function(file_path, binary_matrix, comparison_name) {

  metadata <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  metadata <- standardise_metadata_names(metadata)

  metadata$profile_id <- trimws(as.character(metadata$profile_id))
  metadata$city <- trimws(as.character(metadata$city))
  metadata$environment <- trimws(as.character(metadata$environment))

  metadata$environment[tolower(metadata$environment) == "sewage"] <- "Sewage"
  metadata$environment[tolower(metadata$environment) == "transit"] <- "Transit"

  if (anyDuplicated(metadata$profile_id) > 0) {
    stop(
      paste0(
        comparison_name,
        ": duplicate profile IDs found in metadata."
      )
    )
  }

  matrix_profiles <- rownames(binary_matrix)
  metadata_profiles <- metadata$profile_id

  missing_from_metadata <- setdiff(matrix_profiles, metadata_profiles)
  missing_from_matrix <- setdiff(metadata_profiles, matrix_profiles)

  if (length(missing_from_metadata) > 0 || length(missing_from_matrix) > 0) {
    stop(
      paste0(
        comparison_name,
        ": matrix and metadata profile IDs do not match.\n\n",
        "Missing from metadata:\n",
        paste(missing_from_metadata, collapse = ", "),
        "\n\nMissing from matrix:\n",
        paste(missing_from_matrix, collapse = ", ")
      )
    )
  }

  metadata <- metadata[
    match(matrix_profiles, metadata$profile_id),
    ,
    drop = FALSE
  ]

  if (!identical(metadata$profile_id, matrix_profiles)) {
    stop(
      paste0(
        comparison_name,
        ": metadata could not be aligned to the binary matrix."
      )
    )
  }

  rownames(metadata) <- metadata$profile_id

  metadata$environment <- factor(
    metadata$environment,
    levels = c("Sewage", "Transit")
  )

  metadata$city <- factor(metadata$city)

  if (anyNA(metadata$environment)) {
    stop(
      paste0(
        comparison_name,
        ": environment labels must be Sewage or Transit."
      )
    )
  }

  metadata
}


verify_city_pairing <- function(metadata, comparison_name) {

  pairing_table <- table(metadata$city, metadata$environment)

  if (!all(c("Sewage", "Transit") %in% colnames(pairing_table))) {
    stop(
      paste0(
        comparison_name,
        ": both Sewage and Transit profiles must be present."
      )
    )
  }

  valid_pairing <- all(
    pairing_table[, "Sewage"] == 1 &
      pairing_table[, "Transit"] == 1
  )

  if (!valid_pairing) {
    print(pairing_table)

    stop(
      paste0(
        comparison_name,
        ": each city must contain exactly one sewage and one transit profile."
      )
    )
  }

  pairing_table
}


identify_environment_row <- function(permanova_result, comparison_name) {

  result_rows <- rownames(permanova_result)

  environment_row <- which(tolower(result_rows) == "environment")

  if (length(environment_row) == 0) {
    environment_row <- grep("environment", result_rows, ignore.case = TRUE)
  }

  if (length(environment_row) != 1) {
    print(permanova_result)

    stop(
      paste0(
        comparison_name,
        ": could not uniquely identify the environment row in PERMANOVA output."
      )
    )
  }

  environment_row
}


calculate_eigenvalue_diagnostics <- function(eigenvalues) {

  positive_eigenvalues <- eigenvalues[eigenvalues > 0]
  negative_eigenvalues <- eigenvalues[eigenvalues < 0]

  positive_sum <- sum(positive_eigenvalues)
  negative_sum <- sum(negative_eigenvalues)

  if (length(negative_eigenvalues) == 0) {
    most_negative <- 0
  } else {
    most_negative <- min(negative_eigenvalues)
  }

  data.frame(
    number_of_positive_eigenvalues = length(positive_eigenvalues),
    number_of_negative_eigenvalues = length(negative_eigenvalues),
    largest_positive_eigenvalue = ifelse(length(positive_eigenvalues) > 0, max(positive_eigenvalues), NA_real_),
    most_negative_eigenvalue = most_negative,
    sum_positive_eigenvalues = positive_sum,
    sum_negative_eigenvalues = negative_sum,
    absolute_negative_to_positive_ratio = ifelse(
      positive_sum > 0,
      abs(negative_sum) / positive_sum,
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
}


run_one_comparison <- function(comparison_name, comparison_spec) {

  message("")
  message("====================================================")
  message("Running: ", comparison_name)
  message("====================================================")

  matrix_file <- file.path(input_dir, comparison_spec$binary_matrix)
  metadata_file <- file.path(input_dir, comparison_spec$metadata)

  if (!file.exists(matrix_file)) {
    stop(
      paste0(
        "Binary matrix not found:\n",
        matrix_file
      )
    )
  }

  if (!file.exists(metadata_file)) {
    stop(
      paste0(
        "Metadata file not found:\n",
        metadata_file
      )
    )
  }

  binary_matrix <- read_binary_matrix(
    file_path = matrix_file,
    comparison_name = comparison_name
  )

  metadata <- read_and_align_metadata(
    file_path = metadata_file,
    binary_matrix = binary_matrix,
    comparison_name = comparison_name
  )

  pairing_table <- verify_city_pairing(
    metadata = metadata,
    comparison_name = comparison_name
  )

  n_profiles <- nrow(binary_matrix)
  n_cities <- nlevels(metadata$city)
  n_arg_groups <- ncol(binary_matrix)

  message("Profiles: ", n_profiles)
  message("Cities: ", n_cities)
  message("Retained ARG groups: ", n_arg_groups)

  # ----------------------------------------------------------
  # Binary Jaccard dissimilarity
  # ----------------------------------------------------------

  jaccard_distance <- vegan::vegdist(
    binary_matrix,
    method = "jaccard",
    binary = TRUE
  )

  jaccard_matrix <- as.matrix(jaccard_distance)

  jaccard_output <- data.frame(
    profile_id = rownames(jaccard_matrix),
    jaccard_matrix,
    check.names = FALSE
  )

  write.table(
    jaccard_output,
    file = file.path(
      distance_dir,
      paste0(comparison_name, "_ARG_group_binary_Jaccard_distance_matrix.tsv")
    ),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # PCoA and eigenvalue diagnostics
  # ----------------------------------------------------------

  pcoa_result <- stats::cmdscale(
    jaccard_distance,
    k = 2,
    eig = TRUE
  )

  positive_eigenvalue_sum <- sum(pcoa_result$eig[pcoa_result$eig > 0])

  pcoa1_percent <- (pcoa_result$eig[1] / positive_eigenvalue_sum) * 100
  pcoa2_percent <- (pcoa_result$eig[2] / positive_eigenvalue_sum) * 100

  pcoa_coordinates <- data.frame(
    profile_id = metadata$profile_id,
    city = as.character(metadata$city),
    environment = as.character(metadata$environment),
    PCoA1 = pcoa_result$points[, 1],
    PCoA2 = pcoa_result$points[, 2],
    PCoA1_percent = pcoa1_percent,
    PCoA2_percent = pcoa2_percent,
    stringsAsFactors = FALSE
  )

  write.table(
    pcoa_coordinates,
    file = file.path(
      ordination_dir,
      paste0(comparison_name, "_PCoA_coordinates.tsv")
    ),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  eigenvalue_table <- data.frame(
    axis = seq_along(pcoa_result$eig),
    eigenvalue = pcoa_result$eig,
    sign = ifelse(
      pcoa_result$eig > 0,
      "positive",
      ifelse(pcoa_result$eig < 0, "negative", "zero")
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    eigenvalue_table,
    file = file.path(
      diagnostics_dir,
      paste0(comparison_name, "_PCoA_eigenvalues.tsv")
    ),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  eigenvalue_diagnostics <- calculate_eigenvalue_diagnostics(
    eigenvalues = pcoa_result$eig
  )

  eigenvalue_diagnostics$comparison <- comparison_name
  eigenvalue_diagnostics$sewage_depth <- comparison_spec$sewage_depth
  eigenvalue_diagnostics$transit_depth <- comparison_spec$transit_depth

  # ----------------------------------------------------------
  # City-restricted PERMANOVA
  # ----------------------------------------------------------

  set.seed(analysis_seed)

  permanova_result <- vegan::adonis2(
    jaccard_distance ~ environment,
    data = metadata,
    permutations = n_permutations,
    strata = metadata$city,
    by = "terms"
  )

  environment_row <- identify_environment_row(
    permanova_result = permanova_result,
    comparison_name = comparison_name
  )

  permanova_pseudo_F <- as.numeric(permanova_result[environment_row, "F"])
  permanova_R2 <- as.numeric(permanova_result[environment_row, "R2"])
  permanova_p_value <- as.numeric(permanova_result[environment_row, "Pr(>F)"])

  # ----------------------------------------------------------
  # PERMDISP
  # ----------------------------------------------------------

  dispersion_model <- vegan::betadisper(
    jaccard_distance,
    group = metadata$environment,
    type = "median",
    bias.adjust = FALSE,
    sqrt.dist = FALSE
  )

  set.seed(analysis_seed)

  permdisp_result <- vegan::permutest(
    dispersion_model,
    permutations = n_permutations,
    pairwise = FALSE
  )

  permdisp_table <- as.data.frame(permdisp_result$tab)

  permdisp_F <- as.numeric(permdisp_table[1, "F"])
  permdisp_p_value <- as.numeric(permdisp_table[1, "Pr(>F)"])

  mean_distance_to_group_median <- tapply(
    dispersion_model$distances,
    metadata$environment,
    mean,
    na.rm = TRUE
  )

  median_distance_to_group_median <- tapply(
    dispersion_model$distances,
    metadata$environment,
    median,
    na.rm = TRUE
  )

  dispersion_distances <- data.frame(
    profile_id = metadata$profile_id,
    city = as.character(metadata$city),
    environment = as.character(metadata$environment),
    distance_to_group_median = as.numeric(dispersion_model$distances),
    stringsAsFactors = FALSE
  )

  write.table(
    dispersion_distances,
    file = file.path(
      diagnostics_dir,
      paste0(comparison_name, "_PERMDISP_distances_to_group_median.tsv")
    ),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # Save full statistical output
  # ----------------------------------------------------------

  permanova_full_file <- file.path(
    diagnostics_dir,
    paste0(comparison_name, "_city_restricted_PERMANOVA_9999_full.txt")
  )

  permanova_text <- c(
    paste0("Comparison: ", comparison_name),
    paste0("Binary matrix: ", matrix_file),
    paste0("Metadata file: ", metadata_file),
    paste0("Profiles: ", n_profiles),
    paste0("Cities: ", n_cities),
    paste0("Retained ARG groups: ", n_arg_groups),
    "Distance measure: binary Jaccard dissimilarity",
    paste0("PERMANOVA permutations requested: ", n_permutations),
    "PERMANOVA permutation restriction: within city",
    "PERMANOVA implementation: strata = metadata$city",
    paste0("Random seed: ", analysis_seed),
    "",
    "Matched city structure:",
    capture.output(print(pairing_table)),
    "",
    "City-restricted PERMANOVA result:",
    capture.output(print(permanova_result))
  )

  writeLines(
    permanova_text,
    con = permanova_full_file
  )

  permdisp_full_file <- file.path(
    diagnostics_dir,
    paste0(comparison_name, "_PERMDISP_9999_full.txt")
  )

  permdisp_text <- c(
    paste0("Comparison: ", comparison_name),
    paste0("PERMDISP permutations requested: ", n_permutations),
    "PERMDISP permutation restriction: unrestricted",
    "Group-centre estimator: spatial median",
    paste0("Random seed: ", analysis_seed),
    "",
    "PERMDISP result:",
    capture.output(print(permdisp_result)),
    "",
    "Mean distance to group median:",
    capture.output(print(mean_distance_to_group_median)),
    "",
    "Median distance to group median:",
    capture.output(print(median_distance_to_group_median))
  )

  writeLines(
    permdisp_text,
    con = permdisp_full_file
  )

  # ----------------------------------------------------------
  # Summary row
  # ----------------------------------------------------------

  summary_row <- data.frame(
    comparison = comparison_name,
    sewage_depth = comparison_spec$sewage_depth,
    transit_depth = comparison_spec$transit_depth,
    number_of_profiles = n_profiles,
    number_of_cities = n_cities,
    number_of_ARG_groups = n_arg_groups,
    distance_measure = "Binary Jaccard",
    PCoA1_percent = pcoa1_percent,
    PCoA2_percent = pcoa2_percent,
    PERMANOVA_pseudo_F = permanova_pseudo_F,
    PERMANOVA_R2 = permanova_R2,
    PERMANOVA_p_value = permanova_p_value,
    PERMANOVA_permutations = n_permutations,
    PERMANOVA_restriction = "Restricted within city",
    PERMANOVA_seed = analysis_seed,
    PERMDISP_F = permdisp_F,
    PERMDISP_p_value = permdisp_p_value,
    PERMDISP_permutations = n_permutations,
    PERMDISP_restriction = "Unrestricted",
    mean_distance_sewage = as.numeric(mean_distance_to_group_median["Sewage"]),
    mean_distance_transit = as.numeric(mean_distance_to_group_median["Transit"]),
    median_distance_sewage = as.numeric(median_distance_to_group_median["Sewage"]),
    median_distance_transit = as.numeric(median_distance_to_group_median["Transit"]),
    stringsAsFactors = FALSE
  )

  message("")
  message("City-restricted PERMANOVA:")
  print(permanova_result)

  message("")
  message("PERMDISP:")
  print(permdisp_result)

  list(
    summary = summary_row,
    eigenvalue_diagnostics = eigenvalue_diagnostics
  )
}


# ------------------------------------------------------------
# 6. Run all comparisons
# ------------------------------------------------------------

all_results <- lapply(
  names(comparisons),
  function(comparison_name) {
    run_one_comparison(
      comparison_name = comparison_name,
      comparison_spec = comparisons[[comparison_name]]
    )
  }
)

names(all_results) <- names(comparisons)

summary_table <- do.call(
  rbind,
  lapply(all_results, function(x) x$summary)
)

eigenvalue_summary <- do.call(
  rbind,
  lapply(all_results, function(x) x$eigenvalue_diagnostics)
)

rownames(summary_table) <- NULL
rownames(eigenvalue_summary) <- NULL


# ------------------------------------------------------------
# 7. Save combined outputs
# ------------------------------------------------------------

write.table(
  summary_table,
  file = file.path(
    tables_dir,
    "reduced_depth_city_restricted_PERMANOVA_PERMDISP_complete_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.csv(
  summary_table,
  file = file.path(
    tables_dir,
    "reduced_depth_city_restricted_PERMANOVA_PERMDISP_complete_summary.csv"
  ),
  row.names = FALSE
)

thesis_table <- summary_table[
  ,
  c(
    "comparison",
    "number_of_ARG_groups",
    "PCoA1_percent",
    "PCoA2_percent",
    "PERMANOVA_pseudo_F",
    "PERMANOVA_R2",
    "PERMANOVA_p_value",
    "PERMDISP_F",
    "PERMDISP_p_value"
  )
]

thesis_table$PCoA1_percent <- round(thesis_table$PCoA1_percent, 2)
thesis_table$PCoA2_percent <- round(thesis_table$PCoA2_percent, 2)
thesis_table$PERMANOVA_pseudo_F <- round(thesis_table$PERMANOVA_pseudo_F, 4)
thesis_table$PERMANOVA_R2 <- round(thesis_table$PERMANOVA_R2, 4)
thesis_table$PERMANOVA_p_value <- round(thesis_table$PERMANOVA_p_value, 4)
thesis_table$PERMDISP_F <- round(thesis_table$PERMDISP_F, 4)
thesis_table$PERMDISP_p_value <- round(thesis_table$PERMDISP_p_value, 4)

write.table(
  thesis_table,
  file = file.path(
    tables_dir,
    "reduced_depth_thesis_reporting_table.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  eigenvalue_summary,
  file = file.path(
    diagnostics_dir,
    "reduced_depth_PCoA_negative_eigenvalue_diagnostics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.csv(
  eigenvalue_summary,
  file = file.path(
    diagnostics_dir,
    "reduced_depth_PCoA_negative_eigenvalue_diagnostics.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Save session information
# ------------------------------------------------------------

session_file <- file.path(
  diagnostics_dir,
  "reduced_depth_analysis_sessionInfo.txt"
)

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("R version: ", R.version.string),
  paste0("vegan version: ", as.character(packageVersion("vegan"))),
  paste0("PERMANOVA permutations requested: ", n_permutations),
  "PERMANOVA restriction: within city",
  "PERMANOVA implementation: strata = metadata$city",
  paste0("PERMDISP permutations requested: ", n_permutations),
  "PERMDISP restriction: unrestricted",
  paste0("Random seed: ", analysis_seed),
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_text,
  con = session_file
)


# ------------------------------------------------------------
# 9. Print final summary
# ------------------------------------------------------------

message("")
message("====================================================")
message("REDUCED-DEPTH SENSITIVITY ANALYSIS COMPLETED")
message("====================================================")

message("")
message("Thesis reporting table:")
print(thesis_table, row.names = FALSE)

message("")
message("Negative eigenvalue diagnostics:")
print(eigenvalue_summary, row.names = FALSE)

message("")
message("Outputs saved in:")
message(output_dir)

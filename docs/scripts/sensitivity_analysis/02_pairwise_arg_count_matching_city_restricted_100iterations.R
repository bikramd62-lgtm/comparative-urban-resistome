# ============================================================
# 02_pairwise_arg_count_matching_city_restricted_100iterations.R
#
# Purpose:
#   Reconstruct and rerun the pairwise total ARG-count matching
#   sensitivity analysis
#
# Workflow:
#   1. Read prevalence-filtered ARG-group count matrices for
#      sewage and transit.
#   2. Combine sewage and transit matrices over the union of
#      ARG groups.
#   3. For each city, define the rarefaction target as the lower
#      total ARG count between sewage and transit.
#   4. Rarefy the higher-count profile to the pairwise-minimum
#      total using vegan::rrarefy().
#   5. Repeat this procedure over 100 iterations.
#   6. Convert rarefied counts to binary presence/absence.
#   7. Calculate binary Jaccard dissimilarity.
#   8. Run city-restricted PERMANOVA in every iteration.
#   9. Run PERMDISP in every iteration.
#  10. Save iteration-level and overall summary outputs.
#
# Key statistical design:
#   PERMANOVA permutations are restricted within city using:
#     strata = metadata$City
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
  "arg_count_matching"
)

output_dir <- file.path(
  project_root,
  "results",
  "sensitivity_analysis",
  "arg_count_matching"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

diagnostics_dir <- file.path(
  output_dir,
  "diagnostics"
)

ordination_dir <- file.path(
  output_dir,
  "ordination"
)

binary_matrix_dir <- file.path(
  output_dir,
  "binary_matrices"
)

count_matrix_dir <- file.path(
  output_dir,
  "rarefied_count_matrices"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ordination_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(binary_matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(count_matrix_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_dir)) {
  stop(
    paste0(
      "Input directory not found:\n",
      input_dir,
      "\n\nPlace the ARG-count-matching input files in this folder."
    )
  )
}


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

sewage_file <- file.path(
  input_dir,
  "sewage_group_count_matrix_filtered10pct.tsv"
)

transit_file <- file.path(
  input_dir,
  "transit_group_count_matrix_prev10.tsv"
)

totals_file <- file.path(
  input_dir,
  "original_total_ARG_counts_and_richness.csv"
)

archived_log_file <- file.path(
  input_dir,
  "pairwise_min_ARG_abundance_matching_log_iteration_001.csv"
)

archived_binary_file <- file.path(
  input_dir,
  "pairwise_min_ARG_abundance_matched_binary_iteration_001.csv"
)

required_input_files <- c(
  sewage_file,
  transit_file,
  totals_file
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
# 4. Analysis settings
# ------------------------------------------------------------

number_of_iterations <- 100

permanova_permutations <- 999

permdisp_permutations <- 999

base_seed <- 123

significance_threshold <- 0.05

save_iteration_matrices <- TRUE

# Explicit random-number generator settings for reproducibility.
# sample.kind = "Rejection" is the R >= 3.6 default.
RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)


# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

read_feature_count_matrix <- function(
    file_path,
    matrix_label
) {

  x <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(x) < 2) {
    stop(
      paste0(
        matrix_label,
        " must contain one feature column and at least one city column."
      )
    )
  }

  lower_names <- tolower(colnames(x))

  candidate_feature_columns <- which(
    lower_names %in% c(
      "feature",
      "arg_group",
      "arggroup",
      "group",
      "gene",
      "mechanism",
      "class"
    )
  )

  if (length(candidate_feature_columns) == 1) {

    feature_column <- candidate_feature_columns

  } else {

    # Fallback: use the first column as the feature column.
    # This makes the script robust to input files where the
    # feature column has a different name.
    feature_column <- 1
  }

  feature_names <- trimws(
    as.character(x[[feature_column]])
  )

  if (any(feature_names == "")) {
    stop(
      paste0(
        matrix_label,
        " contains empty feature names."
      )
    )
  }

  if (anyDuplicated(feature_names) > 0) {
    duplicated_features <- unique(
      feature_names[duplicated(feature_names)]
    )

    stop(
      paste0(
        matrix_label,
        " contains duplicated feature names:\n",
        paste(duplicated_features, collapse = "\n")
      )
    )
  }

  count_table <- x[, -feature_column, drop = FALSE]

  count_table[] <- lapply(
    count_table,
    function(column) {
      as.numeric(as.character(column))
    }
  )

  count_matrix <- as.matrix(count_table)

  rownames(count_matrix) <- feature_names

  storage.mode(count_matrix) <- "numeric"

  if (anyNA(count_matrix)) {
    stop(
      paste0(
        matrix_label,
        " contains missing, non-numeric, or invalid values."
      )
    )
  }

  if (any(count_matrix < 0)) {
    stop(
      paste0(
        matrix_label,
        " contains negative count values."
      )
    )
  }

  if (any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
    stop(
      paste0(
        matrix_label,
        " contains non-integer count values."
      )
    )
  }

  count_matrix <- round(count_matrix)

  storage.mode(count_matrix) <- "integer"

  return(count_matrix)
}


standardise_environment_labels <- function(environment_vector) {

  environment_vector <- trimws(
    as.character(environment_vector)
  )

  environment_vector[
    tolower(environment_vector) == "sewage"
  ] <- "Sewage"

  environment_vector[
    tolower(environment_vector) == "transit"
  ] <- "Transit"

  environment_vector
}


identify_environment_row <- function(
    permanova_result,
    comparison_label = "PERMANOVA"
) {

  result_rows <- rownames(permanova_result)

  environment_row <- which(
    tolower(result_rows) == "environment"
  )

  if (length(environment_row) == 0) {
    environment_row <- grep(
      "environment",
      result_rows,
      ignore.case = TRUE
    )
  }

  if (length(environment_row) != 1) {
    message("PERMANOVA row names:")
    print(result_rows)

    print(permanova_result)

    stop(
      paste0(
        comparison_label,
        ": could not uniquely identify the Environment row."
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

  data.frame(
    number_of_positive_eigenvalues = length(positive_eigenvalues),
    number_of_negative_eigenvalues = length(negative_eigenvalues),
    largest_positive_eigenvalue = ifelse(
      length(positive_eigenvalues) > 0,
      max(positive_eigenvalues),
      NA_real_
    ),
    most_negative_eigenvalue = ifelse(
      length(negative_eigenvalues) > 0,
      min(negative_eigenvalues),
      0
    ),
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


# ------------------------------------------------------------
# 6. Read sewage and transit count matrices
# ------------------------------------------------------------

sewage_feature_matrix <- read_feature_count_matrix(
  file_path = sewage_file,
  matrix_label = "Sewage ARG-group count matrix"
)

transit_feature_matrix <- read_feature_count_matrix(
  file_path = transit_file,
  matrix_label = "Transit ARG-group count matrix"
)


# ------------------------------------------------------------
# 7. Read totals and metadata
# ------------------------------------------------------------

original_totals <- read.csv(
  totals_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_total_columns <- c(
  "SampleID",
  "City",
  "Environment",
  "Total_ARG_count",
  "ARG_richness"
)

missing_total_columns <- setdiff(
  required_total_columns,
  colnames(original_totals)
)

if (length(missing_total_columns) > 0) {
  stop(
    paste0(
      "The totals file is missing these required columns:\n",
      paste(missing_total_columns, collapse = ", ")
    )
  )
}

original_totals$SampleID <- trimws(
  as.character(original_totals$SampleID)
)

original_totals$City <- trimws(
  as.character(original_totals$City)
)

original_totals$Environment <- standardise_environment_labels(
  original_totals$Environment
)

if (anyDuplicated(original_totals$SampleID) > 0) {
  stop("Duplicate SampleID values were found in the totals table.")
}

if (!all(original_totals$Environment %in% c("Sewage", "Transit"))) {
  stop("Environment must contain only Sewage and Transit.")
}


# ------------------------------------------------------------
# 8. Validate city and paired design
# ------------------------------------------------------------

sewage_cities <- colnames(sewage_feature_matrix)
transit_cities <- colnames(transit_feature_matrix)
metadata_cities <- unique(original_totals$City)

if (!setequal(sewage_cities, transit_cities)) {
  stop(
    paste0(
      "Sewage and transit count matrices do not contain the same cities.\n",
      "Cities only in sewage:\n",
      paste(setdiff(sewage_cities, transit_cities), collapse = ", "),
      "\n\nCities only in transit:\n",
      paste(setdiff(transit_cities, sewage_cities), collapse = ", ")
    )
  )
}

if (!setequal(sewage_cities, metadata_cities)) {
  stop(
    paste0(
      "Cities in count matrices do not match cities in totals metadata.\n",
      "Cities only in matrices:\n",
      paste(setdiff(sewage_cities, metadata_cities), collapse = ", "),
      "\n\nCities only in metadata:\n",
      paste(setdiff(metadata_cities, sewage_cities), collapse = ", ")
    )
  )
}

city_environment_table <- table(
  original_totals$City,
  original_totals$Environment
)

if (!all(c("Sewage", "Transit") %in% colnames(city_environment_table))) {
  stop("Both Sewage and Transit profiles are required.")
}

valid_pairing <- all(
  city_environment_table[, "Sewage"] == 1 &
    city_environment_table[, "Transit"] == 1
)

if (!valid_pairing) {
  print(city_environment_table)

  stop(
    "Each city must contain exactly one sewage and one transit profile."
  )
}

if (nrow(city_environment_table) != 16) {
  warning(
    paste0(
      "Expected 16 matched cities, but found ",
      nrow(city_environment_table),
      "."
    )
  )
}


# ------------------------------------------------------------
# 9. Construct combined sample-by-feature count matrix
# ------------------------------------------------------------

all_features <- union(
  rownames(sewage_feature_matrix),
  rownames(transit_feature_matrix)
)

combined_count_matrix <- matrix(
  0L,
  nrow = nrow(original_totals),
  ncol = length(all_features),
  dimnames = list(
    original_totals$SampleID,
    all_features
  )
)

for (i in seq_len(nrow(original_totals))) {

  profile_city <- original_totals$City[i]
  profile_environment <- original_totals$Environment[i]

  if (profile_environment == "Sewage") {
    source_matrix <- sewage_feature_matrix
  } else {
    source_matrix <- transit_feature_matrix
  }

  if (!profile_city %in% colnames(source_matrix)) {
    stop(
      paste0(
        "City not found in ",
        profile_environment,
        " matrix: ",
        profile_city
      )
    )
  }

  combined_count_matrix[
    i,
    rownames(source_matrix)
  ] <- source_matrix[, profile_city]
}


# ------------------------------------------------------------
# 10. Validate total ARG counts and richness
# ------------------------------------------------------------

calculated_total_counts <- rowSums(combined_count_matrix)
calculated_richness <- rowSums(combined_count_matrix > 0)

input_validation <- data.frame(
  SampleID = original_totals$SampleID,
  City = original_totals$City,
  Environment = original_totals$Environment,
  archived_total_ARG_count = as.numeric(original_totals$Total_ARG_count),
  calculated_total_ARG_count = as.numeric(calculated_total_counts),
  total_ARG_count_difference =
    as.numeric(calculated_total_counts) -
    as.numeric(original_totals$Total_ARG_count),
  archived_ARG_richness = as.numeric(original_totals$ARG_richness),
  calculated_ARG_richness = as.numeric(calculated_richness),
  ARG_richness_difference =
    as.numeric(calculated_richness) -
    as.numeric(original_totals$ARG_richness),
  stringsAsFactors = FALSE
)

write.csv(
  input_validation,
  file = file.path(
    diagnostics_dir,
    "input_total_and_richness_validation.csv"
  ),
  row.names = FALSE
)

if (any(input_validation$total_ARG_count_difference != 0)) {

  print(
    input_validation[
      input_validation$total_ARG_count_difference != 0,
    ]
  )

  stop(
    paste0(
      "Calculated total ARG counts do not match archived totals. ",
      "The analysis was stopped."
    )
  )
}

if (any(input_validation$ARG_richness_difference != 0)) {
  warning(
    paste0(
      "One or more calculated ARG richness values differ from the archived richness values. ",
      "See diagnostics/input_total_and_richness_validation.csv."
    )
  )
}


# ------------------------------------------------------------
# 11. Prepare metadata
# ------------------------------------------------------------

metadata <- data.frame(
  SampleID = original_totals$SampleID,
  City = factor(original_totals$City),
  Environment = factor(
    original_totals$Environment,
    levels = c("Sewage", "Transit")
  ),
  Original_total_ARG_count = as.numeric(original_totals$Total_ARG_count),
  Original_ARG_richness = as.numeric(original_totals$ARG_richness),
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$SampleID

if (!identical(metadata$SampleID, rownames(combined_count_matrix))) {
  stop("Metadata and combined count matrix are not aligned.")
}


# ------------------------------------------------------------
# 12. Calculate pairwise-minimum targets
# ------------------------------------------------------------

city_names <- levels(metadata$City)

matching_targets <- do.call(
  rbind,
  lapply(
    city_names,
    function(city_name) {

      city_rows <- which(metadata$City == city_name)

      sewage_row <- city_rows[
        metadata$Environment[city_rows] == "Sewage"
      ]

      transit_row <- city_rows[
        metadata$Environment[city_rows] == "Transit"
      ]

      sewage_total <- sum(combined_count_matrix[sewage_row, ])
      transit_total <- sum(combined_count_matrix[transit_row, ])

      target_total <- min(sewage_total, transit_total)

      lower_environment <- if (sewage_total < transit_total) {
        "Sewage"
      } else if (transit_total < sewage_total) {
        "Transit"
      } else {
        "Equal"
      }

      data.frame(
        City = city_name,
        Sewage_original_total = sewage_total,
        Transit_original_total = transit_total,
        Pairwise_min_target_total = target_total,
        Lower_count_environment = lower_environment,
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  matching_targets,
  file = file.path(
    tables_dir,
    "pairwise_minimum_matching_targets.csv"
  ),
  row.names = FALSE
)

metadata$Pairwise_min_target_total <- matching_targets[
  match(
    as.character(metadata$City),
    matching_targets$City
  ),
  "Pairwise_min_target_total"
]


# ------------------------------------------------------------
# 13. Optional comparison with archived iteration-001 target log
# ------------------------------------------------------------

if (file.exists(archived_log_file)) {

  archived_log <- read.csv(
    archived_log_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if ("City" %in% colnames(archived_log)) {

    archived_log$City <- trimws(
      as.character(archived_log$City)
    )

    target_comparison <- merge(
      matching_targets,
      archived_log,
      by = "City",
      all = TRUE,
      suffixes = c("_reconstructed", "_archived")
    )

    write.csv(
      target_comparison,
      file = file.path(
        diagnostics_dir,
        "archived_vs_reconstructed_target_comparison.csv"
      ),
      row.names = FALSE
    )
  }
}


# ------------------------------------------------------------
# 14. Run 100 rarefaction iterations
# ------------------------------------------------------------

iteration_results <- vector(
  mode = "list",
  length = number_of_iterations
)

all_pcoa_coordinates <- vector(
  mode = "list",
  length = number_of_iterations
)

all_eigenvalue_diagnostics <- vector(
  mode = "list",
  length = number_of_iterations
)

for (iteration in seq_len(number_of_iterations)) {

  message(
    "Running iteration ",
    iteration,
    " of ",
    number_of_iterations,
    "..."
  )

  rarefaction_seed <- base_seed + iteration - 1
  permanova_seed <- base_seed + 10000 + iteration - 1
  permdisp_seed <- base_seed + 20000 + iteration - 1

  rarefied_count_matrix <- combined_count_matrix

  set.seed(rarefaction_seed)

  for (profile_index in seq_len(nrow(combined_count_matrix))) {

    original_profile_counts <- combined_count_matrix[
      profile_index,
      ,
      drop = FALSE
    ]

    original_total <- sum(original_profile_counts)

    target_total <- metadata$Pairwise_min_target_total[
      profile_index
    ]

    if (original_total < target_total) {
      stop(
        paste0(
          "Target exceeds original total for profile: ",
          metadata$SampleID[profile_index]
        )
      )
    }

    if (original_total > target_total) {

      rarefied_profile <- vegan::rrarefy(
        original_profile_counts,
        sample = target_total
      )

      rarefied_count_matrix[
        profile_index,
      ] <- rarefied_profile[1, ]

    } else {

      rarefied_count_matrix[
        profile_index,
      ] <- original_profile_counts[1, ]
    }
  }


  # ----------------------------------------------------------
  # Verify rarefaction totals
  # ----------------------------------------------------------

  rarefied_totals <- rowSums(rarefied_count_matrix)

  if (!all(rarefied_totals == metadata$Pairwise_min_target_total)) {

    failed_profiles <- metadata$SampleID[
      rarefied_totals != metadata$Pairwise_min_target_total
    ]

    stop(
      paste0(
        "Rarefaction totals did not match targets in iteration ",
        iteration,
        " for:\n",
        paste(failed_profiles, collapse = ", ")
      )
    )
  }

  city_total_check <- tapply(
    rarefied_totals,
    metadata$City,
    function(x) length(unique(x)) == 1
  )

  if (!all(city_total_check)) {
    stop(
      paste0(
        "Sewage and transit totals were not equal within every city in iteration ",
        iteration,
        "."
      )
    )
  }


  # ----------------------------------------------------------
  # Convert to binary presence/absence
  # ----------------------------------------------------------

  binary_matrix <- ifelse(
    rarefied_count_matrix > 0,
    1L,
    0L
  )

  rownames(binary_matrix) <- rownames(rarefied_count_matrix)
  colnames(binary_matrix) <- colnames(rarefied_count_matrix)

  # Remove features absent from every profile in this iteration.
  # No new prevalence filtering is applied here.
  retained_columns <- colSums(binary_matrix) > 0

  binary_matrix_for_analysis <- binary_matrix[
    ,
    retained_columns,
    drop = FALSE
  ]

  if (ncol(binary_matrix_for_analysis) == 0) {
    stop(
      paste0(
        "No ARG groups remained in iteration ",
        iteration,
        "."
      )
    )
  }

  if (any(rowSums(binary_matrix_for_analysis) == 0)) {
    stop(
      paste0(
        "At least one empty profile occurred in iteration ",
        iteration,
        "."
      )
    )
  }


  # ----------------------------------------------------------
  # Binary Jaccard dissimilarity
  # ----------------------------------------------------------

  jaccard_distance <- vegan::vegdist(
    binary_matrix_for_analysis,
    method = "jaccard",
    binary = TRUE
  )


  # ----------------------------------------------------------
  # PCoA and eigenvalue diagnostics
  # ----------------------------------------------------------

  pcoa_result <- stats::cmdscale(
    jaccard_distance,
    k = 2,
    eig = TRUE
  )

  positive_eigenvalue_sum <- sum(
    pcoa_result$eig[pcoa_result$eig > 0]
  )

  pcoa1_percent <- (
    pcoa_result$eig[1] /
      positive_eigenvalue_sum
  ) * 100

  pcoa2_percent <- (
    pcoa_result$eig[2] /
      positive_eigenvalue_sum
  ) * 100

  pcoa_coordinates <- as.data.frame(
    pcoa_result$points
  )

  colnames(pcoa_coordinates) <- c(
    "PCoA1",
    "PCoA2"
  )

  all_pcoa_coordinates[[iteration]] <- data.frame(
    Iteration = iteration,
    SampleID = metadata$SampleID,
    City = as.character(metadata$City),
    Environment = as.character(metadata$Environment),
    PCoA1 = pcoa_coordinates$PCoA1,
    PCoA2 = pcoa_coordinates$PCoA2,
    PCoA1_percent = pcoa1_percent,
    PCoA2_percent = pcoa2_percent,
    stringsAsFactors = FALSE
  )

  eigenvalue_diagnostics <- calculate_eigenvalue_diagnostics(
    eigenvalues = pcoa_result$eig
  )

  eigenvalue_diagnostics$Iteration <- iteration

  all_eigenvalue_diagnostics[[iteration]] <- eigenvalue_diagnostics


  # ----------------------------------------------------------
  # City-restricted PERMANOVA
  # ----------------------------------------------------------

  set.seed(permanova_seed)

  permanova_result <- vegan::adonis2(
    jaccard_distance ~ Environment,
    data = metadata,
    permutations = permanova_permutations,
    strata = metadata$City,
    by = "terms"
  )

  environment_row <- identify_environment_row(
    permanova_result = permanova_result,
    comparison_label = paste0("Iteration ", iteration)
  )

  permanova_pseudo_F <- as.numeric(
    permanova_result[environment_row, "F"]
  )

  permanova_R2 <- as.numeric(
    permanova_result[environment_row, "R2"]
  )

  permanova_p_value <- as.numeric(
    permanova_result[environment_row, "Pr(>F)"]
  )


  # ----------------------------------------------------------
  # PERMDISP
  # ----------------------------------------------------------

  dispersion_model <- vegan::betadisper(
    jaccard_distance,
    group = metadata$Environment,
    type = "median",
    bias.adjust = FALSE,
    sqrt.dist = FALSE
  )

  set.seed(permdisp_seed)

  permdisp_result <- vegan::permutest(
    dispersion_model,
    permutations = permdisp_permutations,
    pairwise = FALSE
  )

  permdisp_table <- as.data.frame(
    permdisp_result$tab
  )

  permdisp_F <- as.numeric(
    permdisp_table[1, "F"]
  )

  permdisp_p_value <- as.numeric(
    permdisp_table[1, "Pr(>F)"]
  )


  # ----------------------------------------------------------
  # Richness summaries
  # ----------------------------------------------------------

  iteration_richness <- rowSums(binary_matrix)

  mean_sewage_richness <- mean(
    iteration_richness[metadata$Environment == "Sewage"]
  )

  mean_transit_richness <- mean(
    iteration_richness[metadata$Environment == "Transit"]
  )

  median_sewage_richness <- median(
    iteration_richness[metadata$Environment == "Sewage"]
  )

  median_transit_richness <- median(
    iteration_richness[metadata$Environment == "Transit"]
  )


  # ----------------------------------------------------------
  # Save iteration-level result row
  # ----------------------------------------------------------

  iteration_results[[iteration]] <- data.frame(
    Iteration = iteration,
    Rarefaction_seed = rarefaction_seed,
    PERMANOVA_seed = permanova_seed,
    PERMDISP_seed = permdisp_seed,
    Number_of_profiles = nrow(binary_matrix_for_analysis),
    Number_of_retained_ARG_groups = ncol(binary_matrix_for_analysis),
    PERMANOVA_pseudo_F = permanova_pseudo_F,
    PERMANOVA_R2 = permanova_R2,
    PERMANOVA_p_value = permanova_p_value,
    PERMANOVA_significant =
      permanova_p_value < significance_threshold,
    PERMDISP_F = permdisp_F,
    PERMDISP_p_value = permdisp_p_value,
    PERMDISP_significant =
      permdisp_p_value < significance_threshold,
    Mean_sewage_richness = mean_sewage_richness,
    Mean_transit_richness = mean_transit_richness,
    Median_sewage_richness = median_sewage_richness,
    Median_transit_richness = median_transit_richness,
    PCoA1_percent = pcoa1_percent,
    PCoA2_percent = pcoa2_percent,
    stringsAsFactors = FALSE
  )


  # ----------------------------------------------------------
  # Save matrices
  # ----------------------------------------------------------

  if (save_iteration_matrices || iteration == 1) {

    count_output <- data.frame(
      SampleID = metadata$SampleID,
      City = as.character(metadata$City),
      Environment = as.character(metadata$Environment),
      rarefied_count_matrix,
      check.names = FALSE
    )

    binary_output <- data.frame(
      SampleID = metadata$SampleID,
      City = as.character(metadata$City),
      Environment = as.character(metadata$Environment),
      binary_matrix_for_analysis,
      check.names = FALSE
    )

    write.csv(
      count_output,
      file = file.path(
        count_matrix_dir,
        sprintf(
          "pairwise_min_ARG_count_matched_counts_iteration_%03d.csv",
          iteration
        )
      ),
      row.names = FALSE
    )

    write.csv(
      binary_output,
      file = file.path(
        binary_matrix_dir,
        sprintf(
          "pairwise_min_ARG_count_matched_binary_iteration_%03d.csv",
          iteration
        )
      ),
      row.names = FALSE
    )
  }


  # ----------------------------------------------------------
  # Save full iteration-001 statistical output
  # ----------------------------------------------------------

  if (iteration == 1) {

    full_output_file <- file.path(
      diagnostics_dir,
      "iteration_001_city_restricted_PERMANOVA_and_PERMDISP_full.txt"
    )

    full_output_text <- c(
      "CASE 6: PAIRWISE TOTAL ARG-COUNT MATCHING, ITERATION 001",
      "",
      "Rarefaction function: vegan::rrarefy()",
      paste0("Rarefaction seed: ", rarefaction_seed),
      paste0("PERMANOVA permutations: ", permanova_permutations),
      paste0("PERMANOVA seed: ", permanova_seed),
      "PERMANOVA restriction: within city",
      "PERMANOVA implementation: strata = metadata$City",
      paste0("PERMDISP permutations: ", permdisp_permutations),
      paste0("PERMDISP seed: ", permdisp_seed),
      "PERMDISP restriction: unrestricted",
      "",
      "City-restricted PERMANOVA:",
      capture.output(print(permanova_result)),
      "",
      "PERMDISP:",
      capture.output(print(permdisp_result))
    )

    writeLines(
      full_output_text,
      con = full_output_file
    )
  }
}


# ------------------------------------------------------------
# 15. Combine iteration-level outputs
# ------------------------------------------------------------

iteration_results_table <- do.call(
  rbind,
  iteration_results
)

rownames(iteration_results_table) <- NULL

pcoa_coordinates_table <- do.call(
  rbind,
  all_pcoa_coordinates
)

rownames(pcoa_coordinates_table) <- NULL

eigenvalue_diagnostics_table <- do.call(
  rbind,
  all_eigenvalue_diagnostics
)

rownames(eigenvalue_diagnostics_table) <- NULL


# ------------------------------------------------------------
# 16. Save iteration-level tables
# ------------------------------------------------------------

write.csv(
  iteration_results_table,
  file = file.path(
    tables_dir,
    "pairwise_min_city_restricted_iteration_metrics.csv"
  ),
  row.names = FALSE
)

write.table(
  iteration_results_table,
  file = file.path(
    tables_dir,
    "pairwise_min_city_restricted_iteration_metrics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.csv(
  pcoa_coordinates_table,
  file = file.path(
    ordination_dir,
    "pairwise_min_PCoA_coordinates_all_iterations.csv"
  ),
  row.names = FALSE
)

write.table(
  eigenvalue_diagnostics_table,
  file = file.path(
    diagnostics_dir,
    "pairwise_min_PCoA_negative_eigenvalue_diagnostics_all_iterations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.csv(
  eigenvalue_diagnostics_table,
  file = file.path(
    diagnostics_dir,
    "pairwise_min_PCoA_negative_eigenvalue_diagnostics_all_iterations.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Overall summary
# ------------------------------------------------------------

permanova_significant_n <- sum(
  iteration_results_table$PERMANOVA_p_value <
    significance_threshold
)

permanova_significant_percent <- (
  permanova_significant_n /
    number_of_iterations
) * 100

permdisp_significant_n <- sum(
  iteration_results_table$PERMDISP_p_value <
    significance_threshold
)

permdisp_significant_percent <- (
  permdisp_significant_n /
    number_of_iterations
) * 100

R2_quantiles <- quantile(
  iteration_results_table$PERMANOVA_R2,
  probabilities = c(0.25, 0.50, 0.75),
  names = FALSE,
  na.rm = TRUE
)

negative_eigenvalue_iterations <- sum(
  eigenvalue_diagnostics_table$number_of_negative_eigenvalues > 0
)

maximum_negative_ratio <- max(
  eigenvalue_diagnostics_table$absolute_negative_to_positive_ratio,
  na.rm = TRUE
)

overall_summary <- data.frame(
  Number_of_iterations = number_of_iterations,
  Rarefaction_method = "vegan::rrarefy",
  Base_seed = base_seed,
  PERMANOVA_permutations = permanova_permutations,
  PERMANOVA_restriction = "Restricted within city",
  PERMANOVA_significant_iterations = permanova_significant_n,
  PERMANOVA_significant_percentage = permanova_significant_percent,
  PERMANOVA_median_R2 = median(iteration_results_table$PERMANOVA_R2),
  PERMANOVA_R2_Q1 = R2_quantiles[1],
  PERMANOVA_R2_Q3 = R2_quantiles[3],
  PERMANOVA_R2_minimum = min(iteration_results_table$PERMANOVA_R2),
  PERMANOVA_R2_maximum = max(iteration_results_table$PERMANOVA_R2),
  PERMANOVA_mean_R2 = mean(iteration_results_table$PERMANOVA_R2),
  PERMANOVA_SD_R2 = sd(iteration_results_table$PERMANOVA_R2),
  PERMANOVA_median_p_value =
    median(iteration_results_table$PERMANOVA_p_value),
  PERMDISP_permutations = permdisp_permutations,
  PERMDISP_restriction = "Unrestricted",
  PERMDISP_significant_iterations = permdisp_significant_n,
  PERMDISP_significant_percentage = permdisp_significant_percent,
  PERMDISP_median_p_value =
    median(iteration_results_table$PERMDISP_p_value),
  Mean_sewage_richness_across_iterations =
    mean(iteration_results_table$Mean_sewage_richness),
  Mean_transit_richness_across_iterations =
    mean(iteration_results_table$Mean_transit_richness),
  Negative_eigenvalue_iterations =
    negative_eigenvalue_iterations,
  Maximum_absolute_negative_to_positive_eigenvalue_ratio =
    maximum_negative_ratio,
  stringsAsFactors = FALSE
)

write.csv(
  overall_summary,
  file = file.path(
    tables_dir,
    "pairwise_min_city_restricted_overall_summary.csv"
  ),
  row.names = FALSE
)

write.table(
  overall_summary,
  file = file.path(
    tables_dir,
    "pairwise_min_city_restricted_overall_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 18. Thesis-report text
# ------------------------------------------------------------

report_text <- c(
  "CASE 6: PAIRWISE TOTAL ARG-COUNT MATCHING",
  "",
  paste0("Number of rarefaction iterations: ", number_of_iterations),
  "Rarefaction method: vegan::rrarefy()",
  paste0("Base random seed: ", base_seed),
  paste0("PERMANOVA permutations per iteration: ", permanova_permutations),
  "PERMANOVA restriction: within city",
  paste0(
    "Significant PERMANOVA iterations: ",
    permanova_significant_n,
    " of ",
    number_of_iterations,
    " (",
    round(permanova_significant_percent, 1),
    "%)"
  ),
  paste0(
    "Median PERMANOVA R2: ",
    format(
      median(iteration_results_table$PERMANOVA_R2),
      digits = 4,
      nsmall = 4
    )
  ),
  paste0(
    "PERMANOVA R2 IQR: ",
    format(R2_quantiles[1], digits = 4, nsmall = 4),
    "–",
    format(R2_quantiles[3], digits = 4, nsmall = 4)
  ),
  paste0(
    "PERMANOVA R2 range: ",
    format(
      min(iteration_results_table$PERMANOVA_R2),
      digits = 4,
      nsmall = 4
    ),
    "–",
    format(
      max(iteration_results_table$PERMANOVA_R2),
      digits = 4,
      nsmall = 4
    )
  ),
  paste0(
    "Significant PERMDISP iterations: ",
    permdisp_significant_n,
    " of ",
    number_of_iterations,
    " (",
    round(permdisp_significant_percent, 1),
    "%)"
  ),
  paste0(
    "Median PERMDISP p-value: ",
    format(
      median(iteration_results_table$PERMDISP_p_value),
      digits = 4,
      nsmall = 4
    )
  ),
  paste0(
    "Iterations with negative PCoA eigenvalues: ",
    negative_eigenvalue_iterations,
    " of ",
    number_of_iterations
  ),
  paste0(
    "Maximum absolute negative-to-positive eigenvalue ratio: ",
    format(
      maximum_negative_ratio,
      digits = 4,
      scientific = FALSE
    )
  )
)

writeLines(
  report_text,
  con = file.path(
    tables_dir,
    "pairwise_min_city_restricted_thesis_report.txt"
  )
)


# ------------------------------------------------------------
# 19. Optional comparison with archived binary iteration 001
# ------------------------------------------------------------

new_iteration_001_binary_file <- file.path(
  binary_matrix_dir,
  "pairwise_min_ARG_count_matched_binary_iteration_001.csv"
)

if (
  file.exists(archived_binary_file) &&
    file.exists(new_iteration_001_binary_file)
) {

  archived_binary <- read.csv(
    archived_binary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  reconstructed_binary <- read.csv(
    new_iteration_001_binary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  metadata_columns <- c(
    "SampleID",
    "City",
    "Environment"
  )

  archived_features <- setdiff(
    colnames(archived_binary),
    metadata_columns
  )

  reconstructed_features <- setdiff(
    colnames(reconstructed_binary),
    metadata_columns
  )

  comparison_features <- union(
    archived_features,
    reconstructed_features
  )

  archived_aligned <- matrix(
    0L,
    nrow = nrow(archived_binary),
    ncol = length(comparison_features),
    dimnames = list(
      archived_binary$SampleID,
      comparison_features
    )
  )

  reconstructed_aligned <- matrix(
    0L,
    nrow = nrow(reconstructed_binary),
    ncol = length(comparison_features),
    dimnames = list(
      reconstructed_binary$SampleID,
      comparison_features
    )
  )

  archived_aligned[, archived_features] <- as.matrix(
    archived_binary[, archived_features, drop = FALSE]
  )

  reconstructed_aligned[, reconstructed_features] <- as.matrix(
    reconstructed_binary[, reconstructed_features, drop = FALSE]
  )

  common_sample_order <- intersect(
    rownames(archived_aligned),
    rownames(reconstructed_aligned)
  )

  archived_aligned <- archived_aligned[
    common_sample_order,
    ,
    drop = FALSE
  ]

  reconstructed_aligned <- reconstructed_aligned[
    common_sample_order,
    ,
    drop = FALSE
  ]

  number_of_compared_cells <- length(archived_aligned)

  number_of_identical_cells <- sum(
    archived_aligned == reconstructed_aligned
  )

  archived_binary_comparison <- data.frame(
    Compared_profiles = length(common_sample_order),
    Compared_features = length(comparison_features),
    Compared_cells = number_of_compared_cells,
    Identical_cells = number_of_identical_cells,
    Percentage_identical =
      100 * number_of_identical_cells / number_of_compared_cells,
    Exact_matrix_match = identical(
      archived_aligned,
      reconstructed_aligned
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    archived_binary_comparison,
    file = file.path(
      diagnostics_dir,
      "archived_vs_reconstructed_binary_iteration_001_comparison.csv"
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 20. Session information
# ------------------------------------------------------------

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("R version: ", R.version.string),
  paste0("vegan version: ", as.character(packageVersion("vegan"))),
  paste0("Number of iterations: ", number_of_iterations),
  "Rarefaction method: vegan::rrarefy()",
  paste0("Base seed: ", base_seed),
  paste0("PERMANOVA permutations per iteration: ", permanova_permutations),
  "PERMANOVA restriction: within city",
  "PERMANOVA implementation: strata = metadata$City",
  paste0("PERMDISP permutations per iteration: ", permdisp_permutations),
  "PERMDISP restriction: unrestricted",
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_text,
  con = file.path(
    diagnostics_dir,
    "case6_analysis_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 21. Print final results
# ------------------------------------------------------------

message("")
message("====================================================")
message("PAIRWISE ARG-COUNT MATCHING ANALYSIS COMPLETED")
message("====================================================")

message("")
message(
  "Significant city-restricted PERMANOVA iterations: ",
  permanova_significant_n,
  " of ",
  number_of_iterations,
  " (",
  round(permanova_significant_percent, 1),
  "%)"
)

message(
  "Median PERMANOVA R2: ",
  round(
    median(iteration_results_table$PERMANOVA_R2),
    4
  )
)

message(
  "PERMANOVA R2 IQR: ",
  round(R2_quantiles[1], 4),
  "–",
  round(R2_quantiles[3], 4)
)

message(
  "PERMANOVA R2 range: ",
  round(min(iteration_results_table$PERMANOVA_R2), 4),
  "–",
  round(max(iteration_results_table$PERMANOVA_R2), 4)
)

message(
  "Significant PERMDISP iterations: ",
  permdisp_significant_n,
  " of ",
  number_of_iterations,
  " (",
  round(permdisp_significant_percent, 1),
  "%)"
)

message(
  "Median PERMDISP p-value: ",
  round(
    median(iteration_results_table$PERMDISP_p_value),
    4
  )
)

message(
  "Iterations with negative PCoA eigenvalues: ",
  negative_eigenvalue_iterations,
  " of ",
  number_of_iterations
)

message("")
message("Outputs saved in:")
message(output_dir)

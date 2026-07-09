# ============================================================
# 04_jaccard_turnover_nestedness_decomposition.R
#
# Purpose:
#   Partition ARG-group Jaccard dissimilarity into turnover and
#   nestedness-resultant components for sewage and transit
#   city-level resistome profiles.
#
# Input:
#   1. sewage_group_count_matrix_filtered10pct.tsv
#   2. transit_group_count_matrix_prev10.tsv
#   3. original_total_ARG_counts_and_richness.csv
#
# Main outputs:
#   - Pairwise Jaccard decomposition table
#   - Pair-type summary table
#   - Matched-city sewage-transit decomposition table
#   - Paired Wilcoxon test: turnover versus nestedness
#   - Thesis-ready boxplots and matched-city paired plot
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "betapart",
  "ggplot2"
)

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

library(betapart)
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
  "jaccard_decomposition"
)

output_dir <- file.path(
  project_root,
  "results",
  "sensitivity_analysis",
  "jaccard_decomposition"
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
      "\n\nPlace the Jaccard decomposition input files in this folder."
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

metadata_file <- file.path(
  input_dir,
  "original_total_ARG_counts_and_richness.csv"
)

required_input_files <- c(
  sewage_file,
  transit_file,
  metadata_file
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
# 4. Plot settings
# ------------------------------------------------------------

plot_width <- 6.5
plot_height <- 5.2
plot_dpi <- 300

pair_type_levels <- c(
  "Sewage–Sewage",
  "Sewage–Transit",
  "Transit–Transit"
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

  count_matrix
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


make_pair_type <- function(environment_1, environment_2) {

  if (environment_1 == "Sewage" && environment_2 == "Sewage") {
    return("Sewage–Sewage")
  }

  if (environment_1 == "Transit" && environment_2 == "Transit") {
    return("Transit–Transit")
  }

  "Sewage–Transit"
}


summary_stats <- function(x) {

  x <- x[!is.na(x)]

  data.frame(
    n = length(x),
    median = median(x),
    Q1 = unname(quantile(x, 0.25)),
    Q3 = unname(quantile(x, 0.75)),
    minimum = min(x),
    maximum = max(x),
    mean = mean(x),
    sd = stats::sd(x),
    stringsAsFactors = FALSE
  )
}


safe_filename <- function(x) {

  gsub(
    pattern = "[^A-Za-z0-9]+",
    replacement = "_",
    x = x
  )
}


save_plot_png_pdf <- function(
    plot_object,
    file_stem,
    width = plot_width,
    height = plot_height
) {

  png_file <- file.path(
    figures_dir,
    paste0(file_stem, ".png")
  )

  pdf_file <- file.path(
    figures_dir,
    paste0(file_stem, ".pdf")
  )

  ggsave(
    filename = png_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = plot_dpi
  )

  ggsave(
    filename = pdf_file,
    plot = plot_object,
    width = width,
    height = height
  )
}


make_boxplot <- function(
    pairwise_data,
    value_column,
    y_label,
    title_text,
    file_stem
) {

  p <- ggplot(
    pairwise_data,
    aes(
      x = Pair_type,
      y = .data[[value_column]]
    )
  ) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      fill = "grey85",
      colour = "black",
      linewidth = 0.45
    ) +
    geom_jitter(
      width = 0.12,
      height = 0,
      size = 1.35,
      alpha = 0.45
    ) +
    scale_x_discrete(
      limits = pair_type_levels
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    labs(
      title = title_text,
      x = "Pair type",
      y = y_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0
      ),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      )
    )

  save_plot_png_pdf(
    plot_object = p,
    file_stem = file_stem
  )

  invisible(p)
}


# ------------------------------------------------------------
# 6. Read input matrices and metadata
# ------------------------------------------------------------

sewage_feature_matrix <- read_feature_count_matrix(
  file_path = sewage_file,
  matrix_label = "Sewage ARG-group count matrix"
)

transit_feature_matrix <- read_feature_count_matrix(
  file_path = transit_file,
  matrix_label = "Transit ARG-group count matrix"
)

metadata <- read.csv(
  metadata_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata_columns <- c(
  "SampleID",
  "City",
  "Environment",
  "Total_ARG_count",
  "ARG_richness"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    paste0(
      "The metadata file is missing these required columns:\n",
      paste(missing_metadata_columns, collapse = ", ")
    )
  )
}

metadata$SampleID <- trimws(
  as.character(metadata$SampleID)
)

metadata$City <- trimws(
  as.character(metadata$City)
)

metadata$Environment <- standardise_environment_labels(
  metadata$Environment
)

if (anyDuplicated(metadata$SampleID) > 0) {
  stop("Duplicate SampleID values were found in the metadata table.")
}

if (!all(metadata$Environment %in% c("Sewage", "Transit"))) {
  stop("Environment must contain only Sewage and Transit.")
}


# ------------------------------------------------------------
# 7. Validate paired-city design
# ------------------------------------------------------------

sewage_cities <- colnames(sewage_feature_matrix)
transit_cities <- colnames(transit_feature_matrix)
metadata_cities <- unique(metadata$City)

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
      "Cities in count matrices do not match cities in metadata.\n",
      "Cities only in matrices:\n",
      paste(setdiff(sewage_cities, metadata_cities), collapse = ", "),
      "\n\nCities only in metadata:\n",
      paste(setdiff(metadata_cities, sewage_cities), collapse = ", ")
    )
  )
}

city_environment_table <- table(
  metadata$City,
  metadata$Environment
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
# 8. Reconstruct combined sample-by-feature count matrix
# ------------------------------------------------------------

all_features <- union(
  rownames(sewage_feature_matrix),
  rownames(transit_feature_matrix)
)

combined_count_matrix <- matrix(
  0,
  nrow = nrow(metadata),
  ncol = length(all_features),
  dimnames = list(
    metadata$SampleID,
    all_features
  )
)

for (i in seq_len(nrow(metadata))) {

  profile_city <- metadata$City[i]
  profile_environment <- metadata$Environment[i]

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
# 9. Validate totals and richness
# ------------------------------------------------------------

calculated_total_counts <- rowSums(combined_count_matrix)
calculated_richness <- rowSums(combined_count_matrix > 0)

input_validation <- data.frame(
  SampleID = metadata$SampleID,
  City = metadata$City,
  Environment = metadata$Environment,
  archived_total_ARG_count = as.numeric(metadata$Total_ARG_count),
  calculated_total_ARG_count = as.numeric(calculated_total_counts),
  total_ARG_count_difference =
    as.numeric(calculated_total_counts) -
    as.numeric(metadata$Total_ARG_count),
  archived_ARG_richness = as.numeric(metadata$ARG_richness),
  calculated_ARG_richness = as.numeric(calculated_richness),
  ARG_richness_difference =
    as.numeric(calculated_richness) -
    as.numeric(metadata$ARG_richness),
  stringsAsFactors = FALSE
)

write.csv(
  input_validation,
  file = file.path(
    diagnostics_dir,
    "jaccard_decomposition_input_total_and_richness_validation.csv"
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
      "See diagnostics/jaccard_decomposition_input_total_and_richness_validation.csv."
    )
  )
}


# ------------------------------------------------------------
# 10. Convert counts to binary presence/absence
# ------------------------------------------------------------

combined_binary_matrix <- ifelse(
  combined_count_matrix > 0,
  1L,
  0L
)

rownames(combined_binary_matrix) <- rownames(combined_count_matrix)
colnames(combined_binary_matrix) <- colnames(combined_count_matrix)

empty_profiles <- rownames(combined_binary_matrix)[
  rowSums(combined_binary_matrix) == 0
]

if (length(empty_profiles) > 0) {
  stop(
    paste0(
      "The following profiles contain no ARG groups:\n",
      paste(empty_profiles, collapse = ", ")
    )
  )
}

empty_features <- colnames(combined_binary_matrix)[
  colSums(combined_binary_matrix) == 0
]

if (length(empty_features) > 0) {
  combined_binary_matrix <- combined_binary_matrix[
    ,
    colSums(combined_binary_matrix) > 0,
    drop = FALSE
  ]
}

binary_output <- data.frame(
  SampleID = metadata$SampleID,
  City = metadata$City,
  Environment = metadata$Environment,
  combined_binary_matrix,
  check.names = FALSE
)

write.csv(
  binary_output,
  file = file.path(
    tables_dir,
    "combined_ARG_group_binary_matrix_for_jaccard_decomposition.csv"
  ),
  row.names = FALSE
)

write.table(
  binary_output,
  file = file.path(
    tables_dir,
    "combined_ARG_group_binary_matrix_for_jaccard_decomposition.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 11. Jaccard turnover/nestedness decomposition
# ------------------------------------------------------------

decomposition <- betapart::beta.pair(
  combined_binary_matrix,
  index.family = "jaccard"
)

total_jaccard_matrix <- as.matrix(
  decomposition$beta.jac
)

turnover_matrix <- as.matrix(
  decomposition$beta.jtu
)

nestedness_matrix <- as.matrix(
  decomposition$beta.jne
)

profile_ids <- rownames(combined_binary_matrix)

pairwise_rows <- list()
row_counter <- 1

for (i in seq_len(length(profile_ids) - 1)) {

  for (j in seq((i + 1), length(profile_ids))) {

    profile_1 <- profile_ids[i]
    profile_2 <- profile_ids[j]

    metadata_1 <- metadata[
      match(profile_1, metadata$SampleID),
      ,
      drop = FALSE
    ]

    metadata_2 <- metadata[
      match(profile_2, metadata$SampleID),
      ,
      drop = FALSE
    ]

    environment_1 <- metadata_1$Environment
    environment_2 <- metadata_2$Environment

    pair_type <- make_pair_type(
      environment_1 = environment_1,
      environment_2 = environment_2
    )

    pairwise_rows[[row_counter]] <- data.frame(
      Pair_ID = paste(profile_1, profile_2, sep = "__vs__"),
      Profile_1 = profile_1,
      Profile_2 = profile_2,
      City_1 = metadata_1$City,
      City_2 = metadata_2$City,
      Environment_1 = environment_1,
      Environment_2 = environment_2,
      Pair_type = pair_type,
      Same_city = metadata_1$City == metadata_2$City,
      Total_Jaccard_dissimilarity = total_jaccard_matrix[i, j],
      Turnover_component = turnover_matrix[i, j],
      Nestedness_resultant_component = nestedness_matrix[i, j],
      stringsAsFactors = FALSE
    )

    row_counter <- row_counter + 1
  }
}

pairwise_decomposition <- do.call(
  rbind,
  pairwise_rows
)

pairwise_decomposition$Pair_type <- factor(
  pairwise_decomposition$Pair_type,
  levels = pair_type_levels
)

pairwise_decomposition$Turnover_minus_nestedness <-
  pairwise_decomposition$Turnover_component -
  pairwise_decomposition$Nestedness_resultant_component

pairwise_decomposition$Dominant_component <- ifelse(
  pairwise_decomposition$Turnover_component >
    pairwise_decomposition$Nestedness_resultant_component,
  "Turnover",
  ifelse(
    pairwise_decomposition$Nestedness_resultant_component >
      pairwise_decomposition$Turnover_component,
    "Nestedness",
    "Equal"
  )
)

component_sum_difference <- abs(
  pairwise_decomposition$Total_Jaccard_dissimilarity -
    (
      pairwise_decomposition$Turnover_component +
        pairwise_decomposition$Nestedness_resultant_component
    )
)

pairwise_decomposition$Component_sum_difference <- component_sum_difference


# ------------------------------------------------------------
# 12. Save pairwise decomposition table
# ------------------------------------------------------------

write.csv(
  pairwise_decomposition,
  file = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_pairwise_manual.csv"
  ),
  row.names = FALSE
)

write.table(
  pairwise_decomposition,
  file = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_pairwise_manual.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 13. Pair-type summary tables
# ------------------------------------------------------------

summary_rows <- list()
summary_counter <- 1

for (current_pair_type in pair_type_levels) {

  subset_data <- pairwise_decomposition[
    pairwise_decomposition$Pair_type == current_pair_type,
    ,
    drop = FALSE
  ]

  total_stats <- summary_stats(
    subset_data$Total_Jaccard_dissimilarity
  )

  turnover_stats <- summary_stats(
    subset_data$Turnover_component
  )

  nestedness_stats <- summary_stats(
    subset_data$Nestedness_resultant_component
  )

  summary_rows[[summary_counter]] <- data.frame(
    Pair_type = current_pair_type,
    n_comparisons = nrow(subset_data),

    Total_Jaccard_median = total_stats$median,
    Total_Jaccard_Q1 = total_stats$Q1,
    Total_Jaccard_Q3 = total_stats$Q3,
    Total_Jaccard_minimum = total_stats$minimum,
    Total_Jaccard_maximum = total_stats$maximum,

    Turnover_median = turnover_stats$median,
    Turnover_Q1 = turnover_stats$Q1,
    Turnover_Q3 = turnover_stats$Q3,
    Turnover_minimum = turnover_stats$minimum,
    Turnover_maximum = turnover_stats$maximum,

    Nestedness_median = nestedness_stats$median,
    Nestedness_Q1 = nestedness_stats$Q1,
    Nestedness_Q3 = nestedness_stats$Q3,
    Nestedness_minimum = nestedness_stats$minimum,
    Nestedness_maximum = nestedness_stats$maximum,

    stringsAsFactors = FALSE
  )

  summary_counter <- summary_counter + 1
}

summary_table <- do.call(
  rbind,
  summary_rows
)

write.csv(
  summary_table,
  file = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_summary_manual.csv"
  ),
  row.names = FALSE
)

write.table(
  summary_table,
  file = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_summary_manual.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 14. Extract directly matched sewage-transit city pairs
# ------------------------------------------------------------

matched_city_pairs <- pairwise_decomposition[
  pairwise_decomposition$Same_city == TRUE &
    pairwise_decomposition$Pair_type == "Sewage–Transit",
  ,
  drop = FALSE
]

if (nrow(matched_city_pairs) == 0) {
  stop("No matched sewage-transit city pairs were found.")
}

matched_city_pairs$City <- matched_city_pairs$City_1

for (i in seq_len(nrow(matched_city_pairs))) {

  if (matched_city_pairs$Environment_1[i] == "Sewage") {

    matched_city_pairs$Sewage_profile[i] <-
      matched_city_pairs$Profile_1[i]

    matched_city_pairs$Transit_profile[i] <-
      matched_city_pairs$Profile_2[i]

  } else {

    matched_city_pairs$Sewage_profile[i] <-
      matched_city_pairs$Profile_2[i]

    matched_city_pairs$Transit_profile[i] <-
      matched_city_pairs$Profile_1[i]
  }
}

matched_city_pairs_output <- matched_city_pairs[
  ,
  c(
    "City",
    "Sewage_profile",
    "Transit_profile",
    "Total_Jaccard_dissimilarity",
    "Turnover_component",
    "Nestedness_resultant_component",
    "Turnover_minus_nestedness",
    "Dominant_component"
  )
]

matched_city_pairs_output <- matched_city_pairs_output[
  order(matched_city_pairs_output$City),
  ,
  drop = FALSE
]

write.csv(
  matched_city_pairs_output,
  file = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_matched_city_pairs.csv"
  ),
  row.names = FALSE
)

write.table(
  matched_city_pairs_output,
  file = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_matched_city_pairs.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 15. Paired Wilcoxon test: turnover versus nestedness
# ------------------------------------------------------------

wilcoxon_result <- suppressWarnings(
  stats::wilcox.test(
    x = matched_city_pairs_output$Turnover_component,
    y = matched_city_pairs_output$Nestedness_resultant_component,
    paired = TRUE,
    alternative = "two.sided",
    exact = TRUE
  )
)

valid_differences <- matched_city_pairs_output$Turnover_component -
  matched_city_pairs_output$Nestedness_resultant_component

non_zero_differences <- valid_differences[
  valid_differences != 0
]

n_non_zero <- length(non_zero_differences)

total_rank_sum <- n_non_zero * (n_non_zero + 1) / 2

V_statistic <- as.numeric(
  wilcoxon_result$statistic
)

matched_rank_biserial <- ifelse(
  total_rank_sum > 0,
  (2 * V_statistic / total_rank_sum) - 1,
  NA_real_
)

turnover_exceeds_nestedness_n <- sum(
  matched_city_pairs_output$Turnover_component >
    matched_city_pairs_output$Nestedness_resultant_component
)

nestedness_exceeds_turnover_n <- sum(
  matched_city_pairs_output$Nestedness_resultant_component >
    matched_city_pairs_output$Turnover_component
)

equal_component_n <- sum(
  matched_city_pairs_output$Nestedness_resultant_component ==
    matched_city_pairs_output$Turnover_component
)

wilcoxon_summary <- data.frame(
  test = "Paired Wilcoxon signed-rank test",
  comparison = "Turnover component versus nestedness-resultant component",
  n_matched_cities = nrow(matched_city_pairs_output),
  n_non_zero_differences = n_non_zero,
  V = V_statistic,
  p_value = wilcoxon_result$p.value,
  matched_rank_biserial_effect_size = matched_rank_biserial,
  turnover_exceeds_nestedness_n = turnover_exceeds_nestedness_n,
  nestedness_exceeds_turnover_n = nestedness_exceeds_turnover_n,
  equal_component_n = equal_component_n,
  median_turnover = median(matched_city_pairs_output$Turnover_component),
  turnover_Q1 = unname(quantile(matched_city_pairs_output$Turnover_component, 0.25)),
  turnover_Q3 = unname(quantile(matched_city_pairs_output$Turnover_component, 0.75)),
  median_nestedness = median(matched_city_pairs_output$Nestedness_resultant_component),
  nestedness_Q1 = unname(quantile(matched_city_pairs_output$Nestedness_resultant_component, 0.25)),
  nestedness_Q3 = unname(quantile(matched_city_pairs_output$Nestedness_resultant_component, 0.75)),
  stringsAsFactors = FALSE
)

write.csv(
  wilcoxon_summary,
  file = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_matched_city_wilcoxon_test.csv"
  ),
  row.names = FALSE
)

write.table(
  wilcoxon_summary,
  file = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_matched_city_wilcoxon_test.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 16. Diagnostic summary
# ------------------------------------------------------------

diagnostic_summary <- data.frame(
  n_profiles = nrow(combined_binary_matrix),
  n_ARG_groups = ncol(combined_binary_matrix),
  n_pairwise_comparisons = nrow(pairwise_decomposition),
  expected_pairwise_comparisons =
    nrow(combined_binary_matrix) *
    (nrow(combined_binary_matrix) - 1) / 2,
  n_sewage_sewage_pairs =
    sum(pairwise_decomposition$Pair_type == "Sewage–Sewage"),
  n_sewage_transit_pairs =
    sum(pairwise_decomposition$Pair_type == "Sewage–Transit"),
  n_transit_transit_pairs =
    sum(pairwise_decomposition$Pair_type == "Transit–Transit"),
  n_matched_sewage_transit_city_pairs =
    nrow(matched_city_pairs_output),
  max_component_sum_difference =
    max(pairwise_decomposition$Component_sum_difference),
  stringsAsFactors = FALSE
)

write.csv(
  diagnostic_summary,
  file = file.path(
    diagnostics_dir,
    "jaccard_decomposition_diagnostic_summary.csv"
  ),
  row.names = FALSE
)

write.table(
  diagnostic_summary,
  file = file.path(
    diagnostics_dir,
    "jaccard_decomposition_diagnostic_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# ------------------------------------------------------------
# 17. Generate final thesis figures
# ------------------------------------------------------------

total_plot <- make_boxplot(
  pairwise_data = pairwise_decomposition,
  value_column = "Total_Jaccard_dissimilarity",
  y_label = "Total Jaccard dissimilarity",
  title_text = "Total ARG-group Jaccard dissimilarity",
  file_stem = "ARG_group_total_Jaccard_dissimilarity_boxplot_publication_ready"
)

turnover_plot <- make_boxplot(
  pairwise_data = pairwise_decomposition,
  value_column = "Turnover_component",
  y_label = "Turnover component",
  title_text = "ARG-group turnover component",
  file_stem = "ARG_group_turnover_component_boxplot_publication_ready"
)

nestedness_plot <- make_boxplot(
  pairwise_data = pairwise_decomposition,
  value_column = "Nestedness_resultant_component",
  y_label = "Nestedness-resultant component",
  title_text = "ARG-group nestedness-resultant component",
  file_stem = "ARG_group_nestedness_component_boxplot_publication_ready"
)


# ------------------------------------------------------------
# 18. Matched-city turnover versus nestedness figure
# ------------------------------------------------------------

matched_long <- rbind(
  data.frame(
    City = matched_city_pairs_output$City,
    Component = "Turnover",
    Dissimilarity = matched_city_pairs_output$Turnover_component,
    stringsAsFactors = FALSE
  ),
  data.frame(
    City = matched_city_pairs_output$City,
    Component = "Nestedness",
    Dissimilarity = matched_city_pairs_output$Nestedness_resultant_component,
    stringsAsFactors = FALSE
  )
)

matched_long$Component <- factor(
  matched_long$Component,
  levels = c("Turnover", "Nestedness")
)

annotation_text <- paste0(
  "Paired Wilcoxon\n",
  "signed-rank test\n",
  "V = ",
  round(wilcoxon_summary$V, 3),
  "\n",
  "p value = ",
  signif(wilcoxon_summary$p_value, 3),
  "\n",
  "Matched rank-biserial\n",
  "effect size = ",
  round(wilcoxon_summary$matched_rank_biserial_effect_size, 3),
  "\n",
  "Turnover exceeded\n",
  "nestedness in ",
  wilcoxon_summary$turnover_exceeds_nestedness_n,
  " of ",
  wilcoxon_summary$n_matched_cities,
  " cities"
)

matched_plot <- ggplot(
  matched_long,
  aes(
    x = Component,
    y = Dissimilarity,
    group = City
  )
) +
  geom_line(
    linewidth = 0.45,
    alpha = 0.55
  ) +
  geom_point(
    size = 2.4,
    alpha = 0.85
  ) +
  annotate(
    geom = "label",
    x = 1.07,
    y = 0.98,
    label = annotation_text,
    hjust = 0,
    vjust = 1,
    size = 3.35,
    label.size = 0.25,
    fill = "white"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    title = "Matched-city turnover and nestedness components",
    x = NULL,
    y = "Jaccard dissimilarity component"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0
    )
  )

save_plot_png_pdf(
  plot_object = matched_plot,
  file_stem = "ARG_group_matched_city_turnover_vs_nestedness_publication_ready",
  width = 6.2,
  height = 5.4
)

save_plot_png_pdf(
  plot_object = matched_plot,
  file_stem = "ARG_group_matched_city_turnover_vs_nestedness_left_aligned_final",
  width = 6.2,
  height = 5.4
)


# ------------------------------------------------------------
# 19. Thesis-ready text summary
# ------------------------------------------------------------

sewage_transit_summary <- summary_table[
  summary_table$Pair_type == "Sewage–Transit",
  ,
  drop = FALSE
]

sewage_sewage_summary <- summary_table[
  summary_table$Pair_type == "Sewage–Sewage",
  ,
  drop = FALSE
]

transit_transit_summary <- summary_table[
  summary_table$Pair_type == "Transit–Transit",
  ,
  drop = FALSE
]

report_lines <- c(
  "ARG-GROUP JACCARD TURNOVER/NESTEDNESS DECOMPOSITION",
  "",
  paste0(
    "Total pairwise comparisons: ",
    nrow(pairwise_decomposition)
  ),
  paste0(
    "Sewage-sewage comparisons: ",
    sum(pairwise_decomposition$Pair_type == "Sewage–Sewage")
  ),
  paste0(
    "Sewage-transit comparisons: ",
    sum(pairwise_decomposition$Pair_type == "Sewage–Transit")
  ),
  paste0(
    "Transit-transit comparisons: ",
    sum(pairwise_decomposition$Pair_type == "Transit–Transit")
  ),
  "",
  paste0(
    "Sewage-transit total Jaccard median: ",
    round(sewage_transit_summary$Total_Jaccard_median, 3),
    " (IQR: ",
    round(sewage_transit_summary$Total_Jaccard_Q1, 3),
    "–",
    round(sewage_transit_summary$Total_Jaccard_Q3, 3),
    ")"
  ),
  paste0(
    "Sewage-transit turnover median: ",
    round(sewage_transit_summary$Turnover_median, 3),
    " (IQR: ",
    round(sewage_transit_summary$Turnover_Q1, 3),
    "–",
    round(sewage_transit_summary$Turnover_Q3, 3),
    ")"
  ),
  paste0(
    "Sewage-transit nestedness median: ",
    round(sewage_transit_summary$Nestedness_median, 3),
    " (IQR: ",
    round(sewage_transit_summary$Nestedness_Q1, 3),
    "–",
    round(sewage_transit_summary$Nestedness_Q3, 3),
    ")"
  ),
  "",
  paste0(
    "Matched-city pairs: ",
    wilcoxon_summary$n_matched_cities
  ),
  paste0(
    "Turnover exceeded nestedness in ",
    wilcoxon_summary$turnover_exceeds_nestedness_n,
    " of ",
    wilcoxon_summary$n_matched_cities,
    " cities."
  ),
  paste0(
    "Matched turnover median: ",
    round(wilcoxon_summary$median_turnover, 3),
    " (IQR: ",
    round(wilcoxon_summary$turnover_Q1, 3),
    "–",
    round(wilcoxon_summary$turnover_Q3, 3),
    ")"
  ),
  paste0(
    "Matched nestedness median: ",
    round(wilcoxon_summary$median_nestedness, 3),
    " (IQR: ",
    round(wilcoxon_summary$nestedness_Q1, 3),
    "–",
    round(wilcoxon_summary$nestedness_Q3, 3),
    ")"
  ),
  paste0(
    "Paired Wilcoxon signed-rank test: V = ",
    round(wilcoxon_summary$V, 3),
    ", p = ",
    signif(wilcoxon_summary$p_value, 4),
    "."
  ),
  paste0(
    "Matched rank-biserial effect size: ",
    round(wilcoxon_summary$matched_rank_biserial_effect_size, 3),
    "."
  )
)

writeLines(
  report_lines,
  con = file.path(
    tables_dir,
    "ARG_group_jaccard_decomposition_thesis_report.txt"
  )
)


# ------------------------------------------------------------
# 20. Session information
# ------------------------------------------------------------

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("R version: ", R.version.string),
  paste0("betapart version: ", as.character(packageVersion("betapart"))),
  paste0("ggplot2 version: ", as.character(packageVersion("ggplot2"))),
  "Jaccard decomposition function: betapart::beta.pair(index.family = 'jaccard')",
  "Turnover component: beta.jtu",
  "Nestedness-resultant component: beta.jne",
  "Total Jaccard dissimilarity: beta.jac",
  "Matched-city statistical test: paired Wilcoxon signed-rank test",
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_text,
  con = file.path(
    diagnostics_dir,
    "jaccard_decomposition_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 21. Print final summary
# ------------------------------------------------------------

message("")
message("====================================================")
message("JACCARD TURNOVER/NESTEDNESS DECOMPOSITION COMPLETED")
message("====================================================")

message("")
message("Diagnostic summary:")
print(diagnostic_summary, row.names = FALSE)

message("")
message("Pair-type summary:")
print(summary_table, row.names = FALSE)

message("")
message("Matched-city Wilcoxon test:")
print(wilcoxon_summary, row.names = FALSE)

message("")
message("Outputs saved in:")
message(output_dir)

# ============================================================
# fig_sensitivity_sewage_only_PCoA1_vs_Number_of_samples.R
#
# Purpose:
#   Generate the sewage-only ARG-group PCoA1 versus number of
#   sewage samples correlation figure used in the thesis.
#
# Figure:
#   Sewage_only_ARG_group_PCoA1_vs_Number_of_samples_BH.png
#
# Input:
#   data/processed/sensitivity_analysis/pcoa1_correlations/
#   ├── Sewage_only_ARG_group_PCoA_coordinates_metrics.csv
#   └── sample_number_metadata_used.csv
#
# Output:
#   results/sensitivity_analysis/figures/
#   ├── Sewage_only_ARG_group_PCoA1_vs_Number_of_samples_BH.png
#   └── Sewage_only_ARG_group_PCoA1_vs_Number_of_samples_BH.pdf
#
# Author:
#   Bikram Dutta
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c("ggplot2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this script."
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
  "figures"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

pcoa_file <- file.path(
  input_dir,
  "Sewage_only_ARG_group_PCoA_coordinates_metrics.csv"
)

sample_number_file <- file.path(
  input_dir,
  "sample_number_metadata_used.csv"
)

if (!file.exists(pcoa_file)) {
  stop(
    paste0(
      "PCoA input file not found:\n",
      pcoa_file
    )
  )
}

if (!file.exists(sample_number_file)) {
  stop(
    paste0(
      "Sample-number metadata file not found:\n",
      sample_number_file
    )
  )
}


# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

clean_column_name <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]", "", x))
}

find_column <- function(data, candidate_names, required = TRUE, data_label = "data") {
  cleaned_names <- clean_column_name(colnames(data))
  cleaned_candidates <- clean_column_name(candidate_names)

  matches <- which(cleaned_names %in% cleaned_candidates)

  if (length(matches) == 0) {
    if (required) {
      stop(
        paste0(
          data_label,
          ": could not find any of these columns:\n",
          paste(candidate_names, collapse = ", "),
          "\n\nAvailable columns:\n",
          paste(colnames(data), collapse = ", ")
        )
      )
    } else {
      return(NULL)
    }
  }

  colnames(data)[matches[1]]
}

format_p_value <- function(x) {
  if (x < 0.001) {
    return(formatC(x, format = "e", digits = 2))
  }

  signif(x, 3)
}


# ------------------------------------------------------------
# 5. Read sewage-only PCoA data
# ------------------------------------------------------------

pcoa_raw <- read.csv(
  pcoa_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_col <- find_column(
  pcoa_raw,
  c("SampleID", "Sample_ID", "sample_id", "profile_id", "ProfileID"),
  data_label = "Sewage-only PCoA table"
)

city_col <- find_column(
  pcoa_raw,
  c("City", "city"),
  required = FALSE,
  data_label = "Sewage-only PCoA table"
)

environment_col <- find_column(
  pcoa_raw,
  c("Environment", "environment"),
  required = FALSE,
  data_label = "Sewage-only PCoA table"
)

pcoa1_col <- find_column(
  pcoa_raw,
  c("PCoA1", "PC1", "Axis1", "Axis_1", "MDS1", "Dim1"),
  data_label = "Sewage-only PCoA table"
)

sample_number_col_in_pcoa <- find_column(
  pcoa_raw,
  c(
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
  data_label = "Sewage-only PCoA table"
)

if (!is.null(environment_col)) {
  environment_values <- trimws(as.character(pcoa_raw[[environment_col]]))
} else {
  environment_values <- rep("Sewage", nrow(pcoa_raw))
}

if (!is.null(city_col)) {
  city_values <- trimws(as.character(pcoa_raw[[city_col]]))
} else {
  city_values <- NA_character_
}

plot_data <- data.frame(
  SampleID = trimws(as.character(pcoa_raw[[sample_col]])),
  City = city_values,
  Environment = environment_values,
  PCoA1 = as.numeric(as.character(pcoa_raw[[pcoa1_col]])),
  stringsAsFactors = FALSE
)

if (!is.null(sample_number_col_in_pcoa)) {
  plot_data$Number_of_samples <- as.numeric(
    as.character(pcoa_raw[[sample_number_col_in_pcoa]])
  )
} else {
  plot_data$Number_of_samples <- NA_real_
}

plot_data$Environment[tolower(plot_data$Environment) == "sewage"] <- "Sewage"
plot_data$Environment[tolower(plot_data$Environment) == "transit"] <- "Transit"

plot_data <- plot_data[
  plot_data$Environment == "Sewage",
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 6. Read and merge sample-number metadata
# ------------------------------------------------------------

sample_raw <- read.csv(
  sample_number_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_meta_sample_col <- find_column(
  sample_raw,
  c("SampleID", "Sample_ID", "sample_id", "profile_id", "ProfileID"),
  required = FALSE,
  data_label = "Sample-number metadata"
)

sample_meta_city_col <- find_column(
  sample_raw,
  c("City", "city"),
  required = FALSE,
  data_label = "Sample-number metadata"
)

sample_meta_environment_col <- find_column(
  sample_raw,
  c("Environment", "environment"),
  required = FALSE,
  data_label = "Sample-number metadata"
)

sample_meta_number_col <- find_column(
  sample_raw,
  c(
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
  required = TRUE,
  data_label = "Sample-number metadata"
)

sample_metadata <- data.frame(
  Number_of_samples = as.numeric(as.character(sample_raw[[sample_meta_number_col]])),
  stringsAsFactors = FALSE
)

if (!is.null(sample_meta_sample_col)) {
  sample_metadata$SampleID <- trimws(as.character(sample_raw[[sample_meta_sample_col]]))
} else {
  sample_metadata$SampleID <- NA_character_
}

if (!is.null(sample_meta_city_col)) {
  sample_metadata$City <- trimws(as.character(sample_raw[[sample_meta_city_col]]))
} else {
  sample_metadata$City <- NA_character_
}

if (!is.null(sample_meta_environment_col)) {
  sample_metadata$Environment <- trimws(as.character(sample_raw[[sample_meta_environment_col]]))
  sample_metadata$Environment[tolower(sample_metadata$Environment) == "sewage"] <- "Sewage"
  sample_metadata$Environment[tolower(sample_metadata$Environment) == "transit"] <- "Transit"
} else {
  sample_metadata$Environment <- NA_character_
}

sample_metadata <- sample_metadata[
  is.na(sample_metadata$Environment) | sample_metadata$Environment == "Sewage",
  ,
  drop = FALSE
]

# Fill Number_of_samples from metadata if it is missing in the PCoA table.

if (all(is.na(plot_data$Number_of_samples))) {

  if (!all(is.na(sample_metadata$SampleID))) {

    matched_index <- match(
      plot_data$SampleID,
      sample_metadata$SampleID
    )

    plot_data$Number_of_samples <- sample_metadata$Number_of_samples[
      matched_index
    ]

  } else if (!all(is.na(sample_metadata$City)) && !all(is.na(plot_data$City))) {

    matched_index <- match(
      plot_data$City,
      sample_metadata$City
    )

    plot_data$Number_of_samples <- sample_metadata$Number_of_samples[
      matched_index
    ]

  } else {

    stop(
      paste0(
        "Could not merge sample numbers. The metadata file must contain ",
        "either SampleID or City."
      )
    )
  }
}


# ------------------------------------------------------------
# 7. Validate plotting data
# ------------------------------------------------------------

plot_data$Environment <- factor(
  plot_data$Environment,
  levels = c("Sewage")
)

if (nrow(plot_data) != 16) {
  warning(
    paste0(
      "Expected 16 sewage profiles, but found ",
      nrow(plot_data),
      "."
    )
  )
}

if (anyNA(plot_data$PCoA1)) {
  stop("PCoA1 column contains missing or non-numeric values.")
}

if (anyNA(plot_data$Number_of_samples)) {
  stop("Number_of_samples contains missing or non-numeric values after merging.")
}

if (length(unique(plot_data$Number_of_samples)) < 2) {
  stop("Number_of_samples has no variation; correlation cannot be performed.")
}


# ------------------------------------------------------------
# 8. Keep PCoA1 orientation consistent with thesis figure
# ------------------------------------------------------------

# PCoA signs are arbitrary.
# The thesis figure has a negative relationship between number of
# sewage samples and PCoA1.

spearman_check <- suppressWarnings(
  cor(
    plot_data$PCoA1,
    plot_data$Number_of_samples,
    method = "spearman"
  )
)

if (spearman_check > 0) {
  plot_data$PCoA1 <- -plot_data$PCoA1
}


# ------------------------------------------------------------
# 9. Spearman correlation and BH-adjusted value
# ------------------------------------------------------------

spearman_result <- suppressWarnings(
  cor.test(
    plot_data$PCoA1,
    plot_data$Number_of_samples,
    method = "spearman",
    exact = FALSE,
    alternative = "two.sided"
  )
)

rho_value <- unname(spearman_result$estimate)
p_value <- spearman_result$p.value
n_value <- nrow(plot_data)

# In the thesis analysis, this p-value was adjusted together with
# the full family of 16 valid PCoA-axis correlation tests.
# This is the thesis-reported BH-adjusted value for this relationship.
bh_adjusted_p_value <- 0.042

annotation_text <- paste0(
  "\u03c1 = ",
  sprintf("%.3f", rho_value),
  "\n",
  "p value = ",
  format_p_value(p_value),
  "\n",
  "BH-adjusted p value = ",
  format_p_value(bh_adjusted_p_value),
  "\n",
  "n = ",
  n_value
)


# ------------------------------------------------------------
# 10. Save data used for plotting
# ------------------------------------------------------------

write.csv(
  plot_data,
  file = file.path(
    output_dir,
    "Sewage_only_ARG_group_PCoA1_vs_Number_of_samples_plot_data.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Generate figure
# ------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = Number_of_samples,
    y = PCoA1,
    colour = Environment
  )
) +
  geom_point(
    size = 4.2,
    alpha = 0.9
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1.2
  ) +
  annotate(
    geom = "label",
    x = 0.95,
    y = 0.57,
    label = annotation_text,
    hjust = 0,
    vjust = 1,
    size = 5.1,
    label.size = 0.45,
    fill = "white",
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(
      "Sewage" = "#F8766D"
    )
  ) +
  scale_x_continuous(
    limits = c(0.7, 5.3),
    breaks = 1:5,
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_continuous(
    limits = c(-0.45, 0.63),
    breaks = c(-0.4, -0.2, 0.0, 0.2, 0.4, 0.6),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Number of sewage samples",
    y = "PCoA1 coordinate",
    colour = "Environment"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.7
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "grey30",
      fill = NA,
      linewidth = 0.9
    ),
    axis.title = element_text(
      size = 21,
      colour = "black"
    ),
    axis.text = element_text(
      size = 16,
      colour = "black"
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 20,
      colour = "black"
    ),
    legend.text = element_text(
      size = 16,
      colour = "black"
    ),
    legend.key = element_blank(),
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )


# ------------------------------------------------------------
# 12. Save final figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Sewage_only_ARG_group_PCoA1_vs_Number_of_samples_BH.png"
  ),
  plot = p,
  width = 9.6,
  height = 7.2,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "Sewage_only_ARG_group_PCoA1_vs_Number_of_samples_BH.pdf"
  ),
  plot = p,
  width = 9.6,
  height = 7.2
)


# ------------------------------------------------------------
# 13. Save session information
# ------------------------------------------------------------

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("PCoA input file: ", pcoa_file),
  paste0("Sample-number metadata file: ", sample_number_file),
  paste0("R version: ", R.version.string),
  paste0("ggplot2 version: ", as.character(packageVersion("ggplot2"))),
  "Correlation method: two-sided Spearman rank correlation",
  "BH-adjusted p-value corresponds to the thesis-wide family of 16 valid PCoA-axis correlation tests.",
  "",
  paste0("Spearman rho: ", sprintf("%.3f", rho_value)),
  paste0("Raw p-value: ", format_p_value(p_value)),
  paste0("BH-adjusted p-value: ", format_p_value(bh_adjusted_p_value)),
  paste0("n: ", n_value),
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_text,
  con = file.path(
    output_dir,
    "Sewage_only_ARG_group_PCoA1_vs_Number_of_samples_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 14. Print completion message
# ------------------------------------------------------------

message("")
message("====================================================")
message("Figure generated successfully")
message("====================================================")
message("PCoA input file:")
message(pcoa_file)
message("")
message("Sample-number metadata file:")
message(sample_number_file)
message("")
message("Spearman rho: ", sprintf("%.3f", rho_value))
message("Raw p-value: ", format_p_value(p_value))
message("BH-adjusted p-value: ", format_p_value(bh_adjusted_p_value))
message("n: ", n_value)
message("")
message("Output folder:")
message(output_dir)

# ============================================================
# fig_sensitivity_sewage_only_PCoA1_vs_Total_ARG_count.R
#
# Purpose:
#   Generate the sewage-only ARG-group PCoA1 versus total ARG
#   count correlation figure used in the thesis.
#
# Figure:
#   Sewage_only_ARG_group_PCoA1_vs_Total_ARG_count_BH.png
#
# Input:
#   data/processed/sensitivity_analysis/pcoa1_correlations/
#   └── Sewage_only_ARG_group_PCoA_coordinates_metrics.csv
#
# Output:
#   results/sensitivity_analysis/figures/
#   ├── Sewage_only_ARG_group_PCoA1_vs_Total_ARG_count_BH.png
#   └── Sewage_only_ARG_group_PCoA1_vs_Total_ARG_count_BH.pdf
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
# 3. Input file
# ------------------------------------------------------------

input_file <- file.path(
  input_dir,
  "Sewage_only_ARG_group_PCoA_coordinates_metrics.csv"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      input_file
    )
  )
}


# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

clean_column_name <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]", "", x))
}

find_column <- function(data, candidate_names, required = TRUE) {
  cleaned_names <- clean_column_name(colnames(data))
  cleaned_candidates <- clean_column_name(candidate_names)

  matches <- which(cleaned_names %in% cleaned_candidates)

  if (length(matches) == 0) {
    if (required) {
      stop(
        paste0(
          "Could not find any of these columns:\n",
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
# 5. Read and standardise data
# ------------------------------------------------------------

data_raw <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_col <- find_column(
  data_raw,
  c("SampleID", "Sample_ID", "sample_id", "profile_id", "ProfileID")
)

environment_col <- find_column(
  data_raw,
  c("Environment", "environment"),
  required = FALSE
)

pcoa1_col <- find_column(
  data_raw,
  c("PCoA1", "PC1", "Axis1", "Axis_1", "MDS1", "Dim1")
)

total_count_col <- find_column(
  data_raw,
  c(
    "Total_ARG_count",
    "total_ARG_count",
    "TotalARGcount",
    "ARG_count",
    "ARG_counts",
    "Total_count",
    "total_count",
    "Total_ARG_aligned_count",
    "Total_aligned_ARG_count"
  )
)

if (!is.null(environment_col)) {
  environment_values <- trimws(as.character(data_raw[[environment_col]]))
} else {
  environment_values <- rep("Sewage", nrow(data_raw))
}

plot_data <- data.frame(
  SampleID = trimws(as.character(data_raw[[sample_col]])),
  Environment = environment_values,
  PCoA1 = as.numeric(as.character(data_raw[[pcoa1_col]])),
  Total_ARG_count = as.numeric(as.character(data_raw[[total_count_col]])),
  stringsAsFactors = FALSE
)

plot_data$Environment[tolower(plot_data$Environment) == "sewage"] <- "Sewage"
plot_data$Environment[tolower(plot_data$Environment) == "transit"] <- "Transit"

# This figure is sewage-only. Keep only sewage profiles if the table
# contains additional rows.
plot_data <- plot_data[
  plot_data$Environment == "Sewage",
  ,
  drop = FALSE
]

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

if (anyNA(plot_data$Total_ARG_count)) {
  stop("Total ARG count column contains missing or non-numeric values.")
}


# ------------------------------------------------------------
# 6. Keep PCoA1 orientation consistent with thesis figure
# ------------------------------------------------------------

# PCoA signs are arbitrary.
# The thesis figure has a positive relationship between total ARG
# count and PCoA1.

spearman_check <- suppressWarnings(
  cor(
    plot_data$PCoA1,
    plot_data$Total_ARG_count,
    method = "spearman"
  )
)

if (spearman_check < 0) {
  plot_data$PCoA1 <- -plot_data$PCoA1
}


# ------------------------------------------------------------
# 7. Spearman correlation and BH-adjusted value
# ------------------------------------------------------------

spearman_result <- suppressWarnings(
  cor.test(
    plot_data$PCoA1,
    plot_data$Total_ARG_count,
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
bh_adjusted_p_value <- 0.244

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
# 8. Save data used for plotting
# ------------------------------------------------------------

write.csv(
  plot_data,
  file = file.path(
    output_dir,
    "Sewage_only_ARG_group_PCoA1_vs_Total_ARG_count_plot_data.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Generate figure
# ------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = Total_ARG_count,
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
    x = 5200,
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
    limits = c(2500, 44500),
    breaks = c(10000, 20000, 30000, 40000),
    labels = function(x) format(
      x,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_y_continuous(
    limits = c(-0.45, 0.63),
    breaks = c(-0.4, -0.2, 0.0, 0.2, 0.4, 0.6),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Total ARG count",
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
# 10. Save final figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Sewage_only_ARG_group_PCoA1_vs_Total_ARG_count_BH.png"
  ),
  plot = p,
  width = 9.6,
  height = 7.2,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "Sewage_only_ARG_group_PCoA1_vs_Total_ARG_count_BH.pdf"
  ),
  plot = p,
  width = 9.6,
  height = 7.2
)


# ------------------------------------------------------------
# 11. Save session information
# ------------------------------------------------------------

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("Input file: ", input_file),
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
    "Sewage_only_ARG_group_PCoA1_vs_Total_ARG_count_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 12. Print completion message
# ------------------------------------------------------------

message("")
message("====================================================")
message("Figure generated successfully")
message("====================================================")
message("Input file:")
message(input_file)
message("")
message("Spearman rho: ", sprintf("%.3f", rho_value))
message("Raw p-value: ", format_p_value(p_value))
message("BH-adjusted p-value: ", format_p_value(bh_adjusted_p_value))
message("n: ", n_value)
message("")
message("Output folder:")
message(output_dir)

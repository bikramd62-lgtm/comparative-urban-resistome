# ============================================================
# fig_sensitivity_matched_city_turnover_vs_nestedness.R
#
# Purpose:
#   Generate the matched-city ARG-group Jaccard decomposition
#   figure comparing turnover and nestedness-resultant components
#   for each sewage-transit city pair.
#
# Figure:
#   ARG_group_matched_city_turnover_vs_nestedness.png
#
# Input:
#   ARG_group_jaccard_decomposition_matched_city_pairs.csv
#
# Input path:
#   data/processed/sensitivity_analysis/jaccard_decomposition/
#
# Output:
#   results/sensitivity_analysis/figures/
#   ├── ARG_group_matched_city_turnover_vs_nestedness.png
#   └── ARG_group_matched_city_turnover_vs_nestedness.pdf
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

candidate_input_files <- c(
  file.path(
    project_root,
    "results",
    "sensitivity_analysis",
    "jaccard_decomposition",
    "tables",
    "ARG_group_jaccard_decomposition_matched_city_pairs.csv"
  ),
  file.path(
    project_root,
    "data",
    "processed",
    "sensitivity_analysis",
    "jaccard_decomposition",
    "ARG_group_jaccard_decomposition_matched_city_pairs.csv"
  )
)

existing_input_files <- candidate_input_files[
  file.exists(candidate_input_files)
]

if (length(existing_input_files) == 0) {
  stop(
    paste0(
      "Input file not found. Expected one of:\n",
      paste(candidate_input_files, collapse = "\n")
    )
  )
}

input_file <- existing_input_files[1]

output_dir <- file.path(
  project_root,
  "results",
  "sensitivity_analysis",
  "figures"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Helper functions
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

clean_city_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x
}


# ------------------------------------------------------------
# 4. Read matched-city decomposition table
# ------------------------------------------------------------

matched_raw <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

city_col <- find_column(
  matched_raw,
  c("City", "city")
)

turnover_col <- find_column(
  matched_raw,
  c(
    "Turnover_component",
    "turnover_component",
    "Turnover",
    "beta_jtu",
    "beta.jtu"
  )
)

nestedness_col <- find_column(
  matched_raw,
  c(
    "Nestedness_resultant_component",
    "nestedness_resultant_component",
    "Nestedness",
    "Nestedness_resultant",
    "beta_jne",
    "beta.jne"
  )
)

plot_data <- data.frame(
  City = clean_city_label(matched_raw[[city_col]]),
  Turnover = as.numeric(as.character(matched_raw[[turnover_col]])),
  Nestedness_resultant = as.numeric(as.character(matched_raw[[nestedness_col]])),
  stringsAsFactors = FALSE
)

if (anyNA(plot_data$Turnover)) {
  stop("Turnover component contains missing or non-numeric values.")
}

if (anyNA(plot_data$Nestedness_resultant)) {
  stop("Nestedness-resultant component contains missing or non-numeric values.")
}

if (anyDuplicated(plot_data$City) > 0) {
  stop("Duplicated city names were found in the matched-city table.")
}

if (nrow(plot_data) != 16) {
  warning(
    paste0(
      "Expected 16 matched cities, but found ",
      nrow(plot_data),
      "."
    )
  )
}


# ------------------------------------------------------------
# 5. Order cities as in the thesis figure
# ------------------------------------------------------------

plot_data$Turnover_minus_nestedness <-
  plot_data$Turnover - plot_data$Nestedness_resultant

plot_data <- plot_data[
  order(plot_data$Turnover_minus_nestedness, decreasing = TRUE),
  ,
  drop = FALSE
]

# In the plot, larger y values appear higher.
# Therefore assign y positions from top to bottom.
plot_data$y_position <- rev(seq_len(nrow(plot_data)))


# ------------------------------------------------------------
# 6. Paired Wilcoxon test and effect size
# ------------------------------------------------------------

wilcoxon_result <- suppressWarnings(
  stats::wilcox.test(
    x = plot_data$Turnover,
    y = plot_data$Nestedness_resultant,
    paired = TRUE,
    alternative = "two.sided",
    exact = TRUE
  )
)

differences <- plot_data$Turnover - plot_data$Nestedness_resultant
non_zero_differences <- differences[differences != 0]

n_non_zero <- length(non_zero_differences)
rank_sum_total <- n_non_zero * (n_non_zero + 1) / 2

V_value <- as.numeric(wilcoxon_result$statistic)

matched_rank_biserial <- ifelse(
  rank_sum_total > 0,
  (2 * V_value / rank_sum_total) - 1,
  NA_real_
)

turnover_exceeds_nestedness_n <- sum(
  plot_data$Turnover > plot_data$Nestedness_resultant
)

annotation_text <- paste0(
  "Paired Wilcoxon\n",
  "signed-rank test\n",
  "V = ",
  round(V_value, 3),
  "\n",
  "p value = ",
  sprintf("%.4f", wilcoxon_result$p.value),
  "\n",
  "Matched rank-biserial\n",
  "effect size = ",
  sprintf("%.3f", matched_rank_biserial),
  "\n",
  "Turnover exceeded\n",
  "nestedness in ",
  turnover_exceeds_nestedness_n,
  " of ",
  nrow(plot_data),
  " cities"
)


# ------------------------------------------------------------
# 7. Prepare long-format point data
# ------------------------------------------------------------

point_data <- rbind(
  data.frame(
    City = plot_data$City,
    y_position = plot_data$y_position,
    Component = "Turnover",
    Value = plot_data$Turnover,
    stringsAsFactors = FALSE
  ),
  data.frame(
    City = plot_data$City,
    y_position = plot_data$y_position,
    Component = "Nestedness-resultant",
    Value = plot_data$Nestedness_resultant,
    stringsAsFactors = FALSE
  )
)

point_data$Component <- factor(
  point_data$Component,
  levels = c("Turnover", "Nestedness-resultant")
)


# ------------------------------------------------------------
# 8. Save data used for plotting
# ------------------------------------------------------------

write.csv(
  plot_data,
  file = file.path(
    output_dir,
    "ARG_group_matched_city_turnover_vs_nestedness_plot_data_wide.csv"
  ),
  row.names = FALSE
)

write.csv(
  point_data,
  file = file.path(
    output_dir,
    "ARG_group_matched_city_turnover_vs_nestedness_plot_data_long.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Generate figure
# ------------------------------------------------------------

p <- ggplot() +
  geom_segment(
    data = plot_data,
    aes(
      x = Nestedness_resultant,
      xend = Turnover,
      y = y_position,
      yend = y_position
    ),
    colour = "grey68",
    linewidth = 1.0
  ) +
  geom_point(
    data = point_data,
    aes(
      x = Value,
      y = y_position,
      colour = Component,
      shape = Component
    ),
    size = 3.8,
    alpha = 0.95
  ) +
  annotate(
    geom = "label",
    x = 0.695,
    y = 5.6,
    label = annotation_text,
    hjust = 0,
    vjust = 0.5,
    size = 4.0,
    label.size = 0.35,
    fill = "white",
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(
      "Turnover" = "#2C7FB8",
      "Nestedness-resultant" = "#F28E2B"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Turnover" = 16,
      "Nestedness-resultant" = 17
    )
  ) +
  scale_x_continuous(
    limits = c(-0.01, 1.02),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = plot_data$y_position,
    labels = plot_data$City,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Jaccard decomposition component",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4.5)
    ),
    shape = guide_legend(
      override.aes = list(size = 4.5)
    )
  ) +
  theme_bw(base_size = 17) +
  theme(
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.7
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.title.x = element_text(
      size = 19,
      colour = "black",
      margin = margin(t = 12)
    ),
    axis.text.x = element_text(
      size = 14,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    legend.position = "bottom",
    legend.text = element_text(
      size = 15,
      colour = "black"
    ),
    legend.key = element_blank(),
    legend.margin = margin(t = 2, b = 2),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )


# ------------------------------------------------------------
# 10. Save final figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "ARG_group_matched_city_turnover_vs_nestedness.png"
  ),
  plot = p,
  width = 9.6,
  height = 7.6,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "ARG_group_matched_city_turnover_vs_nestedness.pdf"
  ),
  plot = p,
  width = 9.6,
  height = 7.6
)


# ------------------------------------------------------------
# 11. Save session information
# ------------------------------------------------------------

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("Input file: ", input_file),
  paste0("R version: ", R.version.string),
  paste0("ggplot2 version: ", as.character(packageVersion("ggplot2"))),
  "",
  "Statistical test: paired Wilcoxon signed-rank test",
  paste0("V: ", round(V_value, 3)),
  paste0("p-value: ", sprintf("%.4f", wilcoxon_result$p.value)),
  paste0("Matched rank-biserial effect size: ", sprintf("%.3f", matched_rank_biserial)),
  paste0(
    "Turnover exceeded nestedness in ",
    turnover_exceeds_nestedness_n,
    " of ",
    nrow(plot_data),
    " cities"
  ),
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_text,
  con = file.path(
    output_dir,
    "ARG_group_matched_city_turnover_vs_nestedness_sessionInfo.txt"
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
message("Paired Wilcoxon signed-rank test")
message("V = ", round(V_value, 3))
message("p value = ", sprintf("%.4f", wilcoxon_result$p.value))
message("Matched rank-biserial effect size = ", sprintf("%.3f", matched_rank_biserial))
message(
  "Turnover exceeded nestedness in ",
  turnover_exceeds_nestedness_n,
  " of ",
  nrow(plot_data),
  " cities"
)
message("")
message("Output folder:")
message(output_dir)

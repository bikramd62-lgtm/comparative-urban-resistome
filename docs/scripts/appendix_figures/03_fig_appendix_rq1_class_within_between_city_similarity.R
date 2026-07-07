# ============================================================
# fig_appendix_rq1_class_within_between_city_similarity.R
#
# Purpose:
# Generate the appendix Resistance Class within-city vs
# between-city Jaccard similarity figure for the RQ1 comparative
# resistome analysis.
#
# Input:
# results/rq1/within_between_similarity/rq1_class_within_between_city_similarity_plotdata.tsv
#
# Output:
# results/appendix_figures/appendix_rq1_class_within_between_city_jaccard_similarity.png
# results/appendix_figures/appendix_rq1_class_within_between_city_jaccard_similarity.pdf
# ============================================================

# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "dplyr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    FUN.VALUE = logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following required R packages are missing:\n",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this script."
    )
  )
}

library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# 2. User-adjustable paths
# ------------------------------------------------------------

input_dir <- file.path(
  "results",
  "rq1",
  "within_between_similarity"
)

output_dir <- file.path(
  "results",
  "appendix_figures"
)

input_file <- file.path(
  input_dir,
  "rq1_class_within_between_city_similarity_plotdata.tsv"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      input_file,
      "\nRun 04_rq1_within_between_city_similarity.R before generating this appendix figure."
    )
  )
}

# ------------------------------------------------------------
# 3. Read plot data
# ------------------------------------------------------------

plot_data <- read.delim(
  file = input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Comparison",
  "Jaccard_similarity"
)

missing_columns <- setdiff(
  required_columns,
  colnames(plot_data)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the input file:\n",
      paste(missing_columns, collapse = ", ")
    )
  )
}

plot_data$Jaccard_similarity <- suppressWarnings(
  as.numeric(
    plot_data$Jaccard_similarity
  )
)

if (any(is.na(plot_data$Jaccard_similarity))) {
  stop(
    "Jaccard_similarity contains missing or non-numeric values."
  )
}

# ------------------------------------------------------------
# 4. Standardise comparison labels and order
# ------------------------------------------------------------

plot_data$Comparison <- gsub(
  "\\\\n",
  "\n",
  plot_data$Comparison
)

plot_data$Comparison <- factor(
  plot_data$Comparison,
  levels = c(
    "Between-city\n(Sewage vs Transit)",
    "Within-city\n(Sewage vs Transit)"
  )
)

if (any(is.na(plot_data$Comparison))) {
  stop(
    "Unexpected values found in the Comparison column."
  )
}

# ------------------------------------------------------------
# 5. Summary statistics and sample-size labels
# ------------------------------------------------------------

summary_stats <- plot_data %>%
  group_by(
    Comparison
  ) %>%
  summarise(
    n = n(),
    mean = mean(Jaccard_similarity),
    median = median(Jaccard_similarity),
    Q1 = quantile(
      Jaccard_similarity,
      0.25
    ),
    Q3 = quantile(
      Jaccard_similarity,
      0.75
    ),
    IQR = IQR(
      Jaccard_similarity
    ),
    min = min(
      Jaccard_similarity
    ),
    max = max(
      Jaccard_similarity
    ),
    .groups = "drop"
  )

annotation_df <- summary_stats %>%
  mutate(
    label = paste0(
      "n = ",
      n
    ),
    y_position = 0.96
  )

# ------------------------------------------------------------
# 6. Define colours and y-axis
# ------------------------------------------------------------

box_fill_colours <- c(
  "Between-city\n(Sewage vs Transit)" = "#F4A29A",
  "Within-city\n(Sewage vs Transit)" = "#8DB9E2"
)

box_line_colours <- c(
  "Between-city\n(Sewage vs Transit)" = "#E64B35",
  "Within-city\n(Sewage vs Transit)" = "#2166AC"
)

y_axis_upper <- 1.01

y_axis_breaks <- seq(
  0,
  1.0,
  by = 0.2
)

# ------------------------------------------------------------
# 7. Generate final appendix figure
# ------------------------------------------------------------

set.seed(123)

p <- ggplot(
  plot_data,
  aes(
    x = Comparison,
    y = Jaccard_similarity,
    fill = Comparison,
    colour = Comparison
  )
) +

  geom_boxplot(
    width = 0.52,
    alpha = 0.55,
    linewidth = 1.05,
    outlier.shape = NA
  ) +

  geom_jitter(
    width = 0.105,
    height = 0,
    size = 0.75,
    alpha = 0.34,
    shape = 16,
    stroke = 0
  ) +

  geom_text(
    data = annotation_df,
    aes(
      x = Comparison,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5.1,
    colour = "black"
  ) +

  scale_fill_manual(
    values = box_fill_colours,
    guide = "none"
  ) +

  scale_colour_manual(
    values = box_line_colours,
    guide = "none"
  ) +

  scale_y_continuous(
    limits = c(
      0,
      y_axis_upper
    ),
    breaks = y_axis_breaks,
    labels = function(x) {
      sprintf(
        "%.1f",
        x
      )
    },
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    x = NULL,
    y = "Jaccard similarity"
  ) +

  theme_minimal(
    base_size = 18
  ) +

  theme(
    panel.grid.major.x = element_blank(),

    panel.grid.major.y = element_line(
      colour = "grey84",
      linewidth = 0.8
    ),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1.05
    ),

    axis.title.y = element_text(
      size = 20,
      colour = "black",
      margin = margin(
        r = 14
      )
    ),

    axis.text.y = element_text(
      size = 16,
      colour = "black"
    ),

    axis.text.x = element_text(
      size = 17,
      colour = "black",
      lineheight = 1.08,
      margin = margin(
        t = 10
      )
    ),

    axis.ticks.y = element_line(
      colour = "black",
      linewidth = 0.6
    ),

    axis.ticks.x = element_blank(),

    axis.ticks.length = unit(
      4,
      "pt"
    ),

    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),

    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),

    plot.margin = margin(
      t = 16,
      r = 18,
      b = 16,
      l = 18
    )
  )

# ------------------------------------------------------------
# 8. Display figure
# ------------------------------------------------------------

print(p)

# ------------------------------------------------------------
# 9. Save figure and plot data
# ------------------------------------------------------------

png_output <- file.path(
  output_dir,
  "appendix_rq1_class_within_between_city_jaccard_similarity.png"
)

pdf_output <- file.path(
  output_dir,
  "appendix_rq1_class_within_between_city_jaccard_similarity.pdf"
)

plot_data_output <- file.path(
  output_dir,
  "appendix_rq1_class_within_between_city_jaccard_similarity_plotdata.tsv"
)

summary_output <- file.path(
  output_dir,
  "appendix_rq1_class_within_between_city_jaccard_similarity_summary.tsv"
)

ggsave(
  filename = png_output,
  plot = p,
  width = 10.8,
  height = 7.6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_output,
  plot = p,
  width = 10.8,
  height = 7.6,
  units = "in",
  bg = "white"
)

write.table(
  plot_data,
  file = plot_data_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_stats,
  file = summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "\nAppendix figure saved:"
)

message(
  png_output
)

message(
  pdf_output
)

message(
  "\nPlot data saved:"
)

message(
  plot_data_output
)

message(
  "\nSummary saved:"
)

message(
  summary_output
)

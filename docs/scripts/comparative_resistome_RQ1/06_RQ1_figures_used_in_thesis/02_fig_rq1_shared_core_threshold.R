# ============================================================
# fig_rq1_shared_core_threshold.R
#
# Purpose:
# Generate the RQ1 matched-city shared-core threshold figure
# used in the thesis.
#
# Input:
# results/rq1/overlap_core/rq1_combined_shared_core_threshold_summary.tsv
#
# Output:
# results/rq1/figures/rq1_shared_core_threshold_refined.png
# results/rq1/figures/rq1_shared_core_threshold_refined.pdf
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

# ------------------------------------------------------------
# 2. User-adjustable paths
# ------------------------------------------------------------

input_dir <- file.path(
  "results",
  "rq1",
  "overlap_core"
)

output_dir <- file.path(
  "results",
  "rq1",
  "figures"
)

input_file <- file.path(
  input_dir,
  "rq1_combined_shared_core_threshold_summary.tsv"
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
      "\nRun 02_rq1_overlap_and_shared_core.R before generating this figure."
    )
  )
}

# ------------------------------------------------------------
# 3. Read shared-core threshold summary
# ------------------------------------------------------------

shared_core_summary <- read.delim(
  file = input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Level",
  "Threshold_fraction",
  "Number_of_features"
)

missing_columns <- setdiff(
  required_columns,
  colnames(shared_core_summary)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the input file:\n",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 4. Prepare plot data
# ------------------------------------------------------------

plot_data <- shared_core_summary[
  ,
  required_columns,
  drop = FALSE
]

plot_data$Feature_level <- ifelse(
  plot_data$Level == "Resistance Class",
  "ARG Drug Class",
  ifelse(
    plot_data$Level == "ARG Group",
    "ARG Group",
    plot_data$Level
  )
)

plot_data$Threshold <- ifelse(
  plot_data$Threshold_fraction == 1.00,
  "100%",
  ifelse(
    plot_data$Threshold_fraction == 0.75,
    "\u226575%",
    ifelse(
      plot_data$Threshold_fraction == 0.50,
      "\u226550%",
      paste0(
        "\u2265",
        plot_data$Threshold_fraction * 100,
        "%"
      )
    )
  )
)

plot_data$Threshold <- factor(
  plot_data$Threshold,
  levels = c(
    "100%",
    "\u226575%",
    "\u226550%"
  )
)

plot_data$Feature_level <- factor(
  plot_data$Feature_level,
  levels = c(
    "ARG Drug Class",
    "ARG Group"
  )
)

plot_data <- plot_data[
  order(
    plot_data$Threshold,
    plot_data$Feature_level
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 5. Define colours and y-axis range
# ------------------------------------------------------------

bar_fill_colours <- c(
  "ARG Drug Class" = "#9BBFD1",
  "ARG Group" = "#173F7A"
)

maximum_count <- max(
  plot_data$Number_of_features,
  na.rm = TRUE
)

y_axis_upper <- ceiling(
  (maximum_count + 1.0) / 2.5
) * 2.5

y_axis_breaks <- seq(
  0,
  y_axis_upper,
  by = 2.5
)

# ------------------------------------------------------------
# 6. Generate final thesis figure
# ------------------------------------------------------------

dodge_width <- 0.72

p <- ggplot(
  plot_data,
  aes(
    x = Threshold,
    y = Number_of_features,
    fill = Feature_level
  )
) +

  geom_col(
    position = position_dodge(
      width = dodge_width
    ),
    width = 0.62,
    colour = "black",
    linewidth = 0.65
  ) +

  geom_text(
    aes(
      label = Number_of_features
    ),
    position = position_dodge(
      width = dodge_width
    ),
    vjust = -0.25,
    size = 5.8,
    colour = "black"
  ) +

  scale_fill_manual(
    values = bar_fill_colours,
    name = NULL
  ) +

  scale_y_continuous(
    limits = c(
      0,
      y_axis_upper
    ),
    breaks = y_axis_breaks,
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    x = "Matched-city shared threshold",
    y = "Number of shared features"
  ) +

  theme_minimal(
    base_size = 18
  ) +

  theme(
    legend.position = "top",

    legend.direction = "horizontal",

    legend.justification = "center",

    legend.text = element_text(
      size = 16,
      colour = "black"
    ),

    legend.key.size = unit(
      0.9,
      "cm"
    ),

    legend.margin = margin(
      b = 10
    ),

    panel.grid.major.x = element_blank(),

    panel.grid.major.y = element_line(
      colour = "grey84",
      linewidth = 0.8
    ),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      colour = "grey25",
      fill = NA,
      linewidth = 1.15
    ),

    axis.title.x = element_text(
      size = 19,
      colour = "black",
      margin = margin(
        t = 12
      )
    ),

    axis.title.y = element_text(
      size = 19,
      colour = "black",
      margin = margin(
        r = 14
      )
    ),

    axis.text.x = element_text(
      size = 16,
      colour = "grey20",
      margin = margin(
        t = 8
      )
    ),

    axis.text.y = element_text(
      size = 16,
      colour = "grey20"
    ),

    axis.ticks.x = element_line(
      colour = "grey25",
      linewidth = 0.7
    ),

    axis.ticks.y = element_line(
      colour = "grey25",
      linewidth = 0.7
    ),

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
      t = 18,
      r = 20,
      b = 18,
      l = 20
    )
  )

# ------------------------------------------------------------
# 7. Display figure
# ------------------------------------------------------------

print(p)

# ------------------------------------------------------------
# 8. Save figure and plot data
# ------------------------------------------------------------

png_output <- file.path(
  output_dir,
  "rq1_shared_core_threshold_refined.png"
)

pdf_output <- file.path(
  output_dir,
  "rq1_shared_core_threshold_refined.pdf"
)

plot_data_output <- file.path(
  output_dir,
  "rq1_shared_core_threshold_plotdata.tsv"
)

ggsave(
  filename = png_output,
  plot = p,
  width = 10.5,
  height = 6.6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_output,
  plot = p,
  width = 10.5,
  height = 6.6,
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

message(
  "\nFigure saved:"
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

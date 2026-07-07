# ============================================================
# fig_rq1_group_overlap_summary.R
#
# Purpose:
# Generate the RQ1 ARG Group overlap summary figure used in
# the thesis.
#
# Input:
# results/rq1/overlap_core/rq1_group_global_overlap_summary.tsv
#
# Output:
# results/rq1/figures/rq1_group_overlap_summary_refined.png
# results/rq1/figures/rq1_group_overlap_summary_refined.pdf
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
  "rq1_group_global_overlap_summary.tsv"
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
# 3. Read overlap summary
# ------------------------------------------------------------

overlap_summary <- read.delim(
  file = input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Sewage_detected_features",
  "Transit_detected_features",
  "Shared_features",
  "Sewage_only_features",
  "Transit_only_features"
)

missing_columns <- setdiff(
  required_columns,
  colnames(overlap_summary)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the input file:\n",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (nrow(overlap_summary) != 1) {
  stop(
    "The ARG Group overlap summary file should contain exactly one row."
  )
}

# ------------------------------------------------------------
# 4. Prepare plot data
# ------------------------------------------------------------

plot_data <- data.frame(
  Category = c(
    "Total in\nsewage",
    "Total in\ntransit",
    "Shared",
    "Sewage-only",
    "Transit-only"
  ),
  Count = c(
    overlap_summary$Sewage_detected_features,
    overlap_summary$Transit_detected_features,
    overlap_summary$Shared_features,
    overlap_summary$Sewage_only_features,
    overlap_summary$Transit_only_features
  ),
  stringsAsFactors = FALSE
)

plot_data$Category <- factor(
  plot_data$Category,
  levels = plot_data$Category
)

# ------------------------------------------------------------
# 5. Define y-axis range
# ------------------------------------------------------------

y_axis_upper <- ceiling(
  max(plot_data$Count) / 50
) * 50 + 25

y_axis_breaks <- seq(
  0,
  y_axis_upper,
  by = 50
)

# ------------------------------------------------------------
# 6. Generate final thesis figure
# ------------------------------------------------------------

bar_colour <- "#173F7A"

p <- ggplot(
  plot_data,
  aes(
    x = Category,
    y = Count
  )
) +

  geom_col(
    fill = bar_colour,
    colour = "black",
    linewidth = 0.65,
    width = 0.56
  ) +

  geom_text(
    aes(
      label = Count
    ),
    vjust = -0.28,
    size = 5.8,
    colour = "black"
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
    x = NULL,
    y = "Number of ARG Groups"
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
      colour = "grey25",
      fill = NA,
      linewidth = 1.15
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
      colour = "grey20"
    ),

    axis.text.x = element_text(
      size = 16,
      colour = "grey20",
      lineheight = 1.05,
      margin = margin(
        t = 8
      )
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
  "rq1_group_overlap_summary_refined.png"
)

pdf_output <- file.path(
  output_dir,
  "rq1_group_overlap_summary_refined.pdf"
)

plot_data_output <- file.path(
  output_dir,
  "rq1_group_overlap_summary_plotdata.tsv"
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

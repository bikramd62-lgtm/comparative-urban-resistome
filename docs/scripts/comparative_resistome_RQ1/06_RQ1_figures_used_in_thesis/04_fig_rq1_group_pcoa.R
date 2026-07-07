# ============================================================
# fig_rq1_group_pcoa.R
#
# Purpose:
# Generate the RQ1 ARG Group Jaccard PCoA figure used in the
# thesis.
#
# Input:
# results/rq1/ordination_stats/rq1_group_pcoa_coordinates.tsv
# results/rq1/ordination_stats/rq1_group_pcoa_eigenvalue_summary.tsv
#
# Output:
# results/rq1/figures/rq1_group_pcoa_final_thesis.png
# results/rq1/figures/rq1_group_pcoa_final_thesis.pdf
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
  "ordination_stats"
)

output_dir <- file.path(
  "results",
  "rq1",
  "figures"
)

coordinates_file <- file.path(
  input_dir,
  "rq1_group_pcoa_coordinates.tsv"
)

eigenvalue_summary_file <- file.path(
  input_dir,
  "rq1_group_pcoa_eigenvalue_summary.tsv"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

required_input_files <- c(
  coordinates_file,
  eigenvalue_summary_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "The following required input files are missing:\n",
      paste(missing_input_files, collapse = "\n"),
      "\nRun 03_rq1_jaccard_pcoa_permanova.R before generating this figure."
    )
  )
}

# ------------------------------------------------------------
# 3. Read PCoA coordinates
# ------------------------------------------------------------

pcoa_coordinates <- read.delim(
  file = coordinates_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_coordinate_columns <- c(
  "Profile",
  "City",
  "Environment",
  "PCoA1",
  "PCoA2"
)

missing_coordinate_columns <- setdiff(
  required_coordinate_columns,
  colnames(pcoa_coordinates)
)

if (length(missing_coordinate_columns) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the coordinate file:\n",
      paste(missing_coordinate_columns, collapse = ", ")
    )
  )
}

pcoa_coordinates$Environment <- factor(
  pcoa_coordinates$Environment,
  levels = c(
    "Sewage",
    "Transit"
  )
)

if (any(is.na(pcoa_coordinates$Environment))) {
  stop(
    "Environment column must contain only 'Sewage' and 'Transit'."
  )
}

# ------------------------------------------------------------
# 4. Read axis percentages
# ------------------------------------------------------------

eigenvalue_summary <- read.delim(
  file = eigenvalue_summary_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_eigen_columns <- c(
  "PCoA1_percent_positive_inertia",
  "PCoA2_percent_positive_inertia"
)

missing_eigen_columns <- setdiff(
  required_eigen_columns,
  colnames(eigenvalue_summary)
)

if (length(missing_eigen_columns) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the eigenvalue summary file:\n",
      paste(missing_eigen_columns, collapse = ", ")
    )
  )
}

if (nrow(eigenvalue_summary) != 1) {
  stop(
    "The ARG Group eigenvalue summary file should contain exactly one row."
  )
}

pcoa1_percent <- eigenvalue_summary$PCoA1_percent_positive_inertia[1]
pcoa2_percent <- eigenvalue_summary$PCoA2_percent_positive_inertia[1]

x_axis_label <- paste0(
  "PCoA1 (",
  sprintf("%.1f", pcoa1_percent),
  "%)"
)

y_axis_label <- paste0(
  "PCoA2 (",
  sprintf("%.1f", pcoa2_percent),
  "%)"
)

# ------------------------------------------------------------
# 5. Define colours and shapes
# ------------------------------------------------------------

environment_colours <- c(
  "Sewage" = "#F8766D",
  "Transit" = "#00BFC4"
)

environment_shapes <- c(
  "Sewage" = 16,
  "Transit" = 17
)

# ------------------------------------------------------------
# 6. Generate final thesis figure
# ------------------------------------------------------------

p <- ggplot(
  pcoa_coordinates,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Environment,
    shape = Environment
  )
) +

  geom_point(
    size = 4.8,
    alpha = 0.95
  ) +

  scale_colour_manual(
    values = environment_colours,
    name = "Environment"
  ) +

  scale_shape_manual(
    values = environment_shapes,
    name = "Environment"
  ) +

  labs(
    x = x_axis_label,
    y = y_axis_label
  ) +

  theme_minimal(
    base_size = 18
  ) +

  theme(
    panel.grid.major = element_line(
      colour = "grey84",
      linewidth = 0.85
    ),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      colour = "grey25",
      fill = NA,
      linewidth = 1.05
    ),

    axis.title.x = element_text(
      size = 18,
      colour = "black",
      margin = margin(
        t = 10
      )
    ),

    axis.title.y = element_text(
      size = 18,
      colour = "black",
      margin = margin(
        r = 12
      )
    ),

    axis.text.x = element_text(
      size = 14,
      colour = "grey20"
    ),

    axis.text.y = element_text(
      size = 14,
      colour = "grey20"
    ),

    axis.ticks = element_line(
      colour = "grey25",
      linewidth = 0.6
    ),

    axis.ticks.length = unit(
      4,
      "pt"
    ),

    legend.position = "right",

    legend.title = element_text(
      size = 16,
      colour = "black"
    ),

    legend.text = element_text(
      size = 14,
      colour = "black"
    ),

    legend.key.size = unit(
      0.8,
      "cm"
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
      r = 20,
      b = 16,
      l = 18
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
  "rq1_group_pcoa_final_thesis.png"
)

pdf_output <- file.path(
  output_dir,
  "rq1_group_pcoa_final_thesis.pdf"
)

plot_data_output <- file.path(
  output_dir,
  "rq1_group_pcoa_plotdata.tsv"
)

ggsave(
  filename = png_output,
  plot = p,
  width = 8.8,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_output,
  plot = p,
  width = 8.8,
  height = 7.2,
  units = "in",
  bg = "white"
)

write.table(
  pcoa_coordinates,
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

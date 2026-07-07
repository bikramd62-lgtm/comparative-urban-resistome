# ============================================================
# fig_rq2_group_prevalence_difference_barplot.R
#
# Purpose:
# Generate the RQ2 ARG Group descriptive prevalence-difference
# barplot used in the thesis.
#
# This figure is descriptive. It does not show formal
# statistical significance
#
# Prevalence difference:
# Transit prevalence - Sewage prevalence
#
# Interpretation:
# - Positive bars: higher city-level prevalence in transit
# - Negative bars: higher city-level prevalence in sewage
#
# Input:
# results/rq2/prevalence_difference/rq2_group_prevalence_comparison.tsv
#
# Output:
# results/rq2/figures/rq2_group_prevalence_difference_barplot.png
# results/rq2/figures/rq2_group_prevalence_difference_barplot.pdf
# results/rq2/figures/rq2_group_prevalence_difference_barplot_plotdata.tsv
# ============================================================

# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "dplyr",
  "stringr"
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
library(stringr)

# ------------------------------------------------------------
# 2. User-adjustable settings
# ------------------------------------------------------------

input_file <- file.path(
  "results",
  "rq2",
  "prevalence_difference",
  "rq2_group_prevalence_comparison.tsv"
)

output_dir <- file.path(
  "results",
  "rq2",
  "figures"
)

# The thesis figure shows the 15 ARG Groups with the strongest
# descriptive difference toward sewage and the 15 strongest
# descriptive difference toward transit.
top_n_per_direction <- 15

# ------------------------------------------------------------
# 3. Create output directory
# ------------------------------------------------------------

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# 4. Check input file
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      input_file,
      "\nRun 02_rq2_prevalence_difference_analysis.R before generating this figure."
    )
  )
}

# ------------------------------------------------------------
# 5. Read prevalence-comparison table
# ------------------------------------------------------------

prevalence_table <- read.delim(
  file = input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

required_columns <- c(
  "Feature",
  "Feature_label",
  "Sewage_prevalence",
  "Transit_prevalence",
  "Prevalence_difference",
  "Absolute_prevalence_difference",
  "Prevalence_direction"
)

missing_columns <- setdiff(
  required_columns,
  colnames(prevalence_table)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the input file:\n",
      paste(missing_columns, collapse = ", ")
    )
  )
}

numeric_columns <- c(
  "Sewage_prevalence",
  "Transit_prevalence",
  "Prevalence_difference",
  "Absolute_prevalence_difference"
)

for (column_name in numeric_columns) {
  prevalence_table[[column_name]] <- suppressWarnings(
    as.numeric(
      prevalence_table[[column_name]]
    )
  )

  if (any(is.na(prevalence_table[[column_name]]))) {
    stop(
      paste0(
        "Column contains missing or non-numeric values: ",
        column_name
      )
    )
  }
}

# ------------------------------------------------------------
# 6. Prepare ARG Group labels
# ------------------------------------------------------------

prevalence_table <- prevalence_table %>%
  mutate(
    Feature_label = ifelse(
      is.na(Feature_label) | Feature_label == "",
      Feature,
      Feature_label
    ),
    Feature_label = str_replace_all(
      Feature_label,
      "_",
      " "
    ),
    Feature_label = str_squish(
      Feature_label
    ),
    Feature_label = toupper(
      Feature_label
    )
  )

# ------------------------------------------------------------
# 7. Select strongest prevalence differences
# ------------------------------------------------------------

higher_in_sewage <- prevalence_table %>%
  filter(
    Prevalence_direction == "Higher_in_sewage"
  ) %>%
  arrange(
    desc(Absolute_prevalence_difference),
    desc(Sewage_prevalence),
    Feature_label,
    Feature
  ) %>%
  slice_head(
    n = top_n_per_direction
  )

higher_in_transit <- prevalence_table %>%
  filter(
    Prevalence_direction == "Higher_in_transit"
  ) %>%
  arrange(
    desc(Absolute_prevalence_difference),
    desc(Transit_prevalence),
    Feature_label,
    Feature
  ) %>%
  slice_head(
    n = top_n_per_direction
  )

plot_data <- bind_rows(
  higher_in_sewage,
  higher_in_transit
)

if (nrow(plot_data) == 0) {
  stop(
    "No ARG Groups were selected for plotting. Check the prevalence-difference table."
  )
}

plot_data <- plot_data %>%
  mutate(
    Higher_prevalence = ifelse(
      Prevalence_difference < 0,
      "Sewage",
      "Transit"
    )
  ) %>%
  arrange(
    Prevalence_difference,
    Feature_label,
    Feature
  )

plot_data$Feature_label <- factor(
  plot_data$Feature_label,
  levels = unique(
    plot_data$Feature_label
  )
)

plot_data$Higher_prevalence <- factor(
  plot_data$Higher_prevalence,
  levels = c(
    "Sewage",
    "Transit"
  )
)

# ------------------------------------------------------------
# 8. Define colours and axis formatting
# ------------------------------------------------------------

bar_colours <- c(
  "Sewage" = "#A6CEE3",
  "Transit" = "#F6DA7B"
)

x_axis_limits <- c(
  -1,
  1
)

x_axis_breaks <- seq(
  -1,
  1,
  by = 0.25
)

# ------------------------------------------------------------
# 9. Generate final thesis figure
# ------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = Prevalence_difference,
    y = Feature_label,
    fill = Higher_prevalence
  )
) +

  geom_col(
    width = 0.70,
    colour = "grey35",
    linewidth = 0.35
  ) +

  geom_vline(
    xintercept = 0,
    colour = "black",
    linewidth = 0.8
  ) +

  scale_fill_manual(
    values = bar_colours,
    name = "Higher prevalence"
  ) +

  scale_x_continuous(
    limits = x_axis_limits,
    breaks = x_axis_breaks,
    labels = function(x) {
      sprintf(
        "%.2f",
        x
      )
    },
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +

  labs(
    x = "Prevalence difference (transit proportion - sewage proportion)",
    y = "ARG Group"
  ) +

  guides(
    fill = guide_legend(
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        colour = "grey35",
        linewidth = 0.35
      )
    )
  ) +

  theme_minimal(
    base_size = 18
  ) +

  theme(
    panel.grid.major.y = element_blank(),

    panel.grid.major.x = element_line(
      colour = "grey84",
      linewidth = 0.8
    ),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      colour = "grey25",
      fill = NA,
      linewidth = 1.0
    ),

    axis.title.x = element_text(
      size = 18,
      colour = "black",
      margin = margin(
        t = 12
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
      size = 12,
      colour = "grey20"
    ),

    axis.ticks.x = element_line(
      colour = "grey25",
      linewidth = 0.6
    ),

    axis.ticks.y = element_line(
      colour = "grey25",
      linewidth = 0.6
    ),

    axis.ticks.length = unit(
      4,
      "pt"
    ),

    legend.position = "top",

    legend.title = element_text(
      size = 14,
      colour = "black"
    ),

    legend.text = element_text(
      size = 14,
      colour = "black"
    ),

    legend.key.size = unit(
      0.70,
      "cm"
    ),

    legend.spacing.x = unit(
      0.25,
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
      t = 12,
      r = 18,
      b = 14,
      l = 18
    )
  )

# ------------------------------------------------------------
# 10. Display figure
# ------------------------------------------------------------

print(p)

# ------------------------------------------------------------
# 11. Save figure and plot data
# ------------------------------------------------------------

png_output <- file.path(
  output_dir,
  "rq2_group_prevalence_difference_barplot.png"
)

pdf_output <- file.path(
  output_dir,
  "rq2_group_prevalence_difference_barplot.pdf"
)

plot_data_output <- file.path(
  output_dir,
  "rq2_group_prevalence_difference_barplot_plotdata.tsv"
)

ggsave(
  filename = png_output,
  plot = p,
  width = 12.2,
  height = 9.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_output,
  plot = p,
  width = 12.2,
  height = 9.4,
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

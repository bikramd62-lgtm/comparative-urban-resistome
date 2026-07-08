# ============================================================
# 03_RQ3 Jaccard similarity vs June relative humidity
#
# Figure:
# RQ3 Jaccard similarity vs June relative humidity
#
# Input:
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
#
# Output:
#   results/RQ3_figures/RQ3_jaccard_vs_humidity.png
#   results/RQ3_figures/RQ3_jaccard_vs_humidity.pdf
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_dir  <- "results/RQ3_output_tables"
figure_dir <- "results/RQ3_figures"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

master_file <- file.path(input_dir, "RQ3_master_table_FINAL_corrected.csv")

if (!file.exists(master_file)) {
  stop("Input file not found: ", master_file)
}

# ------------------------------------------------------------
# 2. Load corrected RQ3 master table
# ------------------------------------------------------------

rq3 <- read_csv(master_file, show_col_types = FALSE)

# ------------------------------------------------------------
# 3. Check required columns
# ------------------------------------------------------------

required_columns <- c(
  "city",
  "june_relative_humidity_percent",
  "jaccard_similarity_sewage_transit_ARG_group"
)

missing_columns <- setdiff(required_columns, colnames(rq3))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing from the master table:\n",
    paste(missing_columns, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 4. Prepare plotting data
# ------------------------------------------------------------

plot_df <- rq3 %>%
  transmute(
    city = city,
    city_label = str_replace_all(city, "_", " "),
    humidity_june_percent = june_relative_humidity_percent,
    jaccard_similarity = jaccard_similarity_sewage_transit_ARG_group
  ) %>%
  filter(
    !is.na(humidity_june_percent),
    !is.na(jaccard_similarity)
  )

# ------------------------------------------------------------
# 5. Manual label placement to match thesis figure
# ------------------------------------------------------------

label_offsets <- tibble(
  city = c(
    "Barcelona", "Berlin", "Bogota", "Hanoi", "Hong Kong", "Ilorin",
    "Kuala Lumpur", "Lisbon", "Oslo", "Porto", "Rio de Janeiro",
    "Santiago", "Singapore", "Sofia", "Taipei", "Vienna"
  ),
  dx = c(
    -2.3, -1.7, -0.5, -0.2, -1.1, -0.2,
    -1.5, -1.8, 0.4, -1.6, -1.5,
    -0.7, -1.0, -0.1, -1.7, -0.4
  ),
  dy = c(
    -0.005, -0.005, 0.004, 0.004, -0.005, 0.004,
    0.004, -0.005, -0.005, 0.004, 0.004,
    0.004, 0.004, 0.004, 0.004, 0.004
  )
)

plot_df <- plot_df %>%
  left_join(label_offsets, by = "city") %>%
  mutate(
    dx = if_else(is.na(dx), 0, dx),
    dy = if_else(is.na(dy), 0.004, dy),
    label_x = humidity_june_percent + dx,
    label_y = jaccard_similarity + dy
  )

# ------------------------------------------------------------
# 6. Statistics label
# ------------------------------------------------------------

stat_label <- "\u03c1 = 0.091\np = 0.737\nBH-adjusted p = 0.845"

# ------------------------------------------------------------
# 7. Generate figure
# ------------------------------------------------------------

p_jaccard_humidity <- ggplot(
  plot_df,
  aes(x = humidity_june_percent, y = jaccard_similarity)
) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "#2F66D0",
    linewidth = 1.1
  ) +
  geom_point(
    size = 2.8,
    color = "black"
  ) +
  geom_text(
    aes(x = label_x, y = label_y, label = city_label),
    size = 4.7,
    color = "black"
  ) +
  annotate(
    "label",
    x = 48.4,
    y = 0.272,
    label = stat_label,
    hjust = 0,
    vjust = 1,
    size = 4.8,
    label.size = 0.7,
    fill = "white",
    color = "black"
  ) +
  scale_x_continuous(
    limits = c(46.5, 90.5),
    breaks = c(50, 60, 70, 80, 90)
  ) +
  scale_y_continuous(
    limits = c(0.022, 0.282),
    breaks = c(0.05, 0.10, 0.15, 0.20, 0.25)
  ) +
  labs(
    x = "June Relative Humidity (%)",
    y = "Sewage-Transit Jaccard Similarity"
  ) +
  theme_classic(base_size = 17) +
  theme(
    axis.title = element_text(size = 19, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.22, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    plot.margin = margin(12, 18, 12, 12)
  )

# ------------------------------------------------------------
# 8. Save figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "RQ3_jaccard_vs_humidity.png"),
  plot = p_jaccard_humidity,
  width = 8.8,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(figure_dir, "RQ3_jaccard_vs_humidity.pdf"),
  plot = p_jaccard_humidity,
  width = 8.8,
  height = 6.8
)

message("Saved: RQ3_jaccard_vs_humidity.png and RQ3_jaccard_vs_humidity.pdf")

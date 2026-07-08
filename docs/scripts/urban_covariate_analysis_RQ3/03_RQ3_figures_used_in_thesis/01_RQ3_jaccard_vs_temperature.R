# ============================================================
# 01_RQ3_jaccard_vs_temperature
#
# Figure:
# RQ3 Jaccard similarity vs average June temperature
#
# Input:
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
#
# Output:
#   results/RQ3_figures/RQ3_jaccard_vs_temperature.png
#   results/RQ3_figures/RQ3_jaccard_vs_temperature.pdf
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
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
# 3. Prepare plotting data
# ------------------------------------------------------------

plot_df <- rq3 %>%
  transmute(
    city = city,
    city_label = str_replace_all(city, "_", " "),
    temperature_june_c = city_avg_june_temperature,
    jaccard_similarity = jaccard_similarity_sewage_transit_ARG_group
  ) %>%
  filter(
    !is.na(temperature_june_c),
    !is.na(jaccard_similarity)
  )

# ------------------------------------------------------------
# 4. Manual label placement to match thesis figure
# ------------------------------------------------------------

label_offsets <- tibble(
  city = c(
    "Barcelona", "Berlin", "Bogota", "Hanoi", "Hong Kong", "Ilorin",
    "Kuala Lumpur", "Lisbon", "Oslo", "Porto", "Rio de Janeiro",
    "Santiago", "Singapore", "Sofia", "Taipei", "Vienna"
  ),
  dx = c(
    -0.20, -0.45, -0.30, -0.15, -0.15, -0.10,
    0.25, -0.55, -0.45, -0.90, -0.90,
    -0.40, -0.55, 0.10, -0.10, -0.25
  ),
  dy = c(
    -0.005, -0.005, 0.004, 0.004, -0.005, 0.004,
    0.004, -0.005, -0.005, -0.005, 0.004,
    0.004, 0.004, 0.004, 0.004, 0.004
  )
)

plot_df <- plot_df %>%
  left_join(label_offsets, by = "city") %>%
  mutate(
    dx = if_else(is.na(dx), 0, dx),
    dy = if_else(is.na(dy), 0.004, dy),
    label_x = temperature_june_c + dx,
    label_y = jaccard_similarity + dy
  )

# ------------------------------------------------------------
# 5. Statistics label
# ------------------------------------------------------------

stat_label <- "\u03c1 = 0.050\np = 0.854\nBH-adjusted p = 0.854"

# ------------------------------------------------------------
# 6. Generate plot
# ------------------------------------------------------------

p_jaccard_temperature <- ggplot(
  plot_df,
  aes(x = temperature_june_c, y = jaccard_similarity)
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
    x = 7.7,
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
    limits = c(6.8, 31.6),
    breaks = c(10, 15, 20, 25, 30)
  ) +
  scale_y_continuous(
    limits = c(0.022, 0.282),
    breaks = c(0.05, 0.10, 0.15, 0.20, 0.25)
  ) +
  labs(
    x = expression("Average June Temperature (" * degree * "C)"),
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
# 7. Save figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "RQ3_jaccard_vs_temperature.png"),
  plot = p_jaccard_temperature,
  width = 8.8,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(figure_dir, "RQ3_jaccard_vs_temperature.pdf"),
  plot = p_jaccard_temperature,
  width = 8.8,
  height = 6.8
)

message("Saved: RQ3_jaccard_vs_temperature.png and RQ3_jaccard_vs_temperature.pdf")

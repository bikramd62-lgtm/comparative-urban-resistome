# ============================================================
# 06_RQ3_shared_ARG_group_count_vs_population_density.R
#
# Figure:
# Shared ARG-group count vs population density
#
# Input:
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
#
# Output:
#   results/RQ3_figures/RQ3_shared_ARG_group_count_vs_population_density.png
#   results/RQ3_figures/RQ3_shared_ARG_group_count_vs_population_density.pdf
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(scales)
  library(grid)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_dir  <- "results/RQ3_output_tables"
figure_dir <- "results/RQ3_figures"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

master_file <- file.path(input_dir, "RQ3_master_table_FINAL_corrected.csv")

if (!file.exists(master_file)) {
  stop(
    "Input file not found: ", master_file, "\n",
    "Please check that the corrected master table is present in results/RQ3_output_tables."
  )
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
  "city_population_density",
  "shared_ARG_group_count"
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
    population_density = as.numeric(city_population_density),
    shared_ARG_group_count = as.numeric(shared_ARG_group_count)
  ) %>%
  filter(
    !is.na(population_density),
    !is.na(shared_ARG_group_count)
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
    -850, -650, -200, -100, 0, -70,
    -600, -700, -550, -600, -750,
    -250, -350, -550, -120, -150
  ),
  dy = c(
    -1.2, -1.2, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, -1.2, -1.2, 0.8,
    0.8, 0.8, 0.8, -1.2, 0.8
  )
)

plot_df <- plot_df %>%
  left_join(label_offsets, by = "city") %>%
  mutate(
    dx = if_else(is.na(dx), 0, dx),
    dy = if_else(is.na(dy), 0.8, dy),
    label_x = population_density + dx,
    label_y = shared_ARG_group_count + dy
  )

# ------------------------------------------------------------
# 6. Statistics label
# ------------------------------------------------------------

stat_label <- "\u03c1 = -0.432\np = 0.0945\nBH-adjusted p = 0.283"

# ------------------------------------------------------------
# 7. Generate figure
# ------------------------------------------------------------

p_shared_pop_density <- ggplot(
  plot_df,
  aes(x = population_density, y = shared_ARG_group_count)
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
    x = 16300,
    y = 54.8,
    label = stat_label,
    hjust = 1,
    vjust = 1,
    size = 4.8,
    label.size = 0.7,
    fill = "white",
    color = "black"
  ) +
  scale_x_continuous(
    limits = c(500, 16900),
    breaks = c(4000, 8000, 12000, 16000),
    labels = comma
  ) +
  scale_y_continuous(
    limits = c(2.5, 55.5),
    breaks = c(10, 20, 30, 40, 50)
  ) +
  labs(
    x = expression("Population Density (inhabitants km"^{-2}*")"),
    y = "Shared ARG-Group Count"
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
  filename = file.path(figure_dir, "RQ3_shared_ARG_group_count_vs_population_density.png"),
  plot = p_shared_pop_density,
  width = 8.8,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(figure_dir, "RQ3_shared_ARG_group_count_vs_population_density.pdf"),
  plot = p_shared_pop_density,
  width = 8.8,
  height = 6.8
)

message("Saved:")
message(" - ", file.path(figure_dir, "RQ3_shared_ARG_group_count_vs_population_density.png"))
message(" - ", file.path(figure_dir, "RQ3_shared_ARG_group_count_vs_population_density.pdf"))

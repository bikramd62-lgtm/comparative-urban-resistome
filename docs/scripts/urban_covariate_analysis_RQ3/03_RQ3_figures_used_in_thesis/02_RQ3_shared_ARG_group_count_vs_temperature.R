# ============================================================
# 02_RQ3_shared_ARG_group_count_vs_temperature.R
#
# Figure:
# RQ3 shared ARG-group count vs average June temperature
#
# Input:
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
#
# Output:
#   results/RQ3_figures/RQ3_shared_ARG_group_count_vs_temperature.png
#   results/RQ3_figures/RQ3_shared_ARG_group_count_vs_temperature.pdf
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

# ------------------------------------------------------------
# 1.Messages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  Paths
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

required_columns <- c(
  "city",
  "city_avg_june_temperature",
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
# 3. Prepare plotting data
# ------------------------------------------------------------

plot_df <- rq3 %>%
  transmute(
    city = city,
    city_label = str_replace_all(city, "_", " "),
    temperature_june_c = city_avg_june_temperature,
    shared_ARG_group_count = shared_ARG_group_count
  ) %>%
  filter(
    !is.na(temperature_june_c),
    !is.na(shared_ARG_group_count)
  )

# ------------------------------------------------------------
# 4. Manual label placement to match thesis figure
# ------------------------------------------------------------

label_offsets <- tibble(
  city = c(
    "Barcelona", "Berlin", "Bogota", "Hanoi", "Hong_Kong", "Ilorin",
    "Kuala_Lumpur", "Lisbon", "Oslo", "Porto", "Rio_de_Janeiro",
    "Santiago", "Singapore", "Sofia", "Taipei", "Vienna"
  ),
  dx = c(
    -0.90, -0.90, -0.25, 0.00, 0.15, -0.10,
    -0.90, -0.15, -0.80, -0.90, -0.85,
    -0.40, -0.50, 0.10, -0.90, -0.20
  ),
  dy = c(
    -1.00, -1.00, 0.80, 0.90, 0.80, 0.80,
    0.80, 0.80, 0.80, -1.00, 0.80,
    0.80, 0.80, 0.80, -1.00, 0.80
  )
)

plot_df <- plot_df %>%
  left_join(label_offsets, by = "city") %>%
  mutate(
    dx = if_else(is.na(dx), 0, dx),
    dy = if_else(is.na(dy), 0.8, dy),
    label_x = temperature_june_c + dx,
    label_y = shared_ARG_group_count + dy
  )

# ------------------------------------------------------------
# 5. Statistics label
# ------------------------------------------------------------

stat_label <- "\u03c1 = 0.521\np = 0.0385\nBH-adjusted p = 0.231"

# ------------------------------------------------------------
# 6. Generate figure
# ------------------------------------------------------------

p_shared_temperature <- ggplot(
  plot_df,
  aes(x = temperature_june_c, y = shared_ARG_group_count)
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
    x = 7.8,
    y = 53.6,
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
    limits = c(2.5, 55.5),
    breaks = c(10, 20, 30, 40, 50)
  ) +
  labs(
    x = expression("Average June Temperature (" * degree * "C)"),
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
# 7. Save figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "RQ3_shared_ARG_group_count_vs_temperature.png"),
  plot = p_shared_temperature,
  width = 8.8,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(figure_dir, "RQ3_shared_ARG_group_count_vs_temperature.pdf"),
  plot = p_shared_temperature,
  width = 8.8,
  height = 6.8
)

message("Saved: RQ3_shared_ARG_group_count_vs_temperature.png and RQ3_shared_ARG_group_count_vs_temperature.pdf")

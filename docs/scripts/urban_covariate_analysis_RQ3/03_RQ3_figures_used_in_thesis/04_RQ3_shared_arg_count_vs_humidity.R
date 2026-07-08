# ============================================================
# 04_RQ3_shared_arg_count_vs_humidity.R
#
# Figure:
# Shared ARG-group count vs June relative humidity
#
# Input:
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
#
# Output:
#   results/RQ3_figures/RQ3_shared_arg_count_vs_humidity.png
#   results/RQ3_figures/RQ3_shared_arg_count_vs_humidity.pdf
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
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
# 3. Helper function to find required columns robustly
# ------------------------------------------------------------

find_first_existing_column <- function(df, candidates, label) {
  found <- candidates[candidates %in% colnames(df)]
  if (length(found) == 0) {
    stop(
      "Could not find the required column for: ", label, "\n",
      "Checked these candidates:\n",
      paste(candidates, collapse = "\n")
    )
  }
  found[1]
}

city_col <- find_first_existing_column(
  rq3,
  c("city"),
  "city"
)

humidity_col <- find_first_existing_column(
  rq3,
  c(
    "june_relative_humidity_percent",
    "humidity_june_percent",
    "june_humidity_percent",
    "relative_humidity_june_percent"
  ),
  "June relative humidity"
)

shared_count_col <- find_first_existing_column(
  rq3,
  c(
    "shared_ARG_group_count",
    "shared_arg_group_count",
    "shared_ARG_groups",
    "shared_arg_groups",
    "shared_count"
  ),
  "shared ARG-group count"
)

# ------------------------------------------------------------
# 4. Prepare plotting data
# ------------------------------------------------------------

plot_df <- rq3 %>%
  transmute(
    city = .data[[city_col]],
    city_label = str_replace_all(.data[[city_col]], "_", " "),
    humidity_june_percent = .data[[humidity_col]],
    shared_arg_group_count = .data[[shared_count_col]]
  ) %>%
  filter(
    !is.na(humidity_june_percent),
    !is.na(shared_arg_group_count)
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
    -2.2, -1.4, -0.5, -0.2, 0.2, 0.0,
    -1.5, -1.8, -1.4, -1.6, -2.8,
    -0.7, -0.9, -0.1, -1.6, -0.4
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
    label_x = humidity_june_percent + dx,
    label_y = shared_arg_group_count + dy
  )

# ------------------------------------------------------------
# 6. Statistics label
# ------------------------------------------------------------

stat_label <- "\u03c1 = 0.530\np = 0.0347\nBH-adjusted p = 0.231"

# ------------------------------------------------------------
# 7. Generate figure
# ------------------------------------------------------------

p_shared_humidity <- ggplot(
  plot_df,
  aes(x = humidity_june_percent, y = shared_arg_group_count)
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
    y = 54.6,
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
    limits = c(3, 55.5),
    breaks = c(10, 20, 30, 40, 50)
  ) +
  labs(
    x = "June Relative Humidity (%)",
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
  filename = file.path(figure_dir, "RQ3_shared_arg_count_vs_humidity.png"),
  plot = p_shared_humidity,
  width = 8.8,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(figure_dir, "RQ3_shared_arg_count_vs_humidity.pdf"),
  plot = p_shared_humidity,
  width = 8.8,
  height = 6.8
)

message("Saved:")
message(" - ", file.path(figure_dir, "RQ3_shared_arg_count_vs_humidity.png"))
message(" - ", file.path(figure_dir, "RQ3_shared_arg_count_vs_humidity.pdf"))

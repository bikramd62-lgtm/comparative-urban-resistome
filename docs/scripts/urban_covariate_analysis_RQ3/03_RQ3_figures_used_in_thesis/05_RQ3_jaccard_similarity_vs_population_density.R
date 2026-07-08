# ============================================================
# 05_RQ3_jaccard_similarity_vs_population_density.R
#
# Figure:
# Sewage-Transit Jaccard similarity vs population density
#
# Input:
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
#
# Output:
#   results/RQ3_figures/RQ3_jaccard_similarity_vs_population_density.png
#   results/RQ3_figures/RQ3_jaccard_similarity_vs_population_density.pdf
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
# 3. Helper function to locate columns robustly
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

popdens_col <- find_first_existing_column(
  rq3,
  c(
    "population_density",
    "city_population_density",
    "population_density_inhabitants_km2",
    "population_density_inhabitants_km_2",
    "pop_density",
    "city_population_density_km2"
  ),
  "population density"
)

jaccard_col <- find_first_existing_column(
  rq3,
  c(
    "sewage_transit_jaccard_similarity",
    "jaccard_similarity",
    "Jaccard_similarity",
    "sewage_transit_jaccard",
    "jaccard"
  ),
  "Sewage-Transit Jaccard similarity"
)

# ------------------------------------------------------------
# 4. Prepare plotting data
# ------------------------------------------------------------

plot_df <- rq3 %>%
  transmute(
    city = .data[[city_col]],
    city_label = str_replace_all(.data[[city_col]], "_", " "),
    population_density = as.numeric(.data[[popdens_col]]),
    jaccard_similarity = as.numeric(.data[[jaccard_col]])
  ) %>%
  filter(
    !is.na(population_density),
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
    -900,   -450,  -200,  -650,  -950,   -50,
    -650,   -300,  -500,  -650, -600,
      40,     40,   120,    90,   120
  ),
  dy = c(
    -0.006, -0.006, 0.004, -0.006, -0.006, 0.004,
     0.004, -0.006, -0.005, -0.005, 0.004,
     0.004,  0.004, 0.004,  0.004, 0.004
  )
)

plot_df <- plot_df %>%
  left_join(label_offsets, by = "city") %>%
  mutate(
    dx = if_else(is.na(dx), 0, dx),
    dy = if_else(is.na(dy), 0.004, dy),
    label_x = population_density + dx,
    label_y = jaccard_similarity + dy
  )

# ------------------------------------------------------------
# 6. Statistical annotation
# ------------------------------------------------------------

stat_label <- "\u03c1 = -0.559\np = 0.0244\nBH-adjusted p = 0.231"

# ------------------------------------------------------------
# 7. Generate figure
# ------------------------------------------------------------

p_jaccard_popdens <- ggplot(
  plot_df,
  aes(x = population_density, y = jaccard_similarity)
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
    y = 0.279,
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
    limits = c(0.02, 0.285),
    breaks = c(0.05, 0.10, 0.15, 0.20, 0.25)
  ) +
  labs(
    x = expression("Population Density (inhabitants km"^{-2}*")"),
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
  filename = file.path(figure_dir, "RQ3_jaccard_similarity_vs_population_density.png"),
  plot = p_jaccard_popdens,
  width = 8.8,
  height = 6.8,
  dpi = 600
)

ggsave(
  filename = file.path(figure_dir, "RQ3_jaccard_similarity_vs_population_density.pdf"),
  plot = p_jaccard_popdens,
  width = 8.8,
  height = 6.8
)

message("Saved:")
message(" - ", file.path(figure_dir, "RQ3_jaccard_similarity_vs_population_density.png"))
message(" - ", file.path(figure_dir, "RQ3_jaccard_similarity_vs_population_density.pdf"))

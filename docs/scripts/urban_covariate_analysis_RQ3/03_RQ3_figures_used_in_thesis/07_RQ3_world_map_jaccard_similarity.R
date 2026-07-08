# ============================================================
# 07_RQ3_world_map_jaccard_similarity.R
#
# Figure:
# World map of city-level sewage-transit Jaccard similarity
#
# Input:
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
#
# Output:
#   results/RQ3_figures/RQ3_world_map_corrected_jaccard.png
#   results/RQ3_figures/RQ3_world_map_corrected_jaccard.pdf
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(scales)
  library(maps)
  library(viridis)
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
  c("city", "City", "city_name", "City_name"),
  "city"
)

latitude_col <- find_first_existing_column(
  rq3,
  c("latitude", "Latitude", "lat", "city_latitude", "city_lat"),
  "latitude"
)

longitude_col <- find_first_existing_column(
  rq3,
  c("longitude", "Longitude", "lon", "long", "city_longitude", "city_lon"),
  "longitude"
)

jaccard_col <- find_first_existing_column(
  rq3,
  c(
    "jaccard_similarity_sewage_transit_ARG_group",
    "sewage_transit_jaccard_similarity",
    "jaccard_similarity",
    "Jaccard_similarity",
    "jaccard"
  ),
  "sewage-transit Jaccard similarity"
)

# ------------------------------------------------------------
# 4. Prepare map data
# ------------------------------------------------------------

plot_df <- rq3 %>%
  transmute(
    city = .data[[city_col]],
    city_label = str_replace_all(.data[[city_col]], "_", " "),
    longitude = as.numeric(.data[[longitude_col]]),
    latitude = as.numeric(.data[[latitude_col]]),
    jaccard_similarity = as.numeric(.data[[jaccard_col]])
  ) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    !is.na(jaccard_similarity)
  )

# World map polygon data
world_map <- map_data("world")

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
    3.0, 2.0, -4.0, -8.0, -4.0, -7.0,
    -8.0, -7.0, 2.0, -7.0, -9.0,
    -6.0, 3.0, 3.0, 3.0, 3.0
  ),
  dy = c(
    -5.0, 4.0, -5.0, -9.0, -8.0, -7.0,
    -10.0, -6.0, 4.0, 2.0, -6.0,
    -6.0, -6.0, -5.0, 4.0, 4.0
  )
)

plot_df <- plot_df %>%
  left_join(label_offsets, by = "city") %>%
  mutate(
    dx = if_else(is.na(dx), 2, dx),
    dy = if_else(is.na(dy), 2, dy),
    label_x = longitude + dx,
    label_y = latitude + dy
  )

# ------------------------------------------------------------
# 6. Generate world map
# ------------------------------------------------------------

p_world_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey75",
    linewidth = 0.25
  ) +
  geom_point(
    data = plot_df,
    aes(
      x = longitude,
      y = latitude,
      color = jaccard_similarity
    ),
    size = 2.9
  ) +
  geom_text(
    data = plot_df,
    aes(
      x = label_x,
      y = label_y,
      label = city_label
    ),
    size = 3.4,
    color = "black"
  ) +
  scale_color_viridis(
    option = "plasma",
    direction = -1,
    limits = c(0.04, 0.26),
    breaks = c(0.05, 0.10, 0.15, 0.20, 0.25),
    name = "Sewage-Transit\nJaccard similarity"
  ) +
  scale_x_continuous(
    limits = c(-200, 200),
    breaks = c(-200, -100, 0, 100, 200),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-85, 85),
    breaks = c(-50, 0, 50),
    expand = c(0, 0)
  ) +
  coord_quickmap(
    xlim = c(-200, 200),
    ylim = c(-85, 85)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "grey35", size = 12),
    axis.title = element_blank(),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.key.height = unit(1.4, "cm"),
    legend.key.width = unit(0.45, "cm"),
    plot.margin = margin(10, 20, 10, 10)
  )

# ------------------------------------------------------------
# 7. Save figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "RQ3_world_map_corrected_jaccard.png"),
  plot = p_world_map,
  width = 12.5,
  height = 7.2,
  dpi = 600
)

ggsave(
  filename = file.path(figure_dir, "RQ3_world_map_corrected_jaccard.pdf"),
  plot = p_world_map,
  width = 12.5,
  height = 7.2
)

message("Saved:")
message(" - ", file.path(figure_dir, "RQ3_world_map_corrected_jaccard.png"))
message(" - ", file.path(figure_dir, "RQ3_world_map_corrected_jaccard.pdf"))

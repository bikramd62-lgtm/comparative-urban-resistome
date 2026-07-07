# ============================================================
# fig_rq2_group_specificity_heatmap.R
#
# Purpose:
# Generate the RQ2 ARG Group focused specificity heatmap used
# in the thesis.
#
# The heatmap shows binary detection of ARG Groups selected from
# the strongest descriptive prevalence differences between
# transit and sewage.
#
# This figure is descriptive. It does not show formal
# statistical significance.
#
# Input:
# results/rq1/matrices/rq1_group_combined_binary_matrix.tsv
# results/rq2/prevalence_difference/rq2_group_prevalence_comparison.tsv
#
# Output:
# results/rq2/figures/rq2_group_specificity_heatmap.png
# results/rq2/figures/rq2_group_specificity_heatmap.pdf
# results/rq2/figures/rq2_group_specificity_heatmap_plotdata.tsv
# results/rq2/figures/rq2_group_specificity_heatmap_selected_features.tsv
# ============================================================

# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
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
library(tidyr)
library(stringr)

# ------------------------------------------------------------
# 2. User-adjustable settings
# ------------------------------------------------------------

matrix_file <- file.path(
  "results",
  "rq1",
  "matrices",
  "rq1_group_combined_binary_matrix.tsv"
)

prevalence_file <- file.path(
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
# 4. Check input files
# ------------------------------------------------------------

required_input_files <- c(
  matrix_file,
  prevalence_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "The following required input files are missing:\n",
      paste(missing_input_files, collapse = "\n"),
      "\nRun the RQ1 matrix script and 02_rq2_prevalence_difference_analysis.R first."
    )
  )
}

# ------------------------------------------------------------
# 5. Read ARG Group binary matrix
# ------------------------------------------------------------

raw_matrix <- read.delim(
  file = matrix_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (ncol(raw_matrix) < 2) {
  stop(
    "The binary matrix must contain one profile column and at least one ARG Group column."
  )
}

profile_column <- colnames(raw_matrix)[1]

binary_matrix <- raw_matrix[
  ,
  -1,
  drop = FALSE
]

binary_matrix[] <- lapply(
  binary_matrix,
  function(x) {
    suppressWarnings(
      as.numeric(x)
    )
  }
)

binary_matrix <- as.matrix(
  binary_matrix
)

rownames(binary_matrix) <- raw_matrix[[profile_column]]

if (any(is.na(binary_matrix))) {
  stop(
    "The ARG Group binary matrix contains missing or non-numeric values."
  )
}

binary_matrix <- ifelse(
  binary_matrix > 0,
  1,
  0
)

storage.mode(binary_matrix) <- "numeric"

# ------------------------------------------------------------
# 6. Build profile metadata
# ------------------------------------------------------------

profile_metadata <- data.frame(
  Profile = rownames(binary_matrix),
  stringsAsFactors = FALSE
)

profile_metadata$Environment <- sub(
  "^.*_(Sewage|Transit)$",
  "\\1",
  profile_metadata$Profile
)

profile_metadata$City <- sub(
  "_(Sewage|Transit)$",
  "",
  profile_metadata$Profile
)

if (!all(profile_metadata$Environment %in% c("Sewage", "Transit"))) {
  stop(
    "Environment could not be extracted from some profile names. Expected names such as Berlin_Sewage and Berlin_Transit."
  )
}

profile_metadata$Environment <- factor(
  profile_metadata$Environment,
  levels = c(
    "Sewage",
    "Transit"
  )
)

profile_metadata <- profile_metadata %>%
  arrange(
    City,
    Environment
  )

binary_matrix <- binary_matrix[
  profile_metadata$Profile,
  ,
  drop = FALSE
]

profile_metadata <- profile_metadata %>%
  mutate(
    Profile_label = str_replace(
      Profile,
      "_Sewage$",
      "_S"
    ),
    Profile_label = str_replace(
      Profile_label,
      "_Transit$",
      "_T"
    )
  )

# ------------------------------------------------------------
# 7. Read RQ2 prevalence-difference table
# ------------------------------------------------------------

prevalence_table <- read.delim(
  file = prevalence_file,
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
      "The following required columns are missing from the prevalence-difference file:\n",
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
# 8. Prepare ARG Group labels
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
# 9. Select strongest descriptive prevalence differences
# ------------------------------------------------------------

higher_in_transit <- prevalence_table %>%
  filter(
    Prevalence_direction == "Higher_in_transit"
  ) %>%
  arrange(
    desc(Prevalence_difference),
    desc(Transit_prevalence),
    Feature_label,
    Feature
  ) %>%
  slice_head(
    n = top_n_per_direction
  )

higher_in_sewage <- prevalence_table %>%
  filter(
    Prevalence_direction == "Higher_in_sewage"
  ) %>%
  arrange(
    desc(Prevalence_difference),
    desc(Sewage_prevalence),
    Feature_label,
    Feature
  ) %>%
  slice_head(
    n = top_n_per_direction
  )

selected_features <- bind_rows(
  higher_in_transit,
  higher_in_sewage
) %>%
  filter(
    Feature %in% colnames(binary_matrix)
  ) %>%
  arrange(
    desc(Prevalence_difference),
    Feature_label,
    Feature
  )

if (nrow(selected_features) == 0) {
  stop(
    "No ARG Groups were selected for plotting. Check the prevalence-difference table and binary matrix."
  )
}

if (any(duplicated(selected_features$Feature_label))) {
  selected_features$Feature_label <- make.unique(
    selected_features$Feature_label,
    sep = "_"
  )
}

selected_feature_order <- selected_features$Feature

selected_feature_label_order <- selected_features$Feature_label

names(selected_feature_label_order) <- selected_features$Feature

# ------------------------------------------------------------
# 10. Prepare heatmap data
# ------------------------------------------------------------

heatmap_matrix <- binary_matrix[
  ,
  selected_feature_order,
  drop = FALSE
]

heatmap_data <- as.data.frame(
  heatmap_matrix,
  check.names = FALSE
)

heatmap_data$Profile <- rownames(
  heatmap_matrix
)

heatmap_data <- heatmap_data %>%
  pivot_longer(
    cols = -Profile,
    names_to = "Feature",
    values_to = "Presence"
  ) %>%
  left_join(
    profile_metadata,
    by = "Profile"
  ) %>%
  left_join(
    selected_features %>%
      select(
        Feature,
        Feature_label,
        Sewage_prevalence,
        Transit_prevalence,
        Prevalence_difference,
        Absolute_prevalence_difference,
        Prevalence_direction
      ),
    by = "Feature"
  ) %>%
  mutate(
    Presence = ifelse(
      Presence > 0,
      1,
      0
    ),
    Binary_detection = ifelse(
      Presence == 1,
      "Present = 1",
      "Absent = 0"
    )
  )

# ------------------------------------------------------------
# 11. Factor ordering for plotting
# ------------------------------------------------------------

heatmap_data$Profile_label <- factor(
  heatmap_data$Profile_label,
  levels = profile_metadata$Profile_label
)

heatmap_data$Feature <- factor(
  heatmap_data$Feature,
  levels = rev(
    selected_feature_order
  )
)

heatmap_data$Feature_label <- factor(
  heatmap_data$Feature_label,
  levels = rev(
    selected_feature_label_order
  )
)

heatmap_data$Binary_detection <- factor(
  heatmap_data$Binary_detection,
  levels = c(
    "Absent = 0",
    "Present = 1"
  )
)

# ------------------------------------------------------------
# 12. Define colours
# ------------------------------------------------------------

tile_colours <- c(
  "Absent = 0" = "#D6D6D6",
  "Present = 1" = "#08306B"
)

# ------------------------------------------------------------
# 13. Generate final thesis figure
# ------------------------------------------------------------

p <- ggplot(
  heatmap_data,
  aes(
    x = Profile_label,
    y = Feature_label,
    fill = Binary_detection
  )
) +

  geom_tile(
    colour = "#C8C8C8",
    linewidth = 0.22
  ) +

  scale_fill_manual(
    values = tile_colours,
    name = "Binary detection",
    drop = FALSE
  ) +

  labs(
    x = "City-Environment Profile (S = Sewage, T = Transit)",
    y = "ARG Group"
  ) +

  theme_minimal(
    base_size = 16
  ) +

  theme(
    panel.grid = element_blank(),

    panel.border = element_blank(),

    axis.title.x = element_text(
      size = 15,
      colour = "black",
      margin = margin(
        t = 12
      )
    ),

    axis.title.y = element_text(
      size = 16,
      colour = "black",
      margin = margin(
        r = 12
      )
    ),

    axis.text.x = element_text(
      size = 9,
      colour = "grey25",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),

    axis.text.y = element_text(
      size = 9.5,
      colour = "grey25",
      lineheight = 0.95
    ),

    axis.ticks = element_blank(),

    legend.position = "right",

    legend.title = element_text(
      size = 13,
      colour = "black"
    ),

    legend.text = element_text(
      size = 12,
      colour = "black"
    ),

    legend.key.size = unit(
      0.65,
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
      t = 10,
      r = 18,
      b = 14,
      l = 14
    )
  )

# ------------------------------------------------------------
# 14. Display figure
# ------------------------------------------------------------

print(p)

# ------------------------------------------------------------
# 15. Save figure and plot data
# ------------------------------------------------------------

png_output <- file.path(
  output_dir,
  "rq2_group_specificity_heatmap.png"
)

pdf_output <- file.path(
  output_dir,
  "rq2_group_specificity_heatmap.pdf"
)

plot_data_output <- file.path(
  output_dir,
  "rq2_group_specificity_heatmap_plotdata.tsv"
)

selected_features_output <- file.path(
  output_dir,
  "rq2_group_specificity_heatmap_selected_features.tsv"
)

ggsave(
  filename = png_output,
  plot = p,
  width = 14.8,
  height = 8.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_output,
  plot = p,
  width = 14.8,
  height = 8.8,
  units = "in",
  bg = "white"
)

write.table(
  heatmap_data,
  file = plot_data_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  selected_features,
  file = selected_features_output,
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

message(
  "\nSelected features saved:"
)

message(
  selected_features_output
)

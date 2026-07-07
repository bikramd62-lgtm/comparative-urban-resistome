# ============================================================
# fig_rq1_group_presence_absence_heatmap_top40.R
#
# Purpose:
# Generate the RQ1 ARG Group presence/absence heatmap used in
# the thesis.
#
# Input:
# results/rq1/matrices/rq1_group_combined_binary_matrix.tsv
#
# Output:
# results/rq1/figures/rq1_group_presence_absence_heatmap_top40_navy_darkergrey_border.png
# results/rq1/figures/rq1_group_presence_absence_heatmap_top40_navy_darkergrey_border.pdf
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
# 2. User-adjustable paths and settings
# ------------------------------------------------------------

input_dir <- file.path(
  "results",
  "rq1",
  "matrices"
)

output_dir <- file.path(
  "results",
  "rq1",
  "figures"
)

input_file <- file.path(
  input_dir,
  "rq1_group_combined_binary_matrix.tsv"
)

top_n_features <- 40

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      input_file,
      "\nRun 01_build_rq1_binary_matrices.R before generating this figure."
    )
  )
}

# ------------------------------------------------------------
# 3. Read combined ARG Group binary matrix
# ------------------------------------------------------------

raw_df <- read.delim(
  file = input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (ncol(raw_df) < 2) {
  stop(
    "Input matrix must contain one profile column and at least one ARG Group column."
  )
}

profile_ids <- trimws(
  as.character(raw_df[[1]])
)

binary_matrix <- as.matrix(
  raw_df[
    ,
    -1,
    drop = FALSE
  ]
)

storage.mode(binary_matrix) <- "numeric"

rownames(binary_matrix) <- profile_ids

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

if (any(duplicated(rownames(binary_matrix)))) {
  stop(
    "Duplicated profile names detected in the ARG Group matrix."
  )
}

if (any(duplicated(colnames(binary_matrix)))) {
  stop(
    "Duplicated ARG Group names detected in the ARG Group matrix."
  )
}

# ------------------------------------------------------------
# 4. Build metadata and order profiles by city and environment
# ------------------------------------------------------------

metadata <- data.frame(
  Profile = rownames(binary_matrix),
  stringsAsFactors = FALSE
)

metadata$Environment <- sub(
  "^.*_(Sewage|Transit)$",
  "\\1",
  metadata$Profile
)

metadata$City <- sub(
  "_(Sewage|Transit)$",
  "",
  metadata$Profile
)

if (!all(metadata$Environment %in% c("Sewage", "Transit"))) {
  stop(
    paste0(
      "Environment could not be extracted from some profile names:\n",
      paste(metadata$Profile, collapse = "\n")
    )
  )
}

metadata$Environment <- factor(
  metadata$Environment,
  levels = c(
    "Sewage",
    "Transit"
  )
)

metadata <- metadata[
  order(
    metadata$City,
    metadata$Environment
  ),
  ,
  drop = FALSE
]

binary_matrix <- binary_matrix[
  metadata$Profile,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 5. Select top ARG Groups by total detection frequency
# ------------------------------------------------------------

feature_detection_count <- colSums(
  binary_matrix
)

if (length(feature_detection_count) < top_n_features) {
  warning(
    paste0(
      "Only ",
      length(feature_detection_count),
      " ARG Groups are available. Plotting all features."
    )
  )

  top_n_features <- length(
    feature_detection_count
  )
}

feature_order_table <- data.frame(
  Feature = names(feature_detection_count),
  Detection_count = as.numeric(feature_detection_count),
  stringsAsFactors = FALSE
)

feature_order_table <- feature_order_table[
  order(
    -feature_order_table$Detection_count,
    feature_order_table$Feature
  ),
  ,
  drop = FALSE
]

top_features <- feature_order_table$Feature[
  seq_len(top_n_features)
]

# For the heatmap, display less prevalent of the selected top features
# at the top and the most recurrent ARG Groups at the bottom.
top_feature_counts <- feature_detection_count[
  top_features
]

top_features_ordered_for_plot <- names(
  sort(
    top_feature_counts,
    decreasing = FALSE
  )
)

# ------------------------------------------------------------
# 6. Helper functions for display labels
# ------------------------------------------------------------

make_arg_group_label <- function(feature_name) {

  label <- feature_name

  # If MEGARes-style labels contain comma-separated hierarchy,
  # use the final element as the short ARG Group label.
  label <- sub(
    "^.*,",
    "",
    label
  )

  label <- trimws(
    label
  )

  return(label)
}

make_profile_label <- function(profile_name) {

  label <- profile_name

  label <- sub(
    "_Sewage$",
    "_S",
    label
  )

  label <- sub(
    "_Transit$",
    "_T",
    label
  )

  return(label)
}

feature_label_table <- data.frame(
  Feature = top_features_ordered_for_plot,
  ARG_Group_label = vapply(
    top_features_ordered_for_plot,
    make_arg_group_label,
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

if (any(duplicated(feature_label_table$ARG_Group_label))) {
  feature_label_table$ARG_Group_label <- make.unique(
    feature_label_table$ARG_Group_label,
    sep = "_"
  )
}

profile_label_table <- data.frame(
  Profile = metadata$Profile,
  Profile_label = vapply(
    metadata$Profile,
    make_profile_label,
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 7. Convert matrix to long format for ggplot heatmap
# ------------------------------------------------------------

heatmap_matrix <- binary_matrix[
  ,
  top_features_ordered_for_plot,
  drop = FALSE
]

heatmap_df <- expand.grid(
  Profile = rownames(heatmap_matrix),
  Feature = colnames(heatmap_matrix),
  stringsAsFactors = FALSE
)

heatmap_df$Presence <- mapply(
  function(profile, feature) {
    heatmap_matrix[
      profile,
      feature
    ]
  },
  heatmap_df$Profile,
  heatmap_df$Feature
)

heatmap_df <- merge(
  heatmap_df,
  profile_label_table,
  by = "Profile",
  all.x = TRUE,
  sort = FALSE
)

heatmap_df <- merge(
  heatmap_df,
  feature_label_table,
  by = "Feature",
  all.x = TRUE,
  sort = FALSE
)

heatmap_df$Profile_label <- factor(
  heatmap_df$Profile_label,
  levels = profile_label_table$Profile_label
)

heatmap_df$ARG_Group_label <- factor(
  heatmap_df$ARG_Group_label,
  levels = feature_label_table$ARG_Group_label
)

heatmap_df$Binary_detection <- factor(
  heatmap_df$Presence,
  levels = c(
    0,
    1
  ),
  labels = c(
    "Absent = 0",
    "Present = 1"
  )
)

# ------------------------------------------------------------
# 8. Generate final thesis figure
# ------------------------------------------------------------

detection_colours <- c(
  "Absent = 0" = "#D6D6D6",
  "Present = 1" = "#08306B"
)

p <- ggplot(
  heatmap_df,
  aes(
    x = Profile_label,
    y = ARG_Group_label,
    fill = Binary_detection
  )
) +

  geom_tile(
    colour = "#C8C8C8",
    linewidth = 0.25
  ) +

  scale_fill_manual(
    values = detection_colours,
    name = "Binary detection",
    drop = FALSE
  ) +

  labs(
    x = "City-Environment Profile (S = Sewage, T = Transit)",
    y = "ARG Group"
  ) +

  theme_minimal(
    base_size = 13
  ) +

  theme(
    panel.grid = element_blank(),

    panel.border = element_blank(),

    axis.title.x = element_text(
      size = 13,
      colour = "black",
      margin = margin(
        t = 12
      )
    ),

    axis.title.y = element_text(
      size = 13,
      colour = "black",
      margin = margin(
        r = 10
      )
    ),

    axis.text.x = element_text(
      size = 8.5,
      colour = "grey25",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),

    axis.text.y = element_text(
      size = 8.5,
      colour = "grey25"
    ),

    axis.ticks = element_blank(),

    legend.position = "right",

    legend.title = element_text(
      size = 12,
      colour = "black"
    ),

    legend.text = element_text(
      size = 10,
      colour = "black"
    ),

    legend.key.size = unit(
      0.55,
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
      t = 8,
      r = 18,
      b = 8,
      l = 8
    )
  )

# ------------------------------------------------------------
# 9. Display figure
# ------------------------------------------------------------

print(p)

# ------------------------------------------------------------
# 10. Save figure and plot data
# ------------------------------------------------------------

png_output <- file.path(
  output_dir,
  "rq1_group_presence_absence_heatmap_top40_navy_darkergrey_border.png"
)

pdf_output <- file.path(
  output_dir,
  "rq1_group_presence_absence_heatmap_top40_navy_darkergrey_border.pdf"
)

plot_data_output <- file.path(
  output_dir,
  "rq1_group_presence_absence_heatmap_top40_plotdata.tsv"
)

selected_features_output <- file.path(
  output_dir,
  "rq1_group_presence_absence_heatmap_top40_selected_features.tsv"
)

ggsave(
  filename = png_output,
  plot = p,
  width = 13.8,
  height = 9.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_output,
  plot = p,
  width = 13.8,
  height = 9.2,
  units = "in",
  bg = "white"
)

write.table(
  heatmap_df,
  file = plot_data_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  feature_order_table[
    feature_order_table$Feature %in% top_features_ordered_for_plot,
    ,
    drop = FALSE
  ],
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
  "\nSelected top features saved:"
)

message(
  selected_features_output
)

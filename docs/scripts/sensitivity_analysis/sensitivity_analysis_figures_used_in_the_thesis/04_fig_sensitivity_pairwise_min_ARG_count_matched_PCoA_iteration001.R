# ============================================================
# fig_sensitivity_pairwise_min_ARG_count_matched_PCoA_iteration001.R
#
# Purpose:
#   Generate the PCoA figure for the pairwise total ARG-count
#   matched sensitivity analysis
#
# Figure:
#   Pairwise-minimum ARG-count-matched PCoA
#
# Input:
#   processed/sensitivit_analysis/arg_count matching/
#   └── pairwise_min_ARG_count_matched_binary_iteration_001.csv
#
# Output:
#   results/sensitivity_analysis/figures/
#   ├── pairwise_min_ARG_count_matched_PCoA_iteration_001.png
#   └── pairwise_min_ARG_count_matched_PCoA_iteration_001.pdf
#
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "vegan",
  "ggplot2",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this script."
    )
  )
}

library(vegan)
library(ggplot2)
library(ggrepel)


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_root <- Sys.getenv(
  "PROJECT_ROOT",
  unset = getwd()
)

input_dir <- file.path(
  project_root,
  "processed",
  "sensitivit_analysis",
  "arg_count matching"
)

output_dir <- file.path(
  project_root,
  "results",
  "sensitivity_analysis",
  "figures"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Input file
# ------------------------------------------------------------

binary_file <- file.path(
  input_dir,
  "pairwise_min_ARG_count_matched_binary_iteration_001.csv"
)

if (!file.exists(binary_file)) {
  stop(
    paste0(
      "Iteration-001 binary matrix was not found:\n",
      binary_file,
      "\n\nCheck that the folder and filename match exactly."
    )
  )
}


# ------------------------------------------------------------
# 4. Read binary matrix
# ------------------------------------------------------------

binary_table <- read.csv(
  binary_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata_columns <- c(
  "SampleID",
  "City",
  "Environment"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(binary_table)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    paste0(
      "The binary matrix is missing required metadata columns:\n",
      paste(missing_metadata_columns, collapse = ", ")
    )
  )
}

metadata <- binary_table[
  ,
  required_metadata_columns,
  drop = FALSE
]

metadata$SampleID <- trimws(as.character(metadata$SampleID))
metadata$City <- trimws(as.character(metadata$City))
metadata$Environment <- trimws(as.character(metadata$Environment))

metadata$Environment[tolower(metadata$Environment) == "sewage"] <- "Sewage"
metadata$Environment[tolower(metadata$Environment) == "transit"] <- "Transit"

metadata$Environment <- factor(
  metadata$Environment,
  levels = c("Sewage", "Transit")
)

if (anyNA(metadata$Environment)) {
  stop("Environment column must contain only Sewage and Transit.")
}

if (anyDuplicated(metadata$SampleID) > 0) {
  stop("Duplicate SampleID values were found.")
}


# ------------------------------------------------------------
# 5. Extract ARG-group binary feature matrix
# ------------------------------------------------------------

feature_columns <- setdiff(
  colnames(binary_table),
  required_metadata_columns
)

if (length(feature_columns) == 0) {
  stop("No ARG-group feature columns were found in the binary matrix.")
}

binary_matrix <- binary_table[
  ,
  feature_columns,
  drop = FALSE
]

binary_matrix[] <- lapply(
  binary_matrix,
  function(x) as.numeric(as.character(x))
)

binary_matrix <- as.matrix(binary_matrix)
rownames(binary_matrix) <- metadata$SampleID
storage.mode(binary_matrix) <- "numeric"

if (anyNA(binary_matrix)) {
  stop("The binary matrix contains missing or non-numeric values.")
}

observed_values <- sort(unique(as.vector(binary_matrix)))

if (!all(observed_values %in% c(0, 1))) {
  stop(
    "The matrix contains values other than 0 and 1. Observed values: ",
    paste(observed_values, collapse = ", ")
  )
}

if (any(rowSums(binary_matrix) == 0)) {
  empty_profiles <- rownames(binary_matrix)[rowSums(binary_matrix) == 0]

  stop(
    paste0(
      "One or more profiles contain no detected ARG groups:\n",
      paste(empty_profiles, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 6. Calculate ARG-group richness
# ------------------------------------------------------------

metadata$ARG_group_richness <- rowSums(binary_matrix > 0)


# ------------------------------------------------------------
# 7. Binary Jaccard distance and PCoA
# ------------------------------------------------------------

jaccard_distance <- vegan::vegdist(
  binary_matrix,
  method = "jaccard",
  binary = TRUE
)

pcoa_result <- stats::cmdscale(
  jaccard_distance,
  k = 2,
  eig = TRUE
)

positive_eigenvalue_sum <- sum(
  pcoa_result$eig[pcoa_result$eig > 0]
)

pcoa1_percent <- pcoa_result$eig[1] / positive_eigenvalue_sum * 100
pcoa2_percent <- pcoa_result$eig[2] / positive_eigenvalue_sum * 100

plot_data <- data.frame(
  SampleID = metadata$SampleID,
  City = metadata$City,
  Environment = metadata$Environment,
  ARG_group_richness = metadata$ARG_group_richness,
  PCoA1 = pcoa_result$points[, 1],
  PCoA2 = pcoa_result$points[, 2],
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 8. Keep PCoA orientation consistent with thesis figure
# ------------------------------------------------------------

# PCoA axis signs are arbitrary.
# This keeps sewage mostly on the left and transit mostly on the right.

if (
  mean(plot_data$PCoA1[plot_data$Environment == "Sewage"]) >
    mean(plot_data$PCoA1[plot_data$Environment == "Transit"])
) {
  plot_data$PCoA1 <- -plot_data$PCoA1
}

city_standardised <- gsub("_", " ", plot_data$City)

sewage_upper_reference <- plot_data$PCoA2[
  plot_data$Environment == "Sewage" &
    city_standardised %in% c(
      "Bogota",
      "Hanoi",
      "Hong Kong",
      "Singapore",
      "Taipei",
      "Kuala Lumpur"
    )
]

sewage_lower_reference <- plot_data$PCoA2[
  plot_data$Environment == "Sewage" &
    city_standardised %in% c(
      "Barcelona",
      "Berlin",
      "Lisbon",
      "Oslo",
      "Porto",
      "Santiago",
      "Sofia",
      "Vienna"
    )
]

if (
  length(sewage_upper_reference) > 0 &&
    length(sewage_lower_reference) > 0 &&
    mean(sewage_upper_reference) < mean(sewage_lower_reference)
) {
  plot_data$PCoA2 <- -plot_data$PCoA2
}

# Use clean city labels in the plot.

plot_data$City_label <- gsub("_", " ", plot_data$City)


# ------------------------------------------------------------
# 9. Save PCoA coordinates
# ------------------------------------------------------------

coordinate_output <- data.frame(
  plot_data,
  PCoA1_percent = pcoa1_percent,
  PCoA2_percent = pcoa2_percent,
  stringsAsFactors = FALSE
)

write.csv(
  coordinate_output,
  file = file.path(
    output_dir,
    "pairwise_min_ARG_count_matched_PCoA_iteration_001_coordinates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Generate figure
# ------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = PCoA1,
    y = PCoA2
  )
) +
  geom_point(
    aes(
      colour = Environment,
      size = ARG_group_richness
    ),
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    aes(
      label = City_label,
      colour = Environment
    ),
    size = 4.2,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.size = 0.25,
    segment.alpha = 0.6,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "Sewage" = "#F8766D",
      "Transit" = "#00BFC4"
    )
  ) +
  scale_size_continuous(
    name = "ARG group richness",
    range = c(3.5, 10),
    breaks = c(50, 100, 150, 200)
  ) +
  labs(
    x = paste0("PCoA1 (", sprintf("%.1f", pcoa1_percent), "%)"),
    y = paste0("PCoA2 (", sprintf("%.1f", pcoa2_percent), "%)"),
    colour = "Environment"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 5)
    ),
    size = guide_legend()
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.7
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "grey30",
      fill = NA,
      linewidth = 0.9
    ),
    axis.title = element_text(
      size = 20,
      colour = "black"
    ),
    axis.text = element_text(
      size = 15,
      colour = "black"
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 17,
      colour = "black"
    ),
    legend.text = element_text(
      size = 14,
      colour = "black"
    ),
    legend.key = element_blank(),
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )


# ------------------------------------------------------------
# 11. Save final figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "pairwise_min_ARG_count_matched_PCoA_iteration_001.png"
  ),
  plot = p,
  width = 9.6,
  height = 7.2,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "pairwise_min_ARG_count_matched_PCoA_iteration_001.pdf"
  ),
  plot = p,
  width = 9.6,
  height = 7.2
)


# ------------------------------------------------------------
# 12. Save session information
# ------------------------------------------------------------

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("Input file: ", binary_file),
  paste0("R version: ", R.version.string),
  paste0("vegan version: ", as.character(packageVersion("vegan"))),
  paste0("ggplot2 version: ", as.character(packageVersion("ggplot2"))),
  paste0("ggrepel version: ", as.character(packageVersion("ggrepel"))),
  "Distance measure: binary Jaccard dissimilarity",
  "PCoA function: stats::cmdscale(k = 2, eig = TRUE)",
  "PCoA percentages calculated from the sum of positive eigenvalues",
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_text,
  con = file.path(
    output_dir,
    "pairwise_min_ARG_count_matched_PCoA_iteration_001_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 13. Print completion message
# ------------------------------------------------------------

message("")
message("====================================================")
message("Figure generated successfully")
message("====================================================")
message("Input file:")
message(binary_file)
message("")
message("Profiles: ", nrow(binary_matrix))
message("ARG groups: ", ncol(binary_matrix))
message("PCoA1 variation: ", sprintf("%.1f", pcoa1_percent), "%")
message("PCoA2 variation: ", sprintf("%.1f", pcoa2_percent), "%")
message("")
message("Output folder:")
message(output_dir)

# ============================================================
# fig_sensitivity_pcoa_sewage1M_transit5M.R
#
# Purpose:
#   Generate the PCoA figure for the reduced-depth sensitivity
#   analysis: sewage 1M paired reads versus transit 5M paired reads.
#
# Input:
#   data/processed/sensitivity_analysis/reduced_depth/
#   ├── combined_arg_group_binary_matrix_sewage1M_transit5M_t80_prev10.tsv
#   └── metadata_sewage1M_transit5M_arg_group.tsv
#
# Output:
#   results/sensitivity_analysis/figures/
#   ├── pcoa_arg_group_sewage1M_transit5M_t80_prev10_clean.png
#   └── pcoa_arg_group_sewage1M_transit5M_t80_prev10_clean.pdf
#
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "vegan",
  "ggplot2"
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


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_root <- Sys.getenv(
  "PROJECT_ROOT",
  unset = getwd()
)

input_dir <- file.path(
  project_root,
  "data",
  "processed",
  "sensitivity_analysis",
  "reduced_depth"
)

output_dir <- file.path(
  project_root,
  "results",
  "sensitivity_analysis",
  "figures"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

matrix_file <- file.path(
  input_dir,
  "combined_arg_group_binary_matrix_sewage1M_transit5M_t80_prev10.tsv"
)

metadata_file <- file.path(
  input_dir,
  "metadata_sewage1M_transit5M_arg_group.tsv"
)

if (!file.exists(matrix_file)) {
  stop("Binary matrix file not found:\n", matrix_file)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file not found:\n", metadata_file)
}


# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

standardise_metadata <- function(metadata) {

  lower_names <- tolower(colnames(metadata))

  profile_col <- which(lower_names %in% c("profile_id", "sampleid", "sample_id"))
  city_col <- which(lower_names == "city")
  environment_col <- which(lower_names == "environment")

  if (length(profile_col) != 1) {
    stop("Metadata must contain one profile ID column: profile_id, SampleID, or sample_id.")
  }

  if (length(city_col) != 1) {
    stop("Metadata must contain one city column.")
  }

  if (length(environment_col) != 1) {
    stop("Metadata must contain one environment column.")
  }

  colnames(metadata)[profile_col] <- "profile_id"
  colnames(metadata)[city_col] <- "city"
  colnames(metadata)[environment_col] <- "environment"

  metadata$profile_id <- trimws(as.character(metadata$profile_id))
  metadata$city <- trimws(as.character(metadata$city))
  metadata$environment <- trimws(as.character(metadata$environment))

  metadata$environment[tolower(metadata$environment) == "sewage"] <- "Sewage"
  metadata$environment[tolower(metadata$environment) == "transit"] <- "Transit"

  metadata$environment <- factor(
    metadata$environment,
    levels = c("Sewage", "Transit")
  )

  metadata
}


read_binary_matrix <- function(file_path) {

  x <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  lower_names <- tolower(colnames(x))
  profile_col <- which(lower_names %in% c("profile_id", "sampleid", "sample_id"))

  if (length(profile_col) != 1) {
    stop("Binary matrix must contain one profile ID column: profile_id, SampleID, or sample_id.")
  }

  colnames(x)[profile_col] <- "profile_id"

  profile_ids <- trimws(as.character(x$profile_id))
  x$profile_id <- NULL

  x[] <- lapply(
    x,
    function(column) as.numeric(as.character(column))
  )

  binary_matrix <- as.matrix(x)
  rownames(binary_matrix) <- profile_ids
  storage.mode(binary_matrix) <- "numeric"

  if (anyNA(binary_matrix)) {
    stop("Binary matrix contains missing or non-numeric values.")
  }

  observed_values <- sort(unique(as.vector(binary_matrix)))

  if (!all(observed_values %in% c(0, 1))) {
    stop(
      "Binary matrix contains values other than 0 and 1. Observed values: ",
      paste(observed_values, collapse = ", ")
    )
  }

  if (any(rowSums(binary_matrix) == 0)) {
    stop("One or more profiles contain no detected ARG groups.")
  }

  binary_matrix
}


# ------------------------------------------------------------
# 5. Read and validate data
# ------------------------------------------------------------

binary_matrix <- read_binary_matrix(matrix_file)

metadata <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

metadata <- standardise_metadata(metadata)

missing_from_metadata <- setdiff(rownames(binary_matrix), metadata$profile_id)
missing_from_matrix <- setdiff(metadata$profile_id, rownames(binary_matrix))

if (length(missing_from_metadata) > 0 || length(missing_from_matrix) > 0) {
  stop(
    "Profile IDs do not match between matrix and metadata.\n",
    "Missing from metadata: ", paste(missing_from_metadata, collapse = ", "), "\n",
    "Missing from matrix: ", paste(missing_from_matrix, collapse = ", ")
  )
}

metadata <- metadata[
  match(rownames(binary_matrix), metadata$profile_id),
  ,
  drop = FALSE
]

if (!identical(metadata$profile_id, rownames(binary_matrix))) {
  stop("Metadata could not be aligned to matrix rows.")
}


# ------------------------------------------------------------
# 6. Calculate binary Jaccard distance and PCoA
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

positive_eigenvalue_sum <- sum(pcoa_result$eig[pcoa_result$eig > 0])

pcoa1_percent <- pcoa_result$eig[1] / positive_eigenvalue_sum * 100
pcoa2_percent <- pcoa_result$eig[2] / positive_eigenvalue_sum * 100

plot_data <- data.frame(
  profile_id = metadata$profile_id,
  city = metadata$city,
  Environment = metadata$environment,
  PCoA1 = pcoa_result$points[, 1],
  PCoA2 = pcoa_result$points[, 2],
  stringsAsFactors = FALSE
)

# PCoA axis signs are arbitrary.
# sewage on the left and transit on the right.
if (
  mean(plot_data$PCoA1[plot_data$Environment == "Sewage"]) >
    mean(plot_data$PCoA1[plot_data$Environment == "Transit"])
) {
  plot_data$PCoA1 <- -plot_data$PCoA1
}


# ------------------------------------------------------------
# 7. Save PCoA coordinate table
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
    "pcoa_arg_group_sewage1M_transit5M_t80_prev10_coordinates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Generate figure
# ------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Environment,
    shape = Environment
  )
) +
  geom_point(
    size = 4.2,
    alpha = 0.9
  ) +
  scale_colour_manual(
    values = c(
      "Sewage" = "#F8766D",
      "Transit" = "#00BFC4"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Sewage" = 16,
      "Transit" = 17
    )
  ) +
  labs(
    x = paste0("PCoA1 (", sprintf("%.2f", pcoa1_percent), "%)"),
    y = paste0("PCoA2 (", sprintf("%.2f", pcoa2_percent), "%)"),
    colour = "Environment",
    shape = "Environment"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.7
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "grey30",
      fill = NA,
      linewidth = 0.8
    ),
    axis.title = element_text(
      size = 20,
      colour = "black"
    ),
    axis.text = element_text(
      size = 15,
      colour = "grey30"
    ),
    legend.title = element_text(
      size = 20,
      colour = "black"
    ),
    legend.text = element_text(
      size = 16,
      colour = "black"
    ),
    legend.position = "right",
    legend.key = element_blank(),
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )


# ------------------------------------------------------------
# 9. Save final figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "pcoa_arg_group_sewage1M_transit5M_t80_prev10_clean.png"
  ),
  plot = p,
  width = 9.6,
  height = 7.2,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "pcoa_arg_group_sewage1M_transit5M_t80_prev10_clean.pdf"
  ),
  plot = p,
  width = 9.6,
  height = 7.2
)


# ------------------------------------------------------------
# 10. Print completion message
# ------------------------------------------------------------

message("")
message("====================================================")
message("Figure generated successfully")
message("====================================================")
message("PCoA1 variation: ", sprintf("%.2f", pcoa1_percent), "%")
message("PCoA2 variation: ", sprintf("%.2f", pcoa2_percent), "%")
message("Output folder:")
message(output_dir)

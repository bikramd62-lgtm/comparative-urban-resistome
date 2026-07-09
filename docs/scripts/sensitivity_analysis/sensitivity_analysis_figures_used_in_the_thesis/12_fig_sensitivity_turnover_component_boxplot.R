# ============================================================
# fig_sensitivity_turnover_component_boxplot.R
#
# Purpose:
#   Generate the ARG-group turnover component boxplot for the
#   Jaccard decomposition sensitivity analysis.
#
# Figure:
#   ARG_group_turnover_component_boxplot.png
#
# Input:
#   ARG_group_jaccard_decomposition_pairwise_manual.csv
#
# Input path:
#   data/processed/sensitivity_analysis/jaccard_decomposition/
#
# Output:
#   results/sensitivity_analysis/figures/
#   ├── ARG_group_turnover_component_boxplot.png
#   └── ARG_group_turnover_component_boxplot.pdf
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c("ggplot2")

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

library(ggplot2)


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_root <- Sys.getenv(
  "PROJECT_ROOT",
  unset = getwd()
)

candidate_input_files <- c(
  file.path(
    project_root,
    "results",
    "sensitivity_analysis",
    "jaccard_decomposition",
    "tables",
    "ARG_group_jaccard_decomposition_pairwise_manual.csv"
  ),
  file.path(
    project_root,
    "data",
    "processed",
    "sensitivity_analysis",
    "jaccard_decomposition",
    "ARG_group_jaccard_decomposition_pairwise_manual.csv"
  )
)

existing_input_files <- candidate_input_files[
  file.exists(candidate_input_files)
]

if (length(existing_input_files) == 0) {
  stop(
    paste0(
      "Input file not found. Expected one of:\n",
      paste(candidate_input_files, collapse = "\n")
    )
  )
}

input_file <- existing_input_files[1]

output_dir <- file.path(
  project_root,
  "results",
  "sensitivity_analysis",
  "figures"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

clean_column_name <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]", "", x))
}

find_column <- function(data, candidate_names, required = TRUE) {
  cleaned_names <- clean_column_name(colnames(data))
  cleaned_candidates <- clean_column_name(candidate_names)

  matches <- which(cleaned_names %in% cleaned_candidates)

  if (length(matches) == 0) {
    if (required) {
      stop(
        paste0(
          "Could not find any of these columns:\n",
          paste(candidate_names, collapse = ", "),
          "\n\nAvailable columns:\n",
          paste(colnames(data), collapse = ", ")
        )
      )
    } else {
      return(NULL)
    }
  }

  colnames(data)[matches[1]]
}

standardise_pair_type <- function(x) {
  x <- trimws(as.character(x))

  x <- gsub("-", "–", x)
  x <- gsub("Sewage[– ]+sewage", "Sewage–sewage", x, ignore.case = TRUE)
  x <- gsub("Sewage[– ]+transit", "Sewage–transit", x, ignore.case = TRUE)
  x <- gsub("Transit[– ]+transit", "Transit–transit", x, ignore.case = TRUE)

  x
}


# ------------------------------------------------------------
# 4. Read pairwise Jaccard decomposition table
# ------------------------------------------------------------

decomposition_raw <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pair_type_col <- find_column(
  decomposition_raw,
  c(
    "Pair_type",
    "pair_type",
    "PairType",
    "comparison_type"
  )
)

turnover_col <- find_column(
  decomposition_raw,
  c(
    "Turnover_component",
    "turnover_component",
    "Turnover",
    "beta_jtu",
    "beta.jtu"
  )
)

plot_data <- data.frame(
  Pair_type = standardise_pair_type(decomposition_raw[[pair_type_col]]),
  Turnover_component =
    as.numeric(as.character(decomposition_raw[[turnover_col]])),
  stringsAsFactors = FALSE
)

pair_type_levels <- c(
  "Sewage–sewage",
  "Sewage–transit",
  "Transit–transit"
)

plot_data <- plot_data[
  plot_data$Pair_type %in% pair_type_levels,
  ,
  drop = FALSE
]

plot_data$Pair_type <- factor(
  plot_data$Pair_type,
  levels = pair_type_levels
)

if (anyNA(plot_data$Turnover_component)) {
  stop("Turnover component contains missing or non-numeric values.")
}

if (nrow(plot_data) == 0) {
  stop("No valid pairwise decomposition rows were available for plotting.")
}


# ------------------------------------------------------------
# 5. Save plotting data and summary
# ------------------------------------------------------------

write.csv(
  plot_data,
  file = file.path(
    output_dir,
    "ARG_group_turnover_component_boxplot_plot_data.csv"
  ),
  row.names = FALSE
)

summary_table <- do.call(
  rbind,
  lapply(
    pair_type_levels,
    function(current_pair_type) {
      x <- plot_data$Turnover_component[
        plot_data$Pair_type == current_pair_type
      ]

      data.frame(
        Pair_type = current_pair_type,
        n = length(x),
        median = median(x),
        Q1 = unname(quantile(x, 0.25)),
        Q3 = unname(quantile(x, 0.75)),
        minimum = min(x),
        maximum = max(x),
        mean = mean(x),
        sd = sd(x),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  summary_table,
  file = file.path(
    output_dir,
    "ARG_group_turnover_component_boxplot_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Generate figure
# ------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = Pair_type,
    y = Turnover_component,
    fill = Pair_type,
    colour = Pair_type
  )
) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.78,
    linewidth = 0.9,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.55,
    alpha = 0.32
  ) +
  scale_fill_manual(
    values = c(
      "Sewage–sewage" = "#E66B63",
      "Sewage–transit" = "#7B6DB8",
      "Transit–transit" = "#4C78A8"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Sewage–sewage" = "#E66B63",
      "Sewage–transit" = "#7B6DB8",
      "Transit–transit" = "#4C78A8"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Pair type",
    y = "Turnover component"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_line(
      colour = "grey86",
      linewidth = 0.7
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.9
    ),
    axis.title.x = element_text(
      size = 21,
      colour = "black",
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      size = 21,
      colour = "black",
      margin = margin(r = 12)
    ),
    axis.text.x = element_text(
      size = 17,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 16,
      colour = "black"
    ),
    legend.position = "none",
    plot.margin = margin(
      t = 10,
      r = 12,
      b = 10,
      l = 10
    )
  )


# ------------------------------------------------------------
# 7. Save final figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "ARG_group_turnover_component_boxplot.png"
  ),
  plot = p,
  width = 9.6,
  height = 7.2,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "ARG_group_turnover_component_boxplot.pdf"
  ),
  plot = p,
  width = 9.6,
  height = 7.2
)


# ------------------------------------------------------------
# 8. Save session information
# ------------------------------------------------------------

session_text <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("Input file: ", input_file),
  paste0("R version: ", R.version.string),
  paste0("ggplot2 version: ", as.character(packageVersion("ggplot2"))),
  "",
  "Figure: turnover component boxplot by pair type",
  "Pair types: sewage-sewage, sewage-transit, transit-transit",
  "",
  capture.output(sessionInfo())
)

writeLines(
  session_text,
  con = file.path(
    output_dir,
    "ARG_group_turnover_component_boxplot_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 9. Print completion message
# ------------------------------------------------------------

message("")
message("====================================================")
message("Figure generated successfully")
message("====================================================")
message("Input file:")
message(input_file)
message("")
message("Pair-type summary:")
print(summary_table, row.names = FALSE)
message("")
message("Output folder:")
message(output_dir)

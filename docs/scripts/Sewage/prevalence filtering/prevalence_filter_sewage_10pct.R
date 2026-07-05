#!/usr/bin/env Rscript

# 01_prevalence_filter_sewage_10pct.R
#
# Purpose:
# Build sewage ARG count and binary matrices from ResistomeAnalyzer outputs,
# then apply a 10% prevalence filter.
#
# Workflow:
# 1. Read ResistomeAnalyzer output tables for gene, group, mechanism, and class.
# 2. Collapse each level into city-by-feature count tables.
# 3. Convert count matrices to binary presence/absence matrices.
# 4. Calculate feature prevalence across sewage cities.
# 5. Retain features present in at least 10% of sewage cities.
#
# Important:
# This script writes both unfiltered and prevalence-filtered matrices.

#
# Example run:
#   Rscript scripts/09_prevalence_filtering/01_prevalence_filter_sewage_10pct.R \
#     --input-dir results/amr_tables/sewage/megares_v3/resistome_analyzer \
#     --output-dir results/processed_tables/sewage/prevalence_10pct \
#     --prevalence 0.10

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("\n")
  cat("Usage:\n")
  cat("  Rscript scripts/09_prevalence_filtering/01_prevalence_filter_sewage_10pct.R \\\n")
  cat("    --input-dir <resistome_analyzer_output_dir> \\\n")
  cat("    --output-dir <output_dir> \\\n")
  cat("    --prevalence 0.10\n")
  cat("\n")
  cat("Optional:\n")
  cat("  --input-dir     Directory containing ResistomeAnalyzer city folders\n")
  cat("                  Default: results/amr_tables/sewage/megares_v3/resistome_analyzer\n")
  cat("  --output-dir    Output directory for matrices and prevalence-filtered tables\n")
  cat("                  Default: results/processed_tables/sewage/prevalence_10pct\n")
  cat("  --prevalence    Prevalence threshold as a fraction\n")
  cat("                  Default: 0.10\n")
  cat("  --min-cities    Optional fixed minimum city count. If not provided,\n")
  cat("                  ceiling(prevalence * number_of_cities) is used.\n")
  cat("  --levels        Comma-separated levels to process\n")
  cat("                  Default: gene,group,mechanism,class\n")
  cat("\n")
}

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    return(default)
  }
  if (idx == length(args)) {
    stop(paste("Missing value after", flag), call. = FALSE)
  }
  args[idx + 1]
}

if ("--help" %in% args || "-h" %in% args) {
  usage()
  quit(status = 0)
}

INPUT_DIR <- get_arg(
  "--input-dir",
  "results/amr_tables/sewage/megares_v3/resistome_analyzer"
)

OUTPUT_DIR <- get_arg(
  "--output-dir",
  "results/processed_tables/sewage/prevalence_10pct"
)

PREVALENCE <- as.numeric(get_arg("--prevalence", "0.10"))

MIN_CITIES_ARG <- get_arg("--min-cities", NA)

LEVELS <- strsplit(
  get_arg("--levels", "gene,group,mechanism,class"),
  ",",
  fixed = TRUE
)[[1]]

LEVELS <- trimws(LEVELS)

if (!dir.exists(INPUT_DIR)) {
  stop(paste("Input directory not found:", INPUT_DIR), call. = FALSE)
}

if (is.na(PREVALENCE) || PREVALENCE <= 0 || PREVALENCE > 1) {
  stop("--prevalence must be a number between 0 and 1.", call. = FALSE)
}

if (!is.na(MIN_CITIES_ARG)) {
  MIN_CITIES_ARG <- as.integer(MIN_CITIES_ARG)
  if (is.na(MIN_CITIES_ARG) || MIN_CITIES_ARG < 1) {
    stop("--min-cities must be a positive integer.", call. = FALSE)
  }
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

UNFILTERED_DIR <- file.path(OUTPUT_DIR, "unfiltered")
FILTERED_DIR <- file.path(OUTPUT_DIR, "filtered_10pct")
FEATURE_DIR <- file.path(OUTPUT_DIR, "retained_features")
LOG_DIR <- file.path(OUTPUT_DIR, "logs")

dir.create(UNFILTERED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FILTERED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FEATURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

session_file <- file.path(LOG_DIR, "R_session_info_prevalence_filtering.txt")
sink(session_file)
cat("Prevalence filtering script run date:", as.character(Sys.time()), "\n\n")
cat("Input directory:", INPUT_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Prevalence threshold:", PREVALENCE, "\n")
cat("Levels:", paste(LEVELS, collapse = ", "), "\n\n")
cat("R session info:\n")
print(sessionInfo())
sink()

clean_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

find_column_exact <- function(df, candidates) {
  col_names <- names(df)
  col_clean <- clean_name(col_names)
  candidate_clean <- clean_name(candidates)

  idx <- match(candidate_clean, col_clean, nomatch = 0)
  idx <- idx[idx > 0]

  if (length(idx) == 0) {
    return(NA_character_)
  }

  col_names[idx[1]]
}

feature_candidates_for_level <- function(level) {
  if (level == "gene") {
    return(c(
      "Gene",
      "ARG",
      "Feature",
      "Resistance Gene",
      "MEGARes Gene"
    ))
  }

  if (level == "group") {
    return(c(
      "Group",
      "ARG Group",
      "Resistance Group",
      "Feature",
      "MEGARes Group"
    ))
  }

  if (level == "mechanism") {
    return(c(
      "Mechanism",
      "Resistance Mechanism",
      "ARG Mechanism",
      "Feature",
      "MEGARes Mechanism"
    ))
  }

  if (level == "class") {
    return(c(
      "Class",
      "Drug Class",
      "Resistance Class",
      "ARG Class",
      "Feature",
      "MEGARes Class"
    ))
  }

  stop(paste("Unknown annotation level:", level), call. = FALSE)
}

find_feature_column <- function(df, level) {
  candidate_col <- find_column_exact(df, feature_candidates_for_level(level))

  if (!is.na(candidate_col)) {
    return(candidate_col)
  }

  stop(
    paste(
      "Could not identify feature column for level:",
      level,
      "\nAvailable columns:",
      paste(names(df), collapse = ", ")
    ),
    call. = FALSE
  )
}

find_count_column <- function(df, feature_col) {
  exact_candidates <- c(
    "Hits",
    "Hit",
    "Count",
    "Counts",
    "Read Count",
    "Reads",
    "Mapped Reads",
    "Mapped_Reads",
    "Total Hits"
  )

  exact_col <- find_column_exact(df, exact_candidates)

  if (!is.na(exact_col)) {
    return(exact_col)
  }

  col_names <- names(df)
  col_clean <- clean_name(col_names)

  feature_idx <- which(col_names == feature_col)
  fraction_like <- grepl("fraction|percent|percentage|coverage|length", col_clean)
  count_like <- grepl("hit|count|read", col_clean)

  candidate_idx <- which(count_like & !fraction_like)
  candidate_idx <- setdiff(candidate_idx, feature_idx)

  if (length(candidate_idx) > 0) {
    for (idx in candidate_idx) {
      values <- suppressWarnings(as.numeric(df[[idx]]))
      if (sum(!is.na(values)) > 0) {
        return(col_names[idx])
      }
    }
  }

  numeric_idx <- integer(0)

  for (i in seq_along(col_names)) {
    if (i %in% feature_idx) {
      next
    }

    if (fraction_like[i]) {
      next
    }

    values <- suppressWarnings(as.numeric(df[[i]]))

    if (sum(!is.na(values)) > 0) {
      numeric_idx <- c(numeric_idx, i)
    }
  }

  if (length(numeric_idx) > 0) {
    return(col_names[numeric_idx[1]])
  }

  stop(
    paste(
      "Could not identify count column.\nAvailable columns:",
      paste(names(df), collapse = ", ")
    ),
    call. = FALSE
  )
}

get_city_from_file <- function(file_path, input_dir, level) {
  parent_dir <- dirname(file_path)

  same_as_input <- normalizePath(parent_dir, mustWork = FALSE) ==
    normalizePath(input_dir, mustWork = FALSE)

  if (!same_as_input) {
    return(basename(parent_dir))
  }

  file_name <- basename(file_path)
  city <- sub(
    paste0("_", level, "\\.tsv$"),
    "",
    file_name,
    ignore.case = TRUE
  )

  city <- sub("\\.tsv$", "", city, ignore.case = TRUE)
  city
}

read_resistome_table <- function(file_path, level, input_dir) {
  df <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    fill = TRUE
  )

  if (nrow(df) == 0) {
    return(NULL)
  }

  names(df) <- trimws(names(df))

  feature_col <- find_feature_column(df, level)
  count_col <- find_count_column(df, feature_col)

  city <- get_city_from_file(file_path, input_dir, level)

  features <- trimws(as.character(df[[feature_col]]))
  counts <- suppressWarnings(as.numeric(df[[count_col]]))

  out <- data.frame(
    City = city,
    Feature = features,
    Count = counts,
    Source_file = file_path,
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$Feature) & out$Feature != "", , drop = FALSE]
  out <- out[!is.na(out$Count), , drop = FALSE]
  out <- out[out$Count > 0, , drop = FALSE]

  if (nrow(out) == 0) {
    return(NULL)
  }

  out
}

write_table <- function(df, path) {
  write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

write_matrix_table <- function(mat, path) {
  df <- data.frame(
    Feature = rownames(mat),
    as.data.frame(mat, check.names = FALSE),
    check.names = FALSE
  )

  write_table(df, path)
}

make_prevalence_table <- function(binary_matrix) {
  if (nrow(binary_matrix) == 0) {
    return(data.frame(
      Feature = character(),
      Cities_present = integer(),
      Prevalence_fraction = numeric(),
      Present_in_cities = character(),
      stringsAsFactors = FALSE
    ))
  }

  cities_present <- rowSums(binary_matrix > 0)

  present_in_cities <- apply(binary_matrix, 1, function(x) {
    paste(names(x)[x > 0], collapse = ";")
  })

  data.frame(
    Feature = rownames(binary_matrix),
    Cities_present = as.integer(cities_present),
    Prevalence_fraction = cities_present / ncol(binary_matrix),
    Present_in_cities = present_in_cities,
    stringsAsFactors = FALSE
  )
}

make_city_richness_table <- function(binary_matrix) {
  data.frame(
    City = colnames(binary_matrix),
    ARG_feature_richness = as.integer(colSums(binary_matrix > 0)),
    stringsAsFactors = FALSE
  )
}

process_level <- function(level) {
  cat("\n------------------------------------------------------------\n")
  cat("Processing level:", level, "\n")

  files <- list.files(
    INPUT_DIR,
    pattern = paste0("_", level, "\\.tsv$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0) {
    warning(paste("No files found for level:", level))
    return(data.frame(
      Level = level,
      N_cities = 0,
      N_features_unfiltered = 0,
      Prevalence_threshold = PREVALENCE,
      Minimum_cities_required = NA_integer_,
      N_features_retained = 0,
      Status = "no_files_found",
      stringsAsFactors = FALSE
    ))
  }

  cat("Files found:", length(files), "\n")

  long_list <- list()

  for (file_path in files) {
    cat("Reading:", file_path, "\n")

    one <- tryCatch(
      read_resistome_table(file_path, level, INPUT_DIR),
      error = function(e) {
        warning(paste("Failed to read", file_path, ":", e$message))
        NULL
      }
    )

    if (!is.null(one)) {
      long_list[[length(long_list) + 1]] <- one
    }
  }

  if (length(long_list) == 0) {
    warning(paste("No valid non-zero records for level:", level))
    return(data.frame(
      Level = level,
      N_cities = 0,
      N_features_unfiltered = 0,
      Prevalence_threshold = PREVALENCE,
      Minimum_cities_required = NA_integer_,
      N_features_retained = 0,
      Status = "no_valid_records",
      stringsAsFactors = FALSE
    ))
  }

  long_df <- do.call(rbind, long_list)

  collapsed <- aggregate(
    Count ~ City + Feature,
    data = long_df,
    FUN = sum
  )

  cities <- sort(unique(collapsed$City))
  n_cities <- length(cities)

  if (is.na(MIN_CITIES_ARG)) {
    min_cities <- ceiling(PREVALENCE * n_cities)
    min_cities <- max(min_cities, 1)
  } else {
    min_cities <- MIN_CITIES_ARG
  }

  count_matrix <- xtabs(
    Count ~ Feature + City,
    data = collapsed
  )

  count_matrix <- as.matrix(count_matrix)
  count_matrix <- count_matrix[, cities, drop = FALSE]

  binary_matrix <- (count_matrix > 0) * 1L

  prevalence_table <- make_prevalence_table(binary_matrix)
  prevalence_table <- prevalence_table[
    order(
      -prevalence_table$Cities_present,
      prevalence_table$Feature
    ),
    ,
    drop = FALSE
  ]

  retained_features <- prevalence_table$Feature[
    prevalence_table$Cities_present >= min_cities
  ]

  filtered_count_matrix <- count_matrix[retained_features, , drop = FALSE]
  filtered_binary_matrix <- binary_matrix[retained_features, , drop = FALSE]

  filtered_prevalence_table <- make_prevalence_table(filtered_binary_matrix)
  filtered_prevalence_table <- filtered_prevalence_table[
    order(
      -filtered_prevalence_table$Cities_present,
      filtered_prevalence_table$Feature
    ),
    ,
    drop = FALSE
  ]

  unfiltered_richness <- make_city_richness_table(binary_matrix)
  filtered_richness <- make_city_richness_table(filtered_binary_matrix)

  prev_label <- paste0(round(PREVALENCE * 100), "pct")

  long_path <- file.path(
    UNFILTERED_DIR,
    paste0("sewage_", level, "_long.tsv")
  )

  count_path <- file.path(
    UNFILTERED_DIR,
    paste0("sewage_", level, "_count_matrix.tsv")
  )

  binary_path <- file.path(
    UNFILTERED_DIR,
    paste0("sewage_", level, "_binary_matrix.tsv")
  )

  prevalence_path <- file.path(
    UNFILTERED_DIR,
    paste0("sewage_", level, "_feature_prevalence.tsv")
  )

  richness_path <- file.path(
    UNFILTERED_DIR,
    paste0("sewage_", level, "_city_richness.tsv")
  )

  filtered_count_path <- file.path(
    FILTERED_DIR,
    paste0("sewage_", level, "_count_matrix_prev", prev_label, ".tsv")
  )

  filtered_binary_path <- file.path(
    FILTERED_DIR,
    paste0("sewage_", level, "_binary_matrix_prev", prev_label, ".tsv")
  )

  filtered_prevalence_path <- file.path(
    FILTERED_DIR,
    paste0("sewage_", level, "_feature_prevalence_prev", prev_label, ".tsv")
  )

  filtered_richness_path <- file.path(
    FILTERED_DIR,
    paste0("sewage_", level, "_city_richness_prev", prev_label, ".tsv")
  )

  retained_features_path <- file.path(
    FEATURE_DIR,
    paste0("sewage_", level, "_retained_features_prev", prev_label, ".txt")
  )

  write_table(long_df, long_path)
  write_matrix_table(count_matrix, count_path)
  write_matrix_table(binary_matrix, binary_path)
  write_table(prevalence_table, prevalence_path)
  write_table(unfiltered_richness, richness_path)

  write_matrix_table(filtered_count_matrix, filtered_count_path)
  write_matrix_table(filtered_binary_matrix, filtered_binary_path)
  write_table(filtered_prevalence_table, filtered_prevalence_path)
  write_table(filtered_richness, filtered_richness_path)

  writeLines(retained_features, retained_features_path)

  cat("Cities detected:", n_cities, "\n")
  cat("Unfiltered features:", nrow(binary_matrix), "\n")
  cat("Prevalence threshold:", PREVALENCE, "\n")
  cat("Minimum cities required:", min_cities, "\n")
  cat("Retained features:", length(retained_features), "\n")

  data.frame(
    Level = level,
    N_cities = n_cities,
    N_features_unfiltered = nrow(binary_matrix),
    Prevalence_threshold = PREVALENCE,
    Minimum_cities_required = min_cities,
    N_features_retained = length(retained_features),
    Status = "completed",
    stringsAsFactors = FALSE
  )
}

summary_list <- lapply(LEVELS, process_level)
summary_df <- do.call(rbind, summary_list)

summary_path <- file.path(
  OUTPUT_DIR,
  "sewage_prevalence_filtering_summary.tsv"
)

write_table(summary_df, summary_path)

cat("\n============================================================\n")
cat("Sewage prevalence filtering complete.\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Summary table:", summary_path, "\n")
cat("Session info:", session_file, "\n")
cat("============================================================\n")

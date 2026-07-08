# ============================================================
# 02_RQ3_urban_covariate_statistical_analysis.R
#
# Purpose:
# Run the RQ3 urban covariate analysis:
#   1. Spearman correlations for continuous covariates
#   2. Kruskal-Wallis / Wilcoxon tests for categorical covariates
#   3. PERMANOVA for shared ARG-group composition
#   4. PERMDISP for categorical predictors
#   5. Benjamini-Hochberg FDR correction
#   6. Leave-one-city-out robustness checks
#
# Inputs:
#   results/RQ3_output_tables/RQ3_master_table_FINAL_corrected.csv
#   results/RQ3_output_tables/RQ3_shared_ARG_group_matrix.csv
#
# Outputs:
#   results/RQ3_output_tables/RQ3_spearman_continuous_covariates.csv
#   results/RQ3_output_tables/RQ3_categorical_overlap_tests.csv
#   results/RQ3_output_tables/RQ3_shared_composition_PERMANOVA.csv
#   results/RQ3_output_tables/RQ3_categorical_PERMDISP.csv
#   results/RQ3_output_tables/RQ3_leave_one_city_out_spearman.csv
#   results/RQ3_output_tables/RQ3_leave_one_city_out_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(vegan)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_dir  <- "results/RQ3_output_tables"
output_dir <- "results/RQ3_output_tables"

master_file <- file.path(input_dir, "RQ3_master_table_FINAL_corrected.csv")
shared_matrix_file <- file.path(input_dir, "RQ3_shared_ARG_group_matrix.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(master_file, shared_matrix_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required input files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

standardize_city <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[[:space:]-]+", "_", x)
  x
}

pick_col <- function(df, candidates, variable_name) {
  hit <- candidates[candidates %in% names(df)]

  if (length(hit) == 0) {
    stop(
      "Could not find column for: ", variable_name, "\n",
      "Tried these names:\n",
      paste(candidates, collapse = ", "), "\n\n",
      "Available columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }

  hit[1]
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

clean_categorical <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "na", "n/a", "N/A")] <- NA
  factor(x)
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

# ------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------

master_raw <- read_csv(master_file, show_col_types = FALSE, name_repair = "minimal")

city_col <- pick_col(
  master_raw,
  candidates = c("city", "City", "CITY", "city_name", "City_name", "city_clean", "City_clean"),
  variable_name = "city"
)

master <- master_raw %>%
  rename(city = all_of(city_col)) %>%
  mutate(city = standardize_city(city))

if (any(duplicated(master$city))) {
  stop(
    "Duplicated city names detected in master table:\n",
    paste(unique(master$city[duplicated(master$city)]), collapse = ", ")
  )
}

shared_raw <- read_csv(shared_matrix_file, show_col_types = FALSE, name_repair = "minimal")

shared_city_col <- pick_col(
  shared_raw,
  candidates = c("city", "City", "CITY", "city_name", "City_name", "city_clean", "City_clean"),
  variable_name = "city in shared ARG matrix"
)

shared_matrix_df <- shared_raw %>%
  rename(city = all_of(shared_city_col)) %>%
  mutate(city = standardize_city(city))

if (any(duplicated(shared_matrix_df$city))) {
  stop(
    "Duplicated city names detected in shared ARG matrix:\n",
    paste(unique(shared_matrix_df$city[duplicated(shared_matrix_df$city)]), collapse = ", ")
  )
}

common_cities <- intersect(master$city, shared_matrix_df$city)

if (length(common_cities) == 0) {
  stop("No common cities found between master table and shared ARG matrix.")
}

master <- master %>%
  filter(city %in% common_cities) %>%
  arrange(match(city, common_cities))

shared_matrix_df <- shared_matrix_df %>%
  filter(city %in% common_cities) %>%
  arrange(match(city, master$city))

shared_mat <- shared_matrix_df %>%
  select(-city) %>%
  mutate(across(everything(), ~ as.integer(safe_numeric(.x) > 0))) %>%
  as.matrix()

rownames(shared_mat) <- shared_matrix_df$city

# Remove features absent from all cities, if any
shared_mat <- shared_mat[, colSums(shared_mat, na.rm = TRUE) > 0, drop = FALSE]

message("Loaded RQ3 master table with ", nrow(master), " cities.")
message("Loaded shared ARG matrix with ", nrow(shared_mat), " cities and ", ncol(shared_mat), " shared ARG groups.")

# ------------------------------------------------------------
# 4. Define response and covariate columns
# ------------------------------------------------------------

jaccard_col <- pick_col(
  master,
  candidates = c(
    "jaccard_similarity_sewage_transit_ARG_group",
    "jaccard_similarity",
    "sewage_transit_jaccard_similarity",
    "Jaccard_similarity",
    "Jaccard"
  ),
  variable_name = "Jaccard similarity"
)

shared_count_col <- pick_col(
  master,
  candidates = c(
    "shared_ARG_group_count",
    "shared_arg_group_count",
    "Shared_ARG_group_count",
    "shared_ARG_count",
    "shared_count"
  ),
  variable_name = "shared ARG-group count"
)

temperature_col <- pick_col(
  master,
  candidates = c(
    "city_avg_june_temperature",
    "city_ave_june_temp_c",
    "average_june_temperature",
    "june_temperature_c",
    "temperature_june_c",
    "june_temp_c"
  ),
  variable_name = "average June temperature"
)

humidity_col <- pick_col(
  master,
  candidates = c(
    "june_relative_humidity_percent",
    "relative_humidity_june_percent",
    "humidity_percent",
    "june_humidity_percent",
    "humidity"
  ),
  variable_name = "June relative humidity"
)

precipitation_col <- pick_col(
  master,
  candidates = c(
    "june_precipitation_mm_total",
    "june_precipitation_total_mm",
    "precipitation_june_mm",
    "june_precipitation_mm",
    "precipitation_mm",
    "precipitation"
  ),
  variable_name = "June precipitation"
)

gdp_col <- pick_col(
  master,
  candidates = c(
    "GDP_per_capita",
    "gdp_per_capita",
    "gdp_pc",
    "gdp_per_capita_usd",
    "GDP_per_capita_current_USD",
    "gdp_per_capita_current_usd"
  ),
  variable_name = "GDP per capita"
)

sanitation_col <- pick_col(
  master,
  candidates = c(
    "sanitation_coverage_percent",
    "safely_managed_sanitation_percent",
    "sanitation_percent",
    "sanitation",
    "safe_sanitation"
  ),
  variable_name = "sanitation coverage"
)

wastewater_col <- pick_col(
  master,
  candidates = c(
    "wastewater_treatment_coverage_percent",
    "treated_wastewater_percent",
    "safely_treated_domestic_wastewater_percent",
    "wastewater_percent",
    "wastewater",
    "wastewater_treatment"
  ),
  variable_name = "treated wastewater"
)

population_col <- pick_col(
  master,
  candidates = c(
    "city_population",
    "city_total_population",
    "population",
    "Population",
    "total_population"
  ),
  variable_name = "city population"
)

population_density_col <- pick_col(
  master,
  candidates = c(
    "city_population_density",
    "population_density",
    "pop_density",
    "Population_density",
    "population_density_per_km2"
  ),
  variable_name = "population density"
)

elevation_col <- pick_col(
  master,
  candidates = c(
    "elevation",
    "city_elevation_meters",
    "elevation_m",
    "Elevation_m",
    "city_elevation_m"
  ),
  variable_name = "elevation"
)

region_col <- pick_col(
  master,
  candidates = c(
    "region",
    "Region",
    "continent",
    "Continent",
    "world_region"
  ),
  variable_name = "region"
)

climate_col <- pick_col(
  master,
  candidates = c(
    "city_climate",
    "city_koppen_climate",
    "koppen_climate",
    "climate_class",
    "Climate_class",
    "climate"
  ),
  variable_name = "climate class"
)

coastal_col <- pick_col(
  master,
  candidates = c(
    "coastal",
    "coastal_city",
    "Coastal",
    "is_coastal"
  ),
  variable_name = "coastal status"
)

message("\nColumn mapping used:")
print(tibble(
  variable = c(
    "Jaccard similarity",
    "Shared ARG-group count",
    "Temperature",
    "Humidity",
    "Precipitation",
    "GDP per capita",
    "Sanitation",
    "Wastewater",
    "Population",
    "Population density",
    "Elevation",
    "Region",
    "Climate class",
    "Coastal status"
  ),
  column = c(
    jaccard_col,
    shared_count_col,
    temperature_col,
    humidity_col,
    precipitation_col,
    gdp_col,
    sanitation_col,
    wastewater_col,
    population_col,
    population_density_col,
    elevation_col,
    region_col,
    climate_col,
    coastal_col
  )
))

analysis_df <- master %>%
  transmute(
    city = city,

    jaccard_similarity = safe_numeric(.data[[jaccard_col]]),
    shared_ARG_group_count = safe_numeric(.data[[shared_count_col]]),

    temperature_june_c = safe_numeric(.data[[temperature_col]]),
    humidity_june_percent = safe_numeric(.data[[humidity_col]]),
    precipitation_june_mm = safe_numeric(.data[[precipitation_col]]),
    gdp_per_capita = safe_numeric(.data[[gdp_col]]),
    sanitation_percent = safe_numeric(.data[[sanitation_col]]),
    wastewater_percent = safe_numeric(.data[[wastewater_col]]),
    population = safe_numeric(.data[[population_col]]),
    population_density = safe_numeric(.data[[population_density_col]]),
    elevation_m = safe_numeric(.data[[elevation_col]]),

    region = clean_categorical(.data[[region_col]]),
    climate_class = clean_categorical(.data[[climate_col]]),
    coastal = clean_categorical(.data[[coastal_col]])
  )

continuous_covariates <- c(
  "temperature_june_c",
  "humidity_june_percent",
  "precipitation_june_mm",
  "gdp_per_capita",
  "sanitation_percent",
  "wastewater_percent",
  "population",
  "population_density",
  "elevation_m"
)

categorical_covariates <- c(
  "region",
  "climate_class",
  "coastal"
)

responses <- c(
  "jaccard_similarity",
  "shared_ARG_group_count"
)

# ------------------------------------------------------------
# 5. Spearman correlations for continuous covariates
# ------------------------------------------------------------

run_spearman <- function(data, response, covariate) {
  test_df <- data %>%
    select(city, response = all_of(response), covariate = all_of(covariate)) %>%
    filter(!is.na(response), !is.na(covariate))

  if (nrow(test_df) < 4) {
    return(tibble(
      response = response,
      covariate = covariate,
      n = nrow(test_df),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }

  ct <- suppressWarnings(
    cor.test(
      test_df$response,
      test_df$covariate,
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    response = response,
    covariate = covariate,
    n = nrow(test_df),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

spearman_results <- expand_grid(
  response = responses,
  covariate = continuous_covariates
) %>%
  pmap_dfr(~ run_spearman(
    data = analysis_df,
    response = ..1,
    covariate = ..2
  )) %>%
  mutate(
    p_FDR = p.adjust(p_value, method = "BH"),
    p_value_formatted = format_p(p_value),
    p_FDR_formatted = format_p(p_FDR)
  ) %>%
  arrange(response, p_FDR, p_value)

write_csv(
  spearman_results,
  file.path(output_dir, "RQ3_spearman_continuous_covariates.csv")
)

message("\nSpearman analysis completed.")
print(spearman_results)

# ------------------------------------------------------------
# 6. Categorical tests for overlap magnitude
# ------------------------------------------------------------

run_categorical_test <- function(data, response, covariate) {
  test_df <- data %>%
    select(city, response = all_of(response), group = all_of(covariate)) %>%
    filter(!is.na(response), !is.na(group)) %>%
    mutate(group = droplevels(group))

  n_groups <- nlevels(test_df$group)

  if (nrow(test_df) < 4 || n_groups < 2) {
    return(tibble(
      response = response,
      covariate = covariate,
      test = NA_character_,
      n = nrow(test_df),
      groups = n_groups,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }

  if (n_groups == 2) {
    wt <- suppressWarnings(
      wilcox.test(response ~ group, data = test_df, exact = FALSE)
    )

    tibble(
      response = response,
      covariate = covariate,
      test = "Wilcoxon rank-sum",
      n = nrow(test_df),
      groups = n_groups,
      statistic = unname(wt$statistic),
      p_value = wt$p.value
    )

  } else {
    kt <- kruskal.test(response ~ group, data = test_df)

    tibble(
      response = response,
      covariate = covariate,
      test = "Kruskal-Wallis",
      n = nrow(test_df),
      groups = n_groups,
      statistic = unname(kt$statistic),
      p_value = kt$p.value
    )
  }
}

categorical_results <- expand_grid(
  response = responses,
  covariate = categorical_covariates
) %>%
  pmap_dfr(~ run_categorical_test(
    data = analysis_df,
    response = ..1,
    covariate = ..2
  )) %>%
  mutate(
    p_FDR = p.adjust(p_value, method = "BH"),
    p_value_formatted = format_p(p_value),
    p_FDR_formatted = format_p(p_FDR)
  ) %>%
  arrange(response, p_FDR, p_value)

write_csv(
  categorical_results,
  file.path(output_dir, "RQ3_categorical_overlap_tests.csv")
)

message("\nCategorical overlap tests completed.")
print(categorical_results)

# ------------------------------------------------------------
# 7. PERMANOVA for shared ARG-group composition
# ------------------------------------------------------------

run_permanova <- function(shared_mat, metadata, covariate) {
  model_df <- metadata %>%
    select(city, predictor = all_of(covariate)) %>%
    filter(!is.na(predictor))

  model_cities <- intersect(model_df$city, rownames(shared_mat))

  model_df <- model_df %>%
    filter(city %in% model_cities) %>%
    arrange(match(city, model_cities))

  mat_sub <- shared_mat[model_df$city, , drop = FALSE]

  if (nrow(model_df) < 4) {
    return(tibble(
      covariate = covariate,
      n = nrow(model_df),
      df = NA_real_,
      sum_of_squares = NA_real_,
      R2 = NA_real_,
      pseudo_F = NA_real_,
      p_value = NA_real_
    ))
  }

  if (is.factor(model_df$predictor) || is.character(model_df$predictor)) {
    model_df$predictor <- droplevels(factor(model_df$predictor))
    if (nlevels(model_df$predictor) < 2) {
      return(tibble(
        covariate = covariate,
        n = nrow(model_df),
        df = NA_real_,
        sum_of_squares = NA_real_,
        R2 = NA_real_,
        pseudo_F = NA_real_,
        p_value = NA_real_
      ))
    }
  } else {
    model_df$predictor <- safe_numeric(model_df$predictor)
  }

  dist_sub <- vegdist(mat_sub, method = "jaccard", binary = TRUE)

  set.seed(123)

  ad <- adonis2(
    dist_sub ~ predictor,
    data = model_df,
    permutations = 999
  )

  tibble(
    covariate = covariate,
    n = nrow(model_df),
    df = ad$Df[1],
    sum_of_squares = ad$SumOfSqs[1],
    R2 = ad$R2[1],
    pseudo_F = ad$F[1],
    p_value = ad$`Pr(>F)`[1]
  )
}

permanova_covariates <- c(continuous_covariates, categorical_covariates)

permanova_results <- map_dfr(
  permanova_covariates,
  ~ run_permanova(
    shared_mat = shared_mat,
    metadata = analysis_df,
    covariate = .x
  )
) %>%
  mutate(
    p_FDR = p.adjust(p_value, method = "BH"),
    p_value_formatted = format_p(p_value),
    p_FDR_formatted = format_p(p_FDR)
  ) %>%
  arrange(p_FDR, p_value)

write_csv(
  permanova_results,
  file.path(output_dir, "RQ3_shared_composition_PERMANOVA.csv")
)

message("\nPERMANOVA analysis completed.")
print(permanova_results)

# ------------------------------------------------------------
# 8. PERMDISP for categorical covariates
# ------------------------------------------------------------

run_permdisp <- function(shared_mat, metadata, covariate) {
  model_df <- metadata %>%
    select(city, group = all_of(covariate)) %>%
    filter(!is.na(group)) %>%
    mutate(group = droplevels(factor(group)))

  model_cities <- intersect(model_df$city, rownames(shared_mat))

  model_df <- model_df %>%
    filter(city %in% model_cities) %>%
    arrange(match(city, model_cities)) %>%
    mutate(group = droplevels(group))

  mat_sub <- shared_mat[model_df$city, , drop = FALSE]

  if (nrow(model_df) < 4 || nlevels(model_df$group) < 2) {
    return(tibble(
      covariate = covariate,
      n = nrow(model_df),
      groups = nlevels(model_df$group),
      F_value = NA_real_,
      p_value = NA_real_
    ))
  }

  dist_sub <- vegdist(mat_sub, method = "jaccard", binary = TRUE)

  bd <- betadisper(dist_sub, model_df$group)

  set.seed(123)

  pt <- permutest(bd, permutations = 999)

  tibble(
    covariate = covariate,
    n = nrow(model_df),
    groups = nlevels(model_df$group),
    F_value = unname(pt$tab$F[1]),
    p_value = unname(pt$tab$`Pr(>F)`[1])
  )
}

permdisp_results <- map_dfr(
  categorical_covariates,
  ~ run_permdisp(
    shared_mat = shared_mat,
    metadata = analysis_df,
    covariate = .x
  )
) %>%
  mutate(
    p_FDR = p.adjust(p_value, method = "BH"),
    p_value_formatted = format_p(p_value),
    p_FDR_formatted = format_p(p_FDR)
  ) %>%
  arrange(p_FDR, p_value)

write_csv(
  permdisp_results,
  file.path(output_dir, "RQ3_categorical_PERMDISP.csv")
)

message("\nPERMDISP analysis completed.")
print(permdisp_results)

# ------------------------------------------------------------
# 9. Leave-one-city-out robustness checks for Spearman tests
# ------------------------------------------------------------

run_loco_spearman <- function(data, response, covariate) {
  cities <- data$city

  map_dfr(cities, function(excluded_city) {
    test_df <- data %>%
      filter(city != excluded_city) %>%
      select(city, response = all_of(response), covariate = all_of(covariate)) %>%
      filter(!is.na(response), !is.na(covariate))

    if (nrow(test_df) < 4) {
      return(tibble(
        response = response,
        covariate = covariate,
        excluded_city = excluded_city,
        n = nrow(test_df),
        rho = NA_real_,
        p_value = NA_real_
      ))
    }

    ct <- suppressWarnings(
      cor.test(
        test_df$response,
        test_df$covariate,
        method = "spearman",
        exact = FALSE
      )
    )

    tibble(
      response = response,
      covariate = covariate,
      excluded_city = excluded_city,
      n = nrow(test_df),
      rho = unname(ct$estimate),
      p_value = ct$p.value
    )
  })
}

loco_results <- expand_grid(
  response = responses,
  covariate = continuous_covariates
) %>%
  pmap_dfr(~ run_loco_spearman(
    data = analysis_df,
    response = ..1,
    covariate = ..2
  ))

write_csv(
  loco_results,
  file.path(output_dir, "RQ3_leave_one_city_out_spearman.csv")
)

loco_summary <- loco_results %>%
  group_by(response, covariate) %>%
  summarise(
    n_exclusions = sum(!is.na(rho)),
    min_rho = min(rho, na.rm = TRUE),
    max_rho = max(rho, na.rm = TRUE),
    mean_rho = mean(rho, na.rm = TRUE),
    sd_rho = sd(rho, na.rm = TRUE),
    same_sign_as_full = {
      full_rho <- spearman_results %>%
        filter(response == cur_group()$response, covariate == cur_group()$covariate) %>%
        pull(rho)

      if (length(full_rho) == 0 || is.na(full_rho)) {
        NA_integer_
      } else {
        sum(sign(rho) == sign(full_rho), na.rm = TRUE)
      }
    },
    most_negative_exclusion = excluded_city[which.min(rho)],
    most_positive_exclusion = excluded_city[which.max(rho)],
    .groups = "drop"
  )

write_csv(
  loco_summary,
  file.path(output_dir, "RQ3_leave_one_city_out_summary.csv")
)

message("\nLeave-one-city-out analysis completed.")
print(loco_summary)

# ------------------------------------------------------------
# 10. Compact key-result summary
# ------------------------------------------------------------

key_loco <- loco_summary %>%
  filter(
    (response == "jaccard_similarity" & covariate == "population_density") |
      (response == "shared_ARG_group_count" & covariate == "temperature_june_c") |
      (response == "shared_ARG_group_count" & covariate == "humidity_june_percent")
  )

write_csv(
  key_loco,
  file.path(output_dir, "RQ3_leave_one_city_out_key_relationships.csv")
)

message("\nKey leave-one-city-out relationships:")
print(key_loco)

# ------------------------------------------------------------
# 11. Final console message
# ------------------------------------------------------------

message("\nRQ3 urban covariate statistical analysis completed successfully.")
message("Output files written to: ", output_dir)
message(" - RQ3_spearman_continuous_covariates.csv")
message(" - RQ3_categorical_overlap_tests.csv")
message(" - RQ3_shared_composition_PERMANOVA.csv")
message(" - RQ3_categorical_PERMDISP.csv")
message(" - RQ3_leave_one_city_out_spearman.csv")
message(" - RQ3_leave_one_city_out_summary.csv")
message(" - RQ3_leave_one_city_out_key_relationships.csv")

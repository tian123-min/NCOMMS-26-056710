# Bootstrap and Monte Carlo uncertainty estimates.
# CSV-oriented implementation following the R3 MXT Monte Carlo workflow while
# retaining the 01-04 model-training data and RF parameters used in this upload.

source("1.Load packages.R")

set.seed(123)

n_bootstrap <- as.integer(Sys.getenv("N_BOOTSTRAP", "100"))
n_mc <- as.integer(Sys.getenv("N_MC", "1000"))
var_noise <- as.numeric(Sys.getenv("VAR_NOISE", "0.05"))
num_trees <- as.integer(Sys.getenv("NUM_TREES", "500"))

output_root <- Sys.getenv("OUTPUT_DIR", "outputs")
uncertainty_dir <- file.path(output_root, "uncertainty")
model_dir <- file.path(output_root, "models")

dir.create(uncertainty_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

dir_01_2 <- file.path(uncertainty_dir, "01-2_monte_carlo_4models_5scenarios")
dir_01_3 <- file.path(uncertainty_dir, "01-3_monte_carlo_with_id")
dir_01_7 <- file.path(uncertainty_dir, "01-7_uncertainty_cv_by_ci")
for (out_dir in c(dir_01_2, dir_01_3, dir_01_7)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

model_specs <- data.frame(
  indicator = c("CS", "yield", "CO2", "N2O"),
  target = c("CUE", "yield", "CO2", "N2O"),
  file = c(
    "model_training/C_seq.csv",
    "model_training/yield.csv",
    "model_training/CO2.csv",
    "model_training/N2O.csv"
  ),
  mtry = c(8, 7, 5, 7),
  stringsAsFactors = FALSE
)

scenario_names <- c("Wheat_irr", "Wheat_rain", "Maize_irr", "Maize_rain", "Wheat_maize")
model_predictors <- c(
  "MAT", "MAP", "SOC", "pH", "Clay", "Duration", "C_input", "N_input",
  "Planting_system", "Return", "Tillage", "Irrigation", "Topdressing"
)
strata_cols <- c("Planting_system", "Tillage", "Irrigation", "Topdressing", "Return")
indicators <- c("CS", "yield", "CO2", "N2O")

write_utf8_csv <- function(dat, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(dat, file, row.names = FALSE, fileEncoding = "UTF-8")
}

check_required_columns <- function(dat, required, file) {
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop(file, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

read_training_data <- function(file, target_var) {
  dat <- read.csv(file, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, check.names = FALSE)
  status_col <- intToUtf8(89)
  if (status_col %in% names(dat)) {
    stop(file, " still contains the raw status-filter column. Use processed model_training CSV files.")
  }
  check_required_columns(dat, c(model_predictors, target_var), file)
  dat <- dat[, c(model_predictors, target_var), drop = FALSE]
  dat %>% mutate(across(where(is.character), as.factor))
}

read_prediction_input <- function(file) {
  dat <- read.csv(file, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, check.names = FALSE)
  check_required_columns(dat, model_predictors, file)
  dat %>% mutate(across(where(is.character), as.factor))
}

train_bootstrap_model <- function(data_file, target_var, mtry_val) {
  dat <- read_training_data(data_file, target_var)
  strata <- interaction(dat[, strata_cols, drop = FALSE], sep = "_", drop = TRUE)

  boot_models <- vector("list", n_bootstrap)
  p <- length(model_predictors)
  for (b in seq_len(n_bootstrap)) {
    idx <- caret::createDataPartition(strata, p = 0.8, list = FALSE)
    boot_dat <- dat[idx, , drop = FALSE]
    boot_models[[b]] <- ranger(
      dependent.variable.name = target_var,
      data = boot_dat,
      mtry = min(mtry_val, p),
      min.node.size = 5,
      splitrule = "extratrees",
      respect.unordered.factors = TRUE,
      num.trees = num_trees
    )
  }
  cat(sprintf("  Trained %s bootstrap models (mtry=%d)\n", target_var, mtry_val))
  boot_models
}

read_current_inputs <- function() {
  result <- list()
  for (scenario in scenario_names) {
    file <- file.path("current", paste0(scenario, ".csv"))
    dat <- read_prediction_input(file)
    dat$ID <- seq_len(nrow(dat))
    dat$crop_type <- scenario
    result[[scenario]] <- dat
  }
  result
}

run_monte_carlo <- function(mc_df, boot_models) {
  pred_store <- lapply(names(boot_models), function(x) list())
  names(pred_store) <- names(boot_models)

  for (mc in seq_len(n_mc)) {
    cat("  MC iteration: ", mc, "/", n_mc, "\n", sep = "")
    noisy_df <- mc_df
    num_cols <- sapply(noisy_df, is.numeric) & !names(noisy_df) %in% c("ID")
    if (sum(num_cols) > 0) {
      noisy_df[, num_cols] <- noisy_df[, num_cols] * rnorm(nrow(noisy_df), 1, var_noise)
    }
    noisy_df <- noisy_df %>% mutate(across(where(is.character), as.factor))

    for (ind in names(boot_models)) {
      for (b in seq_along(boot_models[[ind]])) {
        pred_store[[ind]][[length(pred_store[[ind]]) + 1]] <-
          predict(boot_models[[ind]][[b]], noisy_df)$predictions
      }
    }
  }

  res <- mc_df %>% select(ID, crop_type)
  for (ind in names(pred_store)) {
    mat <- do.call(cbind, pred_store[[ind]])
    res[[paste0(ind, "_mean")]] <- rowMeans(mat, na.rm = TRUE)
    res[[paste0(ind, "_2.5")]] <- rowQuantiles(mat, probs = 0.025, na.rm = TRUE)
    res[[paste0(ind, "_97.5")]] <- rowQuantiles(mat, probs = 0.975, na.rm = TRUE)
  }
  res
}

read_standard_ids <- function() {
  list(
    Wheat_irr = read.csv("id/wheat.csv", fileEncoding = "UTF-8-BOM")$ID,
    Wheat_rain = read.csv("id/wheat.csv", fileEncoding = "UTF-8-BOM")$ID,
    Maize_irr = read.csv("id/maize.csv", fileEncoding = "UTF-8-BOM")$ID,
    Maize_rain = read.csv("id/maize.csv", fileEncoding = "UTF-8-BOM")$ID,
    Wheat_maize = read.csv("id/wheat_maize.csv", fileEncoding = "UTF-8-BOM")$ID
  )
}

bind_standard_id <- function(result_list) {
  id_map <- read_standard_ids()
  for (scenario in names(result_list)) {
    if (length(id_map[[scenario]]) != nrow(result_list[[scenario]])) {
      stop("Standard ID row count does not match scenario: ", scenario)
    }
    result_list[[scenario]]$ID <- id_map[[scenario]]
  }
  result_list
}

weighted_average <- function(irr, rain, area, prefix) {
  out <- data.frame(ID = area$ID)
  area_irr <- area[[paste0(prefix, "_irr_area")]]
  area_rain <- area[[paste0(prefix, "_rain_area")]]
  area_total <- area[[paste0(prefix, "_total_area")]]
  mc_cols <- grep("_(mean|2.5|97.5)$", names(irr), value = TRUE)

  for (col in mc_cols) {
    irr_vec <- irr[[col]][match(area$ID, irr$ID)]
    rain_vec <- rain[[col]][match(area$ID, rain$ID)]
    out[[col]] <- ifelse(
      is.na(area_total) | area_total == 0,
      NA_real_,
      (irr_vec * area_irr + rain_vec * area_rain) / area_total
    )
  }
  out
}

combine_crop_systems <- function(scenario_results) {
  area <- read.csv("area_weights.csv", fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  check_required_columns(
    area,
    c("ID", "wheat_irr_area", "wheat_rain_area", "maize_irr_area", "maize_rain_area",
      "wheat_total_area", "maize_total_area"),
    "area_weights.csv"
  )

  wheat <- weighted_average(scenario_results$Wheat_irr, scenario_results$Wheat_rain, area, "wheat")
  maize <- weighted_average(scenario_results$Maize_irr, scenario_results$Maize_rain, area, "maize")
  wheat_maize <- area %>%
    select(ID) %>%
    left_join(scenario_results$Wheat_maize %>% select(-crop_type), by = "ID")

  wheat$crop_system <- "Wheat"
  maize$crop_system <- "Maize"
  wheat_maize$crop_system <- "Wheat_maize"

  bind_rows(wheat, maize, wheat_maize) %>% relocate(ID, crop_system)
}

add_sd_cv <- function(df) {
  for (ind in indicators) {
    mean_col <- paste0(ind, "_mean")
    lower_col <- paste0(ind, "_2.5")
    upper_col <- paste0(ind, "_97.5")
    sd_col <- paste0(ind, "_sd")
    cv_col <- paste0(ind, "_cv")

    df[[sd_col]] <- (df[[upper_col]] - df[[lower_col]]) / 3.92
    df[[cv_col]] <- df[[sd_col]] / df[[mean_col]]
  }
  df
}

calc_cv_by_ci <- function(df, indicator) {
  mean_val <- df[[paste0(indicator, "_mean")]]
  lower <- df[[paste0(indicator, "_2.5")]]
  upper <- df[[paste0(indicator, "_97.5")]]
  cv <- ((upper - lower) / 2) / abs(mean_val) * 100
  cv[!is.finite(cv)] <- NA_real_
  cv
}

write_cv_by_ci_tables <- function(final_results) {
  for (ind in indicators) {
    temp <- final_results
    temp$cv <- calc_cv_by_ci(temp, ind)
    cv_df <- temp %>%
      select(ID, crop_system, cv) %>%
      pivot_wider(names_from = crop_system, values_from = cv)

    output <- data.frame(ID = cv_df$ID)
    for (system in c("Wheat", "Maize", "Wheat_maize")) {
      output[[paste0(system, "_CV")]] <- cv_df[[system]]
    }
    write_utf8_csv(output, file.path(dir_01_7, paste0(ind, "_cv_by_ci.csv")))
  }
}

cat("\n===== Training bootstrap models =====\n")
boot_models <- list()
for (i in seq_len(nrow(model_specs))) {
  spec <- model_specs[i, ]
  cat("Training ", spec$indicator, "\n", sep = "")
  boot_models[[spec$indicator]] <- train_bootstrap_model(spec$file, spec$target, spec$mtry)
  saveRDS(boot_models[[spec$indicator]], file.path(model_dir, paste0(spec$indicator, "_bootstrap_models.rds")))
}

cat("\n===== Running Monte Carlo simulations =====\n")
current_inputs <- read_current_inputs()
scenario_results <- list()
for (scenario in names(current_inputs)) {
  cat("\nScenario: ", scenario, "\n", sep = "")
  scenario_results[[scenario]] <- run_monte_carlo(current_inputs[[scenario]], boot_models)
  write_utf8_csv(scenario_results[[scenario]], file.path(dir_01_2, paste0(scenario, ".csv")))
}

scenario_results_with_id <- bind_standard_id(scenario_results)
for (scenario in names(scenario_results_with_id)) {
  write_utf8_csv(scenario_results_with_id[[scenario]], file.path(dir_01_3, paste0(scenario, ".csv")))
}

final_results <- combine_crop_systems(scenario_results_with_id)
write_cv_by_ci_tables(final_results)

final_results_with_sd_cv <- final_results %>% add_sd_cv()
write_utf8_csv(final_results_with_sd_cv, file.path(uncertainty_dir, "01-8_final_results_with_sd_cv.csv"))
write_utf8_csv(final_results_with_sd_cv, file.path(uncertainty_dir, "final_results_with_sd_cv.csv"))

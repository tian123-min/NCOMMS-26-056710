# County-level predictions for the four indicators.

source("1.Load packages.R")

set.seed(123)
quick_run <- tolower(Sys.getenv("QUICK_RUN", "false")) %in% c("1", "true", "yes")
num_trees <- if (quick_run) 100 else 500
output_root <- Sys.getenv("OUTPUT_DIR", "outputs")
prediction_dir <- file.path(output_root, "county_predictions")
model_dir <- file.path(output_root, "models")

dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

indicator_specs <- data.frame(
  indicator = c("C_seq", "yield", "CO2", "N2O"),
  target = c("CUE", "yield", "CO2", "N2O"),
  file = c(
    "data/model_training/C_seq.csv",
    "data/model_training/yield.csv",
    "data/model_training/CO2.csv",
    "data/model_training/N2O.csv"
  ),
  mtry = c(8, 7, 5, 7),
  stringsAsFactors = FALSE
)

read_training_data <- function(file) {
  dat <- read.csv(file, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, check.names = FALSE)
  status_col <- intToUtf8(89)
  if (status_col %in% names(dat)) {
    stop(
      file,
      " still contains the raw status-filter column. Re-run tools/prepare_upload_data.py ",
      "and use the processed data/model_training CSV files."
    )
  }
  dat <- dat %>% mutate(across(where(is.character), as.factor))
  dat
}

fit_prediction_model <- function(spec) {
  dat <- read_training_data(spec$file)
  p <- ncol(dat) - 1
  ranger(
    dependent.variable.name = spec$target,
    data = dat,
    mtry = min(spec$mtry, p),
    min.node.size = 5,
    splitrule = "extratrees",
    respect.unordered.factors = TRUE,
    importance = "permutation",
    num.trees = num_trees
  )
}

predict_scenario_dir <- function(model, indicator, system_name, input_dir) {
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)

  combined <- NULL
  for (file in files) {
    dat <- read.csv(file, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, check.names = FALSE)
    dat <- dat %>% mutate(across(where(is.character), as.factor))
    pred <- predict(model, data = dat)$predictions
    scenario <- tools::file_path_sans_ext(basename(file))
    scenario <- sub("^Table_", "", scenario)

    if (is.null(combined)) {
      id_cols <- intersect(c("ID", "MAT", "MAP", "SOC", "pH", "Clay", "Duration", "C_input", "N_input", "Planting_system"), names(dat))
      combined <- dat[, id_cols, drop = FALSE]
    }
    combined[[scenario]] <- pred
  }

  out_file <- file.path(prediction_dir, paste(indicator, system_name, "scenarios.csv", sep = "_"))
  write.csv(combined, out_file, row.names = FALSE, fileEncoding = "UTF-8")
  out_file
}

predict_current_inputs <- function(model, indicator) {
  input_dir <- "data/county_inputs/current"
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  rows <- list()

  for (file in files) {
    dat <- read.csv(file, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE, check.names = FALSE)
    dat <- dat %>% mutate(across(where(is.character), as.factor))
    pred <- predict(model, data = dat)$predictions
    rows[[basename(file)]] <- data.frame(
      row_id = seq_len(nrow(dat)),
      crop_system = tools::file_path_sans_ext(basename(file)),
      prediction = pred
    )
  }

  result <- bind_rows(rows)
  write.csv(result, file.path(prediction_dir, paste0(indicator, "_current_baseline.csv")), row.names = FALSE, fileEncoding = "UTF-8")
}

scenario_dirs <- c(
  wheat = "data/county_inputs/scenarios/wheat",
  maize = "data/county_inputs/scenarios/maize",
  wheat_maize = "data/county_inputs/scenarios/wheat_maize"
)

for (i in seq_len(nrow(indicator_specs))) {
  spec <- indicator_specs[i, ]
  cat("\n===== ", spec$indicator, " county predictions =====\n", sep = "")
  model <- fit_prediction_model(spec)
  saveRDS(model, file.path(model_dir, paste0(spec$indicator, "_county_RF.rds")))
  predict_current_inputs(model, spec$indicator)
  for (system_name in names(scenario_dirs)) {
    predict_scenario_dir(model, spec$indicator, system_name, scenario_dirs[[system_name]])
  }
}

# Machine-learning models for yield, soil carbon sequestration, CO2, and N2O.
# This script preserves the original per-indicator modelling logic and only
# changes paths, encodings, and upload-friendly outputs.

source("1.Load packages.R")

set.seed(123)

output_root <- Sys.getenv("OUTPUT_DIR", "outputs")
performance_dir <- file.path(output_root, "model_performance")
model_dir <- file.path(output_root, "models")
dir.create(performance_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

indicator_specs <- list(
  C_seq = list(
    label = "C_seq",
    target = "CUE",
    data_file = "data/model_training/C_seq.csv",
    gbm_grid_id = "GBM_CUE_grid",
    relevel = TRUE
  ),
  yield = list(
    label = "yield",
    target = "yield",
    data_file = "data/model_training/yield.csv",
    gbm_grid_id = "GBM_yield_grid",
    relevel = FALSE
  ),
  CO2 = list(
    label = "CO2",
    target = "CO2",
    data_file = "data/model_training/CO2.csv",
    gbm_grid_id = "GBM_CO2_grid",
    relevel = FALSE
  ),
  N2O = list(
    label = "N2O",
    target = "N2O",
    data_file = "data/model_training/N2O.csv",
    gbm_grid_id = "GBM_N2O_grid",
    relevel = FALSE
  )
)

num_vars <- c("MAT", "MAP", "SOC", "Clay", "pH", "Duration", "C_input", "N_input")
fac_vars <- c("Planting_system", "Return", "Tillage", "Irrigation", "Topdressing")
model_predictors <- c("MAT", "MAP", "SOC", "pH", "Clay", "Planting_system",
                      "Duration", "N_input", "C_input", "Tillage", "Return",
                      "Irrigation", "Topdressing")

read_indicator_data <- function(spec) {
  dat <- read.csv(spec$data_file, fileEncoding = "UTF-8-BOM",
                  stringsAsFactors = FALSE, check.names = FALSE)
  status_col <- intToUtf8(89)
  if (status_col %in% names(dat)) {
    stop(
      spec$data_file,
      " still contains the raw status-filter column. Re-run tools/prepare_upload_data.py ",
      "and use the processed data/model_training CSV files."
    )
  }
  required_cols <- c(model_predictors, spec$target)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(spec$data_file, " is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  dat <- dat %>% mutate(across(all_of(fac_vars), as.factor))
  if (isTRUE(spec$relevel)) {
    dat$Planting_system <- relevel(dat$Planting_system, ref = "Wheat")
    dat$Tillage <- relevel(dat$Tillage, ref = "CT")
    dat$Irrigation <- relevel(dat$Irrigation, ref = "NI")
    dat$Topdressing <- relevel(dat$Topdressing, ref = "No")
    dat$Return <- relevel(dat$Return, ref = "SI")
  }
  dat
}

split_train_test <- function(dat) {
  strata_group <- interaction(dat$Planting_system, dat$Tillage,
                              dat$Irrigation, dat$Topdressing,
                              dat$Return,
                              sep = "_", drop = TRUE)
  train_index <- createDataPartition(strata_group, p = 0.8, list = FALSE)
  list(
    train = dat[train_index, , drop = FALSE],
    test = dat[-train_index, , drop = FALSE]
  )
}

metric_fun <- function(obs, pred, name) {
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae <- mean(abs(obs - pred), na.rm = TRUE)
  r2 <- cor(obs, pred, use = "complete.obs")^2
  rng <- diff(range(obs, na.rm = TRUE))
  data.frame(
    Dataset = name,
    RMSE = rmse,
    nRMSE_range_percent = rmse / rng * 100,
    MAE = mae,
    R2 = r2
  )
}

write_result_bundle <- function(result_list, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(result_list, path = output_file)

  csv_dir <- file.path(dirname(output_file), tools::file_path_sans_ext(basename(output_file)))
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  for (sheet_name in names(result_list)) {
    safe_name <- gsub("[^A-Za-z0-9_\\-]+", "_", sheet_name)
    write.csv(result_list[[sheet_name]],
              file.path(csv_dir, paste0(safe_name, ".csv")),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
}

model_output_file <- function(indicator, model_name) {
  file.path(
    performance_dir,
    indicator,
    paste0(model_name, "_", indicator, "_model_results.xlsx")
  )
}

save_model <- function(model, indicator, model_name) {
  saveRDS(model, file.path(model_dir, paste(indicator, model_name, "rds", sep = ".")))
}

make_formula <- function(target) {
  as.formula(paste(target, "~", paste(model_predictors, collapse = " + ")))
}

with_caret_parallel <- function(expr) {
  detected_cores <- parallel::detectCores()
  num_cores <- if (is.na(detected_cores)) 1 else max(1, detected_cores - 1)
  if (num_cores <= 1) {
    return(eval.parent(substitute(expr)))
  }

  cl <- parallel::makePSOCKcluster(num_cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }, add = TRUE)

  eval.parent(substitute(expr))
}

run_rf <- function(spec, dat, split) {
  target <- spec$target
  train_dat <- split$train
  test_dat <- split$test

  control <- trainControl(method = "cv", number = 10, savePredictions = "final")
  grid <- expand.grid(
    mtry = seq(2, 8, by = 1),
    splitrule = c("extratrees"),
    min.node.size = c(5, 10, 20, 30)
  )

  set.seed(123)
  rf_model_cv <- with_caret_parallel(caret::train(
    as.formula(paste(target, "~ .")),
    data = train_dat,
    method = "ranger",
    trControl = control,
    tuneGrid = grid,
    metric = "RMSE",
    importance = "permutation",
    respect.unordered.factors = TRUE,
    num.trees = 500
  ))
  save_model(rf_model_cv, spec$label, "RF")

  pred_fun_rf <- function(model, input_dat, name) {
    pred <- as.numeric(predict(model, newdata = input_dat))
    pred_result <- input_dat %>%
      mutate(
        prediction = pred,
        residual = .data[[target]] - prediction
      )
    list(
      pred = pred_result,
      metrics = metric_fun(pred_result[[target]], pred_result$prediction, name)
    )
  }

  test_out <- pred_fun_rf(rf_model_cv, test_dat, "Test")
  train_out <- pred_fun_rf(rf_model_cv, train_dat, "Train")
  all_out <- pred_fun_rf(rf_model_cv, dat, "All data")

  train_range <- diff(range(train_dat[[target]], na.rm = TRUE))
  cv_results <- rf_model_cv$results %>%
    mutate(nRMSE_range_percent = RMSE / train_range * 100) %>%
    arrange(RMSE)

  best_cv_result <- cv_results %>%
    semi_join(rf_model_cv$bestTune, by = c("mtry", "splitrule", "min.node.size")) %>%
    transmute(
      Dataset = "10-fold CV best",
      mtry,
      splitrule,
      min.node.size,
      RMSE,
      nRMSE_range_percent,
      MAE,
      R2 = Rsquared
    )

  performance_summary <- bind_rows(
    train_out$metrics,
    test_out$metrics,
    all_out$metrics,
    best_cv_result %>% select(Dataset, RMSE, nRMSE_range_percent, MAE, R2)
  )

  var_importance <- varImp(rf_model_cv)$importance %>%
    rownames_to_column("variable") %>%
    arrange(desc(Overall))

  output_file <- model_output_file(spec$label, "RF")
  write_result_bundle(
    list(
      "01_Test_prediction" = test_out$pred,
      "02_Train_prediction" = train_out$pred,
      "03_All_prediction" = all_out$pred,
      "04_CV_all_results" = cv_results,
      "05_Best_parameter" = best_cv_result,
      "06_Performance_summary" = performance_summary,
      "07_Variable_importance" = var_importance
    ),
    output_file
  )
}

run_ann <- function(spec, dat, split) {
  target <- spec$target
  train_dat <- split$train
  test_dat <- split$test

  train_control <- trainControl(method = "cv", number = 10, savePredictions = "final")
  tune_grid <- expand.grid(
    size = seq(2, 20, by = 2),
    decay = c(0, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2)
  )

  set.seed(123)
  nn_model_cv <- with_caret_parallel(caret::train(
    make_formula(target),
    data = train_dat,
    method = "nnet",
    trControl = train_control,
    tuneGrid = tune_grid,
    linout = TRUE,
    preProcess = c("center", "scale"),
    maxit = 1000,
    trace = FALSE,
    MaxNWts = 10000
  ))
  save_model(nn_model_cv, spec$label, "ANN")

  pred_fun <- function(input_dat, name) {
    pred <- as.numeric(predict(nn_model_cv, newdata = input_dat))
    pred_result <- input_dat %>%
      mutate(
        data_pred = pred,
        residual = .data[[target]] - data_pred
      )
    list(
      pred = pred_result,
      metrics = metric_fun(pred_result[[target]], pred_result$data_pred, name)
    )
  }

  test_out <- pred_fun(test_dat, "Test")
  train_out <- pred_fun(train_dat, "Train")
  all_out <- pred_fun(dat, "All data")

  train_range <- diff(range(train_dat[[target]], na.rm = TRUE))
  cv_results <- nn_model_cv$results %>%
    mutate(nRMSE_range_percent = RMSE / train_range * 100) %>%
    arrange(RMSE)

  best_cv_result <- cv_results %>%
    semi_join(nn_model_cv$bestTune, by = c("size", "decay")) %>%
    transmute(
      Dataset = "10-fold CV best",
      size,
      decay,
      RMSE,
      nRMSE_range_percent,
      MAE,
      R2 = Rsquared
    )

  performance_summary <- bind_rows(
    train_out$metrics,
    test_out$metrics,
    all_out$metrics,
    best_cv_result %>% select(Dataset, RMSE, nRMSE_range_percent, MAE, R2)
  )

  cv_pred_result <- if (!is.null(nn_model_cv$pred) && nrow(nn_model_cv$pred) > 0) {
    nn_model_cv$pred %>% mutate(residual = obs - pred)
  } else {
    data.frame(Note = "No CV predictions saved. Use savePredictions = 'final' in trainControl.")
  }

  output_file <- model_output_file(spec$label, "ANN")
  write_result_bundle(
    list(
      "01_Test_prediction" = test_out$pred,
      "02_Train_prediction" = train_out$pred,
      "03_All_prediction" = all_out$pred,
      "04_CV_all_results" = cv_results,
      "05_Best_parameter" = best_cv_result,
      "06_Performance_summary" = performance_summary,
      "07_CV_prediction" = cv_pred_result
    ),
    output_file
  )
}

run_svm <- function(spec, dat, split) {
  target <- spec$target
  train_dat <- split$train
  test_dat <- split$test

  ctrl <- trainControl(method = "cv", number = 10, savePredictions = "final")
  svm_grid <- expand.grid(
    sigma = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1),
    C = c(0.1, 0.3, 1, 3, 10, 30, 100)
  )

  set.seed(123)
  svm_tune <- with_caret_parallel(caret::train(
    make_formula(target),
    data = train_dat,
    method = "svmRadial",
    trControl = ctrl,
    tuneGrid = svm_grid,
    preProcess = c("center", "scale"),
    metric = "RMSE"
  ))
  save_model(svm_tune, spec$label, "SVM")

  pred_fun <- function(input_dat, name) {
    pred <- as.numeric(predict(svm_tune, newdata = input_dat))
    pred_result <- input_dat %>%
      mutate(
        data_pred = pred,
        residual = .data[[target]] - data_pred
      )
    list(
      pred = pred_result,
      metrics = metric_fun(pred_result[[target]], pred_result$data_pred, name)
    )
  }

  test_out <- pred_fun(test_dat, "Test")
  train_out <- pred_fun(train_dat, "Train")
  all_out <- pred_fun(dat, "All data")

  train_range <- diff(range(train_dat[[target]], na.rm = TRUE))
  cv_results <- svm_tune$results %>%
    mutate(nRMSE_range_percent = RMSE / train_range * 100) %>%
    arrange(RMSE)

  best_cv_result <- cv_results %>%
    semi_join(svm_tune$bestTune, by = c("sigma", "C")) %>%
    transmute(
      Dataset = "10-fold CV best",
      sigma,
      C,
      RMSE,
      nRMSE_range_percent,
      MAE,
      R2 = Rsquared
    )

  num_support_vectors <- tryCatch(
    length(svm_tune$finalModel@SVindex),
    error = function(e) sum(sapply(svm_tune$finalModel@alphaindex, length))
  )

  support_vector_info <- data.frame(
    Model = "SVM radial",
    Train_n = nrow(train_dat),
    Support_vectors = num_support_vectors,
    Support_vector_percent = num_support_vectors / nrow(train_dat) * 100
  )

  performance_summary <- bind_rows(
    train_out$metrics,
    test_out$metrics,
    all_out$metrics,
    best_cv_result %>% select(Dataset, RMSE, nRMSE_range_percent, MAE, R2)
  )

  cv_pred_result <- if (!is.null(svm_tune$pred) && nrow(svm_tune$pred) > 0) {
    svm_tune$pred %>% mutate(residual = obs - pred)
  } else {
    data.frame(Note = "No CV predictions saved. Use savePredictions = 'final' in trainControl.")
  }

  output_file <- model_output_file(spec$label, "SVM")
  write_result_bundle(
    list(
      "01_Test_prediction" = test_out$pred,
      "02_Train_prediction" = train_out$pred,
      "03_All_prediction" = all_out$pred,
      "04_CV_all_results" = cv_results,
      "05_Best_parameter" = best_cv_result,
      "06_Performance_summary" = performance_summary,
      "07_CV_prediction" = cv_pred_result,
      "08_Support_vectors" = support_vector_info
    ),
    output_file
  )
}

run_gbm <- function(spec, dat, split) {
  target <- spec$target
  train_dat <- split$train
  test_dat <- split$test

  h2o.init()
  x_vars <- c(num_vars, fac_vars)
  y_var <- target

  train_h2o <- as.h2o(train_dat)
  test_h2o <- as.h2o(test_dat)
  all_h2o <- as.h2o(dat)

  hyper_params <- list(
    learn_rate = c(0.01, 0.03, 0.05, 0.1),
    max_depth = c(2, 3, 4, 5),
    ntrees = c(100, 200, 500),
    min_rows = c(5, 10, 20),
    sample_rate = c(0.7, 0.8, 1.0),
    col_sample_rate = c(0.7, 0.8, 1.0)
  )

  search_criteria <- list(
    strategy = "RandomDiscrete",
    max_models = 80,
    seed = 123,
    stopping_metric = "RMSE",
    stopping_tolerance = 0.001,
    stopping_rounds = 5
  )

  set.seed(123)
  gbm_grid <- h2o.grid(
    algorithm = "gbm",
    grid_id = spec$gbm_grid_id,
    x = x_vars,
    y = y_var,
    training_frame = train_h2o,
    distribution = "gaussian",
    hyper_params = hyper_params,
    search_criteria = search_criteria,
    nfolds = 10,
    fold_assignment = "Random",
    seed = 123,
    keep_cross_validation_predictions = TRUE,
    keep_cross_validation_models = FALSE
  )

  gbm_grid_sorted <- h2o.getGrid(
    grid_id = spec$gbm_grid_id,
    sort_by = "rmse",
    decreasing = FALSE
  )

  best_model_id <- gbm_grid_sorted@model_ids[[1]]
  gbm_model <- h2o.getModel(best_model_id)
  h2o.saveModel(gbm_model, path = normalizePath(model_dir, mustWork = FALSE), force = TRUE)

  pred_fun_gbm <- function(model, h2o_dat, raw_dat, name) {
    pred <- as.data.frame(h2o.predict(model, h2o_dat))$predict
    pred_result <- raw_dat %>%
      mutate(
        data_pred = as.numeric(pred),
        residual = .data[[target]] - data_pred
      )
    list(
      pred = pred_result,
      metrics = metric_fun(pred_result[[target]], pred_result$data_pred, name)
    )
  }

  train_out <- pred_fun_gbm(gbm_model, train_h2o, train_dat, "Train")
  test_out <- pred_fun_gbm(gbm_model, test_h2o, test_dat, "Test")
  all_out <- pred_fun_gbm(gbm_model, all_h2o, dat, "All data")

  cv_perf <- h2o.performance(gbm_model, xval = TRUE)
  cv_rmse <- h2o.rmse(cv_perf)
  cv_mae <- h2o.mae(cv_perf)
  cv_r2 <- h2o.r2(cv_perf)

  train_range <- diff(range(train_dat[[target]], na.rm = TRUE))
  cv_best_result <- data.frame(
    Dataset = "10-fold CV best",
    RMSE = cv_rmse,
    nRMSE_range_percent = cv_rmse / train_range * 100,
    MAE = cv_mae,
    R2 = cv_r2
  )

  cv_results <- as.data.frame(gbm_grid_sorted@summary_table)
  num_cols <- c("rmse", "mse", "mae", "r2", "mean_residual_deviance")
  for (cc in intersect(num_cols, names(cv_results))) {
    cv_results[[cc]] <- as.numeric(cv_results[[cc]])
  }
  if ("rmse" %in% names(cv_results)) {
    cv_results <- cv_results %>%
      mutate(nRMSE_range_percent = rmse / train_range * 100) %>%
      arrange(rmse)
  }

  best_params <- data.frame(
    model_id = gbm_model@model_id,
    learn_rate = gbm_model@allparameters$learn_rate,
    max_depth = gbm_model@allparameters$max_depth,
    ntrees = gbm_model@allparameters$ntrees,
    min_rows = gbm_model@allparameters$min_rows,
    sample_rate = gbm_model@allparameters$sample_rate,
    col_sample_rate = gbm_model@allparameters$col_sample_rate
  )

  best_parameter_result <- cbind(
    best_params,
    cv_best_result %>% select(RMSE, nRMSE_range_percent, MAE, R2)
  )

  var_importance <- as.data.frame(h2o.varimp(gbm_model))
  performance_summary <- bind_rows(
    train_out$metrics,
    test_out$metrics,
    all_out$metrics,
    cv_best_result
  )

  output_file <- model_output_file(spec$label, "GBM")
  write_result_bundle(
    list(
      "01_Test_prediction" = test_out$pred,
      "02_Train_prediction" = train_out$pred,
      "03_All_prediction" = all_out$pred,
      "04_CV_all_results" = cv_results,
      "05_Best_parameter" = best_parameter_result,
      "06_Performance_summary" = performance_summary,
      "07_Variable_importance" = var_importance
    ),
    output_file
  )
}

run_indicator <- function(spec) {
  cat("\n===== ", spec$label, " =====\n", sep = "")
  dat <- read_indicator_data(spec)
  split <- split_train_test(dat)

  run_rf(spec, dat, split)
  run_ann(spec, dat, split)
  run_svm(spec, dat, split)
  run_gbm(spec, dat, split)
}

for (spec in indicator_specs) {
  run_indicator(spec)
}

try(h2o.shutdown(prompt = FALSE), silent = TRUE)

# ==========================================
# Dependencies
# ==========================================

library(tidyverse)
library(patchwork)
library(pander)
library(DT)
library(tidymodels)
library(glmnet)
library(vip)
library(doParallel)
library(broom)
library(knitr)
library(kableExtra)
library(themis)



# ==========================================
# Functions
# ==========================================

# Helper to clean and format datasets uniformly
clean_dataset <- function(df) {
  df |> 
    mutate(
      is_gxe = as.factor(is_gxe),
      across(where(is.character) | where(is.logical), as.factor)
    ) |> 
    mutate(is_gxe = factor(is_gxe, levels = c(TRUE, FALSE)))
}

# HELPER: Run LASSO Pipeline
run_lasso_pipeline <- function(data, folds, dataset_name = "") {
  recipe <- recipe(is_gxe ~ ., data = data) |> 
    step_dummy(all_nominal_predictors(), -all_outcomes()) |> 
    step_corr(all_numeric_predictors(), threshold = 0.85) |> 
    step_normalize(all_numeric_predictors()) |>
    step_downsample(is_gxe)
  
  spec <- logistic_reg(penalty = tune(), mixture = 1) |> set_engine("glmnet")
  wf <- workflow() |> add_recipe(recipe) |> add_model(spec)
  
  tune_res <- tune_grid(wf, resamples = folds, grid = grid_regular(penalty(), levels = 50), control = control_grid(save_pred = TRUE))
  best_param <- select_best(tune_res, metric = "roc_auc")
  final_fit <- wf |> finalize_workflow(best_param) |> fit(data = data)
  
  p <- final_fit |> extract_fit_parsnip() |> vip(num_features = 15, geom = "col", fill = "steelblue") +
    theme_minimal() + labs(title = paste("LASSO"), x = "Features", y = "Coefficient Magnitude")
  
  return(list(data=tune_res,plot=p))
}

# HELPER: Run SVM Pipeline (COMPLETELY SILENCED)
run_svm_pipeline <- function(data, folds, dataset_name = "") {
  recipe <- recipe(is_gxe ~ ., data = data) |> 
    step_dummy(all_nominal_predictors(), -all_outcomes()) |> 
    step_corr(all_numeric_predictors(), threshold = 0.85) |> 
    step_normalize(all_numeric_predictors()) |>
    step_downsample(is_gxe)
  
  spec <- svm_rbf(cost = tune(), rbf_sigma = tune()) |> 
    set_engine("kernlab", prob.model = TRUE, trace = FALSE) |> 
    set_mode("classification")
  
  wf <- workflow() |> add_recipe(recipe) |> add_model(spec)
  
  # CRITICAL SILENCE FIX: capture.output() intercepts C-level console output streams
  ignored_output <- capture.output({
    tune_res <- suppressMessages(suppressWarnings(
      tune_grid(wf, resamples = folds, grid = grid_regular(cost(), rbf_sigma(), levels = 3), control = control_grid(save_pred = TRUE))
    ))
  })
  
  best_param <- select_best(tune_res, metric = "roc_auc")
  
  ignored_output2 <- capture.output({
    final_fit <- wf |> finalize_workflow(best_param) |> fit(data = data)
  })
  
  processed_data <- bake(prep(recipe), new_data = data)
  svm_pred_wrapper <- function(object, newdata) { kernlab::predict(object, newdata = newdata, type = "probabilities")[, "TRUE"] }
  
  p <- vip::vi_permute(object = extract_fit_engine(final_fit), train = processed_data |> select(-is_gxe), target = processed_data$is_gxe, metric = "roc_auc", event_level = "first", pred_wrapper = svm_pred_wrapper, reference_class = "TRUE") |> 
    vip(num_features = 15, geom = "col", fill = "darkred") + theme_minimal() + labs(title = paste("SVM"), x = "Features", y = "Drop in ROC AUC")
  
  return(list(data=tune_res,plot=p))
}

# HELPER: Run Random Forest Pipeline
run_rf_pipeline <- function(data, folds, dataset_name = "") {
  recipe <- recipe(is_gxe ~ ., data = data) |> 
    step_dummy(all_nominal_predictors(), -all_outcomes()) |>
    step_downsample(is_gxe)
  spec <- rand_forest(mtry = tune(), trees = 500, min_n = tune()) |> set_engine("ranger", importance = "permutation", num.threads = 1) |> set_mode("classification")
  wf <- workflow() |> add_recipe(recipe) |> add_model(spec)
  
  all_cores <- parallel::detectCores(logical = FALSE)
  cl <- makePSOCKcluster(all_cores - 1)
  registerDoParallel(cl)
  
  tune_res <- tune_grid(wf, resamples = folds, grid = grid_regular(mtry(range = c(2, min(10, ncol(data)-1))), min_n(), levels = 3), control = control_grid(save_pred = TRUE, parallel_over = "resamples"))
  
  stopCluster(cl)
  registerDoSEQ()
  
  best_param <- select_best(tune_res, metric = "roc_auc")
  final_fit <- wf |> finalize_workflow(best_param) |> fit(data = data)
  
  p <- final_fit |> extract_fit_parsnip() |> vip(num_features = 15, geom = "col", fill = "darkgreen") +
    theme_minimal() + labs(title = paste("Random Forest"), x = "Features", y = "Importance Score")
  
  
  return(list(data = tune_res, plot = p))
}



# HELPER: Run Standard Logistic Regression Pipeline (With Fixed Prediction Storage)
run_logistic_pipeline <- function(data, folds, dataset_name = "") {
  
  recipe <- recipe(is_gxe ~ ., data = data) |> 
    step_dummy(all_nominal_predictors(), -all_outcomes()) |> 
    step_corr(all_numeric_predictors(), threshold = 0.85) |> 
    step_normalize(all_numeric_predictors()) |>
    step_downsample(is_gxe)
  
  spec <- logistic_reg() |> set_engine("glm")
  wf   <- workflow() |> add_recipe(recipe) |> add_model(spec)
  
  # 1. CRITICAL: Define the control object to explicitly force prediction capture
  eval_control <- control_resamples(save_pred = TRUE, event_level = "first")
  
  # 2. Evaluate across folds while tracking the explicit control object
  resample_res <- fit_resamples(
    wf, 
    resamples = folds, 
    metrics = metric_set(roc_auc),
    control = eval_control # CRITICAL FIX: Make sure this parameter is applied
  )
  
  # 3. EXTRA SAFEGUARD: Check if predictions were extracted successfully
  preds <- collect_predictions(resample_res)
  
  if (nrow(preds) > 0) {
    # If standard collection works, append the dummy penalty column directly to the dataframe layout
    preds <- preds |> mutate(penalty = 0)
    attr(resample_res, "predictions") <- preds
  } else {
    stop("Predictions were not captured. Double-check that your cross-validation folds match your training split.")
  }
  
  # 4. Generate VIP plot data 
  final_fit <- wf |> fit(data = data)
  p <- final_fit |> 
    extract_fit_parsnip() |> 
    vip(num_features = 15, geom = "col", fill = "seagreen") +
    theme_minimal() + 
    labs(title = paste("Logistic"), x = "Features", y = "Absolute Significance")
  
  # 5. Return the list setup safely
  return(list(data=resample_res,  plot = p))
}

render_bold_summary <- function(model, table_title) {
  # 1. Tidy the model data first
  tidy_df <- broom::tidy(model)
  
  # 2. Identify significant rows, but exclude row 1 (the Intercept)
  significant_rows <- which(tidy_df$p.value < 0.05)
  significant_rows <- setdiff(significant_rows, 1) # Force-remove the intercept row
  
  # 3. Format and render the table safely
  tidy_df |> 
    mutate(
      estimate  = round(estimate, 4),
      std.error = round(std.error, 4),
      statistic = round(statistic, 4),
      p.value   = ifelse(p.value < 0.0001, "< 0.0001", as.character(round(p.value, 4)))
    ) |> 
    knitr::kable(
      format = "html", 
      escape = FALSE, 
      caption = table_title,
      col.names = c("Term (Variable)", "Estimate (Slope)", "Std. Error", "z-statistic", "p-value")
    ) |> 
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) |> 
    # Highlight only the predictor rows we flagged
    kableExtra::row_spec(
      row = significant_rows, 
      bold = TRUE, 
      background = "grey90" 
    )
}

# HELPER: Plot Model Comparisons with Safe HTML Table Output
plot_model_comparison <- function(lasso_tune, svm_tune, rf_tune, log_tune, dataset_name = "") {
  best_pen <- select_best(lasso_tune, metric = "roc_auc")
  best_svm <- select_best(svm_tune, metric = "roc_auc")
  best_rf  <- select_best(rf_tune, metric = "roc_auc")
  
  lasso_preds <- collect_predictions(lasso_tune) |> filter(penalty == best_pen$penalty) |> mutate(model = "LASSO")
  svm_preds   <- collect_predictions(svm_tune) |> inner_join(best_svm, by = c("cost", "rbf_sigma")) |> mutate(model = "SVM")
  rf_preds    <- collect_predictions(rf_tune) |> inner_join(best_rf, by = c("mtry", "min_n")) |> mutate(model = "Random Forest")
  log_preds   <- collect_predictions(log_tune) |> mutate(model = "Logistic Regression")
  
  combined_preds <- bind_rows(lasso_preds, svm_preds, rf_preds, log_preds)
  
  auc_table <- combined_preds |> 
    group_by(model) |> 
    roc_auc(truth = is_gxe, .pred_TRUE, event_level = "first") |> 
    ungroup() |> 
    select(Model = model, Metric = .metric, Estimate = .estimate) |> 
    mutate(Estimate = round(Estimate, 4)) # Fix: No leading dot here anymore!
  
  roc_plot <- combined_preds |> 
    group_by(model) |> 
    roc_curve(truth = is_gxe, .pred_TRUE, event_level = "first") |> 
    autoplot() +
    labs(
      title = paste("ROC Curve Comparison:", dataset_name), 
      subtitle = "Evaluated using stratified 5-fold cross-validation", 
      x = "False Positive Rate", 
      y = "True Positive Rate"
    ) + 
    theme_minimal()
  
  print(roc_plot)
  
  # FIX: Swap pander out for an explicit HTML kable print statement
  
  auc_table |> 
    knitr::kable(format = "html", caption = paste("Comparative ROC AUC Performance:", dataset_name)) |> 
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) |> 
    print()
}

# Wrapper execution block
run_all_mods <- function(data, titles = "", Group_title=""){
  data_sv <- clean_dataset(data)
  set.seed(123)
  folds_sv <- vfold_cv(data_sv, v = 5, strata = is_gxe)
  
  sv_log   <- run_logistic_pipeline(data_sv, folds_sv, titles)
  sv_lasso <- run_lasso_pipeline(data_sv, folds_sv, titles)
  sv_svm   <- run_svm_pipeline(data_sv, folds_sv, titles)
  sv_rf    <- run_rf_pipeline(data_sv, folds_sv, titles)
  
  # Render Patchwork Feature Importance Map
  print((sv_lasso$plot + sv_svm$plot) / (sv_rf$plot + sv_log$plot) +
          plot_annotation(
            title = paste0("Feature Importance: ", Group_title),
            theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
          ))
  
  # Compare Performance (Plot & Table 1)
  plot_model_comparison(sv_lasso$data, sv_svm$data, sv_rf$data, sv_log$data, titles)
  
  # Fit standard summary GLM
  sv_log1 <- glm(data = data, is_gxe ~., family = "binomial")
  
  # FIX: Force print() around your bold summary table so it displays in the HTML layout
  
  print(render_bold_summary(sv_log1, paste0("Logistic: ", Group_title)))
  
  return(list(logistic = sv_log, lasso = sv_lasso, svm = sv_svm, rf = sv_rf))
}

# Helper function to compute standalone univariate metrics for every feature
calculate_predictive_metrics <- function(data, target_var = "is_gxe") {
  
  predictors <- setdiff(names(data), target_var)
  results <- list()
  
  # CRITICAL: Explicitly order levels so TRUE is the case of interest (level 2)
  # pROC treats the second level as the "positive" event by default when direction="auto"
  y <- factor(data[[target_var]], levels = c("FALSE", "TRUE"))
  
  for (pred in predictors) {
    x <- data[[pred]]
    
    # 1. Calculate Standalone Univariate ROC AUC
    if (is.numeric(x)) {
      # Let pROC automatically find if a high or low value predicts TRUE
      # This keeps the AUC honest relative to finding your rare TRUEs
      raw_auc <- suppressMessages(
        pROC::auc(y, x, quiet = TRUE, direction = "auto")
      )
      auc_score <- as.numeric(raw_auc)
      
    } else {
      # For categorical data, a quick GLM handles the probabilities
      fit <- glm(y ~ x, family = "binomial")
      preds <- predict(fit, type = "response")
      raw_auc <- pROC::auc(y, preds, quiet = TRUE, direction = "auto")
      auc_score <- as.numeric(raw_auc)
    }
    
    # 2. Calculate Mutual Information
    x_discrete <- if(is.numeric(x)) cut(x, breaks = min(10, length(unique(x))), labels = FALSE, duplicates = "drop") else x
    
    tab <- table(x_discrete, y)
    row_sums <- rowSums(tab)
    col_sums <- colSums(tab)
    total <- sum(tab)
    
    mi_score <- 0
    for(i in 1:nrow(tab)) {
      for(j in 1:ncol(tab)) {
        if(tab[i,j] > 0) {
          mi_score <- mi_score + (tab[i,j]/total) * log2((tab[i,j]/total) / ((row_sums[i]/total) * (col_sums[j]/total)))
        }
      }
    }
    
    # Track origin safely using exists() check
    origin <- case_when(
      exists("sv1") && pred %in% colnames(sv1)  ~ "SV Only",
      exists("dap1") && pred %in% colnames(dap1) ~ "DAP Only",
      TRUE                                      ~ "Combined"
    )
    
    results[[pred]] <- tibble(
      Variable = pred,
      Source = origin,
      Univariate_AUC = round(auc_score, 4),
      Mutual_Information = round(mi_score, 4)
    )
  }
  
  # Combine results and rank
  # We look for AUCs furthest from 0.50 (either close to 1.0 or close to 0.0)
  final_table <- bind_rows(results) |> 
    mutate(Predictive_Magnitude = abs(Univariate_AUC - 0.50)) |> 
    arrange(desc(Predictive_Magnitude), desc(Mutual_Information)) |> 
    select(-Predictive_Magnitude)
  
  return(final_table)
}

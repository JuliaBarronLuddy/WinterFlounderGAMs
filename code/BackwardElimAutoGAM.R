library(data.table)
library(mgcv)
library(dplyr)

# === Backward Elimination Function ===
backward_elimination_gam <- function(data, target_var, predictors, k = 5, familyXYZ = gaussian(), p_threshold = 0.05) {
  current_predictors <- predictors
  best_model <- NULL
  improved <- TRUE
  
  while (improved && length(current_predictors) > 0) {
    smooth_terms <- paste0("s(", current_predictors, ", k=", k, ")")
    formula_str <- paste(target_var, "~", paste(smooth_terms, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    
    # Fit model
    model <- tryCatch({
      gam(formula_obj, family = familyXYZ, method = "REML", data = data)
    }, error = function(e) {
      message("Model fitting failed: ", e$message)
      return(NULL)
    })
    
    if (is.null(model)) break
    
    model_summary <- tryCatch({
      summary(model)
    }, error = function(e) {
      message("Model summary failed: ", e$message)
      return(NULL)
    })
    
    # Bail out if summary failed
    if (is.null(model_summary)) break
    
    # Get p-values of smooth terms
    p_vals <- model_summary$s.pv
    if (any(is.na(p_vals))) break  # Stop if p-values are NA
    
    max_p <- max(p_vals)
    max_p_index <- which.max(p_vals)
    
    # Remove worst predictor if above threshold
    if (max_p > p_threshold) {
      message(sprintf("Removing predictor '%s' with p-value %.3f", current_predictors[max_p_index], max_p))
      current_predictors <- current_predictors[-max_p_index]
    } else {
      improved <- FALSE
      best_model <- model
    }
  }
  
  if (is.null(best_model)) best_model <- model  # Fallback
  
  return(list(model = best_model, predictors = current_predictors))
}


# === Example: Your predictor structure ===
# Define the predictor structure from the table
response_predictors <- list(
  "Age1_NMFS_spring_BTS" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO"),
    GBK = c("SST", "BT", "GSI", "AMO", "NAO", "WCR"),
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO")
  ),
  "Age1_NMFS_fall_BTS" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO"),
    GBK = c("SST", "BT", "GSI", "AMO", "NAO", "WCR"),
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO")
  ),
  "Age1_MADMF_spring_BTS" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO"),
    GBK = NULL,
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO")
  ),
  "rssb_Fall" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO"),
    GBK = c("SST", "BT", "GSI", "AMO", "NAO", "WCR"),
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO")
  ),
  "rssb_Spring" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO"),
    GBK = c("SST", "BT", "GSI", "AMO", "NAO", "WCR"),
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO")
  ),
  "Age1" = list(
    GOM = c("BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "Age2" = list(
    GOM = c("BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "Age3" = list(
    GOM = c("BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "Age4" = list(
    GOM = c("BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "Age5" = list(
    GOM = c("BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "Age6" = list(
    GOM = c("BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "Age7" = list(
    GOM = c("BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "Condition" = list(
    GOM = c("BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "spring_cog_depth" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("SST", "BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "fall_cog_depth" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("SST", "BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "spring_cog_lat" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("SST", "BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO", "SSB")
  ),
  "fall_cog_lat" = list(
    GOM = c("SST", "BT", "GSI", "AMO", "NAO", "SSB"),
    GBK = c("SST", "BT", "GSI", "AMO", "NAO", "WCR", "SSB"),
    SNEMA = c("SST", "BT", "GSI", "AMO", "NAO", "SSB")
  )
)

run_gams_for_winterflounder <- function(data, target_list, stocks, k = 5, familyXYZ = mgcv::tw()) {
  for (target_name in names(target_list)) {
    for (stock_name in stocks) {
      base_vars <- target_list[[target_name]][[stock_name]]
      if (is.null(base_vars)) next
      
      # Create predictors, including lags
      predictors <- unlist(lapply(base_vars, function(var) {
        if (var == "SSB") return("SSB")
        c(paste0(var, "_lag_6mo"), paste0(var, "_lag_12mo"))
      }))
      
      stock_data <- data %>% filter(stock == stock_name)
      
      message(sprintf("Running backward elimination for target: %s, stock: %s", target_name, stock_name))
      
      # Run backward elimination
      be_result <- backward_elimination_gam(
        data = stock_data,
        target_var = target_name,
        predictors = predictors,
        k = k,
        familyXYZ = familyXYZ
      )
      
      final_model <- be_result$model
      selected_predictors <- be_result$predictors
      
      # Safely summarize final model
      model_summary <- tryCatch({
        summary(final_model)
      }, error = function(e) {
        message("Model summary failed: ", e$message)
        return(NULL)
      })
      
      if (is.null(model_summary)) next  # Skip if summary failed
      
      # Safely extract p-values
      pvals_str <- tryCatch({
        if (!is.null(model_summary$s.pv)) {
          paste(round(model_summary$s.pv, 3), collapse = ", ")
        } else {
          NA_character_
        }
      }, error = function(e) {
        message("Failed to extract p-values: ", e$message)
        NA_character_
      })
      
      # AIC and deviance explained
      aic_val <- tryCatch(round(final_model$aic, 3), error = function(e) NA_real_)
      dev_expl <- tryCatch(round(model_summary$dev.expl, 3), error = function(e) NA_real_)
      
      message(sprintf("Finished %s - %s: AIC=%.3f, Deviance explained=%.3f, Predictors: %s, p-values: %s",
                      target_name, stock_name, aic_val, dev_expl,
                      paste(selected_predictors, collapse = ", "), pvals_str))
      
      # Save model and summary
      save_dir <- file.path("data/trial_results", paste0("winter_flounder_", stock_name, "_", target_name))
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      
      saveRDS(final_model, file = file.path(save_dir, paste0("final_model_", target_name, ".RDS")))
      
      summary_dt <- data.table(
        Target = target_name,
        Stock = stock_name,
        AIC = aic_val,
        Deviance_Explained = dev_expl,
        Predictors = paste(selected_predictors, collapse = ", "),
        Smooth_p_values = pvals_str
      )
      
      fwrite(summary_dt, file = file.path(save_dir, paste0("model_summary_", target_name, ".csv")))
    }
  }
}

run_gams_for_winterflounder(data = Edata, target_list = response_predictors, stocks = c("GOM", "GBK", "SNEMA"), k = 5)


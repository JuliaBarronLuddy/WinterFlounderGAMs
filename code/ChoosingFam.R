library(data.table)
library(mgcv)
library(dplyr)

# === Backward Elimination Function (unchanged) ===
backward_elimination_gam <- function(data, target_var, predictors, k = 5, familyXYZ = gaussian(), p_threshold = 0.05) {
  current_predictors <- predictors
  best_model <- NULL
  improved <- TRUE
  
  while (improved && length(current_predictors) > 0) {
    smooth_terms <- paste0("s(", current_predictors, ", k=", k, ")")
    formula_str <- paste(target_var, "~", paste(smooth_terms, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    
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
    
    if (is.null(model_summary)) break
    
    p_vals <- model_summary$s.pv
    if (any(is.na(p_vals))) break
    
    max_p <- max(p_vals)
    max_p_index <- which.max(p_vals)
    
    if (max_p > p_threshold) {
      message(sprintf("Removing predictor '%s' with p-value %.3f", current_predictors[max_p_index], max_p))
      current_predictors <- current_predictors[-max_p_index]
    } else {
      improved <- FALSE
      best_model <- model
    }
  }
  
  if (is.null(best_model)) best_model <- model
  return(list(model = best_model, predictors = current_predictors))
}


# === Helper: choose GAM family automatically ===
choose_family_for_target <- function(target_name, data_subset) {
  fam <- gaussian()  # default
  
  if (target_name %in% c("Age1", "Age1_NMFS_spring_BTS", "Age1_NMFS_fall_BTS",
                         "Age1_MADMF_spring_BTS", "rssb_Fall", "rssb_Spring")) {
    fam <- gaussian()
    
  } else if (target_name %in% c("Condition")) {
    fam <- gaussian()
    
  } else if (target_name %in% c("spring_cog_depth", "fall_cog_depth")) {
    y <- abs(data_subset[[target_name]])
    if (mean(y == 0, na.rm = TRUE) > 0.05) {
      fam <- mgcv::tw()   # Tweedie if zeros present
    } else if (mean(y, na.rm = TRUE) > 0 && sd(y, na.rm = TRUE)/mean(y, na.rm = TRUE) > 1.5) {
      fam <- Gamma(link = "log")   # highly skewed
    } else fam <- gaussian()
    
  } else if (grepl("WAA|Age2|Age3|Age4|Age5|Age6|Age7", target_name)) {
    fam <- mgcv::scat()  # robust scaled t, for bimodal/heavy tails
  }
  
  return(fam)
}


# === MAIN FUNCTION: now includes transformation + family + diagnostics ===
run_gams_for_winterflounder <- function(data, target_list, stocks, k = 5) {
  for (target_name in names(target_list)) {
    for (stock_name in stocks) {
      base_vars <- target_list[[target_name]][[stock_name]]
      if (is.null(base_vars)) next
      
      # === TRANSFORMATION SECTION ===
      data <- data %>%
        mutate("{target_name}" := case_when(
          target_name %in% c("Age1_NMFS_spring_BTS", "Age1_NMFS_fall_BTS", "Age1_MADMF_spring_BTS") ~ log(get(target_name)),
          target_name %in% c("rssb_Fall", "rssb_Spring") ~ log(get(target_name)),
          target_name %in% c("spring_cog_depth", "fall_cog_depth") ~ abs(get(target_name)),
          TRUE ~ get(target_name)
        ))
      # === END TRANSFORMATION SECTION ===
      
      # Define predictors with lag structure
      predictors <- unlist(lapply(base_vars, function(var) {
        if (var == "SSB") return("SSB")
        c(paste0(var, "_lag_6mo"), paste0(var, "_lag_12mo"))
      }))
      
      # Subset by stock
      stock_data <- data %>% filter(stock == stock_name)
      
      # === CHOOSE FAMILY BASED ON TARGET ===
      family_to_use <- choose_family_for_target(target_name, stock_data)
      
      message(sprintf("Running backward elimination for target: %s, stock: %s", target_name, stock_name))
      
      be_result <- backward_elimination_gam(
        data = stock_data,
        target_var = target_name,
        predictors = predictors,
        k = k,
        familyXYZ = family_to_use
      )
      
      final_model <- be_result$model
      selected_predictors <- be_result$predictors
      
      # === DIAGNOSTICS ===
      message("Running diagnostics...")
      try({
        gam.check(final_model)
      }, silent = TRUE)
      
      if (requireNamespace("DHARMa", quietly = TRUE)) {
        try({
          DHARMa::simulateResiduals(final_model, plot = TRUE)
        }, silent = TRUE)
      }
      
      # === SAFE SUMMARY EXTRACTION ===
      model_summary <- tryCatch({
        summary(final_model)
      }, error = function(e) {
        message("Model summary failed: ", e$message)
        return(NULL)
      })
      
      if (is.null(model_summary)) next
      
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
      
      aic_val <- tryCatch(round(final_model$aic, 3), error = function(e) NA_real_)
      dev_expl <- tryCatch(round(model_summary$dev.expl, 3), error = function(e) NA_real_)
      
      message(sprintf("Finished %s - %s: AIC=%.3f, Deviance explained=%.3f, Predictors: %s, p-values: %s",
                      target_name, stock_name, aic_val, dev_expl,
                      paste(selected_predictors, collapse = ", "), pvals_str))
      
      # === SAVE OUTPUTS ===
      save_dir <- file.path("data/famgam_results", paste0("winter_flounder_", stock_name, "_", target_name))
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      
      saveRDS(final_model, file = file.path(save_dir, paste0("final_model_", target_name, ".RDS")))
      
      summary_dt <- data.table(
        Target = target_name,
        Stock = stock_name,
        Family = as.character(family_to_use$family[1]),
        AIC = aic_val,
        Deviance_Explained = dev_expl,
        Predictors = paste(selected_predictors, collapse = ", "),
        Smooth_p_values = pvals_str
      )
      
      fwrite(summary_dt, file = file.path(save_dir, paste0("model_summary_", target_name, ".csv")))
    }
  }
  
  message("✅ All GAMs completed successfully.")
}

run_gams_for_winterflounder(Edata, response_predictors, c("GOM", "GBK", "SNEMA"), k = 5)



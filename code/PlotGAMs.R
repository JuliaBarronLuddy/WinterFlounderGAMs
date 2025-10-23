library(mgcv)
library(ggplot2)
library(gratia)
library(fs)

plot_all_gams_single_folder <- function(base_dir = "data/famgam_results",
                                        plot_dir = "data/famgam_results/plots_all",
                                        p_threshold = 0.05) {
  # Ensure output folder exists
  dir_create(plot_dir)
  
  # List all saved GAM RDS files
  model_files <- dir(base_dir, pattern = "final_model_.*\\.RDS$", recursive = TRUE, full.names = TRUE)
  if (length(model_files) == 0) {
    stop("No model files found in ", base_dir)
  }
  
  for (file_path in model_files) {
    # Read model safely
    model <- tryCatch(readRDS(file_path), error = function(e) NULL)
    if (is.null(model) || !"gam" %in% class(model)) next
    
    folder_name <- basename(dirname(file_path))
    target_name <- gsub("final_model_", "", gsub(".RDS", "", basename(file_path)))
    
    # Extract summary safely
    model_summary <- tryCatch(summary(model), error = function(e) NULL)
    if (is.null(model_summary) || !is.list(model_summary)) next
    
    # Deviance explained
    dev_expl <- tryCatch(round(model_summary$dev.expl * 100, 1), error = function(e) NA_real_)
    
    # Smooth terms and p-values
    if (!is.null(model_summary$s.pv) && length(model_summary$s.pv) > 0) {
      smooth_pvals <- model_summary$s.pv
      smooth_names <- rownames(model_summary$s.table)
      
      # Skip if all smooths are non-significant
      if (all(is.na(smooth_pvals)) || all(smooth_pvals > p_threshold)) {
        message(sprintf("Skipping %s: no significant smooths", target_name))
        next
      }
      
      pval_text <- paste0(smooth_names, ": p=", signif(smooth_pvals, 3), collapse = "; ")
    } else {
      message(sprintf("Skipping %s: no smooth terms", target_name))
      next
    }
    
    # Draw plot with gratia
    g <- tryCatch({
      draw(model, residuals = FALSE, rug = FALSE) +
        ggtitle(paste0(target_name, " (", folder_name, ")")) +
        labs(subtitle = paste0("Deviance explained: ", dev_expl, "%; Smooth p-values: ", pval_text)) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5)
        )
    }, error = function(e) NULL)
    
    if (!is.null(g)) {
      # Save plot
      out_file <- file.path(plot_dir, paste0("smooths_", target_name, "_clean.png"))
      ggsave(out_file, g, width = 10, height = 6, dpi = 300)
      message(sprintf("Saved plot for %s", target_name))
    }
  }
  
  message("All significant GAM plots processed and saved in ", plot_dir)
}

# Run the function
plot_all_gams_single_folder("data/famgam_results")






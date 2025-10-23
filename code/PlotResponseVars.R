library(ggplot2)
library(data.table)
library(dplyr)
library(here)

# Convert Edata to data.table
dt <- as.data.table(Edata)

# Define which combinations to plot
target_list <- response_predictors
stocks <- c("GOM", "GBK", "SNEMA")

# Set output folder path using `here`
output_folder <- here("figures/response_variable_distributions")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# Loop through each response variable and stock
for (response in names(target_list)) {
  for (stock in stocks) {
    
    # Skip if this stock has no predictors for the response variable
    predictors <- target_list[[response]][[stock]]
    if (is.null(predictors)) next
    
    # Skip if response variable doesn't exist
    if (!(response %in% names(dt))) {
      message(paste("Skipping", response, "– column not in data"))
      next
    }
    
    # Filter data for correct Stock and non-missing response
    subset_data <- dt[get("stock") == stock & !is.na(get(response))]
    
    # Skip if not enough data
    if (nrow(subset_data) < 5) {
      message(paste("Skipping", response, "for", stock, "– too few data points"))
      next
    }
    
    # Basic histogram
    p <- ggplot(subset_data, aes_string(x = response)) +
      geom_histogram(bins = 30, fill = "#2C77B3", color = "white") +
      theme_minimal() +
      labs(
        title = paste0("Distribution of ", response, " (", stock, ")"),
        x = response,
        y = "Count"
      )
    
    # Save histogram
    ggsave(
      filename = file.path(output_folder, paste0("hist_", response, "_", stock, ".png")),
      plot = p, width = 6, height = 4
    )
    
    # Optional: log-scale histogram — only if all values are > 0
    if (all(subset_data[[response]] > 0, na.rm = TRUE)) {
      p_log <- ggplot(subset_data, aes_string(x = response)) +
        geom_histogram(bins = 30, fill = "#D95F02", color = "white") +
        theme_minimal() +
        scale_x_log10() +
        labs(
          title = paste0("Log-scale Distribution of ", response, " (", stock, ")"),
          x = paste0("log10(", response, ")"),
          y = "Count"
        )
      
      ggsave(
        filename = file.path(output_folder, paste0("hist_log_", response, "_", stock, ".png")),
        plot = p_log, width = 6, height = 4
      )
    }
  }
}


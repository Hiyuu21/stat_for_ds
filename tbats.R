# ===================================================================
# Libraries
# ===================================================================
library(tidyverse)
library(forecast)
library(ggplot2)

# ===================================================================
# IEEE Plot Settings
# ===================================================================
IEEE_W_SINGLE <- 3.5
IEEE_W_DOUBLE <- 7.16
IEEE_H        <- 2.8
IEEE_DPI      <- 300

theme_ieee <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(family = "serif", size = base_size),
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size),
      strip.text       = element_text(size = base_size, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

save_ieee <- function(fname, plot, width = IEEE_W_DOUBLE, height = IEEE_H) {
  ggsave(
    filename = paste0(fname, ".png"),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = IEEE_DPI,
    units    = "in"
  )
  cat("Saved:", paste0(fname, ".png"), "\n")
}

# ===================================================================
# TBATS Function (train/test validation + diagnostics + future forecast)
# ===================================================================
tbats_analysis <- function(cat_name, data_list, future_h=12) {
  cat("\n======================================================\n")
  cat("ANALYZING CATEGORY:", cat_name, "\n")
  cat("======================================================\n")
  
  # Split dataset
  train_ts <- data_list[[cat_name]]$train
  test_ts  <- data_list[[cat_name]]$test
  full_ts  <- data_list[[cat_name]]$full
  
  # ---- 1. Fit TBATS model ----
  cat("\n1. Fitting TBATS Model\n")
  model <- tbats(train_ts)
  print(model)
  
  # ---- 2. Decomposition ----
  cat("\n2. Inspecting Decomposition\n")
  decomp <- tbats.components(model)
  p_decomp <- autoplot(decomp) +
    theme_ieee() +
    labs(title = paste("TBATS Decomposition:", cat_name))
  print(p_decomp)
  save_ieee(paste0("tbats_decomp_", gsub(" ", "_", tolower(cat_name))), p_decomp)
  
  # ---- 3. Residual Diagnostics ----
  cat("\n3. Residual Diagnostics\n")
  checkresiduals(model, plot = TRUE)
  
  # ---- 4. Forecast on Test Set ----
  cat("\n4. Forecasting on Test Set\n")
  fc <- forecast(model, h = length(test_ts))
  eval <- accuracy(fc, test_ts)
  print(eval)
  
  p_test <- autoplot(fc) +
    autolayer(test_ts, series = "Actual Test Data", linewidth = 0.7) +
    theme_ieee() +
    labs(title = paste("TBATS Forecast vs Actual (Test Set):", cat_name),
         x = "Year", y = "Sales Quantity (units)")
  print(p_test)
  save_ieee(paste0("tbats_testforecast_", gsub(" ", "_", tolower(cat_name))), p_test)
  
  # ---- 5. Refit on Full Data + Future Forecast ----
  cat("\n5. Future Forecast\n")
  final_model <- tbats(full_ts)
  future_fc   <- forecast(final_model, h = future_h)
  print(future_fc)
  
  p_future <- autoplot(future_fc) +
    theme_ieee() +
    labs(title = paste("TBATS Future Forecast (next", future_h, "months):", cat_name),
         x = "Year", y = "Sales Quantity (units)")
  print(p_future)
  save_ieee(paste0("tbats_futureforecast_", gsub(" ", "_", tolower(cat_name))), p_future)
  
  return(list(
    Model = model,
    FinalModel = final_model,
    TestForecast = fc,
    FutureForecast = future_fc,
    Accuracy = eval,
    TestPlot = p_test,
    FuturePlot = p_future
  ))
}

# ===================================================================
# Function Calls
# ===================================================================
results_2w <- tbats_analysis("2-Wheelers", ts_data)
results_3w <- tbats_analysis("3-Wheelers", ts_data)
results_4w <- tbats_analysis("4-Wheelers", ts_data)

# Cross-category comparison
comparison <- bind_rows(
  as.data.frame(t(results_2w$Accuracy["Test set", c("RMSE","MAE","MASE","Theil's U","ME")])) %>%
    mutate(Category = "2-Wheelers"),
  as.data.frame(t(results_3w$Accuracy["Test set", c("RMSE","MAE","MASE","Theil's U","ME")])) %>%
    mutate(Category = "3-Wheelers"),
  as.data.frame(t(results_4w$Accuracy["Test set", c("RMSE","MAE","MASE","Theil's U","ME")])) %>%
    mutate(Category = "4-Wheelers")
) %>% relocate(Category)

cat("\n=================================================\nCross-Category TBATS Model Performance Comparison\n=================================================\n")
print(comparison)


# ============================================================================= 
# 03_tbats.R
# =============================================================================
source("preprocess.R")

tbats_analysis <- function(cat_name, data_list, future_h = 12) {
  cat("\n======================================================\n")
  cat("ANALYZING CATEGORY:", cat_name, "\n")
  cat("======================================================\n")
  
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
  
  fname_decomp <- paste0("tbats_01_decomposition_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_decomp, p_decomp)
  
  # ---- 3. Residual Diagnostics ----
  cat("\n3. Residual Diagnostics\n")
  print(checkresiduals(model, plot = TRUE))
  
  fname_resid <- paste0("tbats_02_residuals_",
                        gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_resid, expression({
    checkresiduals(model, plot = TRUE)
  }), height = 3.5)
  
  # ---- 4. Forecast on Test Set ----
  cat("\n4. Forecasting on Test Set\n")
  fc   <- forecast(model, h = length(test_ts))
  eval <- compute_metrics(
    train_y   = as.numeric(train_ts),
    actual    = as.numeric(test_ts),
    predicted = as.numeric(fc$mean)
  )
  print(eval)
  
  p_test <- autoplot(fc) +
    autolayer(test_ts, series = "Actual Test Data", linewidth = 0.7) +
    theme_ieee() +
    labs(title = paste("TBATS Forecast vs Actual (Test Set):", cat_name),
         x = "Year", y = "Sales Quantity (units)")
  print(p_test)
  
  fname_test <- paste0("tbats_03_testforecast_",
                       gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_test, p_test)
  
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
  
  fname_future <- paste0("tbats_04_futureforecast_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_future, p_future)
  
  return(list(
    Model = model, FinalModel = final_model, TestForecast = fc,
    FutureForecast = future_fc, Accuracy = eval,
    TestPlot = p_test, FuturePlot = p_future
  ))
}

sink("TBATS_Detailed_Log.txt", split = TRUE)

tbats_2w <- tbats_analysis("2-Wheelers", ts_data)
tbats_3w <- tbats_analysis("3-Wheelers", ts_data)
tbats_4w <- tbats_analysis("4-Wheelers", ts_data)

comparison_tbats <- bind_rows(
  tbats_2w$Accuracy |> mutate(Category = "2-Wheelers"),
  tbats_3w$Accuracy |> mutate(Category = "3-Wheelers"),
  tbats_4w$Accuracy |> mutate(Category = "4-Wheelers")
) |> relocate(Category) |> mutate(Model = "TBATS")

cat("\n=================================================\nCross-Category TBATS Model Performance Comparison\n=================================================\n")
print(comparison_tbats)

saveRDS(comparison_tbats, "comparison_tbats.rds")

sink()
cat("\n--- TBATS details saved to 'TBATS_Detailed_Log.txt' ---\n")
# ============================================================================= 
# 02_stlets.R
# =============================================================================
source("preprocess.R")

stl_ets <- function(cat_name, data_list, future_h = 12) {
  cat(paste0("\n======================================================\n"))
  cat(paste0("ANALYZING CATEGORY: ", cat_name, "\n"))
  cat(paste0("======================================================\n"))
  
  train_ts <- data_list[[cat_name]]$train
  test_ts  <- data_list[[cat_name]]$test
  full_ts  <- data_list[[cat_name]]$full
  
  # ---- 1. STL Decomposition ----
  cat("\n1. Performing STL Decomposition on Training Data...\n")
  stl_decomp <- stl(train_ts, s.window = "periodic", robust = TRUE)
  print(summary(stl_decomp))
  
  fname_stl <- paste0("stlets_01_decomposition_",
                      gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_stl, expression({
    plot(stl_decomp, main = paste("STL Decomposition -", cat_name))
  }), height = 4.0)
  
  # ---- 2. Model Fitting (stlm with ETS) ----
  cat("\n2. Fitting STL + ETS Model...\n")
  model <- stlm(train_ts, s.window = "periodic", method = "ets",
                robust = TRUE, damped = TRUE)
  
  cat("\nFitted Model Summary:\n")
  if (!is.null(model$model) && inherits(model$model, "ets")) {
    print(model$model)
  } else {
    print(model)
  }
  
  # ---- 3. Residual Diagnostics ----
  cat("\n3. Residual Diagnostics (Ljung-Box test):\n")
  print(checkresiduals(model, plot = TRUE))
  
  fname_resid <- paste0("stlets_02_residuals_",
                        gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_resid, expression({
    checkresiduals(model, plot = TRUE)
  }), height = 3.5)
  
  # ---- 4. Forecast on held-out test set + evaluation ----
  fc   <- forecast(model, h = length(test_ts))
  eval <- compute_metrics(
    train_y   = as.numeric(train_ts),
    actual    = as.numeric(test_ts),
    predicted = as.numeric(fc$mean)
  )
  
  cat("\n4. Forecasting and Evaluation on Test Set\n")
  cat("\n--- Full Accuracy Matrix ---\n")
  print(eval)
  
  cat("\n--- Target Test Metrics (RMSE, MAE, MASE, Theil's U, ME) ---\n")
  print(eval)
  
  p_test <- autoplot(fc) +
    autolayer(test_ts, series = "Actual Test Data", linewidth = 0.7) +
    theme_ieee() +
    labs(title = paste("STL+ETS Forecast vs Actual (Test Set):", cat_name),
         x = "Year", y = "Sales Quantity (units)")
  print(p_test)
  
  fname_test <- paste0("stlets_03_testforecast_",
                       gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  # ---- 5. Refit on FULL data + future forecast ----
  cat("\n5. Refitting STL+ETS on Full Dataset and Forecasting Future...\n")
  final_model <- stlm(full_ts, s.window = "periodic", method = "ets", robust = TRUE)
  future_fc <- forecast(final_model, h = future_h)
  
  cat("\nFuture Forecast (next", future_h, "months beyond dataset):\n")
  print(future_fc)
  
  p_future <- autoplot(future_fc) +
    theme_ieee() +
    labs(title = paste("STL+ETS Future Forecast (next", future_h, "months):", cat_name),
         x = "Year", y = "Sales Quantity (units)")
  print(p_future)
  
  fname_future <- paste0("stlets_04_futureforecast_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  return(list(
    Model = model, FinalModel = final_model, TestForecast = fc,
    FutureForecast = future_fc, Accuracy = eval,
    TestPlot = p_test, FuturePlot = p_future
  ))
}

sink("STLETS_Detailed_Log.txt", split = TRUE)

ets_2w <- stl_ets("2-Wheelers", ts_data)
ets_3w <- stl_ets("3-Wheelers", ts_data)
ets_4w <- stl_ets("4-Wheelers", ts_data)

comparison_ets <- bind_rows(
  ets_2w$Accuracy |> mutate(Category = "2-Wheelers"),
  ets_3w$Accuracy |> mutate(Category = "3-Wheelers"),
  ets_4w$Accuracy |> mutate(Category = "4-Wheelers")
) |> relocate(Category) |> mutate(Model = "STL+ETS")

cat("\nCross-Category Model Performance Comparison (STL + ETS)\n")
print(comparison_ets)

saveRDS(comparison_ets, "comparison_ets.rds")

sink()
cat("\n--- STL+ETS details saved to 'STLETS_Detailed_Log.txt' ---\n")
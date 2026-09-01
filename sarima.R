# ============================================================================= 
# 01_sarima.R
# =============================================================================

source("preprocess.R")

sarima <- function(cat_name, data_list, future_h = 12) {
  cat(paste0("\n======================================================\n"))
  cat(paste0("ANALYZING CATEGORY: ", cat_name, "\n"))
  cat(paste0("======================================================\n"))
  
  train_ts <- data_list[[cat_name]]$train
  test_ts  <- data_list[[cat_name]]$test
  full_ts  <- data_list[[cat_name]]$full
  
  # ---- 1a. Stationarity checks ----
  adf_res  <- adf.test(train_ts)
  kpss_res <- kpss.test(train_ts)
  
  cat("\n1. Stationarity Tests on Raw Training Data:\n")
  cat("   ADF p-value :", round(adf_res$p.value, 4),
      "-> ", ifelse(adf_res$p.value < 0.05, "Stationary", "Non-stationary"), "\n")
  cat("   KPSS p-value:", round(kpss_res$p.value, 4),
      "-> ", ifelse(kpss_res$p.value < 0.05, "Non-stationary", "Stationary"), "\n")
  
  d_order <- ndiffs(train_ts)
  D_order <- nsdiffs(train_ts)
  cat("   Suggested non-seasonal differencing (d):", d_order, "\n")
  cat("   Suggested seasonal differencing (D)    :", D_order, "\n")
  
  # ---- 1b. ACF/PACF on RAW series ----
  par(mfrow = c(1, 2))
  acf(train_ts,  main = paste("ACF (raw) -",  cat_name))
  pacf(train_ts, main = paste("PACF (raw) -", cat_name))
  par(mfrow = c(1, 1))
  
  fname_acfpacf_raw <- paste0("sarima_00_acfpacf_raw_",
                              gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_acfpacf_raw, expression({
    par(mfrow = c(1, 2))
    acf(train_ts,  main = paste("ACF (raw) -",  cat_name))
    pacf(train_ts, main = paste("PACF (raw) -", cat_name))
    par(mfrow = c(1, 1))
  }))
  
  # ---- 2. ACF/PACF on differenced series ----
  diff_ts <- train_ts
  if (D_order > 0) diff_ts <- diff(diff_ts, lag = 12, differences = D_order)
  if (d_order > 0) diff_ts <- diff(diff_ts, differences = d_order)
  
  par(mfrow = c(1, 2))
  acf(diff_ts,  main = paste("ACF (differenced) -",  cat_name))
  pacf(diff_ts, main = paste("PACF (differenced) -", cat_name))
  par(mfrow = c(1, 1))
  
  fname_acfpacf <- paste0("sarima_01_acfpacf_",
                          gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_acfpacf, expression({
    par(mfrow = c(1, 2))
    acf(diff_ts,  main = paste("ACF (differenced) -",  cat_name))
    pacf(diff_ts, main = paste("PACF (differenced) -", cat_name))
    par(mfrow = c(1, 1))
  }))
  
  # ---- 3. Model fitting (auto.arima, AICc-driven search) ----
  cat("\n2. Fitting SARIMA Model...\n")
  model <- auto.arima(train_ts, seasonal = TRUE, stepwise = FALSE,
                      approximation = FALSE, trace = TRUE, ic = "aicc")
  cat("\nSelected order:", paste(arimaorder(model), collapse = ","), "\n")
  print(summary(model))
  
  # ---- 4. Residual diagnostics ----
  cat("\n3. Residual Diagnostics (Ljung-Box test):\n")
  print(checkresiduals(model, plot = TRUE))
  
  fname_resid <- paste0("sarima_02_residuals_",
                        gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_resid, expression({
    checkresiduals(model, plot = TRUE)
  }), height = 3.5)
  
  # ---- 5. Forecast on held-out test set + evaluation ----
  fc   <- forecast(model, h = length(test_ts))
  eval <- compute_metrics(
    train_y   = as.numeric(train_ts),
    actual    = as.numeric(test_ts),
    predicted = as.numeric(fc$mean)
  )
  
  cat("\n4. Forecasting and Evaluation on Test Set\n")
  cat("\n--- Full Accuracy Matrix ---\n")
  print(eval)
  
  p_test <- autoplot(fc) +
    autolayer(test_ts, series = "Actual Test Data", linewidth = 0.7) +
    theme_ieee() +
    labs(title = paste("SARIMA Forecast vs Actual (Test Set):", cat_name),
         x = "Year", y = "Sales Quantity (units)")
  print(p_test)
  
  fname_test <- paste0("sarima_03_testforecast_",
                       gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  # ---- 6. Refit on FULL data + forecast forward ----
  ord <- arimaorder(model)
  final_model <- Arima(full_ts, order = ord[1:3],
                       seasonal = list(order = ord[4:6], period = 12))
  future_fc <- forecast(final_model, h = future_h)
  
  cat("\n5. Future Forecast (next", future_h, "months beyond dataset):\n")
  print(future_fc)
  
  p_future <- autoplot(future_fc) +
    theme_ieee() +
    labs(title = paste("SARIMA Future Forecast (next", future_h, "months):", cat_name),
         x = "Year", y = "Sales Quantity (units)")
  print(p_future)
  
  fname_future <- paste0("sarima_04_futureforecast_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  return(list(
    Model = model, FinalModel = final_model, TestForecast = fc,
    FutureForecast = future_fc, Accuracy = eval,
    TestPlot = p_test, FuturePlot = p_future
  ))
}

sink("SARIMA_Detailed_Log.txt", split = TRUE)

sarima_2w <- sarima("2-Wheelers", ts_data)
sarima_3w <- sarima("3-Wheelers", ts_data)
sarima_4w <- sarima("4-Wheelers", ts_data)

comparison_sarima <- bind_rows(
  sarima_2w$Accuracy |> mutate(Category = "2-Wheelers"),
  sarima_3w$Accuracy |> mutate(Category = "3-Wheelers"),
  sarima_4w$Accuracy |> mutate(Category = "4-Wheelers")
) |> relocate(Category) |> mutate(Model = "SARIMA")

cat("\nCross-Category Model Performance Comparison (SARIMA)\n")
print(comparison_sarima)

saveRDS(comparison_sarima, "comparison_sarima.rds")

sink()
cat("\n--- SARIMA details saved to 'SARIMA_Detailed_Log.txt' ---\n")
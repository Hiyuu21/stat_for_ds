library(tidyverse)
library(tseries)
library(forecast)

# -----------------------------
# Data Loading & Preparation
# -----------------------------
# Load the preprocessed time-series objects
ts_data <- readRDS("ts_data.rds")

# -----------------------------------------------------------------------
# SARIMA Function (train/test validation + diagnostics + future forecast)
# -----------------------------------------------------------------------

sarima <- function(cat_name, data_list, future_h=12) {
  cat(paste0("\n======================================================\n"))
  cat(paste0("ANALYZING CATEGORY: ", cat_name, "\n"))
  cat(paste0("======================================================\n"))
  
  # filter and split
  train_ts <- data_list[[cat_name]]$train
  test_ts  <- data_list[[cat_name]]$test
  full_ts  <- data_list[[cat_name]]$full
  
  # ---- 1. Stationarity checks ----
  # ADF: null hypothesis = non-stationary (unit root present)
  # KPSS: null hypothesis = stationary (opposite direction)
  # Using both protects against relying on a single test's assumptions.
  adf_res <- adf.test(train_ts)
  kpss_res <- kpss.test(train_ts)
  
  cat("\n1. Stationarity Tests on Raw Training Data:\n")
  cat("   ADF p-value :", round(adf_res$p.value, 4),
      "-> ", ifelse(adf_res$p.value < 0.05, "Stationary", "Non-stationary"), "\n")
  cat("   KPSS p-value:", round(kpss_res$p.value, 4),
      "-> ", ifelse(kpss_res$p.value < 0.05, "Non-stationary", "Stationary"), "\n")
  
  d_order  <- ndiffs(train_ts)
  D_order  <- nsdiffs(train_ts)
  cat("   Suggested non-seasonal differencing (d):", d_order, "\n")
  cat("   Suggested seasonal differencing (D)    :", D_order, "\n")
  
  
  # ---- 2. ACF/PACF for order identification ----
  # Plotted on the (seasonally) differenced series
  # so the correlogram is interpretable 
  diff_ts <- train_ts
  if (D_order > 0) diff_ts <- diff(diff_ts, lag = 12, differences = D_order)
  if (d_order > 0) diff_ts <- diff(diff_ts, differences = d_order)
  
  par(mfrow = c(1,2))
  acf(diff_ts, main = paste("ACF (differenced) -", cat_name))
  pacf(diff_ts, main = paste("PACF (differenced) -", cat_name))
  par(mfrow = c(1,1))
  
  # ---- 3. Model fitting (auto.arima, AICc-driven search) ----
  cat("\n2. Fitting SARIMA Model...\n")
  model <- auto.arima(
    train_ts, seasonal = TRUE, 
    stepwise = FALSE, 
    approximation = FALSE,
    trace = TRUE,
    ic = "aicc"
  )
  cat("\nSelected order:", paste(arimaorder(model), collapse = ","), "\n")
  print(summary(model))
  
  
  # ---- 4. Residual diagnostics ----
  # Ljung-Box test H0: residuals are white noise (no leftover autocorrelation).
  # A high p-value (>0.05) means the model has captured the structure in
  # the data and what's left is just noise - this is what "good fit" means
  # beyond just low error numbers.
  
  cat("\n3. Residual Diagnostics (Ljung-Box test):\n")
  print(checkresiduals(model, plot = TRUE))
  
  
  # ---- 5. Forecast on held-out test set + evaluation ----  
  fc <-  forecast(model, h=length(test_ts))
  eval <- accuracy(fc, test_ts)
  
  cat("\n4. Forecasting and Evaluation on Test Set\n")
  cat("\n--- Full Accuracy Matrix ---\n")
  print(eval)
  
  cat("\n--- Target Test Metrics (MAE, RMSE, MAPE, Theil's U) ---\n")
  if ("Test set" %in% rownames(eval)) {
    print(eval["Test set", c("RMSE", "MAE", "MAPE","Theil's U")])
  } else {
    cat("WARNING: 'Test set' row is missing from the accuracy matrix!\n")
  }
  
  p_test <- autoplot(fc) +
    autolayer(test_ts, series = "Actual Test Data", linewidth = 1)+
    theme_minimal()+
    labs(
      title=paste("SARIMA Forecast vs Actual (Test Set):",cat_name),
      x="Year",
      y="Sales Quantity"
    ) + 
    theme(legend.position = "bottom")
  
  print(p_test)
  
  
  # ---- 6. Refit on FULL data + forecast forward ----
  # The train/test split above exists purely to validate the model
  # honestly (test data it never saw). Once validated, we refit the
  # *same order* on ALL available data (train+test) so the forward
  # forecast uses every observation we have, then project forward
  
  ord <- arimaorder(model)
  final_model <- Arima(
    full_ts,
    order = ord[1:3],
    seasonal = list(order = ord[4:6], period = 12)
  )
  
  future_fc <- forecast(final_model, h = future_h)
  
  cat("\n5. Future Forecast (next", future_h, "months beyond dataset):\n")
  print(future_fc)
  
  p_future <- autoplot(future_fc) +
    theme_minimal() + 
    labs(
      title = paste("SARIMA Future Forecast (next", future_h, "months):", cat_name),
      x="Year", y="Sales Quantity"
    ) +
    theme(legend.position="bottom")
  print(p_future)
  
  
  
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

# -----------------------------
# Function Calling
# -----------------------------

# record results into the log file
sink("SARIMA_Detailed_Log.txt", split = TRUE)

results_2w <- sarima("2-Wheelers", ts_data)
results_3w <- sarima("3-Wheelers", ts_data)
results_4w <- sarima("4-Wheelers", ts_data)

# -------------------------------
# Cross-category comparison table
# -------------------------------
comparison <- bind_rows(
  as.data.frame(t(results_2w$Accuracy["Test set", c("RMSE", "MAE", "MAPE", "Theil's U")])) |>
    mutate(Category = "2-Wheelers"),
  as.data.frame(t(results_3w$Accuracy["Test set", c("RMSE", "MAE", "MAPE", "Theil's U")])) |>
    mutate(Category = "3-Wheelers"),
  as.data.frame(t(results_4w$Accuracy["Test set", c("RMSE", "MAE", "MAPE", "Theil's U")])) |>
    mutate(Category = "4-Wheelers")
) |> relocate(Category)

cat("\nCross-Category Model Performance Comparison\n")
print(comparison)

sink()

cat("\n--- All details have been saved to 'SARIMA_Detailed_Log.txt' ---\n")

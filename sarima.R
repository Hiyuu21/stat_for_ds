library(tidyverse)
library(tseries)
library(forecast)

# =============================================================================
# IEEE PLOT SETTINGS
# =============================================================================
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

# For ggplot objects
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

# For base R graphics (ACF/PACF, checkresiduals)
save_ieee_base <- function(fname, expr, width = IEEE_W_DOUBLE, height = IEEE_H,
                           env = parent.frame()) {
  png(
    filename  = paste0(fname, ".png"),
    width     = width,
    height    = height,
    units     = "in",
    res       = IEEE_DPI,
    pointsize = 7       # matches base_size = 7 in theme_ieee()
  )
  eval(expr, envir = env)  # evaluate in caller's env so local vars are visible
  dev.off()
  cat("Saved:", paste0(fname, ".png"), "\n")
}

# -----------------------------
# Data Loading & Preparation
# -----------------------------
# Load the preprocessed time-series objects
#ts_data <- readRDS("ts_data.rds") use this if github unavailable

url <- gzcon(url("https://raw.githubusercontent.com/Hiyuu21/stat_for_ds/main/ts_data.rds"))
ts_data <- readRDS(url)

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
  
  # Save ACF/PACF (base R graphic)
  fname_acfpacf <- paste0(
    "sarima_01_acfpacf_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee_base(fname_acfpacf, expression({
    par(mfrow = c(1, 2))
    acf(diff_ts,  main = paste("ACF (differenced) -",  cat_name))
    pacf(diff_ts, main = paste("PACF (differenced) -", cat_name))
    par(mfrow = c(1, 1))
  }))
  
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
  
  # Save residual diagnostics (base R graphic — checkresiduals uses base plots)
  fname_resid <- paste0(
    "sarima_02_residuals_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee_base(fname_resid, expression({
    checkresiduals(model, plot = TRUE)
  }), height = 3.5)  # slightly taller to fit the 3-panel layout
  
  
  # ---- 5. Forecast on held-out test set + evaluation ----  
  fc <-  forecast(model, h=length(test_ts))
  eval <- accuracy(fc, test_ts)
  
  cat("\n4. Forecasting and Evaluation on Test Set\n")
  cat("\n--- Full Accuracy Matrix ---\n")
  print(eval)
  
  cat("\n--- Target Test Metrics (MAE, RMSE, MASE, Theil's U, ME) ---\n")
  if ("Test set" %in% rownames(eval)) {
    print(eval["Test set", c("RMSE", "MAE", "MASE","Theil's U","ME")])
  } else {
    cat("WARNING: 'Test set' row is missing from the accuracy matrix!\n")
  }
  
  p_test <- autoplot(fc) +
    autolayer(test_ts, series = "Actual Test Data", linewidth = 0.7) +
    theme_ieee() +
    labs(
      title = paste("SARIMA Forecast vs Actual (Test Set):", cat_name),
      x     = "Year",
      y     = "Sales Quantity (units)"
    )
  
  print(p_test)
  
  fname_test <- paste0(
    "sarima_03_testforecast_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  
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
    theme_ieee() +
    labs(
      title = paste("SARIMA Future Forecast (next", future_h, "months):", cat_name),
      x     = "Year",
      y     = "Sales Quantity (units)"
    )
  
  print(p_future)
  
  fname_future <- paste0(
    "sarima_04_futureforecast_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  
  
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
  as.data.frame(t(results_2w$Accuracy["Test set", c("RMSE", "MAE", "MASE", "Theil's U","ME")])) |>
    mutate(Category = "2-Wheelers"),
  as.data.frame(t(results_3w$Accuracy["Test set", c("RMSE", "MAE", "MASE", "Theil's U","ME")])) |>
    mutate(Category = "3-Wheelers"),
  as.data.frame(t(results_4w$Accuracy["Test set", c("RMSE", "MAE", "MASE", "Theil's U","ME")])) |>
    mutate(Category = "4-Wheelers")
) |> relocate(Category)

cat("\nCross-Category Model Performance Comparison\n")
print(comparison)

sink()

cat("\n--- All details have been saved to 'SARIMA_Detailed_Log.txt' ---\n")

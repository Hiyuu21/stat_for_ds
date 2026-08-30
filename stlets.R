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

# For base R graphics (Decomposition, ACF/PACF, checkresiduals)
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
# ts_data <- readRDS("ts_data.rds") # use this if github unavailable

url <- gzcon(url("https://raw.githubusercontent.com/Hiyuu21/stat_for_ds/main/ts_data.rds"))
ts_data <- readRDS(url)

# -----------------------------------------------------------------------
# STL + ETS Function (train/test validation + diagnostics + future forecast)
# -----------------------------------------------------------------------

stl_ets <- function(cat_name, data_list, future_h = 12) {
  cat(paste0("\n======================================================\n"))
  cat(paste0("ANALYZING CATEGORY: ", cat_name, "\n"))
  cat(paste0("======================================================\n"))
  
  # Filter and split
  train_ts <- data_list[[cat_name]]$train
  test_ts  <- data_list[[cat_name]]$test
  full_ts  <- data_list[[cat_name]]$full
  
  # ---- 1. STL Decomposition ----
  cat("\n1. Performing STL Decomposition on Training Data...\n")
  stl_decomp <- stl(train_ts, s.window = "periodic", robust = TRUE)
  print(summary(stl_decomp))
  
  # Plot & Save STL Decomposition
  fname_stl <- paste0(
    "stlets_01_decomposition_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee_base(fname_stl, expression({
    plot(stl_decomp, main = paste("STL Decomposition -", cat_name))
  }), height = 4.0)
  
  # ---- 2. Model Fitting (stlm with ETS) ----
  cat("\n2. Fitting STL + ETS Model...\n")
  model <- stlm(
    train_ts,
    s.window = "periodic",
    method   = "ets",
    robust   = TRUE, 
    damped   = TRUE
  )
  
  cat("\nFitted Model Summary:\n")
  # Safely print the ETS structure without triggering atomic vector errors
  if (!is.null(model$model) && inherits(model$model, "ets")) {
    print(model$model)
  } else {
    print(model)
  }
  
  # ---- 3. Residual Diagnostics ----
  cat("\n3. Residual Diagnostics (Ljung-Box test):\n")
  print(checkresiduals(model, plot = TRUE))
  
  fname_resid <- paste0(
    "stlets_02_residuals_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee_base(fname_resid, expression({
    checkresiduals(model, plot = TRUE)
  }), height = 3.5)
  
  # ---- 4. Forecast on held-out test set + evaluation ----  
  fc <- forecast(model, h = length(test_ts))
  eval <- accuracy(fc, test_ts)
  
  cat("\n4. Forecasting and Evaluation on Test Set\n")
  cat("\n--- Full Accuracy Matrix ---\n")
  print(eval)
  
  cat("\n--- Target Test Metrics (MAE, RMSE, MASE, Theil's U, ME) ---\n")
  if ("Test set" %in% rownames(eval)) {
    # Guard against missing column names in edge cases
    available_cols <- intersect(c("RMSE", "MAE", "MASE", "Theil's U", "ME"), colnames(eval))
    print(eval["Test set", available_cols])
  } else {
    cat("WARNING: 'Test set' row is missing from the accuracy matrix!\n")
  }
  
  p_test <- autoplot(fc) +
    autolayer(test_ts, series = "Actual Test Data", linewidth = 0.7) +
    theme_ieee() +
    labs(
      title = paste("STL+ETS Forecast vs Actual (Test Set):", cat_name),
      x     = "Year",
      y     = "Sales Quantity (units)"
    )
  
  print(p_test)
  
  fname_test <- paste0(
    "stlets_03_testforecast_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  # ---- 5. Refit on FULL data + forecast forward ----
  cat("\n5. Refitting STL+ETS on Full Dataset and Forecasting Future...\n")
  final_model <- stlm(
    full_ts,
    s.window = "periodic",
    method   = "ets",
    robust   = TRUE
  )
  
  future_fc <- forecast(final_model, h = future_h)
  
  cat("\nFuture Forecast (next", future_h, "months beyond dataset):\n")
  print(future_fc)
  
  p_future <- autoplot(future_fc) +
    theme_ieee() +
    labs(
      title = paste("STL+ETS Future Forecast (next", future_h, "months):", cat_name),
      x     = "Year",
      y     = "Sales Quantity (units)"
    )
  
  print(p_future)
  
  fname_future <- paste0(
    "stlets_04_futureforecast_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  return(list(
    Model          = model,
    FinalModel     = final_model,
    TestForecast   = fc,
    FutureForecast = future_fc,
    Accuracy       = eval,
    TestPlot       = p_test,
    FuturePlot     = p_future
  ))
}

# -----------------------------
# Function Calling
# -----------------------------

sink("STLETS_Detailed_Log.txt", split = TRUE)

results_2w <- stl_ets("2-Wheelers", ts_data)
results_3w <- stl_ets("3-Wheelers", ts_data)
results_4w <- stl_ets("4-Wheelers", ts_data)

# -------------------------------
# Cross-category comparison table
# -------------------------------
extract_metrics <- function(res_obj, cat_name) {
  acc <- res_obj$Accuracy
  if ("Test set" %in% rownames(acc)) {
    cols <- intersect(c("RMSE", "MAE", "MASE", "Theil's U", "ME"), colnames(acc))
    df <- as.data.frame(t(acc["Test set", cols]))
    df$Category <- cat_name
    return(df)
  }
  return(NULL)
}

comparison <- bind_rows(
  extract_metrics(results_2w, "2-Wheelers"),
  extract_metrics(results_3w, "3-Wheelers"),
  extract_metrics(results_4w, "4-Wheelers")
) |> relocate(Category)

cat("\nCross-Category Model Performance Comparison (STL + ETS)\n")
print(comparison)

sink()

cat("\n--- All details have been saved to 'STLETS_Detailed_Log.txt' ---\n")
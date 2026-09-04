library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(prophet)
library(ggplot2)

# =============================================================================
# --- SECTION 1: IEEE PLOT SETTINGS & SHARED HELPERS ---
# =============================================================================
# These standardise the output format of all plots for your final report.
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
  ggsave(filename = paste0(fname, ".png"), plot = plot, width = width, height = height, dpi = IEEE_DPI, units = "in")
  cat("Saved:", paste0(fname, ".png"), "\n")
}

save_ieee_base <- function(fname, expr, width = IEEE_W_DOUBLE, height = IEEE_H, env = parent.frame()) {
  png(filename = paste0(fname, ".png"), width = width, height = height, units = "in", res = IEEE_DPI, pointsize = 7)
  eval(expr, envir = env)
  dev.off()
  cat("Saved:", paste0(fname, ".png"), "\n")
}

categories <- c("2-Wheelers", "3-Wheelers", "4-Wheelers")

compute_metrics <- function(train_y, actual, predicted){
  e <- actual-predicted
  me <- mean(e)
  mae <- mean(abs(e))
  rmse <- sqrt(mean(e^2))
  
  n_train <- length(train_y)
  insample_naive_errors <- abs(diff(train_y, lag=12))
  scale <- mean(insample_naive_errors)
  mase <- mae/scale
  
  n_test <- length(actual)
  naive_seed <- tail(train_y,12)
  naive_pred <- rep(naive_seed, length.out=n_test)
  naive_rmse <- sqrt(mean((actual - naive_pred) ^2))
  theils_u <- rmse/naive_rmse
  
  data.frame(RMSE = round(rmse, 4), MAE = round(mae,  4), MASE = round(mase, 4), Theils_U = round(theils_u, 4), ME = round(me, 4))
}

# =============================================================================
# --- SECTION 2: DATA LOADING & CLEANING ---
# =============================================================================
url <- "https://raw.githubusercontent.com/Hiyuu21/stat_for_ds/main/EV_Dataset.csv"
raw_data <- read.csv(url)

clean_data <- raw_data |>
  mutate(Date = mdy(Date)) |>
  filter(Vehicle_Category %in% c("2-Wheelers", "3-Wheelers", "4-Wheelers")) |>
  group_by(Date, Vehicle_Category) |>
  summarise(Total_Sales = sum(EV_Sales_Quantity, na.rm = TRUE), .groups = "drop") |>
  arrange(Date)

# =============================================================================
# --- SECTION 3: EDA ---
# =============================================================================
summary_table <- clean_data |>
  group_by(Vehicle_Category) |>
  summarise(Total_Months = n(), Mean_Sales = round(mean(Total_Sales), 2), Median_Sales = median(Total_Sales), Std_Dev = round(sd(Total_Sales), 2), Min_Sales = min(Total_Sales), Max_Sales = max(Total_Sales))

p_ts <- ggplot(clean_data, aes(x = Date, y = Total_Sales, colour = Vehicle_Category)) +
  geom_line(linewidth = 0.5) + facet_wrap(~ Vehicle_Category, scales = "free_y", ncol = 1) +
  theme_ieee() + labs(title = "Monthly EV Sales by Vehicle Category (India, 2014–2024)", x = "Date", y = "Total Sales (units)") + theme(legend.position = "none")
save_ieee("eda_01_timeseries_all", p_ts, width = IEEE_W_DOUBLE, height = 5.5)

# ---- STL Decomposition per Category (raw, undifferenced, full series) ----
strength_table <- data.frame()

for (cat_name in categories) {
  
  raw_vec <- clean_data |>
    filter(Vehicle_Category == cat_name) |>
    pull(Total_Sales)
  
  detected_S_raw <- findfrequency(raw_vec)
  
  cat_ts <- ts(raw_vec, start = c(2014, 1), frequency = 12)
  
  stl_fit <- stl(cat_ts, s.window = "periodic")
  
  # ---- Base-R STL plot (matches original eda.R style) ----
  fname <- paste0("eda_03_stl_", gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname, expression({
    plot(stl_fit, main = paste("STL Decomposition:", cat_name))
  }), width = IEEE_W_SINGLE, height = 4.5)
  
  # ---- Robustness check: detect S on the DETRENDED series ----
  detrended <- cat_ts - stl_fit$time.series[, "trend"]
  detected_S_detrended <- findfrequency(detrended)
  
  cat("\n", cat_name,
      "- Detected S (raw):", detected_S_raw,
      "| Detected S (detrended):", detected_S_detrended, "\n")
  
  # ---- Seasonal & Trend Strength (Hyndman & Athanasopoulos, FPP3 §3.6) ----
  ts_comp <- stl_fit$time.series
  remainder <- ts_comp[, "remainder"]
  seasonal  <- ts_comp[, "seasonal"]
  trend     <- ts_comp[, "trend"]
  
  Fs <- max(0, 1 - var(remainder) / var(seasonal + remainder))
  Ft <- max(0, 1 - var(remainder) / var(trend + remainder))
  
  seasonal_cycle <- as.numeric(ts_comp[1:12, "seasonal"])
  names(seasonal_cycle) <- month.abb
  
  cat("\n--- Seasonal Indices -", cat_name, "---\n")
  print(round(seasonal_cycle, 1))
  
  strength_table <- rbind(strength_table, data.frame(
    Vehicle_Category  = cat_name,
    Detected_S_Raw     = detected_S_raw,
    Detected_S_Detrend = detected_S_detrended,
    Seasonal_Strength  = round(Fs, 4),
    Trend_Strength     = round(Ft, 4)
  ))
}

cat("\n--- Seasonal Period, Strength & Trend (raw series, EDA stage) ---\n")
print(strength_table)

# =============================================================================
# --- SECTION 4 & 5: TRAIN / TEST SPLIT, STATIONARITY & DIFFERENCING ---
# =============================================================================
# This section ensures that every model receives exactly the same differenced 
# dataset and dynamically finds the best d and D orders based on statistical tests.
train_n <- 96

split_category <- function(dataset, cat_name, train_n = 96) {
  cat_data <- dataset |> filter(Vehicle_Category == cat_name)
  list(full = cat_data, train = cat_data |> slice(1:train_n), test = cat_data |> slice((train_n + 1):n()))
}

splits <- setNames(lapply(categories, split_category, dataset = clean_data, train_n = train_n), categories)

make_ts <- function(df, start_year = 2014, start_month = 1, freq = 12) {
  ts(df$Total_Sales, start = c(start_year, start_month), frequency = freq)
}

# Iterate over category names to trigger console prints and save plots dynamically
ts_data <- lapply(categories, function(cat_name) {
  s <- splits[[cat_name]]
  full_ts  <- make_ts(s$full)
  train_ts <- make_ts(s$train)
  test_ts  <- make_ts(s$test, start_year=2022, start_month=1)

  cat(paste0("\n======================================================\n"))
  cat(paste0("STATIONARITY & DIFFERENCING: ", cat_name, "\n"))
  cat(paste0("======================================================\n"))
  
  # Run ADF and KPSS tests to assess initial stationarity
  adf_res  <- adf.test(train_ts)  # H0: data is not stationary
  kpss_res <- kpss.test(train_ts)
  
  cat("Stationarity Tests on Raw Training Data:\n")
  cat("   ADF p-value :", round(adf_res$p.value, 4), "-> ", ifelse(adf_res$p.value < 0.05, "Stationary", "Non-stationary"), "\n")
  cat("   KPSS p-value:", round(kpss_res$p.value, 4), "-> ", ifelse(kpss_res$p.value < 0.05, "Non-stationary", "Stationary"), "\n")
  
  # Auto-calculate optimal differencing
  # nsdiffs() finds D. ndiffs() evaluates the series AFTER applying D to find d.
  D_order <- nsdiffs(train_ts)
  d_order <- ndiffs(if(D_order > 0) diff(train_ts, lag=12, differences=D_order) else train_ts)
  
  cat("   Universal seasonal differencing (D)    :", D_order, "\n")
  cat("   Universal non-seasonal differencing (d):", d_order, "\n")

  # 1. Plot RAW ACF/PACF (before differencing)
  fname_acfpacf_raw <- paste0("00_acfpacf_raw_", gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_acfpacf_raw, expression({
    par(mfrow = c(1, 2))
    acf(train_ts,  main = paste("ACF (raw) -",  cat_name))
    pacf(train_ts, main = paste("PACF (raw) -", cat_name))
    par(mfrow = c(1, 1))
  }))

  # 2. Apply the calculated Differencing (d and D) to the Training and Full sets
  diff_ts <- train_ts 
  if (D_order > 0) diff_ts <- diff(diff_ts, lag=12, differences = D_order)
  if (d_order > 0) diff_ts <- diff(diff_ts, differences = d_order)

  full_diff_ts <- full_ts
  if (D_order > 0) full_diff_ts <- diff(full_diff_ts, lag=12, differences = D_order)
  if (d_order > 0) full_diff_ts <- diff(full_diff_ts, differences = d_order)

  # 3. Plot DIFFERENCED ACF/PACF (to verify stationarity has been achieved)
  fname_acfpacf <- paste0("01_acfpacf_diff_", gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_acfpacf, expression({
    par(mfrow = c(1, 2))
    acf(diff_ts,  main = paste("ACF (diff) -",  cat_name))
    pacf(diff_ts, main = paste("PACF (diff) -", cat_name))
    par(mfrow = c(1, 1))
  }))

  list(
    full = full_ts, 
    train = train_ts, 
    test = test_ts, 
    train_diff = diff_ts,      # Stored so models can train on differenced data
    full_diff = full_diff_ts,  # Stored so models can run future forecasts on differenced data
    d = d_order, 
    D = D_order
  )
})
names(ts_data) <- categories # Maps the data list to the category names

# Helper for Prophet dataframe formatting
make_prophet_df <- function(df) {
  df |> mutate(ds = floor_date(Date, "month")) |> transmute(ds, y = Total_Sales)
}
prophet_data <- lapply(splits, function(s) {
  list(full = make_prophet_df(s$full), train = make_prophet_df(s$train), test = make_prophet_df(s$test))
})

# =============================================================================
# --- SECTION 6: INVERSE-DIFFERENCING HELPER FUNCTION ---
# =============================================================================
# Because the models are fed differenced data, their predictions will be near 0.
# This mathematically restores the forecasted values back to the actual sales scale 
# so we can calculate MAE, RMSE, etc. fairly against the actual test set.
inverse_difference <- function(forecast_diff, original_train, d, D){
    # Reconstruct non-seasonal differencing using cumulative sum (cumsum)
    if (d > 0) {
        last_val <- tail(original_train, 1)
        if (D > 0) last_val <- tail(diff(original_train, lag = 12, differences = D), 1)
        forecasted_diff <- cumsum(c(last_val, forecast_diff))[-1]
    }
    
    # Reconstruct seasonal differencing by adding the historical lag-12 values
    if (D > 0) {
        hist_seasonal <- tail(original_train, 12)
        reconstructed <- numeric(length(forecasted_diff))
        for (i in seq_along(forecasted_diff)) {
        prev_val <- if (i <= 12) hist_seasonal[i] else reconstructed[i - 12]
        reconstructed[i] <- prev_val + forecasted_diff[i]
        }
        forecasted_diff <- reconstructed
    }
    
    return(forecasted_diff) # Returns un-differenced predictions
}

# Helper to inverse-difference interval matrices (lower/upper bounds)
inverse_difference_matrix <- function(matrix_diff, original_train, d, D) {
  apply(matrix_diff, 2, function(col) {
    inverse_difference(col, original_train, d, D)
  })
}

# =============================================================================
# --- SARIMA MODEL ---
# =============================================================================
sarima <- function(cat_name, data_list, future_h = 12) {
  cat(paste0("\n======================================================\n"))
  cat(paste0("ANALYZING CATEGORY (SARIMA): ", cat_name, "\n"))
  cat(paste0("======================================================\n"))
  
  train_ts   <- data_list[[cat_name]]$train
  train_diff <- data_list[[cat_name]]$train_diff
  test_ts    <- data_list[[cat_name]]$test
  full_ts    <- data_list[[cat_name]]$full
  
  # ---- 1. Model fitting (auto.arima) ----
  cat("\n1. Fitting SARIMA Model...\n")
  # We MUST force d=0 and D=0 because the data is ALREADY differenced globally.
  # Otherwise, auto.arima will difference it again, destroying the signal.
  model <- auto.arima(train_diff, d=0, D=0, seasonal = TRUE, stepwise = FALSE,
                      approximation = FALSE, trace = TRUE, ic = "aicc")
  cat("\nSelected order:", paste(arimaorder(model), collapse = ","), "\n")
  print(summary(model))
  
  # ---- 2. Residual diagnostics ----
  cat("\n2. Residual Diagnostics (Ljung-Box test):\n")
  checkresiduals(model, plot = TRUE)

  fname_resid <- paste0("sarima_02_residuals_",
                        gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_resid, expression({
    checkresiduals(model, plot = TRUE)
  }), height = 3.5)
  
  # ---- 3. Forecast on held-out test set + evaluation ----
  fc <- forecast(model, h = length(test_ts), level = c(80, 95))
  
  undiff_pred  <- inverse_difference(as.numeric(fc$mean), as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower <- inverse_difference_matrix(fc$lower, as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_upper <- inverse_difference_matrix(fc$upper, as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  
  test_dates <- seq(as.Date("2022-01-01"), by = "month", length.out = length(test_ts))
  plot_df <- data.frame(
    ds    = test_dates,
    y     = as.numeric(test_ts),
    yhat  = undiff_pred,
    lower80 = pmax(0, undiff_lower[, 1]),
    upper80 = undiff_upper[, 1],
    lower95 = pmax(0, undiff_lower[, 2]),
    upper95 = undiff_upper[, 2]
  )
  
  eval <- compute_metrics(
    train_y   = as.numeric(train_ts),
    actual    = as.numeric(test_ts),
    predicted = undiff_pred
  )
  
  cat("\n3. Forecasting and Evaluation on Test Set\n")
  cat("\n--- Full Accuracy Matrix ---\n")
  print(eval)
  
  p_test <- ggplot(plot_df, aes(x = ds)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = "95% PI"), alpha = 0.45) +
    geom_ribbon(aes(ymin = lower80, ymax = upper80, fill = "80% PI"), alpha = 0.55) +
    geom_line(aes(y = y, colour = "Actual Test Data"), linewidth = 0.7) +
    geom_line(aes(y = yhat, colour = "SARIMA Forecast"), linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = c("95% PI" = "#9ecae1", "80% PI" = "#4292c6")) +
    scale_colour_manual(values = c("Actual Test Data" = "#08519c", "SARIMA Forecast" = "#cb181d")) +
    theme_ieee() +
    labs(title = paste("SARIMA Forecast vs Actual (Test Set):", cat_name),
         x = "Date", y = "Sales Quantity (units)", fill = NULL, colour = NULL)
  print(p_test)
  
  fname_test <- paste0("sarima_03_testforecast_",
                       gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  # ---- 4. Refit on FULL data + forecast forward ----
  cat("\n4. Future Forecast (next", future_h, "months beyond dataset):\n")
  
  model_full <- auto.arima(data_list[[cat_name]]$full_diff, d = 0, D = 0, 
                           seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
  
  future_fc_diff <- forecast(model_full, h = future_h, level = c(80, 95))
  
  # Inverse difference mean and interval matrices
  undiff_future <- inverse_difference(as.numeric(future_fc_diff$mean), as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower  <- inverse_difference_matrix(future_fc_diff$lower, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_upper  <- inverse_difference_matrix(future_fc_diff$upper, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  
  # Historical dates + future dates
  hist_dates   <- seq(as.Date("2014-01-01"), by = "month", length.out = length(full_ts))
  future_dates <- seq(tail(hist_dates, 1) %m+% months(1), by = "month", length.out = future_h)
  
  hist_df <- data.frame(ds = hist_dates, y = as.numeric(full_ts))
  future_df <- data.frame(
    ds      = future_dates,
    yhat    = undiff_future,
    lower80 = pmax(0, undiff_lower[, 1]),
    upper80 = undiff_upper[, 1],
    lower95 = pmax(0, undiff_lower[, 2]),
    upper95 = undiff_upper[, 2]
  )
  
  p_future <- ggplot() +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower95, ymax = upper95, fill = "95% PI"), alpha = 0.45) +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower80, ymax = upper80, fill = "80% PI"), alpha = 0.55) +
    geom_line(data = hist_df, aes(x = ds, y = y, colour = "Historical"), linewidth = 0.6) +
    geom_line(data = future_df, aes(x = ds, y = yhat, colour = "SARIMA Forecast"), linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = c("95% PI" = "#9ecae1", "80% PI" = "#4292c6")) +
    scale_colour_manual(values = c("Historical" = "#252525", "SARIMA Forecast" = "#cb181d")) +
    theme_ieee() +
    labs(title = paste("SARIMA Future Forecast (Next", future_h, "Months):", cat_name),
         x = "Date", y = "Sales Quantity (units)", fill = NULL, colour = NULL)
  print(p_future)
  
  fname_future <- paste0("sarima_04_futureforecast_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  return(list(
    Model = model, FinalModel = model_full, TestForecast = fc,
    FutureForecast = future_fc_diff, Accuracy = eval,
    TestPlot = p_test, FuturePlot = p_future
  ))
}


# =============================================================================
# --- STL-ETS MODEL ---
# =============================================================================
stl_ets <- function(cat_name, data_list, future_h = 12) {
  cat(paste0("\n======================================================\n"))
  cat(paste0("ANALYZING CATEGORY (STL-ETS): ", cat_name, "\n"))
  cat(paste0("======================================================\n"))
  
  train_ts   <- data_list[[cat_name]]$train
  train_diff <- data_list[[cat_name]]$train_diff
  test_ts    <- data_list[[cat_name]]$test
  full_ts    <- data_list[[cat_name]]$full
  
  # ---- 1. Hyperparameter Grid Search & Optimization ----
  cat("\n1. Running Joint Grid Search for STL Parameters (s.window & robust) & ETS Architecture...\n")
  
  s_windows      <- list(7, 9, 11, 13, 15, 21, "periodic")
  robust_options <- c(TRUE, FALSE)
  grid_results   <- data.frame()
  
  for (sw in s_windows) {
    for (rob in robust_options) {
      
      stl_cand <- tryCatch({
        stl(train_diff, s.window = sw, robust = rob)
      }, error = function(e) NULL)
      
      if (!is.null(stl_cand)) {
        sa_cand <- train_diff - stl_cand$time.series[, "seasonal"]
        
        ets_cand <- tryCatch({
          ets(sa_cand, model = "ZZN", ic = "aicc")
        }, error = function(e) NULL)
        
        if (!is.null(ets_cand)) {
          fitted_cand <- ets_cand$fitted + stl_cand$time.series[, "seasonal"]
          rmse_cand   <- sqrt(mean((train_diff - fitted_cand)^2, na.rm = TRUE))
          mae_cand    <- mean(abs(train_diff - fitted_cand), na.rm = TRUE)
          
          grid_results <- rbind(grid_results, data.frame(
            s_window  = as.character(sw),
            robust    = rob,
            ets_model = ets_cand$method,
            AICc      = ets_cand$aicc,
            RMSE      = rmse_cand,
            MAE       = mae_cand
          ))
        }
      }
    }
  }
  
  grid_results <- grid_results[order(grid_results$AICc), ]
  
  cat("   --- Grid Search Results ---\n")
  print(grid_results)
  
  best_sw  <- grid_results$s_window[1]
  if (best_sw != "periodic") best_sw <- as.numeric(best_sw)
  best_rob <- grid_results$robust[1]
  
  cat(paste0("\n   [OPTIMAL PARAMETERS SELECTED]: s.window = ", best_sw, ", robust = ", best_rob, "\n"))
  
  # ---- 2. Perform Optimal STL Decomposition ----
  cat("\n2. Fitting Optimal STL Decomposition on Differenced Data...\n")
  stl_decomp <- stl(train_diff, s.window = best_sw, robust = best_rob)
  
  fname_stl <- paste0("stlets_01_decomposition_",
                      gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_stl, expression({
    plot(stl_decomp, main = paste("Optimal STL Decomposition -", cat_name))
  }), height = 4.0)
  
  # ---- 3. Fit Final Dynamic STL-ETS Model using stlm ----
  cat("\n3. Fitting Final Dynamic STL-ETS Model...\n")
  
  # etsmodel="ZZN" forces the underlying ets function to search for the best trend dynamically
  model <- stlm(
    train_diff, 
    s.window      = best_sw, 
    robust        = best_rob, 
    modelfunction = ets,
    etsmodel      = "ZZN", 
    ic            = "aicc"
  )
  
  ets_fit <- model$model
  
  cat("   ---------------- ETS MODEL FITTING SUMMARY ----------------\n")
  cat("   Model Type          :", ets_fit$method, "\n")
  cat("   Selected s.window   :", as.character(best_sw), "\n")
  cat("   Selected Robust Flag:", best_rob, "\n")
  cat("   Alpha (Level)       :", round(ets_fit$par["alpha"], 5), "\n")
  if ("beta" %in% names(ets_fit$par)) {
    cat("   Beta (Trend)        :", round(ets_fit$par["beta"], 5), "\n")
  }
  if ("phi" %in% names(ets_fit$par)) {
    cat("   Phi (Damping factor):", round(ets_fit$par["phi"], 5), "\n")
  }
  cat("   Sigma^2 (Variance)  :", round(ets_fit$sigma2, 4), "\n")
  cat("   AICc                :", round(ets_fit$aicc, 4), "\n")
  cat("   -----------------------------------------------------------\n\n")
  
  # ---- 4. Residual Diagnostics ----
  cat("4. Residual Diagnostics (Ljung-Box test):\n")
  print(checkresiduals(model, plot = TRUE))
  
  fname_resid <- paste0("stlets_02_residuals_",
                        gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_resid, expression({
    checkresiduals(model, plot = TRUE)
  }), height = 3.5)
  
  # ---- 5. Out-of-Sample Test Evaluation ----
  cat("\n5. Forecasting and Evaluation on Test Set...\n")
  fc <- forecast(model, h = length(test_ts), level = c(80, 95))
  
  undiff_pred  <- inverse_difference(as.numeric(fc$mean), as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower <- inverse_difference_matrix(fc$lower, as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_upper <- inverse_difference_matrix(fc$upper, as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  
  test_dates <- seq(as.Date("2022-01-01"), by = "month", length.out = length(test_ts))
  plot_df <- data.frame(
    ds      = test_dates,
    y       = as.numeric(test_ts),
    yhat    = undiff_pred,
    lower80 = pmax(0, undiff_lower[, 1]),
    upper80 = undiff_upper[, 1],
    lower95 = pmax(0, undiff_lower[, 2]),
    upper95 = undiff_upper[, 2]
  )
  
  eval <- compute_metrics(
    train_y   = as.numeric(train_ts),
    actual    = as.numeric(test_ts),
    predicted = undiff_pred
  )
  
  cat("--- Out-of-Sample Accuracy Matrix ---\n")
  print(eval)
  
  p_test <- ggplot(plot_df, aes(x = ds)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = "95% PI"), alpha = 0.45) +
    geom_ribbon(aes(ymin = lower80, ymax = upper80, fill = "80% PI"), alpha = 0.55) +
    geom_line(aes(y = y, colour = "Actual Test Data"), linewidth = 0.7) +
    geom_line(aes(y = yhat, colour = "STL+ETS Forecast"), linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = c("95% PI" = "#9ecae1", "80% PI" = "#4292c6")) +
    scale_colour_manual(values = c("Actual Test Data" = "#08519c", "STL+ETS Forecast" = "#cb181d")) +
    theme_ieee() +
    labs(title = paste("Optimized STL+ETS Forecast vs Actual:", cat_name),
         x = "Date", y = "Sales Quantity (units)", fill = NULL, colour = NULL)
  print(p_test)
  
  fname_test <- paste0("stlets_03_testforecast_",
                       gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  # ---- 6. Refit on FULL Data with Dynamic Parameters ----
  cat("\n6. Refitting Full Dataset with Dynamic Selection & Forecasting Future...\n")
  
  final_model <- stlm(
    data_list[[cat_name]]$full_diff, 
    s.window      = best_sw, 
    robust        = best_rob, 
    modelfunction = ets,
    etsmodel      = "ZZN",
    ic            = "aicc"
  )
  
  future_fc_diff <- forecast(final_model, h = future_h, level = c(80, 95))
  
  undiff_future <- inverse_difference(as.numeric(future_fc_diff$mean), as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower  <- inverse_difference_matrix(future_fc_diff$lower, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_upper  <- inverse_difference_matrix(future_fc_diff$upper, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  
  hist_dates   <- seq(as.Date("2014-01-01"), by = "month", length.out = length(full_ts))
  future_dates <- seq(tail(hist_dates, 1) %m+% months(1), by = "month", length.out = future_h)
  
  hist_df   <- data.frame(ds = hist_dates, y = as.numeric(full_ts))
  future_df <- data.frame(
    ds      = future_dates,
    yhat    = undiff_future,
    lower80 = pmax(0, undiff_lower[, 1]),
    upper80 = undiff_upper[, 1],
    lower95 = pmax(0, undiff_lower[, 2]),
    upper95 = undiff_upper[, 2]
  )
  
  p_future <- ggplot() +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower95, ymax = upper95, fill = "95% PI"), alpha = 0.45) +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower80, ymax = upper80, fill = "80% PI"), alpha = 0.55) +
    geom_line(data = hist_df, aes(x = ds, y = y, colour = "Historical"), linewidth = 0.6) +
    geom_line(data = future_df, aes(x = ds, y = yhat, colour = "STL+ETS Forecast"), linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = c("95% PI" = "#9ecae1", "80% PI" = "#4292c6")) +
    scale_colour_manual(values = c("Historical" = "#252525", "STL+ETS Forecast" = "#cb181d")) +
    theme_ieee() +
    labs(title = paste("STL+ETS Future Forecast (Next", future_h, "Months):", cat_name),
         x = "Date", y = "Sales Quantity (units)", fill = NULL, colour = NULL)
  print(p_future)
  
  fname_future <- paste0("stlets_04_futureforecast_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  return(list(
    GridSearch = grid_results, OptimalSW = best_sw, OptimalRobust = best_rob,
    Model = model, FinalModel = final_model, TestForecast = fc,
    FutureForecast = future_fc_diff, Accuracy = eval,
    TestPlot = p_test, FuturePlot = p_future
  ))
}


# =============================================================================
# --- TBATS MODEL ---
# =============================================================================
tbats_analysis <- function(cat_name, data_list, future_h = 12) {
  cat("\n======================================================\n")
  cat("ANALYZING CATEGORY (TBATS):", cat_name, "\n")
  cat("======================================================\n")
  
  train_ts   <- data_list[[cat_name]]$train
  train_diff <- data_list[[cat_name]]$train_diff
  test_ts    <- data_list[[cat_name]]$test
  full_ts    <- data_list[[cat_name]]$full
  
  # ---- 1. Fit TBATS model ----
  cat("\n1. Fitting TBATS Model\n")
  model <- tbats(train_diff, use.box.cox = FALSE)
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
  checkresiduals(model, plot = TRUE)

  fname_resid <- paste0("tbats_02_residuals_",
                        gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee_base(fname_resid, expression({
    checkresiduals(model, plot = TRUE)
  }), height = 3.5)
  
  # ---- 4. Forecast on Test Set ----
  cat("\n4. Forecasting on Test Set\n")
  fc <- forecast(model, h = length(test_ts), level = c(80, 95))
  
  undiff_pred  <- inverse_difference(as.numeric(fc$mean), as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower <- inverse_difference_matrix(fc$lower, as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_upper <- inverse_difference_matrix(fc$upper, as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  
  test_dates <- seq(as.Date("2022-01-01"), by = "month", length.out = length(test_ts))
  plot_df <- data.frame(
    ds    = test_dates,
    y     = as.numeric(test_ts),
    yhat  = undiff_pred,
    lower80 = pmax(0, undiff_lower[, 1]),
    upper80 = undiff_upper[, 1],
    lower95 = pmax(0, undiff_lower[, 2]),
    upper95 = undiff_upper[, 2]
  )
  
  eval <- compute_metrics(
    train_y   = as.numeric(train_ts),
    actual    = as.numeric(test_ts),
    predicted = undiff_pred
  )
  
  cat("\n5. Forecasting and Evaluation on Test Set\n")
  cat("\n--- Full Accuracy Matrix ---\n")
  print(eval)
  
  p_test <- ggplot(plot_df, aes(x = ds)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = "95% PI"), alpha = 0.45) +
    geom_ribbon(aes(ymin = lower80, ymax = upper80, fill = "80% PI"), alpha = 0.55) +
    geom_line(aes(y = y, colour = "Actual Test Data"), linewidth = 0.7) +
    geom_line(aes(y = yhat, colour = "TBATS Forecast"), linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = c("95% PI" = "#9ecae1", "80% PI" = "#4292c6")) +
    scale_colour_manual(values = c("Actual Test Data" = "#08519c", "TBATS Forecast" = "#cb181d")) +
    theme_ieee() +
    labs(title = paste("TBATS Forecast vs Actual (Test Set):", cat_name),
         x = "Date", y = "Sales Quantity (units)", fill = NULL, colour = NULL)
  print(p_test)
  
  fname_test <- paste0("tbats_03_testforecast_",
                       gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_test, p_test)
  
  # ---- 6. Refit on Full Data + Future Forecast ----
  cat("\n6. Refitting TBATS on Full Dataset and Forecasting Future...\n")
  
  final_model <- tbats(data_list[[cat_name]]$full_diff, use.box.cox = FALSE)
  future_fc_diff <- forecast(final_model, h = future_h, level = c(80, 95))
  
  undiff_future <- inverse_difference(as.numeric(future_fc_diff$mean), as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower  <- inverse_difference_matrix(future_fc_diff$lower, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_upper  <- inverse_difference_matrix(future_fc_diff$upper, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  
  hist_dates   <- seq(as.Date("2014-01-01"), by = "month", length.out = length(full_ts))
  future_dates <- seq(tail(hist_dates, 1) %m+% months(1), by = "month", length.out = future_h)
  
  hist_df <- data.frame(ds = hist_dates, y = as.numeric(full_ts))
  future_df <- data.frame(
    ds      = future_dates,
    yhat    = undiff_future,
    lower80 = pmax(0, undiff_lower[, 1]),
    upper80 = undiff_upper[, 1],
    lower95 = pmax(0, undiff_lower[, 2]),
    upper95 = undiff_upper[, 2]
  )
  
  p_future <- ggplot() +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower95, ymax = upper95, fill = "95% PI"), alpha = 0.45) +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower80, ymax = upper80, fill = "80% PI"), alpha = 0.55) +
    geom_line(data = hist_df, aes(x = ds, y = y, colour = "Historical"), linewidth = 0.6) +
    geom_line(data = future_df, aes(x = ds, y = yhat, colour = "TBATS Forecast"), linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = c("95% PI" = "#9ecae1", "80% PI" = "#4292c6")) +
    scale_colour_manual(values = c("Historical" = "#252525", "TBATS Forecast" = "#cb181d")) +
    theme_ieee() +
    labs(title = paste("TBATS Future Forecast (Next", future_h, "Months):", cat_name),
         x = "Date", y = "Sales Quantity (units)", fill = NULL, colour = NULL)
  print(p_future)
  
  fname_future <- paste0("tbats_04_futureforecast_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_future, p_future)
  
  return(list(
    Model = model, FinalModel = final_model, TestForecast = fc,
    FutureForecast = future_fc_diff, Accuracy = eval,
    TestPlot = p_test, FuturePlot = p_future
  ))
}


# =============================================================================
# --- PROPHET MODEL ---
# =============================================================================
CPS <- 0.5 # Changepoint Prior Scale configuration

run_prophet <- function(cat_name, data_list, future_h = 12) {
  cat(paste0("\n=========================================\n"))
  cat(paste0("ANALYZING CATEGORY (PROPHET): ", cat_name, "\n"))
  cat(paste0("=========================================\n"))
  
  train_ts   <- data_list[[cat_name]]$train
  train_diff <- data_list[[cat_name]]$train_diff
  
  # ---- Generate matching dates for Prophet ----
  # Differencing drops initial rows (e.g., lag 12 drops the first year).
  # We must align Prophet's dates directly to the remaining values.
  train_dates <- seq(as.Date("2014-01-01"), by = "month", length.out = length(train_ts))
  valid_dates <- tail(train_dates, length(train_diff))
  
  train_df <- data.frame(ds = valid_dates, y = as.numeric(train_diff))
  test_df  <- prophet_data[[cat_name]]$test
  
  # ---- 1. Fit model on differenced training data (80% + 95% PI) ----
  set.seed(123)
  m80 <- prophet(
    yearly.seasonality      = TRUE,
    weekly.seasonality      = FALSE,
    daily.seasonality       = FALSE,
    interval.width          = 0.8,
    changepoint.prior.scale = CPS
  )
  m80 <- fit.prophet(m80, train_df)
  
  set.seed(123)
  m95 <- prophet(
    yearly.seasonality      = TRUE,
    weekly.seasonality      = FALSE,
    daily.seasonality       = FALSE,
    interval.width          = 0.95,
    changepoint.prior.scale = CPS
  )
  m95 <- fit.prophet(m95, train_df)
  
  # ---- 2. Predict on test period (both widths) ----
  future <- make_future_dataframe(m80, periods = nrow(test_df), freq = "month")
  fc80   <- predict(m80, future)
  fc95   <- predict(m95, future)
  test_pred <- fc80$yhat[(nrow(train_df) + 1):nrow(fc80)]   # yhat identical across both fits
  
  # Inverse difference the Prophet predictions
  undiff_pred <- inverse_difference(
    forecast_diff  = test_pred,
    original_train = as.numeric(train_ts),
    d = data_list[[cat_name]]$d,
    D = data_list[[cat_name]]$D
  )
  
  metrics <- compute_metrics(
    train_y   = as.numeric(train_ts),
    actual    = test_df$y,
    predicted = undiff_pred
  )
  
  cat("\n--- Test Set Metrics ---\n")
  print(metrics)
  
  # ---- 4. Test forecast plot (80% + 95% PI) ----
  idx <- (nrow(train_df) + 1):nrow(fc80)
  undiff_lower80 <- pmax(0, inverse_difference(fc80$yhat_lower[idx], as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D))
  undiff_upper80 <- inverse_difference(fc80$yhat_upper[idx], as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower95 <- pmax(0, inverse_difference(fc95$yhat_lower[idx], as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D))
  undiff_upper95 <- inverse_difference(fc95$yhat_upper[idx], as.numeric(train_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  
  plot_df <- data.frame(
    ds = test_df$ds, y = test_df$y, yhat = undiff_pred,
    lower80 = undiff_lower80, upper80 = undiff_upper80,
    lower95 = undiff_lower95, upper95 = undiff_upper95
  )
  
  p_test <- ggplot(plot_df, aes(x = ds)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = "95% PI"), alpha = 0.45) +
    geom_ribbon(aes(ymin = lower80, ymax = upper80, fill = "80% PI"), alpha = 0.55) +
    geom_line(aes(y = y, colour = "Actual Test Data"), linewidth = 0.7) +
    geom_line(aes(y = yhat, colour = "Prophet Forecast"), linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = c("95% PI" = "#9ecae1", "80% PI" = "#4292c6")) +
    scale_colour_manual(values = c("Actual Test Data" = "#08519c", "Prophet Forecast" = "#cb181d")) +
    theme_ieee() +
    labs(title = paste("Prophet Forecast vs Actual (Test Set):", cat_name),
         x = "Date", y = "Sales Quantity (units)", fill = NULL, colour = NULL)
  print(p_test)
  
  fname_test <- paste0("prophet_01_testforecast_",
                       gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  # ---- 5. Refit on full data + future forecast ----
  cat("\n5. Refitting Prophet on Full Dataset and Forecasting Future...\n")
  
  full_ts   <- data_list[[cat_name]]$full
  full_diff <- data_list[[cat_name]]$full_diff
  
  # Dynamically map dates to the differenced full dataset
  full_dates <- seq(as.Date("2014-01-01"), by = "month", length.out = length(full_ts))
  valid_full_dates <- tail(full_dates, length(full_diff))
  full_diff_df <- data.frame(ds = valid_full_dates, y = as.numeric(full_diff))
  
  set.seed(123)
  m_full80 <- prophet(
    yearly.seasonality      = TRUE,
    weekly.seasonality      = FALSE,
    daily.seasonality       = FALSE,
    interval.width          = 0.8,
    changepoint.prior.scale = CPS
  )
  m_full80 <- fit.prophet(m_full80, full_diff_df)
  
  set.seed(123)
  m_full95 <- prophet(
    yearly.seasonality      = TRUE,
    weekly.seasonality      = FALSE,
    daily.seasonality       = FALSE,
    interval.width          = 0.95,
    changepoint.prior.scale = CPS
  )
  m_full95 <- fit.prophet(m_full95, full_diff_df)
  
  future_full <- make_future_dataframe(m_full80, periods = future_h, freq = "month")
  forecast_full80 <- predict(m_full80, future_full)
  forecast_full95 <- predict(m_full95, future_full)
  
  fidx <- (nrow(full_diff_df) + 1):nrow(forecast_full80)
  future_pred_diff   <- forecast_full80$yhat[fidx]          # identical across both fits
  future_lower80_diff <- forecast_full80$yhat_lower[fidx]
  future_upper80_diff <- forecast_full80$yhat_upper[fidx]
  future_lower95_diff <- forecast_full95$yhat_lower[fidx]
  future_upper95_diff <- forecast_full95$yhat_upper[fidx]
  
  undiff_future  <- inverse_difference(future_pred_diff, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower80 <- pmax(0, inverse_difference(future_lower80_diff, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D))
  undiff_upper80 <- inverse_difference(future_upper80_diff, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  undiff_lower95 <- pmax(0, inverse_difference(future_lower95_diff, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D))
  undiff_upper95 <- inverse_difference(future_upper95_diff, as.numeric(full_ts), data_list[[cat_name]]$d, data_list[[cat_name]]$D)
  
  future_df <- data.frame(
    ds = tail(forecast_full80$ds, future_h),
    yhat = undiff_future,
    lower80 = undiff_lower80, upper80 = undiff_upper80,
    lower95 = undiff_lower95, upper95 = undiff_upper95
  )
  
  hist_df <- data.frame(ds = full_dates, y = as.numeric(full_ts))
  
  p_future <- ggplot() +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower95, ymax = upper95, fill = "95% PI"), alpha = 0.45) +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower80, ymax = upper80, fill = "80% PI"), alpha = 0.55) +
    geom_line(data = hist_df, aes(x = ds, y = y, colour = "Historical"), linewidth = 0.6) +
    geom_line(data = future_df, aes(x = ds, y = yhat, colour = "Prophet Forecast"), linewidth = 0.7, linetype = "dashed") +
    scale_fill_manual(values = c("95% PI" = "#9ecae1", "80% PI" = "#4292c6")) +
    scale_colour_manual(values = c("Historical" = "#252525", "Prophet Forecast" = "#cb181d")) +
    theme_ieee() +
    labs(title = paste("Prophet Future Forecast (Next", future_h, "Months):", cat_name),
         x = "Date", y = "Sales Quantity (units)", fill = NULL, colour = NULL)
  print(p_future)
  
  fname_future <- paste0("prophet_02_futureforecast_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  return(list(
    Model = m80, FinalModel = m_full80, TestForecast = fc80,
    FutureForecast = forecast_full80, Metrics = metrics,
    TestPlot = p_test, FuturePlot = p_future
  ))
}

# =============================================================================
# --- FINAL EXECUTION PIPELINE WITH DETAILED LOGS ---
# =============================================================================

# 1. SARIMA Logs
sink("SARIMA_Detailed_Log.txt", split = TRUE)
sarima_2w <- sarima("2-Wheelers", ts_data)
sarima_3w <- sarima("3-Wheelers", ts_data)
sarima_4w <- sarima("4-Wheelers", ts_data)

comparison_sarima <- bind_rows(
  sarima_2w$Accuracy |> mutate(Category = "2-Wheelers"),
  sarima_3w$Accuracy |> mutate(Category = "3-Wheelers"),
  sarima_4w$Accuracy |> mutate(Category = "4-Wheelers")
) |> relocate(Category) |> mutate(Model = "SARIMA")
print(comparison_sarima)
saveRDS(comparison_sarima, "comparison_sarima.rds")
sink()
cat("\n--- SARIMA details saved to 'SARIMA_Detailed_Log.txt' ---\n")


# 2. STL-ETS Logs
sink("STLETS_Detailed_Log.txt", split = TRUE)
ets_2w <- stl_ets("2-Wheelers", ts_data)
ets_3w <- stl_ets("3-Wheelers", ts_data)
ets_4w <- stl_ets("4-Wheelers", ts_data)

comparison_ets <- bind_rows(
  ets_2w$Accuracy |> mutate(Category = "2-Wheelers"),
  ets_3w$Accuracy |> mutate(Category = "3-Wheelers"),
  ets_4w$Accuracy |> mutate(Category = "4-Wheelers")
) |> relocate(Category) |> mutate(Model = "STL+ETS")
print(comparison_ets)
saveRDS(comparison_ets, "comparison_ets.rds")
sink()
cat("\n--- STL+ETS details saved to 'STLETS_Detailed_Log.txt' ---\n")


# 3. TBATS Logs
sink("TBATS_Detailed_Log.txt", split = TRUE)
tbats_2w <- tbats_analysis("2-Wheelers", ts_data)
tbats_3w <- tbats_analysis("3-Wheelers", ts_data)
tbats_4w <- tbats_analysis("4-Wheelers", ts_data)

comparison_tbats <- bind_rows(
  tbats_2w$Accuracy |> mutate(Category = "2-Wheelers"),
  tbats_3w$Accuracy |> mutate(Category = "3-Wheelers"),
  tbats_4w$Accuracy |> mutate(Category = "4-Wheelers")
) |> relocate(Category) |> mutate(Model = "TBATS")
print(comparison_tbats)
saveRDS(comparison_tbats, "comparison_tbats.rds")
sink()
cat("\n--- TBATS details saved to 'TBATS_Detailed_Log.txt' ---\n")


# 4. Prophet Logs
sink("Prophet_Detailed_Log.txt", split = TRUE)
prophet_2w <- run_prophet("2-Wheelers", ts_data)
prophet_3w <- run_prophet("3-Wheelers", ts_data)
prophet_4w <- run_prophet("4-Wheelers", ts_data)

comparison_prophet <- bind_rows(
  prophet_2w$Metrics |> mutate(Category = "2-Wheelers"),
  prophet_3w$Metrics |> mutate(Category = "3-Wheelers"),
  prophet_4w$Metrics |> mutate(Category = "4-Wheelers")
) |> relocate(Category) |> mutate(Model = "Prophet")
print(comparison_prophet)
saveRDS(comparison_prophet, "comparison_prophet.rds")
sink()
cat("\n--- Prophet details saved to 'Prophet_Detailed_Log.txt' ---\n")


library(patchwork)

# =============================================================================
# --- COMBINED 4-MODEL FORECAST VISUALIZATION ---
# =============================================================================

plot_combined_forecast <- function(cat_name, sarima_res, ets_res, tbats_res, prophet_res, data_list, future_h = 12) {
  
  full_ts <- data_list[[cat_name]]$full
  d_order <- data_list[[cat_name]]$d
  D_order <- data_list[[cat_name]]$D
  
  # 1. Extract and inverse-difference the future point forecasts
  sarima_fc <- inverse_difference(as.numeric(sarima_res$FutureForecast$mean), as.numeric(full_ts), d_order, D_order)
  ets_fc    <- inverse_difference(as.numeric(ets_res$FutureForecast$mean), as.numeric(full_ts), d_order, D_order)
  tbats_fc  <- inverse_difference(as.numeric(tbats_res$FutureForecast$mean), as.numeric(full_ts), d_order, D_order)
  
  prophet_full <- prophet_res$FutureForecast
  prophet_future_diff <- tail(prophet_full$yhat, future_h)
  prophet_fc <- inverse_difference(prophet_future_diff, as.numeric(full_ts), d_order, D_order)
  
  # 2. Align dates for plotting
  hist_dates <- seq(as.Date("2014-01-01"), by = "month", length.out = length(full_ts))
  future_dates <- seq(tail(hist_dates, 1) %m+% months(1), by = "month", length.out = future_h)
  
  df_hist <- data.frame(Date = hist_dates, Sales = as.numeric(full_ts), Model = "Historical") |> 
    filter(Date >= as.Date("2022-01-01"))
  
  df_sarima  <- data.frame(Date = future_dates, Sales = sarima_fc,  Model = "SARIMA")
  df_ets     <- data.frame(Date = future_dates, Sales = ets_fc,     Model = "STL+ETS")
  df_tbats   <- data.frame(Date = future_dates, Sales = tbats_fc,   Model = "TBATS")
  df_prophet <- data.frame(Date = future_dates, Sales = prophet_fc, Model = "Prophet")
  
  df_combined <- bind_rows(df_hist, df_sarima, df_ets, df_tbats, df_prophet)
  
  # 3. Build the plot using highly distinct linetypes and a slightly thicker line
  p <- ggplot(df_combined, aes(x = Date, y = Sales, color = Model, linetype = Model)) +
    geom_line(linewidth = 0.9) + 
    scale_color_manual(values = c("Historical" = "#252525", 
                                  "SARIMA"     = "#377eb8", 
                                  "STL+ETS"    = "#4daf4a", 
                                  "TBATS"      = "#984ea3", 
                                  "Prophet"    = "#e41a1c")) +
    scale_linetype_manual(values = c("Historical" = "solid", 
                                     "SARIMA"     = "longdash", 
                                     "STL+ETS"    = "dotdash", 
                                     "TBATS"      = "twodash", 
                                     "Prophet"    = "dashed")) +
    theme_ieee() +
    labs(title = cat_name, x = NULL, y = "Sales") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  return(p)
}

# --- Generate Individual Plots ---
cat("\nGenerating Combined Forecast Plots...\n")
p_2w <- plot_combined_forecast("2-Wheelers", sarima_2w, ets_2w, tbats_2w, prophet_2w, ts_data)
p_3w <- plot_combined_forecast("3-Wheelers", sarima_3w, ets_3w, tbats_3w, prophet_3w, ts_data)
p_4w <- plot_combined_forecast("4-Wheelers", sarima_4w, ets_4w, tbats_4w, prophet_4w, ts_data)

# --- Merge and Save as a Single Image ---
# Use the '/' operator from patchwork to stack vertically
# plot_layout(guides = 'collect') merges the legends into one shared legend at the bottom
merged_plot <- (p_2w / p_3w / p_4w) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")

# Add an overall title to the merged layout
merged_plot <- merged_plot + plot_annotation(
  title = "Future Forecasts Across All Models and Categories",
  theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5, family = "serif"))
)

print(merged_plot)

# Save the merged plot with a taller height to accommodate all three graphs
save_ieee("combined_forecast_all_categories", merged_plot, width = IEEE_W_DOUBLE, height = 8.5)
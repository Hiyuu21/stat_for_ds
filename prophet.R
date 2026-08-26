library(prophet)
library(dplyr)
library(ggplot2)

# =============================================================================
# 1. DATA LOADING
# =============================================================================

data_list <- readRDS("prophet_data.rds")

# =============================================================================
# 2. IEEE PLOT SETTINGS
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


# =============================================================================
# 3. PROPHET MODELLING & EVALUATION FUNCTION
# =============================================================================

run_prophet <- function(cat_name, future_h = 12) {
  cat(paste0("\n=========================================\n"))
  cat(paste0("ANALYZING CATEGORY: ", cat_name, "\n"))
  cat(paste0("=========================================\n"))
  
  train_df <- data_list[[cat_name]]$train
  test_df  <- data_list[[cat_name]]$test
  full_df  <- data_list[[cat_name]]$full
  
  # ---- 1. Fit model on training data ----
  m <- prophet(
    yearly.seasonality  = TRUE,
    weekly.seasonality  = FALSE,
    daily.seasonality   = FALSE
  )
  m <- fit.prophet(m, train_df)
  
  # ---- 2. Predict on test period ----
  # make_future_dataframe includes the training period rows first,
  # so test predictions start at index (nrow(train_df) + 1)
  future   <- make_future_dataframe(m, periods = nrow(test_df), freq = "month")
  fc       <- predict(m, future)
  test_pred <- fc$yhat[(nrow(train_df) + 1):nrow(fc)]
  actual    <- test_df$y
  n_test    <- length(actual)
  
  # ---- 3. Metrics ----
  mae  <- mean(abs(actual - test_pred))
  rmse <- sqrt(mean((actual - test_pred)^2))
  me <- mean(actual - test_pred)
  
  # Mean Absolute Scaled Error (MASE)
  # The forecast package uses the in-sample seasonal naive MAE as the denominator.
  # Since data is monthly, we use lag = 12.
  insample_naive_mae <- mean(abs(diff(train_df$y, lag = 12)))
  mase <- mae / insample_naive_mae
  
  # Theil's U: model RMSE / seasonal-naive RMSE
  # Seasonal naive benchmark: repeat the last 12 months of training
  # to cover the full test horizon (cycles if test > 12 months).
  # This matches the benchmark used by forecast::accuracy() in SARIMA/ETS/TBATS,
  # keeping cross-model comparisons on equal footing.
  train_y    <- train_df$y
  naive_seed <- tail(train_y, 12)                        # last 12 training months
  naive_pred <- rep(naive_seed, length.out = n_test)     # cycle to cover 25 months
  naive_rmse <- sqrt(mean((actual - naive_pred)^2))
  theils_u   <- rmse / naive_rmse
  
  metrics <- data.frame(
    RMSE     = round(rmse, 4),
    MAE      = round(mae,  4),
    MASE     = round(mase, 4),
    Theils_U = round(theils_u, 4),
    ME       = round(me, 4)
  )
  
  cat("\n--- Test Set Metrics ---\n")
  print(metrics)
  
  # ---- 4. Test forecast plot ----
  p_test <- ggplot() +
    geom_line(
      aes(x = test_df$ds, y = actual, colour = "Actual"),
      linewidth = 0.7
    ) +
    geom_line(
      aes(x = test_df$ds, y = test_pred, colour = "Prophet Forecast"),
      linewidth = 0.7, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Actual" = "#2c7bb6", "Prophet Forecast" = "#d7191c")) +
    theme_ieee() +
    labs(
      title  = paste("Prophet Forecast vs Actual (Test Set):", cat_name),
      x      = "Date",
      y      = "Sales Quantity (units)",
      colour = NULL
    )
  
  print(p_test)
  
  fname_test <- paste0(
    "prophet_01_testforecast_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  # ---- 5. Refit on full data + future forecast ----
  # Same model order, all 121 months of data — gives the forward forecast
  # the benefit of every available observation.
  m_full        <- prophet(
    yearly.seasonality = TRUE,
    weekly.seasonality = FALSE,
    daily.seasonality  = FALSE
  )
  m_full        <- fit.prophet(m_full, full_df)
  future_full   <- make_future_dataframe(m_full, periods = future_h, freq = "month")
  forecast_full <- predict(m_full, future_full)
  
  # Build future forecast plot manually (ggplot) so theme_ieee() applies cleanly
  # prophet's built-in plot() does not accept ggplot theme additions reliably
  hist_df   <- data.frame(ds = full_df$ds,          y    = full_df$y)
  future_df <- data.frame(
    ds    = forecast_full$ds,
    yhat  = forecast_full$yhat,
    lower = forecast_full$yhat_lower,
    upper = forecast_full$yhat_upper
  )
  
  p_future <- ggplot() +
    geom_ribbon(
      data = future_df,
      aes(x = ds, ymin = lower, ymax = upper),
      fill = "#a6cee3", alpha = 0.4
    ) +
    geom_line(
      data = hist_df,
      aes(x = ds, y = y, colour = "Historical"),
      linewidth = 0.6
    ) +
    geom_line(
      data = future_df,
      aes(x = ds, y = yhat, colour = "Prophet Forecast"),
      linewidth = 0.7, linetype = "dashed"
    ) +
    scale_colour_manual(
      values = c("Historical" = "black", "Prophet Forecast" = "#d7191c")
    ) +
    theme_ieee() +
    labs(
      title  = paste("Prophet Future Forecast (Next", future_h, "Months):", cat_name),
      x      = "Date",
      y      = "Sales Quantity (units)",
      colour = NULL
    )
  
  print(p_future)
  
  fname_future <- paste0(
    "prophet_02_futureforecast_",
    gsub("-", "", gsub(" ", "_", tolower(cat_name)))
  )
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  return(list(
    Model          = m,
    FinalModel     = m_full,
    TestForecast   = fc,
    FutureForecast = forecast_full,
    Metrics        = metrics,
    TestPlot       = p_test,
    FuturePlot     = p_future
  ))
}


# =============================================================================
# 4. EXECUTION & CROSS-CATEGORY COMPARISON
# =============================================================================

sink("Prophet_Detailed_Log.txt", split = TRUE)

res_2w <- run_prophet("2-Wheelers")
res_3w <- run_prophet("3-Wheelers")
res_4w <- run_prophet("4-Wheelers")

# Cross-category summary table
comparison <- bind_rows(
  res_2w$Metrics |> mutate(Category = "2-Wheelers"),
  res_3w$Metrics |> mutate(Category = "3-Wheelers"),
  res_4w$Metrics |> mutate(Category = "4-Wheelers")
) |> relocate(Category)

cat("\n--- Cross-Category Model Performance Comparison (Prophet) ---\n")
print(comparison)

sink()

cat("\n--- All details saved to 'Prophet_Detailed_Log.txt' ---\n")

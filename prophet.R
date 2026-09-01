# ============================================================================= 
# 04_prophet.R
# =============================================================================
source("preprocess.R")

# changepoint.prior.scale = 0.5 is FINAL (confirmed over the 0.05 default;
# sensitivity testing showed RMSE reductions of 10-28% across categories).
# The structural ~100% underforecasting bias (ME = MAE) persists at this
# setting regardless -- see report discussion.
CPS <- 0.5

run_prophet <- function(cat_name, data_list, future_h = 12) {
  cat(paste0("\n=========================================\n"))
  cat(paste0("ANALYZING CATEGORY: ", cat_name, "\n"))
  cat(paste0("=========================================\n"))
  
  train_df <- data_list[[cat_name]]$train
  test_df  <- data_list[[cat_name]]$test
  full_df  <- data_list[[cat_name]]$full
  
  # ---- 1. Fit model on training data ----
  set.seed(123)
  m <- prophet(
    yearly.seasonality      = TRUE,
    weekly.seasonality      = FALSE,
    daily.seasonality       = FALSE,
    interval.width           = 0.8,
    changepoint.prior.scale = CPS
  )
  m <- fit.prophet(m, train_df)
  
  # ---- 2. Predict on test period ----
  future    <- make_future_dataframe(m, periods = nrow(test_df), freq = "month")
  fc        <- predict(m, future)
  test_pred <- fc$yhat[(nrow(train_df) + 1):nrow(fc)]
  actual    <- test_df$y
  
  # ---- 3. Metrics ----
  metrics <- compute_metrics(
    train_y   = train_df$y,
    actual    = actual,
    predicted = test_pred
  )
  
  cat("\n--- Test Set Metrics ---\n")
  print(metrics)
  
  # ---- 4. Test forecast plot ----
  p_test <- ggplot() +
    geom_line(aes(x = test_df$ds, y = actual, colour = "Actual"), linewidth = 0.7) +
    geom_line(aes(x = test_df$ds, y = test_pred, colour = "Prophet Forecast"),
              linewidth = 0.7, linetype = "dashed") +
    scale_colour_manual(values = c("Actual" = "#2c7bb6", "Prophet Forecast" = "#d7191c")) +
    theme_ieee() +
    labs(title = paste("Prophet Forecast vs Actual (Test Set):", cat_name),
         x = "Date", y = "Sales Quantity (units)", colour = NULL)
  print(p_test)
  
  fname_test <- paste0("prophet_01_testforecast_",
                       gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_test, p_test, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  # ---- 5. Refit on full data + future forecast ----
  set.seed(123)
  m_full <- prophet(
    yearly.seasonality      = TRUE,
    weekly.seasonality      = FALSE,
    daily.seasonality       = FALSE,
    interval.width           = 0.8,
    changepoint.prior.scale = CPS
  )
  m_full        <- fit.prophet(m_full, full_df)
  future_full   <- make_future_dataframe(m_full, periods = future_h, freq = "month")
  forecast_full <- predict(m_full, future_full)
  
  hist_df   <- data.frame(ds = full_df$ds, y = full_df$y)
  future_df <- data.frame(
    ds    = forecast_full$ds,
    yhat  = forecast_full$yhat,
    lower = forecast_full$yhat_lower,
    upper = forecast_full$yhat_upper
  )
  
  p_future <- ggplot() +
    geom_ribbon(data = future_df, aes(x = ds, ymin = lower, ymax = upper),
                fill = "#a6cee3", alpha = 0.4) +
    geom_line(data = hist_df, aes(x = ds, y = y, colour = "Historical"), linewidth = 0.6) +
    geom_line(data = future_df, aes(x = ds, y = yhat, colour = "Prophet Forecast"),
              linewidth = 0.7, linetype = "dashed") +
    scale_colour_manual(values = c("Historical" = "black", "Prophet Forecast" = "#d7191c")) +
    theme_ieee() +
    labs(title = paste("Prophet Future Forecast (Next", future_h, "Months):", cat_name),
         x = "Date", y = "Sales Quantity (units)", colour = NULL)
  print(p_future)
  
  fname_future <- paste0("prophet_02_futureforecast_",
                         gsub("-", "", gsub(" ", "_", tolower(cat_name))))
  save_ieee(fname_future, p_future, width = IEEE_W_DOUBLE, height = IEEE_H)
  
  return(list(
    Model = m, FinalModel = m_full, TestForecast = fc,
    FutureForecast = forecast_full, Metrics = metrics,
    TestPlot = p_test, FuturePlot = p_future
  ))
}

sink("Prophet_Detailed_Log.txt", split = TRUE)

prophet_2w <- run_prophet("2-Wheelers", prophet_data)
prophet_3w <- run_prophet("3-Wheelers", prophet_data)
prophet_4w <- run_prophet("4-Wheelers", prophet_data)

comparison_prophet <- bind_rows(
  prophet_2w$Metrics |> mutate(Category = "2-Wheelers"),
  prophet_3w$Metrics |> mutate(Category = "3-Wheelers"),
  prophet_4w$Metrics |> mutate(Category = "4-Wheelers")
) |> relocate(Category) |> mutate(Model = "Prophet")

cat("\n--- Cross-Category Model Performance Comparison (Prophet) ---\n")
print(comparison_prophet)

# Save the RDS for the final step
saveRDS(comparison_prophet, "comparison_prophet.rds")

sink()
cat("\n--- Prophet details saved to 'Prophet_Detailed_Log.txt' ---\n")
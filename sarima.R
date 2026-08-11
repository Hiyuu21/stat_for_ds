library(tidyverse)
library(lubridate)
library(tseries)
library(forecast)

# -----------------------------
# Data Loading & Cleaning
# -----------------------------
data <- read.csv("EV_Dataset.csv")

cat("\n--- Data Exploration ---\n")
cat("Total Rows: ", nrow(data), "\n")
cat("Total Columns: ", ncol(data), "\n")
cat("Total Missing Values: ", sum(is.na(data)), "\n")
cat("Total Duplicated Rows: ", sum(duplicated(data)), "\n")
cat("Negative Sales Records: ", sum(data$EV_Sales_Quantity < 0), "\n")

# Zero-sales discovery
zero_sales <- sum(data$EV_Sales_Quantity == 0)
cat("Zero-Sales Records: ", zero_sales, "(", round(zero_sales/nrow(data)*100,2),"% of raw data )\n")


clean_data <- data |>
  mutate(Date = dmy(Date)) |>
  filter(Vehicle_Category %in% c("2-Wheelers","3-Wheelers","4-Wheelers")) |>
  group_by(Date, Vehicle_Category) |>
  summarise(Total_Sales = sum(EV_Sales_Quantity, na.rm = TRUE), .groups = 'drop') |>
  arrange(Date)

# -----------------------------
# Exploratory Data Analysis
# -----------------------------

# --- 1 Time series line plot: shows trend + seasonality visually
#     Justification for using SARIMA

ggplot(clean_data, aes(x=Date, y=Total_Sales, color=Vehicle_Category))+
  geom_line(linewidth=0.8)+
  theme_minimal() +
  labs(
    title = "Monthly EV Sales Trend by Vehicle Category (2014-2024)",
    x="Date", y="Monthly Sales Quantity", color="Category"
  )

# --- 2 Check outliers with boxplot method

ggplot(data = clean_data, aes(x=Vehicle_Category, y=Total_Sales, fill = Vehicle_Category)) +
  geom_boxplot(alpha=0.7, outlier.colour = "red", outlier.size = 2)+
  theme_minimal()+
  labs(
    title="Distribution and Outliers of Monthly EV Sales",
    subtitle="Highlighting extreme market surges (Red Dots)",
    x="Vehicle Category",
    y="Monthly Sales Quantity"
  ) + 
  theme(legend.position = "none")

# --- 3 Seasonal decomposition (STL) 
#     visually separate trend, seasonal, and remainder components per category
#     supporting evidence for why seasonal terms

for(cat_name in c("2-Wheelers","3-Wheelers","4-Wheelers")){
  cat_ts <- clean_data |>
    filter(Vehicle_Category == cat_name) |>
    pull(Total_Sales) |>
    ts(start = c(2014,1), frequency = 12)
  
  decomp <- stl(cat_ts, s.window = "periodic")
  plot(decomp, main=paste("STL Decomposition:", cat_name))
}

# --- 2.4 Descriptive statistics summary 

summary_table <- clean_data |>
  group_by(Vehicle_Category) |>
  summarise(
    Total_Months = n(),
    Mean_Sales = round(mean(Total_Sales), 2),
    Median_Sales = median(Total_Sales),
    Std_Dev = round(sd(Total_Sales), 2),
    Min_Sales = min(Total_Sales),
    Max_Sales = max(Total_Sales)
  )

print(summary_table)

# -----------------------------------------------------------------------
# SARIMA Function (train/test validation + diagnostics + future forecast)
# -----------------------------------------------------------------------

sarima <- function(dataset, cat_name, future_h=12) {
  cat(paste0("\n======================================================\n"))
  cat(paste0("ANALYZING CATEGORY: ", cat_name, "\n"))
  cat(paste0("======================================================\n"))
  
  # filter and split
  cat_data <- dataset |> filter(Vehicle_Category == cat_name)
  train_data <- cat_data |> slice(1:96)
  test_data <- cat_data |> slice(97:n())
  
  full_ts <- ts(cat_data$Total_Sales, start = c(2014,1), frequency = 12)
  train_ts <- ts(train_data$Total_Sales, start = c(2014,1), frequency = 12)
  test_ts <- ts(test_data$Total_Sales, start = c(2022,1), frequency = 12)
  
  
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

results_2w <- sarima(clean_data, "2-Wheelers")
results_3w <- sarima(clean_data, "3-Wheelers")
results_4w <- sarima(clean_data, "4-Wheelers")

# -------------------------------
# Cross-category comparison table
# -------------------------------
comparison <- bind_rows(
  as.data.frame(t(results_2w$Accuracy["Test set", c("RMSE", "MAE", "MAPE", "Theil's U")])) |>
    mutate(Category = "2-Wheelers"),
  as.data.frame(t(results_3w$Accuracy["Test set", c("RMSE", "MAE", "MAPE", "Theil's U")])) |>
    mutate(Category = "3-Wheelers"),
  as.data.frame(t(results_4w$Accuracy["Test set", c("RMSE", "MAE", "MAPE", "Theil's U")])) |>
    mutate(Category = "4-Wheelers"),
) |> relocate(Category)

cat("\nCross-Category Model Performance Comparison\n")
print(comparison)
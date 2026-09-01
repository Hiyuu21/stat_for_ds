library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(prophet)
library(ggplot2)

# --- SECTION 1: IEEE PLOT SETTINGS & SHARED HELPERS ---
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

# --- SECTION 2: DATA LOADING & CLEANING ---
url <- "https://raw.githubusercontent.com/Hiyuu21/stat_for_ds/main/EV_Dataset.csv"
raw_data <- read.csv(url)

clean_data <- raw_data |>
  mutate(Date = mdy(Date)) |>
  filter(Vehicle_Category %in% c("2-Wheelers", "3-Wheelers", "4-Wheelers")) |>
  group_by(Date, Vehicle_Category) |>
  summarise(Total_Sales = sum(EV_Sales_Quantity, na.rm = TRUE), .groups = "drop") |>
  arrange(Date)

# --- SECTION 3: EDA ---
summary_table <- clean_data |>
  group_by(Vehicle_Category) |>
  summarise(Total_Months = n(), Mean_Sales = round(mean(Total_Sales), 2), Median_Sales = median(Total_Sales), Std_Dev = round(sd(Total_Sales), 2), Min_Sales = min(Total_Sales), Max_Sales = max(Total_Sales))

p_ts <- ggplot(clean_data, aes(x = Date, y = Total_Sales, colour = Vehicle_Category)) +
  geom_line(linewidth = 0.5) + facet_wrap(~ Vehicle_Category, scales = "free_y", ncol = 1) +
  theme_ieee() + labs(title = "Monthly EV Sales by Vehicle Category (India, 2014–2024)", x = "Date", y = "Total Sales (units)") + theme(legend.position = "none")
save_ieee("eda_01_timeseries_all", p_ts, width = IEEE_W_DOUBLE, height = 5.5)

# --- SECTION 4 & 5: TRAIN / TEST SPLIT & TS/PROPHET FORMATTING ---
train_n <- 96

split_category <- function(dataset, cat_name, train_n = 96) {
  cat_data <- dataset |> filter(Vehicle_Category == cat_name)
  list(full = cat_data, train = cat_data |> slice(1:train_n), test = cat_data |> slice((train_n + 1):n()))
}

splits <- setNames(lapply(categories, split_category, dataset = clean_data, train_n = train_n), categories)

make_ts <- function(df, start_year = 2014, start_month = 1, freq = 12) {
  ts(df$Total_Sales, start = c(start_year, start_month), frequency = freq)
}

ts_data <- lapply(splits, function(s) {
  list(full = make_ts(s$full), train = make_ts(s$train), test = make_ts(s$test, start_year = 2022, start_month = 1))
})

make_prophet_df <- function(df) {
  df |> mutate(ds = floor_date(Date, "month")) |> transmute(ds, y = Total_Sales)
}

prophet_data <- lapply(splits, function(s) {
  list(full = make_prophet_df(s$full), train = make_prophet_df(s$train), test = make_prophet_df(s$test))
})
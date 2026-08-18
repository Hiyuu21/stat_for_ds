library(tidyverse)
library(lubridate)

# ============================================================
# Shared Preprocessing Script
# ============================================================

# ---- 1. Load & clean ----

url <- "https://raw.githubusercontent.com/Hiyuu21/stat_for_ds/main/EV_Dataset.csv"
data <- read.csv(url)
#data <- read.csv("EV_Dataset.csv")  use this if github unavailable 

clean_data <- data |>
  mutate(Date = dmy(Date)) |>
  filter(Vehicle_Category %in% c("2-Wheelers", "3-Wheelers", "4-Wheelers")) |>
  group_by(Date, Vehicle_Category) |>
  summarise(Total_Sales = sum(EV_Sales_Quantity, na.rm = TRUE), .groups = 'drop') |>
  arrange(Date)


# ---- 2. Per-category train/test split  ----

# Using the same split boundary across all four models 
# to make the model-comparison table (MAE/RMSE/MAPE/Theil's U) valid to compare

split_category <- function(dataset, cat_name, train_n = 96) {
  cat_data <- dataset |> filter(Vehicle_Category == cat_name)
  list(
    full  = cat_data,
    train = cat_data |> slice(1:train_n),
    test  = cat_data |> slice((train_n + 1):n())
  )
}

categories <- c("2-Wheelers", "3-Wheelers", "4-Wheelers")
splits <- setNames(
  lapply(categories, split_category, dataset = clean_data),
  categories
)

# ============================================================
# 3a. ts-object format - for SARIMA, ETS, TBATS
# These three use the same base R ts() object.
# ============================================================
make_ts <- function(df, start_year = 2014, start_month = 1, freq = 12) {
  ts(df$Total_Sales, start = c(start_year, start_month), frequency = freq)
}

ts_data <- lapply(splits, function(s) {
  list(
    full  = make_ts(s$full),
    train = make_ts(s$train),
    # test set here starts Jan 2022 given a 96-month (8-year) train split
    # from a Jan-2014 start - adjust start_year if you change train_n above
    test  = make_ts(s$test, start_year = 2022)
  )
})

# ============================================================
# 3b. Data-frame format - for Prophet
# ============================================================
make_prophet_df <- function(df) {
  df |> transmute(ds = Date, y = Total_Sales)
}

prophet_data <- lapply(splits, function(s) {
  list(
    full  = make_prophet_df(s$full),
    train = make_prophet_df(s$train),
    test  = make_prophet_df(s$test)
  )
})

# ---- 4. Save data ----
saveRDS(ts_data, "ts_data.rds")
saveRDS(prophet_data, "prophet_data.rds")

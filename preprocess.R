library(tidyverse)
library(lubridate)
library(forecast)
library(ggplot2)

# =============================================================================
# 1. LOAD & CLEAN
# =============================================================================

url <- "https://raw.githubusercontent.com/Hiyuu21/stat_for_ds/main/EV_Dataset.csv"
data <- read.csv(url)
# data <- read.csv("EV_Dataset.csv")  # use this if GitHub is unavailable

clean_data <- data |>
  mutate(Date = mdy(Date)) |>
  filter(Vehicle_Category %in% c("2-Wheelers", "3-Wheelers", "4-Wheelers")) |>
  group_by(Date, Vehicle_Category) |>
  summarise(Total_Sales = sum(EV_Sales_Quantity, na.rm = TRUE), .groups = "drop") |>
  arrange(Date)

# Quick sanity check — confirm date range and row count
cat("=== Data Overview ===\n")
cat("Date range :", as.character(min(clean_data$Date)), "to", as.character(max(clean_data$Date)), "\n")
cat("Categories :", paste(unique(clean_data$Vehicle_Category), collapse = ", "), "\n")
cat("Total rows :", nrow(clean_data), "\n\n")

# Confirm each category has the expected number of monthly observations
cat("=== Monthly Observations per Category ===\n")
clean_data |>
  group_by(Vehicle_Category) |>
  summarise(n_months = n(), min_date = min(Date), max_date = max(Date)) |>
  print()


# =============================================================================
# 2. EDA
# =============================================================================

categories <- c("2-Wheelers", "3-Wheelers", "4-Wheelers")

# ---- IEEE plot settings ----
IEEE_W_SINGLE <- 3.5   
IEEE_W_DOUBLE <- 7.16  
IEEE_H        <- 2.8  
IEEE_DPI      <- 300

# Base theme that matches IEEE style
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

# Helper to save a plot in IEEE format
# fname   : output filename (no extension), saved as .png
# plot    : ggplot object
# width   : figure width in inches (use IEEE_W_SINGLE or IEEE_W_DOUBLE)
# height  : figure height in inches
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


# ---- 2a. Overall time series plot (all 3 categories) ----
# Purpose: spot trend direction, growth timing, and scale differences
p_ts <- ggplot(clean_data, aes(x = Date, y = Total_Sales, colour = Vehicle_Category)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ Vehicle_Category, scales = "free_y", ncol = 1) +
  theme_ieee() +
  labs(
    title = "Monthly EV Sales by Vehicle Category (India, 2014–2024)",
    x = "Date", y = "Total Sales (units)"
  ) +
  theme(legend.position = "none")

print(p_ts)
save_ieee("eda_01_timeseries_all", p_ts, width = IEEE_W_DOUBLE, height = 5.5)


# ---- 2b. Seasonal plot per category ----
# Purpose: visualise within-year seasonal patterns across different years
for (cat in categories) {
  cat_ts <- clean_data |>
    filter(Vehicle_Category == cat) |>
    pull(Total_Sales) |>
    ts(start = c(2014, 1), frequency = 12)
  
  # ggseasonplot returns a ggplot — we can add our theme on top
  p <- ggseasonplot(cat_ts, year.labels = TRUE, year.labels.left = TRUE) +
    theme_ieee() +
    labs(
      title = paste("Seasonal Plot:", cat),
      x = "Month", y = "Total Sales (units)"
    )
  
  print(p)
  
  fname <- paste0("eda_02_seasonal_", gsub("-", "", gsub(" ", "_", tolower(cat))))
  save_ieee(fname, p, width = IEEE_W_DOUBLE, height = IEEE_H)
}


# ---- 2c. STL decomposition per category ----
# Purpose: visually separate trend, seasonality, and remainder components;
# confirms that both trend and seasonal components exist before modelling
for (cat in categories) {
  cat_ts <- clean_data |>
    filter(Vehicle_Category == cat) |>
    pull(Total_Sales) |>
    ts(start = c(2014, 1), frequency = 12)
  
  stl_fit <- stl(cat_ts, s.window = "periodic")
  
  p <- autoplot(stl_fit) +
    theme_ieee() +
    labs(title = paste("STL Decomposition:", cat))
  
  print(p)
  
  fname <- paste0("eda_03_stl_", gsub("-", "", gsub(" ", "_", tolower(cat))))
  save_ieee(fname, p, width = IEEE_W_DOUBLE, height = 4.5)
}

# ---- 2d. Summary statistics per category ----
cat("\n=== Summary Statistics per Category ===\n")
clean_data |>
  group_by(Vehicle_Category) |>
  summarise(
    Min    = min(Total_Sales),
    Median = median(Total_Sales),
    Mean   = round(mean(Total_Sales), 1),
    Max    = max(Total_Sales),
    SD     = round(sd(Total_Sales), 1)
  ) |>
  print()


# =============================================================================
# 3. TRAIN / TEST SPLIT
# =============================================================================
# train_n = 96 months = Jan 2014 to Dec 2021 (8 years)
# test    = remaining months = Jan 2022 onwards (~25 months through Jan 2024)
# This boundary is fixed across ALL four models (SARIMA, ETS, TBATS, Prophet)
# so that cross-model metric comparisons (MAE/RMSE/MAPE/Theil's U) are valid.

train_n <- 96

split_category <- function(dataset, cat_name, train_n = 96) {
  cat_data <- dataset |> filter(Vehicle_Category == cat_name)
  list(
    full  = cat_data,
    train = cat_data |> slice(1:train_n),
    test  = cat_data |> slice((train_n + 1):n())
  )
}

splits <- setNames(
  lapply(categories, split_category, dataset = clean_data, train_n = train_n),
  categories
)

# Sanity check split boundaries
cat("\n=== Train/Test Split Verification ===\n")
for (cat in categories) {
  cat(cat, "|\n")
  cat("  Train:", as.character(min(splits[[cat]]$train$Date)),
      "to", as.character(max(splits[[cat]]$train$Date)),
      "(", nrow(splits[[cat]]$train), "months )\n")
  cat("  Test :", as.character(min(splits[[cat]]$test$Date)),
      "to", as.character(max(splits[[cat]]$test$Date)),
      "(", nrow(splits[[cat]]$test), "months )\n")
}


# =============================================================================
# 4a. TS FORMAT — for SARIMA, ETS, TBATS
# =============================================================================

make_ts <- function(df, start_year = 2014, start_month = 1, freq = 12) {
  ts(df$Total_Sales, start = c(start_year, start_month), frequency = freq)
}

ts_data <- lapply(splits, function(s) {
  list(
    full  = make_ts(s$full),
    train = make_ts(s$train),
    test  = make_ts(s$test, start_year = 2022, start_month = 1)
  )
})


# =============================================================================
# 4b. PROPHET FORMAT — ds/y data frames with correct monthly dates
# =============================================================================
# Key fix: ds must be a proper Date (first of each month), NOT a sequential int.
# We extract year and month from the Date column to rebuild it cleanly.

make_prophet_df <- function(df) {
  df |>
    mutate(ds = floor_date(Date, "month")) |>
    transmute(ds, y = Total_Sales)
}

prophet_data <- lapply(splits, function(s) {
  list(
    full  = make_prophet_df(s$full),
    train = make_prophet_df(s$train),
    test  = make_prophet_df(s$test)
  )
})

# Verify Prophet ds column looks correct
cat("\n=== Prophet Data ds Verification (first 3 + last 3 of train, 4-Wheelers) ===\n")
p4w_train <- prophet_data[["4-Wheelers"]]$train
print(rbind(head(p4w_train, 3), tail(p4w_train, 3)))

cat("\n=== Prophet Test Set (4-Wheelers) ===\n")
print(prophet_data[["4-Wheelers"]]$test)


# =============================================================================
# 5. SAVE
# =============================================================================

saveRDS(ts_data,      "ts_data.rds")
saveRDS(prophet_data, "prophet_data.rds")
cat("\n--- Saved: ts_data.rds and prophet_data.rds ---\n")

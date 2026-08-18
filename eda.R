library(tidyverse)
library(lubridate)
library(forecast)

# ============================================================
# 1. Data Loading & Quality Check
# ============================================================

url <- "https://raw.githubusercontent.com/Hiyuu21/stat_for_ds/main/EV_Dataset.csv"
#raw_data <- read.csv(url)
raw_data <- read.csv("EV_Dataset.csv")

cat("\n--- Data Quality Report (For Methodology Section) ---\n")
cat("Total Rows: ", nrow(raw_data), "\n")
cat("Total Missing Values: ", sum(is.na(raw_data)), "\n")
cat("Total Duplicated Rows: ", sum(duplicated(raw_data)), "\n")
cat("Zero-Sales Records: ", sum(raw_data$EV_Sales_Quantity == 0), "\n")

# ============================================================
# 2. Data Cleaning
# ============================================================
clean_data <- raw_data |>
  mutate(Date = dmy(Date)) |>
  filter(Vehicle_Category %in% c("2-Wheelers","3-Wheelers","4-Wheelers")) |>
  group_by(Date, Vehicle_Category) |>
  summarise(Total_Sales = sum(EV_Sales_Quantity, na.rm = TRUE), .groups = 'drop') |>
  arrange(Date)

# ============================================================
# 3. Descriptive Statistics Table 
# ============================================================
cat("\n--- Descriptive Statistics (For Data Analysis Section) ---\n")
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

# ============================================================
# 4. Essential Visualizations
# ============================================================

# Plot A: Overall Time Series Trend
trend_plot <- ggplot(clean_data, aes(x = Date, y = Total_Sales, color = Vehicle_Category)) +
  # 1. Make the line slightly thinner so it doesn't look bulky on a small chart
  geom_line(linewidth = 0.5) + 
  
  # 2. Set the base_size down to 8 (default is usually 11)
  theme_minimal(base_size = 8) + 
  
  labs(
    title = "Monthly EV Sales Trend (2014-2024)", # Shortened slightly to fit better
    x = "Date", 
    y = "Monthly Sales Quantity", 
    color = "Category"
  ) +
  
  # 3. Micro-manage the exact font sizes for the different text elements
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5), # Centered, bold title
    axis.title = element_text(size = 8),                             # Axis labels
    axis.text = element_text(size = 7),                              # The numbers on the axes
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm")                                # Make the legend lines smaller
  )

print(trend_plot)

# 4. Save the plot with a slightly taller height to give the legend room to breathe
ggsave(
  filename = "trend_plot_single.png", 
  plot = trend_plot, 
  width = 3.5,     # IEEE single-column width
  height = 3.0,    # Increased from 2.5 to 3.0 so the bottom legend isn't squished
  units = "in", 
  dpi = 600
)

cat("\nSaved: trend_plot_single.png (Rescaled for IEEE)\n")



# Plot B: STL Decomposition
categories <- c("2-Wheelers", "3-Wheelers", "4-Wheelers")

cat("\n--- Generating and Saving STL Decomposition Plots for All Categories ---\n")

for(cat_name in categories){
  
  # 1. Isolate the time series for the specific category
  cat_ts <- clean_data |>
    filter(Vehicle_Category == cat_name) |>
    pull(Total_Sales) |>
    ts(start = c(2014,1), frequency = 12)
  
  # 2. Perform STL decomposition
  decomp <- stl(cat_ts, s.window = "periodic")
  
  # 3. Create a clean, web-safe filename (e.g., "stl_decomp_2-Wheelers.png")
  file_name <- paste0("stl_decomp_", cat_name, ".png")
  
  # 4. Open the PNG device (Single IEEE column: 3.5 inches, 600 DPI)
  png(
    filename = file_name, 
    width = 3.5, 
    height = 4.5, 
    units = "in", 
    res = 600
  )
  
  # 5. Draw the plot
  plot(decomp, main=paste("STL Decomposition:", cat_name))
  
  # 6. Close and save
  dev.off()
  
  cat("Saved:", file_name, "\n")
}

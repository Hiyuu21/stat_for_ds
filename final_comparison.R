# =============================================================================
# 05_final_comparison.R
# =============================================================================
library(tidyverse)

# Load the saved accuracy metrics from all 4 model scripts
comparison_sarima  <- readRDS("comparison_sarima.rds")
comparison_ets     <- readRDS("comparison_ets.rds")
comparison_tbats   <- readRDS("comparison_tbats.rds")
comparison_prophet <- readRDS("comparison_prophet.rds")

# Bind them together and structure the dataframe
final_comparison <- bind_rows(
  comparison_sarima, comparison_ets, comparison_tbats, comparison_prophet
) |>
  relocate(Model, Category) |>
  arrange(Category, Model)

# Print to console
cat("\n=================================================\n")
cat("FINAL 4-MODEL COMPARISON (SARIMA / STL+ETS / TBATS / Prophet)\n")
cat("=================================================\n")
print(final_comparison)

# Save to CSV
write.csv(final_comparison, "Final_Model_Comparison.csv", row.names = FALSE)
cat("\n--- Saved: Final_Model_Comparison.csv ---\n")
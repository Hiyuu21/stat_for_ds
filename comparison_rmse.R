comparison <- read.csv("Final_Model_Comparison_Differenced.csv")

p_rmse <- ggplot(comparison, aes(x = Category, y = RMSE, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c(
    "Prophet" = "#cb181d", "SARIMA" = "#2c7bb6",
    "STL+ETS" = "#41ab5d", "TBATS" = "#6a51a3"
  )) +
  theme_ieee() +
  labs(
    title = "RMSE Comparison Across Models by Category",
    x = "Category", y = "RMSE (units)", fill = NULL
  )

print(p_rmse)
save_ieee("rmse_comparison_all", p_rmse, width = IEEE_W_DOUBLE, height = IEEE_H)
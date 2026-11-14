# 02_eda.R ----------------------------------------------------------------
# Exploratory analysis: distribution of continuous variables, categorical
# counts, and a baseline characteristics table (Table 1) stratified by survival.

log_step("02_eda: distribution plots + Table 1")

# ---- Continuous variables: boxplots --------------------------------------
box_fun <- function(var, title) {
  ggplot(df_imp, aes(x = "", y = .data[[var]])) +
    geom_boxplot(fill = "#66C2A5", alpha = 0.6, outlier.colour = "red") +
    geom_jitter(width = 0.15, alpha = 0.15) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal()
}

p_boxes <- gridExtra::grid.arrange(
  box_fun("age", "Age (years)"),
  box_fun("pathsize", "Tumor size (cm)"),
  box_fun("lnpos", "Positive lymph nodes"),
  ncol = 3
)
save_fig(p_boxes, "distribution_continuous.png", w = 9, h = 3.6)

# ---- Categorical variables: bar charts -----------------------------------
bar_fun <- function(var, title) {
  ggplot(df_imp, aes(x = .data[[var]])) +
    geom_bar(fill = "#FC8D62", alpha = 0.85, color = "white") +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.4, size = 3) +
    labs(title = title, x = NULL, y = "Count") +
    theme_minimal()
}

p_bars <- gridExtra::grid.arrange(
  bar_fun("histgrad", "Histologic grade"),
  bar_fun("er", "ER status"),
  bar_fun("pr", "PR status"),
  bar_fun("lnpos_cat", "Positive lymph nodes (category)"),
  ncol = 2, nrow = 2
)
save_fig(p_bars, "distribution_categorical.png", w = 8, h = 6)

# ---- Table 1: baseline characteristics by survival status ----------------
table1 <- df_imp %>%
  select(age, pathsize, lnpos_cat, histgrad, er, pr, status) %>%
  tbl_summary(
    by = status,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
    label = list(age ~ "Age (years)",
                 pathsize ~ "Tumor size (cm)",
                 lnpos_cat ~ "Positive lymph nodes",
                 histgrad ~ "Histologic grade",
                 er ~ "ER status",
                 pr ~ "PR status")
  ) %>%
  add_overall() %>%
  add_p()

# Strip gtsummary's markdown markers and clean column names for a tidy CSV.
table1_clean <- as_tibble(table1) %>%
  mutate(across(everything(), ~ gsub("\\*\\*|__", "", .)))
names(table1_clean) <- c("Characteristic", "Overall (N=1,207)",
                         "Alive (N=1,135)", "Dead (N=72)", "p-value")

save_table(table1_clean, "table1_baseline.csv")

log_step("02_eda: done")

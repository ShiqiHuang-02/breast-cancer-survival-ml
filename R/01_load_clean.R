# 01_load_clean.R ----------------------------------------------------------
# Load the raw data, recode variables, diagnose missingness, and prepare two
# analysis branches:
#   (a) `df_imp`  -> MICE-imputed data for statistical inference (KM / Cox).
#   (b) `df`      -> cleaned data (still with NA) for the ML pipeline, where
#                    imputation is performed inside each CV fold to avoid leakage.

log_step("01_load_clean: loading data")

# Primary runtime input is the CSV export of the original SPSS (.sav) file.
raw <- read_csv(out_path("data/raw/breast_cancer_survival.csv"),
                show_col_types = FALSE)

stopifnot(all(c("age", "pathsize", "lnpos", "histgrad", "er", "pr",
                "status", "time") %in% names(raw)))

# ---- Recode & feature engineering ----------------------------------------
df <- raw %>%
  select(-id) %>%
  mutate(
    status     = factor(status, levels = c(0, 1), labels = c("Alive", "Dead")),
    status_num = as.integer(status == "Dead"),            # numeric 0/1 for Surv()
    histgrad   = factor(histgrad, levels = c(1, 2, 3), labels = c("G1", "G2", "G3")),
    er         = factor(er, levels = c(0, 1), labels = c("Negative", "Positive")),
    pr         = factor(pr, levels = c(0, 1), labels = c("Negative", "Positive")),
    lnpos_cat  = cut(lnpos, breaks = c(-1, 0, 3, Inf), labels = c("0", "1-3", ">3"))
  )

log_step(sprintf("01_load_clean: %d patients, %d events (%.1f%%)",
                 nrow(df), sum(df$status_num), 100 * mean(df$status_num)))

# ---- Missing-data diagnostics --------------------------------------------
miss_rate <- df %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_rate") %>%
  filter(missing_rate > 0) %>%
  arrange(desc(missing_rate))

p_missing <- ggplot(miss_rate, aes(x = reorder(variable, missing_rate), y = missing_rate)) +
  geom_col(fill = OKBLUE, alpha = 0.85) +
  geom_text(aes(label = scales::percent(missing_rate, accuracy = 0.1)),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, max(miss_rate$missing_rate) * 1.15)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = OKRED) +
  labs(title = "Missing data by variable",
       subtitle = "Red dashed line = 5% threshold; variables above it are imputed",
       x = NULL, y = "Proportion missing")
save_fig(p_missing, "missing_rate.png", w = 7, h = 4.2)

# ---- MICE multiple imputation (for the inference branch only) ------------
log_step("01_load_clean: MICE multiple imputation (m = 5)")
set.seed(SEED)
imp <- mice(df, m = 5, printFlag = FALSE)
df_imp <- complete(imp, action = 1)

log_step("01_load_clean: done")

# 06_interpret.R ----------------------------------------------------------
# Model interpretability on full-data fits:
#   - RSF variable importance (VIMP)
#   - GBM partial-dependence plots for the top continuous predictors
#   - XGBoost SHAP importance + dependence

log_step("06_interpret: full-data fits for interpretability")

# Impute once for the full-data interpretation fits (no CV here).
df_full <- impute_fit_apply(df_ml, df_ml)$train

# ---- RSF variable importance ---------------------------------------------
set.seed(SEED)
fit_rsf_full <- rfsrc(tree_formula, data = df_full, ntree = 500,
                      importance = "permute", seed = SEED)

vimp_df <- data.frame(
  variable = names(fit_rsf_full$importance),
  importance = as.numeric(fit_rsf_full$importance)
) %>%
  arrange(desc(importance))

p_vimp <- ggplot(vimp_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = OKBLUE, alpha = 0.85) +
  coord_flip() +
  labs(x = NULL, y = "VIMP (permutation importance)",
       title = "Random Survival Forest: variable importance") +
  theme_minimal()
save_fig(p_vimp, "rsf_importance.png", w = 7, h = 4)
save_table(vimp_df, "rsf_importance.csv")

# ---- GBM partial dependence ----------------------------------------------
set.seed(SEED)
fit_gbm_full <- gbm(tree_formula, data = df_full, distribution = "coxph",
                    n.trees = 1000, interaction.depth = 3, shrinkage = 0.01, cv.folds = 5)
best_iter <- gbm.perf(fit_gbm_full, method = "cv", plot.it = FALSE)

pdp_fun <- function(var, label) {
  partial(fit_gbm_full, pred.var = var, n.trees = best_iter,
          plot = FALSE, train = df_full) %>%
    ggplot(aes(x = !!sym(var), y = exp(yhat))) +
    geom_line(color = "#E64B35", linewidth = 1.1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    labs(x = label, y = "Hazard ratio (partial dependence)",
         title = paste("GBM partial dependence:", label)) +
    theme_minimal()
}

p_pdp <- gridExtra::grid.arrange(
  pdp_fun("pathsize", "Tumor size (cm)"),
  pdp_fun("lnpos", "Positive lymph nodes"),
  pdp_fun("age", "Age (years)"),
  ncol = 3
)
save_fig(p_pdp, "gbm_partial_dependence.png", w = 11, h = 3.8)

# ---- XGBoost SHAP --------------------------------------------------------
set.seed(SEED)
x_full <- tree_mm(df_full)
y_full <- ifelse(df_full$status_num == 1, df_full$time, -df_full$time)
fit_xgb_full <- xgb.train(params = list(objective = "survival:cox", eta = 0.05,
                                        max_depth = 3, subsample = 0.8),
                          data = xgb.DMatrix(x_full, label = y_full),
                          nrounds = 300, verbose = 0)

shp <- shapviz(fit_xgb_full, X_pred = x_full)

p_shap_imp <- sv_importance(shp, kind = "bar", show_numbers = TRUE) +
  labs(title = "XGBoost SHAP feature importance")
save_fig(p_shap_imp, "xgb_shap_importance.png", w = 7, h = 4)

p_shap_age <- sv_dependence(shp, v = "age") +
  labs(title = "XGBoost SHAP dependence: age")
save_fig(p_shap_age, "xgb_shap_age.png", w = 7, h = 4.5)

log_step("06_interpret: done")

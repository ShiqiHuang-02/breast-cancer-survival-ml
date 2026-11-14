# 05_evaluation.R ---------------------------------------------------------
# Model evaluation via repeated 5-fold cross-validation:
#   - Harrell's C-index (mean +/- SD over repeats)
#   - Time-dependent AUC (timeROC, averaged over repeats)
# Imputation and model fitting are confined to the training fold (no leakage).

log_step("05_evaluation: repeated CV")

df_ml <- df  # cleaned data, NA retained -> imputed within each fold
n <- nrow(df_ml)

# Collect out-of-fold risk scores for every repeat.
oof_list <- vector("list", N_REPEATS)
cindex_mat <- matrix(NA_real_, N_REPEATS, length(MODELS),
                     dimnames = list(NULL, MODEL_NAMES))

for (r in seq_len(N_REPEATS)) {
  set.seed(SEED + r)
  fold_id <- sample(rep(seq_len(N_FOLDS), length.out = n))
  oof <- matrix(NA_real_, n, length(MODELS), dimnames = list(NULL, MODEL_NAMES))

  for (k in seq_len(N_FOLDS)) {
    tr_idx <- which(fold_id != k)
    te_idx <- which(fold_id == k)
    split  <- impute_fit_apply(df_ml[tr_idx, ], df_ml[te_idx, ])

    for (mi in seq_along(MODELS)) {
      fit  <- MODELS[[mi]]$fit(split$train)
      oof[te_idx, mi] <- MODELS[[mi]]$pred(fit, split$test)
    }
  }

  for (mi in seq_along(MODELS)) {
    # survival::concordance treats higher predictor = LONGER survival, but our
    # scores are risk scores (higher = higher risk), so we negate them.
    risk <- -oof[, mi]
    cindex_mat[r, mi] <- concordance(Surv(df_ml$time, df_ml$status_num) ~ risk)$concordance
  }
  oof_list[[r]] <- oof
  log_step(sprintf("05_evaluation: repeat %d/%d done", r, N_REPEATS))
}

cindex_summary <- data.frame(
  model = MODEL_NAMES,
  c_index_mean = colMeans(cindex_mat),
  c_index_sd = apply(cindex_mat, 2, sd)
) %>%
  arrange(desc(c_index_mean))

save_table(cindex_summary, "cindex_cv.csv")

# ---- C-index comparison plot --------------------------------------------
p_cindex <- ggplot(cindex_summary, aes(x = reorder(model, c_index_mean),
                                       y = c_index_mean, fill = model)) +
  geom_col(alpha = 0.9, show.legend = FALSE) +
  geom_errorbar(aes(ymin = pmax(0, c_index_mean - c_index_sd),
                    ymax = pmin(1, c_index_mean + c_index_sd)),
                width = 0.25) +
  geom_text(aes(label = sprintf("%.3f", c_index_mean)), hjust = -0.15, size = 3.4) +
  coord_flip(ylim = c(0.55, 0.85)) +
  scale_fill_viridis_d(option = "D") +
  labs(x = NULL, y = "C-index (5x5-fold CV, mean +/- SD)",
       title = "Survival model comparison",
       subtitle = "Higher is better; error bars show SD across CV repeats") +
  theme_minimal()
save_fig(p_cindex, "model_comparison_cindex.png", w = 8.5, h = 4.5)

# ---- Time-dependent AUC --------------------------------------------------
tdauc_mean <- matrix(NA_real_, length(MODELS), length(AUC_TIMES),
                     dimnames = list(MODEL_NAMES, paste0("t", AUC_TIMES)))
for (mi in seq_along(MODELS)) {
  auc_per_repeat <- matrix(NA_real_, N_REPEATS, length(AUC_TIMES))
  for (r in seq_len(N_REPEATS)) {
    tr <- tryCatch(
      timeROC(T = df_ml$time, delta = df_ml$status_num, marker = oof_list[[r]][, mi],
              cause = 1, times = AUC_TIMES, iid = FALSE)$AUC,
      error = function(e) rep(NA_real_, length(AUC_TIMES))
    )
    auc_per_repeat[r, ] <- tr
  }
  tdauc_mean[mi, ] <- colMeans(auc_per_repeat, na.rm = TRUE)
}

tdauc_long <- as.data.frame(tdauc_mean) %>%
  mutate(model = rownames(.)) %>%
  pivot_longer(-model, names_to = "time", values_to = "auc") %>%
  mutate(time = as.numeric(gsub("t", "", time)))

p_tdauc <- ggplot(tdauc_long, aes(x = time, y = auc, color = model)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_viridis_d(option = "D") +
  labs(x = "Time (months)", y = "Time-dependent AUC",
       title = "Discrimination over time (timeROC)",
       subtitle = "AUC averaged over 5 CV repeats") +
  theme_minimal() +
  theme(legend.position = "right")
save_fig(p_tdauc, "time_dependent_auc.png", w = 8, h = 4.8)
save_table(tdauc_long, "tdauc.csv")

log_step("05_evaluation: done")

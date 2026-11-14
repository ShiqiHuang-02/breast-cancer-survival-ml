# 03_inference.R ----------------------------------------------------------
# Statistical inference branch:
#   - Kaplan-Meier survival curves by key prognostic groups
#   - Multivariable Cox proportional hazards model (MICE-pooled)
#   - Model diagnostics: collinearity (VIF) and proportional-hazards assumption
#   - Forest plot of hazard ratios

log_step("03_inference: KM curves + Cox model + diagnostics")

# ---- Kaplan-Meier curves -------------------------------------------------
df_km <- df_imp %>%
  mutate(
    pathsize_group = ifelse(pathsize <= median(pathsize, na.rm = TRUE),
                            "Small (<= median)", "Large (> median)"),
    age_group = ifelse(age <= 50, "<= 50 years", "> 50 years")
  )

km_lnpos <- ggsurvplot(
  survfit(Surv(time, status_num) ~ lnpos_cat, data = df_km), data = df_km,
  pval = TRUE, conf.int = TRUE, risk.table = TRUE,
  xlab = "Time (months)", ylab = "Survival probability",
  title = "KM curve by positive lymph nodes",
  legend.title = "LN positive", legend.labs = c("0", "1-3", ">3"),
  ggtheme = theme_minimal(), palette = "jco"
)
km_path <- ggsurvplot(
  survfit(Surv(time, status_num) ~ pathsize_group, data = df_km), data = df_km,
  pval = TRUE, conf.int = TRUE, risk.table = TRUE,
  xlab = "Time (months)", ylab = "Survival probability",
  title = "KM curve by tumor size",
  legend.title = "Tumor size",
  legend.labs = c("Small", "Large"),
  ggtheme = theme_minimal(), palette = "jco"
)
km_age <- ggsurvplot(
  survfit(Surv(time, status_num) ~ age_group, data = df_km), data = df_km,
  pval = TRUE, conf.int = TRUE, risk.table = TRUE,
  xlab = "Time (months)", ylab = "Survival probability",
  title = "KM curve by age group",
  legend.title = "Age", legend.labs = c("<= 50", "> 50"),
  ggtheme = theme_minimal(), palette = "jco"
)

save_fig(gridExtra::grid.arrange(km_lnpos$plot, km_path$plot, km_age$plot, ncol = 3),
         "kaplan_meier.png", w = 12, h = 4.5)

# ---- Multivariable Cox (MICE-pooled) -------------------------------------
cox_formula <- Surv(time, status_num) ~ age + pathsize + lnpos_cat + histgrad + er + pr

fit_full <- coxph(cox_formula, data = df_imp)

# ---- Collinearity (VIF / GVIF) -------------------------------------------
# VIF is computed on an equivalent linear model with the same predictors
# (a standard trick, since coxph has no intercept term).
fit_vif <- lm(as.numeric(time) ~ age + pathsize + lnpos_cat + histgrad + er + pr,
              data = df_imp)
v <- car::vif(fit_vif)
vif_df <- if (is.matrix(v)) {
  data.frame(variable = rownames(v), vif = v[, "GVIF^(1/(2*Df))"])
} else {
  data.frame(variable = names(v), vif = as.numeric(v))
}

p_vif <- ggplot(vif_df, aes(x = reorder(variable, vif), y = vif)) +
  geom_col(fill = OKBLUE, alpha = 0.85) +
  coord_flip() +
  geom_hline(yintercept = 5, linetype = "dashed", color = OKRED) +
  labs(x = NULL, y = "VIF (adjusted GVIF)", title = "Collinearity diagnostics (threshold = 5)") +
  theme_minimal()
save_fig(p_vif, "vif.png", w = 7, h = 4)
save_table(vif_df, "vif.csv")

# ---- Proportional hazards assumption -------------------------------------
test_ph <- cox.zph(fit_full)
ph_plots <- ggcoxzph(test_ph)
if (inherits(ph_plots, "ggplot")) {
  save_fig(ph_plots, "ph_assumption.png", w = 9, h = 6)
} else {
  save_fig(do.call(gridExtra::grid.arrange, c(ph_plots, list(ncol = 2))),
           "ph_assumption.png", w = 11, h = 9)
}

# ---- Pooled Cox across MICE datasets + forest plot -----------------------
# NOTE: the formula is written inline inside with() so its environment points
# to the imputed data (referencing a pre-built formula object fails here).
fit_mice <- with(imp, coxph(Surv(time, status_num) ~ age + pathsize + lnpos_cat +
                              histgrad + er + pr))
pooled <- pool(fit_mice)
res <- summary(pooled, conf.int = TRUE, exponentiate = TRUE)

plot_data <- res %>%
  mutate(
    sig = ifelse(p.value < 0.05, "Significant", "Not significant"),
    term_label = dplyr::recode(term,
      "age" = "Age (per year)",
      "pathsize" = "Tumor size (per cm)",
      "lnpos_cat1-3" = "LN: 1-3 (vs 0)",
      "lnpos_cat>3" = "LN: >3 (vs 0)",
      "histgradG2" = "Grade G2 (vs G1)",
      "histgradG3" = "Grade G3 (vs G1)",
      "erPositive" = "ER positive (vs negative)",
      "prPositive" = "PR positive (vs negative)"
    )
  )

p_forest <- ggplot(plot_data, aes(x = estimate, y = reorder(term_label, estimate), color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), width = 0.25,
                linewidth = 0.8, orientation = "y") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("HR %.2f", estimate)), vjust = -0.8, size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = c("Significant" = OKRED, "Not significant" = "gray60")) +
  labs(x = "Hazard ratio (95% CI)", y = NULL, color = NULL,
       title = "Multivariable Cox regression (pooled MICE)",
       subtitle = "Tumor size and lymph-node involvement drive prognosis") +
  theme_minimal() +
  theme(legend.position = "top")

save_fig(p_forest, "cox_forest.png", w = 8, h = 4.8)
save_table(plot_data, "cox_results.csv")

log_step("03_inference: done")

# 04_ml_models.R ----------------------------------------------------------
# Predictive-modelling branch (the data-science / ML core of the project).
# Defines a unified fit/predict interface for 7 survival models and a
# leak-free imputation step applied inside each CV fold.

log_step("04_ml_models: defining model zoo")

# ---- Leak-free imputation (fit on train, apply to train & test) ----------
impute_fit_apply <- function(train, test) {
  for (nm in c("pathsize", "histgrad", "er", "pr")) {
    if (is.numeric(train[[nm]])) {
      m <- median(train[[nm]], na.rm = TRUE)
      train[[nm]][is.na(train[[nm]])] <- m
      test[[nm]][is.na(test[[nm]])] <- m
    } else {
      m <- names(sort(table(train[[nm]], useNA = "no"), decreasing = TRUE))[1]
      train[[nm]][is.na(train[[nm]])] <- m
      test[[nm]][is.na(test[[nm]])] <- m
    }
  }
  list(train = train, test = test)
}

# ---- Feature matrices ----------------------------------------------------
# Linear models use the clinically-standard categorical cutpoints for lymph
# nodes; tree/boosting models use the raw continuous count to retain full info.
lin_mm  <- function(d) model.matrix(~ age + pathsize + lnpos_cat + histgrad + er + pr, data = d)[, -1]
tree_mm <- function(d) model.matrix(~ age + pathsize + lnpos + histgrad + er + pr, data = d)[, -1]

lin_formula  <- Surv(time, status_num) ~ age + pathsize + lnpos_cat + histgrad + er + pr
tree_formula <- Surv(time, status_num) ~ age + pathsize + lnpos + histgrad + er + pr

# ---- Cox PH --------------------------------------------------------------
fit_cox <- function(d) coxph(lin_formula, data = d)
pred_cox <- function(f, d) predict(f, newdata = d, type = "lp")

# ---- Regularized Cox (glmnet) --------------------------------------------
glm_fit <- function(alpha) {
  function(d) {
    cv.glmnet(lin_mm(d), Surv(d$time, d$status_num),
              family = "cox", alpha = alpha, nfolds = 5, cox.ties = "efron")
  }
}
glm_pred <- function(f, d) as.numeric(predict(f, newx = lin_mm(d), s = "lambda.min", type = "link"))

# ---- Gradient Boosting Machine (Cox) -------------------------------------
fit_gbm <- function(d) {
  gbm(tree_formula, data = d, distribution = "coxph",
      n.trees = GBM_TREES, interaction.depth = 3, shrinkage = 0.01, cv.folds = 0)
}
pred_gbm <- function(f, d) predict(f, newdata = d, n.trees = GBM_TREES, type = "link")

# ---- Random Survival Forest ----------------------------------------------
fit_rsf <- function(d) {
  rfsrc(tree_formula, data = d, ntree = RSF_TREES, seed = SEED, importance = "none")
}
pred_rsf <- function(f, d) predict(f, newdata = d)$predicted

# ---- XGBoost (Cox objective) ---------------------------------------------
fit_xgb <- function(d) {
  y <- ifelse(d$status_num == 1, d$time, -d$time)   # event=+time, censored=-time
  dtrain <- xgb.DMatrix(tree_mm(d), label = y)
  xgb.train(params = list(objective = "survival:cox", eta = 0.1, max_depth = 3,
                          subsample = 0.8, colsample_bytree = 0.8),
            data = dtrain, nrounds = XGB_ROUNDS, verbose = 0)
}
pred_xgb <- function(f, d) predict(f, newdata = tree_mm(d))

# ---- Model registry ------------------------------------------------------
MODELS <- list(
  list(name = "Cox PH",            fit = fit_cox,            pred = pred_cox),
  list(name = "Cox Lasso",         fit = glm_fit(1),         pred = glm_pred),
  list(name = "Cox Ridge",         fit = glm_fit(0),         pred = glm_pred),
  list(name = "Cox ElasticNet",    fit = glm_fit(0.5),       pred = glm_pred),
  list(name = "GBM",               fit = fit_gbm,            pred = pred_gbm),
  list(name = "Random Survival Forest", fit = fit_rsf,       pred = pred_rsf),
  list(name = "XGBoost (Cox)",     fit = fit_xgb,            pred = pred_xgb)
)
MODEL_NAMES <- vapply(MODELS, function(m) m$name, character(1))

log_step("04_ml_models: done")

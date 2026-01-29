get_y0_features <- function(dt) {
  pats <- c(
    "^h_y_l\\d+$",     # lagged outcomes
    "^ybar_\\d+$",      # within-unit means
    "^ydot_\\d+$",      # within-unit slopes
    "^czero_\\d+$",     # zero shares
    "_nbr$",            # spatial neighbors
#    "_coh$",            # cohort averages
    "^ybar_ctrl$"       # contemporaneous control mean
  )

  unique(unlist(lapply(
    pats, \(p) grep(p, names(dt), value = TRUE)
  )))

}


# Add hazard-based propensity scores to a data.table (in place).
# Trains a discrete-time hazard model on the risk set and
# predicts Pr(treated next year) for all unit–year rows.
# Supports ranger or xgboost backends.
get_hazard <- function(dt
  , engine = c("ranger", "xgboost")
  , predictors = NULL
  , y_col = "y"
  , year_col = "year"
  , start_year_col = "start_year"
  , at_risk_col = "n"
  , id_col = "i"
  , risk_min_year = NULL
  # ranger
  , ranger_num_trees = 500
  # xgb
  , train_id_frac = 0.5
  , xgb_params = NULL
  , xgb_nrounds = 200
  , xgb_early_stop = 5
  # shared
  , seed = 42
  , num_threads = 1
  , verbose = TRUE
  , conf = getOption("xpn.conf")
) {

  engine <- match.arg(engine)
  Lset   <- conf$Lset

  # ---- decide predictors here (and only here) ----
  basic_pred <- c(year_col, at_risk_col, get_y0_features(dt))  |> unique()
  preds <- c(basic_pred, NULL)  |> unique()
  if (!is.null(predictors)) preds <- c(preds, predictors) |> unique()
  # avoid obvious leakage / junk even if user passes them accidentally
  preds <- setdiff(preds, c("w", "cpr", "start_year", "i"))

  miss <- setdiff(c(preds), names(dt))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  full_dt <- dt[, ..preds]

  train_dt <- prep_hazard_data(
    full_dt,
    cohort = dt[[start_year_col]],
    Lset = Lset,
    year_col = year_col,
    risk_min_year = risk_min_year
  )

  if (engine == "ranger") {
    fit_model <- fit_hazard_ranger(
      train_dt
      , num_trees = ranger_num_trees, seed = seed, num_threads = num_threads
    )
    cpr <- predict_hazard_ranger(fit_model, full_dt, num_threads)
  }

  if (engine == "xgboost") {
    fit_model <- fit_hazard_xgb(train_dt
      , id_col = id_col
      , train_id_frac = train_id_frac, seed = seed, num_threads = num_threads
      , params = xgb_params, nrounds = xgb_nrounds, early_stop = xgb_early_stop
    )
    cpr <- predict_hazard_xgb(fit_model, full_dt)
  }
  cpr
}

# Prepare person–year hazard training data.
# Builds the at-risk set (year < start_year), defines the event
# w = 1[treated next year], and selects predictor columns.
# Returns the training table plus predictor names.
prep_hazard_data <- function(full_dt
  , cohort
  , Lset = NULL
  , year_col = "year"
  , risk_min_year = NULL
) {
  stopifnot(data.table::is.data.table(full_dt))

  min_y <- min(full_dt[[year_col]], na.rm = TRUE)
  # derive risk_min_year from Lset if not supplied
  if (is.null(risk_min_year)) {
    Lmax <- if (is.null(Lset)) 0L else max(as.integer(Lset), na.rm = TRUE)
    risk_min_year <- min_y + Lmax
  }

  # risk set + outcome
  idx <- full_dt[[year_col]] < cohort & full_dt[[year_col]] >= risk_min_year
  sdt <- full_dt[idx, ]
  w <- as.integer(is.finite(cohort[idx]) &
      sdt[[year_col]] == (cohort[idx] - 1L)
  )
  sdt[, w := factor(w, levels = c(0, 1))]

  sdt
}


# Fit a discrete-time hazard model using ranger.
# Trains a probabilistic random forest on the risk set,
# including calendar year to learn the baseline hazard.
# Returns a fitted ranger model.
fit_hazard_ranger <- function(train_dt
  , num_trees = 500
  , seed = 1
  , num_threads = 1
) {

  if (!is.null(seed)) set.seed(seed)
  ranger::ranger(
    w ~ .,
    data = train_dt,
    probability = TRUE,
    num.trees = num_trees,
    num.threads = num_threads
  )
}



#sdt <- prep_hazard_data(dt)$sdt
#predictors <- prep_hazard_data(dt)$predictors

# Fit a discrete-time hazard model using XGBoost.
# Splits training by unit ID (if available), uses a validation
# watchlist for early stopping, and stores design-matrix metadata.
# Returns the fitted model and encoding information.
#sdt <- train_dt
fit_hazard_xgb <- function(sdt,
                           id_col = "i",
                           train_id_frac = 0.5,
                           seed = 42,
                           num_threads = 1,
                           params = NULL,
                           nrounds = 200,
                           early_stop = 5) {
  id <- id_col

  set.seed(seed)

  # split by ID if possible
  if (!is.null(id) && id %in% names(sdt)) {
    uid       <- unique(sdt[[id]])
    n_train   <- as.integer(train_id_frac * length(uid))
    train_ids <- sample(uid, n_train)
    is_train  <- sdt[[id]] %in% train_ids
  } else {
    is_train <- runif(nrow(sdt)) < train_id_frac
  }


  tr <- sdt[is_train]
  va <- sdt[!is_train]
  dim(tr)
  dim(va)


  dtrain <- xgboost::xgb.DMatrix(tr[, !("w")], label = as.numeric(tr$w) - 1)
  dvalid <- xgboost::xgb.DMatrix(va[, !("w")], label = as.numeric(va$w) - 1)

  if (is.null(params)) {
    params <- list(
      objective = "binary:logistic",
      eta = 0.1,
      max_depth = 6,
      nthread = num_threads,
      tree_method = "hist"
    )
  }

  xgboost::xgb.train(
    params = params,
    data = dtrain,
    evals = list(train = dtrain, eval = dvalid),
    nrounds = nrounds,
    early_stopping_rounds = early_stop,
    verbose = 1
  )

}

# Predict hazard probabilities using a fitted ranger model.
# Applies the model to all rows (not only the risk set),
# leaving enforcement of the treatment process to the simulator.
predict_hazard_ranger <- function(
  mp, full_dt, num_threads = 1
) {
  predict(mp, data = full_dt, num.threads = num_threads)$predictions[, "1"]
}

# Predict hazard probabilities from a fitted XGBoost model.
# Reconstructs the design matrix with training-time columns,
# aligning factors and missing features before prediction.
predict_hazard_xgb <- function(fit, full_dt) {

  x <- xgboost::xgb.DMatrix(full_dt)

  best_iter <- fit$best_iteration
  if (is.null(best_iter) || best_iter < 1) best_iter <- fit$niter

  pred <- predict(fit, x, iteration_range = c(1, best_iter))
  pred
}








if (FALSE) {
  library(data.table)
  library(Matrix)
  library(ggplot2)

  set.seed(1)
  N <- 10
  Tt <- 10
  years <- 2001:2010

  dt <- data.table(
    i = rep(1:N, each = Tt),
    year = rep(years, times = N)
  )

  # unit-specific start year (Inf = never treated)
  start_per_id <- rep(2006:2010, length.out = N)
  start_per_id[start_per_id <= 2008] <- Inf
  dt[, start_year := rep(start_per_id, each = Tt)]

  # covariates
  dt[, x_unit := rnorm(N)[i]]
  dt[, x_time := 0.2 * (year - min(years)) + rnorm(.N, sd = 0.5)]
  dt[, x_bin  := rbinom(.N, 1, plogis(x_unit + 0.1 * x_time))]
  dt[, n  := 20L]

  # ring adjacency G (i matches row/col)
  G <- Matrix(0, nrow = N, ncol = N, sparse = TRUE)
  for (k in 1:N) {
    G[k, ifelse(k == 1, N, k - 1)] <- 1
    G[k, ifelse(k == N, 1, k + 1)] <- 1
  }
  diag(G) <- 0

  # ---- simulate an "untreated outcome" y0 with mild spatial correlation ----
  # baseline linear predictor
  alpha_i <- rnorm(N, sd = 0.4)                 # unit effect
  beta_t  <- 0.05 * (dt$year - min(years))      # common time trend

  lp0 <- 0.3 + alpha_i[dt$i] + beta_t + 0.4 * dt$x_unit + 0.2 * dt$x_time + 0.2 * dt$x_bin

  # add a simple neighbor component: mean neighbor x_unit (time-invariant, easy)
  neigh_xunit <- as.numeric(G %*% alpha_i) / pmax(as.numeric(rowSums(G)), 1)
  lp0 <- lp0 + 0.3 * neigh_xunit[dt$i]

  # Poisson outcome
  dt[, y := rpois(.N, lambda = exp(lp0))]

  # ---- set y0 missing after treatment starts (observed Y(0) only when untreated) ----
  conf <- list(Lset = choose_lset(dt))
  options(xpn.conf = conf)
  dt |> expand_y()

  # ---- feature expansion on y0 ----

  dt |> get_hazard(engine = "ranger") |> print()
  dt |> get_hazard(engine = "xgboost") |> print()

}

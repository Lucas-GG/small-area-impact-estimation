library(xgboost)

pre_xgb <- function(
  .dt
  , shuffle_sample = TRUE
  , wtype =  "exp"
  , shuffle_time = FALSE
  , predictors = NULL
  , year_col = "year"
  , start_year_col = "start_year"
  , at_risk_col = "n"
  , id_col = "i"
) {

  xgb_dt <- copy(.dt)
  # Step 1: simulate new assigment
  A_sim <- xgb_dt |> simulate_assignments()
  xgb_mask <- mask_pre(xgb_dt, A_sim)

  # Step 2: build assigment dependent features
  art <- xgb_dt |> make_y0_context(A_sim = A_sim) |> materialize()
  xgb_dt <- xgb_dt[, names(art) := art]
  xgb_dt[, start_year := factor(start_year)]
  head(xgb_dt)

  # Step 3: bootstrap weights (requieres id_col)
  if (shuffle_sample) {
    wts <- xgb_dt |> cluster_boot(wtype = wtype, shuffle_time = shuffle_time)
    xgb_dt[, wts := wts]
  }

  # Step 4: Clean up columns
  keep_cols <- c(year_col, at_risk_col, start_year_col, "y0", "wts")
  keep_cols <- c(keep_cols, get_y0_features(xgb_dt))  |> unique()
  keep_cols <- c(keep_cols, predictors) |> unique()

  cat("predictors are", setdiff(keep_cols, c("y0", "wts")), "\n")
  xgb_dt <- xgb_dt[, ..keep_cols]
  list(data = xgb_dt, mask = xgb_mask)
}


fit_xgb <- function(
  xgb_dt
  , is_train
  , num_threads = 1
  , params = NULL
  , nrounds = 1000
  , early_stop = 20
  , predictors = NULL
  , ptry = NULL
) {
  k <- ncol(xgb_dt)
  if (is.null(ptry)) ptry <- sqrt(length(k)) / length(k)

  # Step 5: Fit XGBoost Poisson model
  tr <- xgb_dt[is_train]
  va <- xgb_dt[!is_train & !is.na(y0)]

  dtrain <- xgboost::xgb.DMatrix(tr[, !("y0")], label = tr$y0
    , base_margin = log(tr$n), weight = tr$wts
  )

  dvalid <- xgboost::xgb.DMatrix(va[, !("y0")], label = va$y0
    , base_margin = log(va$n), weight = va$wts
  )

  if (is.null(params)) {
    params <- list(
      objective = "count:poisson"
      , eval_metric = "poisson-nloglik"
      , tree_method = "hist"
      , max_depth = 3 # #defautl 6
      , eta = 0.03 #default is .3
      , subsample = 1 # no row subsample
      , nthread = num_threads
      , colsample_bynode = ptry
      , min_child_weight = 5 #minimum effective sample size
      , lambda = 10 #L2 regularization default to 1
      #, gamma = 1 #minimum loss reduction required to make a spli, defualt 0
      #, max_delta_step = 1 #defualt to 0
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

predict_xgb  <- function(fit, xgb_dt) {
  x <- xgboost::xgb.DMatrix(xgb_dt[, !("y0")], base_margin = log(xgb_dt$n))
  best_iter <- fit$best_iteration
  if (is.null(best_iter) || best_iter < 1) best_iter <- fit$niter
  predict(fit, x, iteration_range = c(1, best_iter))
}

single_impute_xgb <- function(.dt
  , ...
  , verbose = TRUE
) {
  # Step 1: Fit
  if (verbose) cat("Building boosted ensemble\n")
  xgb_out <- pre_xgb(.dt, ...)

  xgb_fit <- xgb_out$data |> fit_xgb(is_train = xgb_out$mask, ...)

  # Step 2: Predict
  predict_xgb(xgb_fit, xgb_out$data)
}

#dt |> single_impute_xgb()

multiple_impute_xgb <- function(
  .dt
  , ...
  , verbose = TRUE
  , nboot = 20
  , ecores = 1
) {

  # Start timer for performance tracking
  start_time <- Sys.time()
  if (verbose) cat("Starting impute_xgb with", nboot, "bootstraps\n")

  full_dt <- copy(.dt)
  if (verbose) cat("Adding history...\n")
  full_dt |> add_y0_history()

  # Step 3: compute features that depend on assignment
  if (verbose) cat("Adding context...\n")
  art <- full_dt |> make_y0_context() |> materialize()
  full_dt[, names(art) := art]

  if (verbose) cat("Adding hazard...\n")
  cpr <- get_hazard(full_dt, engine = "xgboost") #add pro includes expandx
  full_dt[, cpr := cpr]

  if (verbose) cat("imputing...\n")
  pred <- parallel::mclapply(seq_len(nboot), \(b) {
    full_dt |> single_impute_xgb(...)
  }, mc.cores = ecores)
  # Combine results
  y0_post_mat <- do.call(cbind, pred)

  # Report timing if verbose
  if (verbose) {
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "secs")
    cat("impute_xgb completed in", round(elapsed, 2), "seconds\n")
  }
  y0_post_mat
}

#dt |> multiple_impute_xgb(nboot = 1000, ecores = 16)


#' Impute untreated outcomes using a forest-based multiple imputation procedure
#'
#' Clears previously created posterior draw columns (\code{y0_post_*}) and
#' summary columns (\code{y0_hat}, \code{y0_sd}, \code{y0_pi}, \code{y0_lower}, \code{y0_upper})
#' to avoid stale values, then refits and imputes using \code{multiple_impute_forest()}.
#'
#' @param .dt A data.table modified by reference.
#' @param ... Arguments forwarded to \code{multiple_impute_forest()}.
#' @param verbose If TRUE, prints status messages.
#'
#' @return The input \code{.dt}, invisibly modified by reference to include
#'   posterior draw columns and summary columns.
impute_xgb <- function(.dt, ..., verbose = TRUE) {
  #' Clears any prior y0_post_* and y0_* summary columns
  #' o avoid stale values, then refits and imputes.
  drop_by_pattern(.dt, "^y0_post_")
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")

  y0_post <-  multiple_impute_xgb(.dt, ..., verbose = verbose)
  inference_xgb(.dt, y0_post)

  .dt
}

inference_xgb <- function(.dt, y0_post) {
  nboot <- ncol(y0_post)
  # point prediction on original rows
  y0_hat <- apply(y0_post, 1, mean)
  .dt[, y0_hat := y0_hat]
  .dt[, (paste0("y0_post_", seq_len(nboot))) := data.frame(y0_post)]
  .dt
}

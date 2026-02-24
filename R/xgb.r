library(xgboost)

pre_xgb <- function(
  .dt
  , shuffle_sample = TRUE
  , wtype =  "exp"
  , shuffle_time = FALSE
  , predictors = NULL
  , cat_preds = NULL
  , year_col = "year"
  , start_year_col = "start_year"
  , at_risk_col = "n"
  , id_col = "i"
  , shift_type = NULL
  , shift = 1
) {

  predictors <- unique(c(predictors, cat_preds))

  xgb_dt <- copy(.dt)
  # Step 1: simulate new assigment
  A_sim <- xgb_dt |> simulate_sticky_assignments(
    shift_type = shift_type
    , shift = shift
  )
  xgb_mask <- mask_pre(xgb_dt, A_sim)

  # Step 2: build assigment dependent features
  art <- xgb_dt |>
    make_y0_context(A_sim = A_sim, cat_preds = cat_preds) |>
    materialize()
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
  , xgb_verbose = 0
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
    verbose = xgb_verbose
  )

}

#a <- dt |> pre_xgb(shift_type = "shift_only", shift = 2)
#m <- a$data |> fit_xgb(a$mask)
#p <- predict_xgb(m, a$data)
#p


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
  dots <- list(...)
  # split args by what each function can take
  pre_args <- dots[names(dots) %in% names(formals(pre_xgb))]
  fit_args <- dots[names(dots) %in% names(formals(fit_xgb))]

  # Step 1: Fit
  if (verbose) cat("Building boosted ensemble\n")
  xgb_out <- do.call(
    pre_xgb, c(list(.dt), pre_args)
  )

  xgb_fit <- do.call(
    fit_xgb, c(list(xgb_out$data, is_train = xgb_out$mask), fit_args)
  )

  # Step 2: Predict
  predict_xgb(xgb_fit, xgb_out$data)
}
#dt |> single_impute_xgb(shift_type = "shift_only", shift = 2)

#dt |> add_y0_history() |> single_impute_xgb()
#colnames(dt)


multiple_impute_xgb <- function(
  .dt
  , ...
  , verbose = TRUE
  , nboot = 20
  , ncores = 1
  , cat_preds = NULL
  , predictors = NULL
) {

  predictors <- unique(c(predictors, cat_preds))

  # Start timer for performance tracking
  start_time <- Sys.time()
  if (verbose) cat("Starting impute_xgb with", nboot, "bootstraps\n")

  full_dt <- copy(.dt)
  if (verbose) cat("Adding history...\n")
  full_dt |> add_y0_history()

  # Step 3: compute features that depend on assignment
  if (verbose) cat("Adding context...\n")
  art <- full_dt |> make_y0_context(cat_preds = cat_preds) |> materialize()
  full_dt[, names(art) := art]

  if (verbose) cat("Adding hazard...\n")
  cpr <- full_dt |>
    get_hazard(engine = "xgboost", predictors = predictors)
  full_dt[, cpr := cpr]
  #print(paste(
  #  nrow(full_dt)
  #  , length(unique(paste0(full_dt$i, full_dt$year)))
  #  , mean(full_dt$cpr)
  #))

  if (verbose) cat("imputing...\n")
  pred <- parallel::mclapply(seq_len(nboot), \(b) {
    full_dt |> single_impute_xgb(
      ..., cat_preds = cat_preds, predictors = predictors
    )
  }, mc.cores = ncores)
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

#dt |> multiple_impute_xgb(nboot = 2, ecores = 1, shift_type = "shift_only", shift = 2)




impute_xgb <- function(.dt, ..., verbose = TRUE) {
  #' Clears any prior y0_post_* and y0_* summary columns
  #' o avoid stale values, then refits and imputes.
  drop_by_pattern(.dt, "^y0_post_")
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")

  y0_post <-  multiple_impute_xgb(.dt, ..., verbose = verbose)
  inference_xgb(.dt, y0_post)

  .dt
}

# merely a summary
inference_xgb <- function(.dt, y0_post) {
  nboot <- ncol(y0_post)
  # point prediction on original rows
  y0_hat <- apply(y0_post, 1, mean)
  .dt[, y0_hat := y0_hat]
  .dt[, (paste0("y0_post_", seq_len(nboot))) := data.frame(y0_post)]
  .dt
}

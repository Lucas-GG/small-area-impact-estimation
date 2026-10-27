library(xgboost)

pre_xgb <- function(
  .dt
  , shuffle_sample = TRUE
  , wtype =  "exp"
  , shuffle_time = FALSE
  , predictors = NULL
  , cat_preds = NULL
  , year_col = "year"
  , cohort_col = "cohort"
  , at_risk_col = "n"
  , id_col = "i"
  , shift_type = NULL
  , shift = 1
  , baseline = c("etwfe", "none")
) {
  baseline <- match.arg(baseline)
  predictors <- unique(c(predictors, cat_preds))

  xgb_dt <- copy(.dt)
  # Step 1: perturb adoption backward from the first masked year
  A_sim <- xgb_dt |> simulate_sticky_assignments(
    shift_type = shift_type
    , shift = shift
  )
  anchor <- pmin(xgb_dt$mask_year, unname(A_sim[as.character(xgb_dt[[id_col]])]))
  xgb_mask <- !is.na(xgb_dt$y0) & xgb_dt[[year_col]] < anchor

  # Step 2: bootstrap weights (requieres id_col); case weights only, never a feature
  wts <- rep(1, nrow(xgb_dt))
  if (shuffle_sample) {
    wts <- xgb_dt |> cluster_boot(wtype = wtype, shuffle_time = shuffle_time)
  }

  select <- function(d) xgb_select(d, predictors = predictors, year_col = year_col
    , cohort_col = cohort_col, at_risk_col = at_risk_col
  )
  if (baseline == "none") {
    # Step 3: features frozen at the member's adoption year, so validation rows
    # (observed y0 after the perturbed year) are built like prediction rows
    xgb_dt |> build_y0_features(anchor = anchor, cohort_col = cohort_col
      , cat_preds = cat_preds
    )
    return(list(data = select(xgb_dt), mask = xgb_mask, wts = wts
      , offset = log(xgb_dt[[at_risk_col]])))
  }

  # Step 3 (baseline = "etwfe"): the member's Poisson ETWFE on its training rows
  # gives the offset; history features are observed/expected ratios against it,
  # frozen at the member's anchor (fit) and at mask_year (prediction).
  lp <- etwfe_baseline(xgb_dt, xgb_mask, wts, predictors = predictors
    , cohort_col = cohort_col, year_col = year_col, at_risk_col = at_risk_col
  )
  merged <- attr(lp, "merged")
  lp <- as.vector(lp)
  expected <- xgb_dt[[at_risk_col]] * exp(lp)
  features <- function(anc) {
    f <- copy(xgb_dt)
    f[, y0 := y0 / expected]
    f |> build_y0_features(anchor = anc, cohort_col = cohort_col, cat_preds = cat_preds)
    f[, ybar_ctrl := NULL]  # composition-driven; the trend is in the offset
    f[, y0 := xgb_dt$y0]
    select(f)
  }
  list(data = features(anchor), mask = xgb_mask, wts = wts
    , offset = log(xgb_dt[[at_risk_col]]) + lp
    , pred_data = features(xgb_dt$mask_year)
    , merged = merged
  )
}

# Linear predictor (without log n) of a Poisson ETWFE fitted on the training rows:
# cohort and year effects plus common coefficients on `predictors`, with log(n)
# as exposure and the bootstrap weights. A year without training rows takes the
# nearest year's effect; a cohort without training rows or events, the mean
# cohort effect.
etwfe_baseline <- function(dt, is_train, wts, predictors = NULL, cohort_col = "cohort"
  , year_col = "year", at_risk_col = "n"
) {
  # Cohort x group cells from the member's training rows
  cell <- if ("ncf" %in% names(dt)) cohort_cell(dt[[cohort_col]], dt$ncf, dt$y0, train = is_train)
    else structure(as.character(dt[[cohort_col]]), merged = character())
  tr <- data.table::data.table(y0 = dt$y0, n = dt[[at_risk_col]], year = dt[[year_col]]
    , coh = as.vector(cell), w = wts)
  if (length(predictors)) tr[, (predictors) := dt[, ..predictors]]
  tr <- tr[is_train]
  rhs <- if (length(predictors)) paste(predictors, collapse = " + ") else "1"
  fit <- fixest::fepois(stats::as.formula(paste("y0 ~", rhs, "| coh + year")), data = tr
    , weights = ~w, offset = ~log(n), notes = FALSE, warn = FALSE)
  fe <- fixest::fixef(fit)
  years <- sort(unique(dt[[year_col]]))
  gam <- unname(fe$year[as.character(years)])
  gam <- data.table::nafill(data.table::nafill(gam, "locf"), "nocb")
  alp <- unname(fe$coh[as.vector(cell)])
  alp[is.na(alp)] <- mean(fe$coh)
  xb <- if (length(predictors)) as.vector(as.matrix(dt[, ..predictors]) %*% stats::coef(fit)[predictors]) else 0
  structure(alp + gam[match(dt[[year_col]], years)] + xb, merged = attr(cell, "merged"))
}

# Booster inputs: year, n, adoption cohort (factor), y0 features and user
# predictors, followed by the label y0. Fails if a forbidden column slips in.
xgb_select <- function(dt
  , predictors = NULL
  , year_col = "year"
  , cohort_col = "cohort"
  , at_risk_col = "n"
  , verbose = FALSE
) {
  feats <- unique(c(year_col, at_risk_col, cohort_col, get_y0_features(dt), predictors))
  bad <- intersect(feats, c("y0", "wts", "cpr", "mask_year", "i", hidden_cols))
  if (length(bad)) stop("Not allowed as GB-Ens predictors: ", paste(bad, collapse = ", "))
  if (verbose) cat("predictors are", feats, "\n")
  out <- dt[, c(feats, "y0"), with = FALSE]
  out[, (cohort_col) := factor(get(cohort_col))]
  out
}


fit_xgb <- function(
  xgb_dt
  , is_train
  , wts = rep(1, nrow(xgb_dt))
  , num_threads = 1
  , params = NULL
  , nrounds = 1000
  , early_stop = 20
  , predictors = NULL
  , ptry = 1 # no column subsampling
  , xgb_verbose = 0
  , offset = log(xgb_dt$n)
) {
  # Step 5: Fit XGBoost Poisson model; the offset is log(n), plus the ETWFE
  # linear predictor with baseline = "etwfe"; n is also a feature
  is_valid <- !is_train & !is.na(xgb_dt$y0)
  x <- xgb_dt[, !("y0")]

  dtrain <- xgboost::xgb.DMatrix(x[is_train], label = xgb_dt$y0[is_train]
    , base_margin = offset[is_train], weight = wts[is_train]
  )

  dvalid <- xgboost::xgb.DMatrix(x[is_valid], label = xgb_dt$y0[is_valid]
    , base_margin = offset[is_valid], weight = wts[is_valid]
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


predict_xgb  <- function(fit, xgb_dt, offset = log(xgb_dt$n)) {
  x <- xgboost::xgb.DMatrix(xgb_dt[, !("y0")], base_margin = offset)
  predict(fit, x, iteration_range = c(1, best_round(fit)))
}

# Early-stopping round (1-based). xgboost >= 2 stores it 0-based as an attribute.
best_round <- function(fit) {
  best <- xgboost::xgb.attr(fit, "best_iteration")
  if (is.null(best)) return(xgboost::xgb.get.num.boosted.rounds(fit))
  as.integer(best) + 1L
}

single_impute_xgb <- function(.dt
  , pred_data
  , ...
  , verbose = TRUE
) {
  dots <- list(...)
  # split args by what each function can take
  pre_args <- dots[names(dots) %in% names(formals(pre_xgb))]
  fit_args <- dots[names(dots) %in% names(formals(fit_xgb))]

  # Step 1: Fit on the member's features
  if (verbose) cat("Building boosted ensemble\n")
  xgb_out <- do.call(
    pre_xgb, c(list(.dt), pre_args)
  )
  # With an ETWFE baseline the prediction features depend on the member's fit.
  if (!is.null(xgb_out$pred_data)) pred_data <- xgb_out$pred_data
  stopifnot(identical(names(xgb_out$data), names(pred_data)))

  xgb_fit <- do.call(
    fit_xgb, c(list(xgb_out$data, is_train = xgb_out$mask, wts = xgb_out$wts
      , offset = xgb_out$offset), fit_args)
  )

  # Step 2: Predict on the outer features, frozen at the first masked year
  list(
    pred = predict_xgb(xgb_fit, pred_data, offset = xgb_out$offset)
    , importance = xgboost::xgb.importance(model = xgb_fit)[, .(Feature, Gain)]
    , best_round = best_round(xgb_fit)
    , features = setdiff(names(pred_data), "y0")
    , cell_merges = xgb_out$merged
  )
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
  , cohort_col = "cohort"
) {

  predictors <- unique(c(predictors, cat_preds))

  # Start timer for performance tracking
  start_time <- Sys.time()
  if (verbose) cat("Starting impute_xgb with", nboot, "bootstraps\n")

  # Step 3: features at the observed adoption (first masked year)
  full_dt <- copy(.dt)
  if (verbose) cat("Adding history and context...\n")
  full_dt |> build_y0_features(cohort_col = cohort_col, cat_preds = cat_preds)

  if (verbose) cat("Adding hazard...\n")
  cpr <- full_dt |>
    get_hazard(engine = "xgboost", predictors = predictors, cohort_col = cohort_col)
  full_dt[, cpr := cpr]
  pred_data <- full_dt |>
    xgb_select(predictors = predictors, cohort_col = cohort_col, verbose = verbose)

  if (verbose) cat("imputing...\n")
  members <- parallel::mclapply(seq_len(nboot), \(b) {
    full_dt |> single_impute_xgb(pred_data
      , ..., cat_preds = cat_preds, predictors = predictors, cohort_col = cohort_col
      , verbose = verbose
    )
  }, mc.cores = ncores)
  failed <- vapply(members, inherits, NA, "try-error")
  if (any(failed)) stop(members[[which(failed)[1]]])

  # Combine results
  y0_post_mat <- do.call(cbind, lapply(members, `[[`, "pred"))
  dots <- list(...)
  attr(y0_post_mat, "provenance") <- list(
    features = members[[1]]$features
    , baseline = if (is.null(dots$baseline)) "etwfe" else dots$baseline
    , importance = rbindlist(lapply(seq_along(members), \(b)
      members[[b]]$importance[, member := b]))
    , best_round = vapply(members, `[[`, 1L, "best_round")
    , cell_merges = lapply(members, `[[`, "cell_merges")
    , git = git_state()
    , session = sessionInfo()
  )

  # Report timing if verbose
  if (verbose) {
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "secs")
    cat("impute_xgb completed in", round(elapsed, 2), "seconds\n")
  }
  y0_post_mat
}

# Commit and dirty flag of the working tree, for provenance.
git_state <- function() {
  run <- function(...) tryCatch(
    suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE)),
    error = function(e) character())
  commit <- run("rev-parse", "HEAD")
  list(
    commit = if (length(commit)) commit else NA_character_
    , dirty = length(run("status", "--porcelain", "--untracked-files=no")) > 0
  )
}

#dt |> multiple_impute_xgb(nboot = 2, ecores = 1, shift_type = "shift_only", shift = 2)




impute_xgb <- function(.dt, ..., verbose = TRUE) {
  #' Clears any prior y0_post_* and y0_* summary columns
  #' o avoid stale values, then refits and imputes.
  drop_by_pattern(.dt, "^y0_post_")
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")

  y0_post <-  multiple_impute_xgb(.dt, ..., verbose = verbose)
  inference_xgb(.dt, y0_post)
  data.table::setattr(.dt, "xgb_provenance", attr(y0_post, "provenance"))

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

# Reset outcomes at and after treatment start.
# Modifies dt by reference:
# - event_time = year - start_year (NA for never-treated if start_year is Inf)
# - sets y0 and y0_25 to NA for year >= start_year (finite start_year only)
reset_y0 <- function(dt) {

  dt_y0 <- data.table(i, year, y0_sim)

  # Preserve observed start year once
  if (!"start_year_obs" %in% names(dt)) dt[, start_year_obs := start_year]

  # Use simulated start year if available
  if ("start_year_sim" %in% names(dt)) dt[, start_year := start_year_sim]
  if ("W_sim" %in% names(dt)) dt[, W := W_sim]

  if ("event_time" %in% names(dt))
    dt[, event_time := fifelse(is.finite(start_year), year - start_year, NA_real_)]

  if ("y0" %in% names(dt))
    dt[is.finite(start_year) & year >= start_year, y0 := NA_real_]

  if ("y0_25" %in% names(dt))
    dt[is.finite(start_year) & year >= start_year, y0_25 := NA_real_]

  dt
}
#Scheme A (my preference)
#fit_*() = returns model object (no mutation)
#inference_*() = given fitted model, attach draws/summaries (mutates dt)
#impute_*() = convenience wrapper: clear old columns, fit, inference (mutates dt)

library(rpart)

#' impute_forest()
#' Fit a small forest of Poisson CARTs on pre-treatment data
#' (after feature expansion) and add counterfactual predictions y0_hat
#' (and optionally y0_post_* draws) to a copy of .dt.
#' Expects .dt to contain at least i, year, start_year, y, n
#' (and whatever add_prop/expand_y0 use).
#' Returns: data.table copy with y0_hat 
#' (and optionally y0_post_1:K). Does NOT modify .dt.
impute_forest <- function(.dt, ntrees = 20
                          , k = 10, min_bucket = 7, max_depth = 30
                          , mtry = NULL, ncores = 1, exvar = TRUE
                          , wtype = "exp", shuffle_time = FALSE
                          , shuffle_lag = 2
                          , interval = FALSE
                          , verbose = TRUE) {

  #' Clears any prior y0_post_* and y0_* summary columns
  #' o avoid stale values, then refits and imputes.
  drop_by_pattern(.dt, "^y0_post_")
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")


  #result <- copy(.dt)
  # Start timer for performance tracking
  start_time <- Sys.time()
  if (verbose) cat("Starting impute_forest_dt with", ntrees, "trees\n")



  # Step 2: Build forest
  if (verbose) cat("Building forest with", ntrees, "trees...\n")
  .forest <- fit_forest(.dt
                    , ntrees = ntrees
                    , k = k
                    , min_bucket = min_bucket
                    , max_depth = max_depth
                    , mtry = mtry
                    , ncores = ncores
                    , exvar = exvar
                    , wtype = wtype
                    , shuffle_time = shuffle_time
                    , shuffle_lag = shuffle_lag)


  inference_forest(.dt, .forest, ncores = ncores)

  # Report timing if verbose
  if (verbose) {
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "secs")
    cat("impute_forest_dt completed in", round(elapsed, 2), "seconds\n")
  }

  .dt
}

inference_forest <- function(.dt, ntrees = 20,
                             k = 10, min_bucket = 7, max_depth = 30,
                             mtry = NULL, nboot = 20, wtype = "exp",
                             ncores = 1, shuffle_time = FALSE, exvar = TRUE,
                             shuffle_lag = 2,
                             verbose = TRUE) {

  # point prediction on original rows
  y0_hat <- predict(m, newdata = .dt)
  .dt[, y0_hat := y0_hat]

  # bootstrap draws: refit on bootstrap sample, predict on original rows
  y0_post <- parallel::mclapply(seq_len(nboot), \(b) {
    cluster_boot(.dt
      , wtype = wtype
      , shuffle_time = shuffle_time
    )
    .forest <- fit_forest(.dt
                    , ntrees = ntrees
                    , k = k
                    , min_bucket = min_bucket
                    , max_depth = max_depth
                    , mtry = mtry
                    , ncores = ncores
                    , exvar = exvar
                    , wtype = wtype
                    , shuffle_time = shuffle_time
                    , shuffle_lag = shuffle_lag)
    predict(.forest, newdata = .dt)
  }, mc.cores = ncores)

  .dt[, (paste0("y0_post_", seq_len(nboot))) := y0_post]
  .dt
}

#' fit_single_tree()
#' Fit one Poisson regression tree (rpart) to model y0 given features,
#' with optional bootstrap weights.
#' Expects columns: y0, n, and predictors;
#' uses wts if present, otherwise sets wts := 1.
#' Returns: fitted rpart object.
#' Modifies .dt by reference only if it needs to create wts.
fit_single_tree <- \(.dt, k, min_bucket, max_depth) {
  #' expecting BB weights
  if (!"wts" %in% names(.dt)) .dt[, wts := 1]

  m0 <- rpart(cbind(n / 10^5, y) ~ .
    , data = .dt
    , method = "poisson"
    , parms = list(shrink = k) # the smaller, the most shrinkage
    , weights = wts
    , control = rpart.control(
      cp = -Inf
      , minbucket = min_bucket
      , maxdepth = max_depth
    )
  )
  m0
}

#' rselect()
#' Randomly select mtry predictor columns
#' while always retaining response columns (yvar).
#' Intended for per-tree feature subsampling;
#' does not copy the underlying data beyond the subset.
#' Returns: data.table with kept columns only; does not modify .dt.
rselect <- function(.dt, mtry = 0, yvar = NULL) {
  # Set default mtry if it's 0
  if (mtry == 0) mtry <- floor(sqrt(ncol(.dt)))

  # Set default yvar
  if (is.null(yvar)) yvar <- c("n", "y0", "wts", "y0_hat")

  # Find which columns to keep as response variables (y)
  y_cols <- intersect(names(.dt), yvar)

  # Find which columns to consider for predictors (x)
  x_cols <- setdiff(names(.dt), y_cols)

  # If mtry is greater than available columns, use all columns
  mtry <- min(mtry, length(x_cols))

  # Sample mtry columns from x
  if (mtry < length(x_cols)) {
    sampled_x_cols <- sample(x_cols, mtry)
  } else {
    sampled_x_cols <- x_cols
  }

  # Combine columns to keep
  keep_cols <- c(y_cols, sampled_x_cols)

  # Select only the columns we need (data.table way)
  return(.dt[, ..keep_cols])
}

#' forest()
#' Build a list of Poisson CARTs (rpart) using
#' simulated treatment assignments and bootstrap weights.
#' Pipeline per tree: simulate_assignments() -> expand_y0()
#' -> keep pre-treatment rows -> (optional) cluster_boot()
#' -> drop response columns y/y25 -> (optional) feature subsample 
#' -> fit_single_tree().
#' Expects .dt to contain at least i, year, start_year, y, n
#' (and hazard cols used by add_prop_hazard()).
#' Returns: list of rpart models.
fit_forest <- \(.dt, ntrees = 20L
  , k = 10,  min_bucket = 7, max_depth = 30
  , mtry = NULL
  , ncores = 1
  , wtype = "exp"
  , shuffle_assignment = TRUE
  , shuffle_sample = TRUE
  , shuffle_time = FALSE
  , shuffle_lag = 2
  , verbose = TRUE
) {

  # Step 1: Expand features and add propensity scores
  if (verbose) cat("Expanding features and calculating propensity scores...\n")
  cpr <- get_hazard(.dt, engine = "xgboost") #add pro includes expandx
  .dt[, cpr := cpr]


  list_of_trees <- lapply(seq_len(ntrees), \(i) {
    if (i %% 10 == 0) message("Processing tree ", i, " of ", ntrees)
    r <- sample(1:10000, 1)
    set.seed(i + r)
    # Step 2: simulate new assigment
    if (shuffle_assignment) {
      A_sim <- .dt |> simulate_assignments()
      mask <- mask_pre(.dt, A_sim)
    }

    # Step 3: compute features that depend on assignment
    art <- .dt |> expand_y0(A_sim) |> materialize()
    tree_dt <- .dt[, names(art) := art][mask]

    # Step 4: bootstrap weights
    if (shuffle_sample) {
      wts <- tree_dt |> cluster_boot(wtype = wtype, shuffle_time = shuffle_time)
      tree_dt[, wts := wts]
    }

    # Step 4: Select variables randomly (bagging)
    keep_cols <- colnames(tree_dt) |> sample_feature()

    # Step 5: Fit tree
    tree_dt[, ..keep_cols] |>
      fit_single_tree(k, min_bucket, max_depth)

  })

  list_of_trees
}

dt |> fit_forest(ntrees = 20)

#' predict_forest()
#' Predict with each tree in the forest and average predictions across trees.
#' Expects .dt to have column n;
#' uses wts only for compatibility (not used in prediction aggregation).
#' Returns: numeric vector of predicted counts
#' on the original scale (pi_hat * n / 1e5).
predict_forest <- function(forest, .dt, ncores = 1) {
  # Add wts column if not present
  if (!"wts" %in% names(.dt)) .dt[, wts := 1]

  # Use mclapply for parallel prediction
  tree_predictions <- mclapply(forest, function(m) {
    predictions <- predict(m, .dt)
    predictions
  }, mc.cores = ncores)

  # Convert list of predictions to matrix for efficiency
  predictions_matrix <- do.call(cbind, tree_predictions)

  # Calculate row means (average predictions across trees)
  pi_hat <- rowMeans(predictions_matrix)

  # Apply scaling formula - note this returns a vector, not a data.table
  pi_hat * .dt$n / 10^5
}

if (FALSE) {
  # Example usage
  dt <- data.table(
    i = rep(1:5, each = 4),
    year = rep(2001:2004, times = 5),
    y = rpois(20, lambda = 5),
    n = rep(10, 20)
  )
  dt[, wts := 1]

  impute_forest(dt, ntrees = 10, ncores = 2, verbose = TRUE)

  print(dt)
  .forest <- fit_forest(dt, ntrees = 10, ncores = 2, verbose = TRUE)

}




sample_feature <- function(cols
  , mtry = floor(sqrt(length(cols)))
  , always = c("n", "y", "wts")
) {
  y_keep <- intersect(cols, always)
  x_cols <- setdiff(cols, y_keep)
  mtry <- min(mtry, length(x_cols))
  x_keep <- if (mtry < length(x_cols)) sample(x_cols, mtry) else x_cols
  c(y_keep, x_keep)
}

#sample_feature(c("a", "b", "c", "d", "e"), 3)

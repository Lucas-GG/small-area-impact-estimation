
library(rpart)

impute_forest <- function(.dt, ntrees = 20,
                          k = 10, min_bucket = 7, max_depth = 30,
                          mtry = NULL, ncores = 1, exvar = FALSE,
                          wtype = "exp", shuffle_time = FALSE,
                          verbose = TRUE) {

  .dt <- copy(.dt)
  # Start timer for performance tracking
  start_time <- Sys.time()
  if (verbose) cat("Starting impute_forest_dt with", ntrees, "trees\n")

  # Step 1: Expand features and add propensity scores
  if (verbose) cat("Expanding features and calculating propensity scores...\n")
  edata <- add_prop(.dt)

  # Step 2: Build forest
  if (verbose) cat("Building forest with", ntrees, "trees...\n")
  forest <- forest(edata,
                   ntrees = ntrees,
                   k = k,
                   min_bucket = min_bucket,
                   max_depth = max_depth,
                   mtry = mtry,
                   ncores = ncores,
                   exvar = exvar,
                   wtype = wtype,
                   shuffle_time = shuffle_time)

  # Step 3: Make predictions
  if (verbose) cat("Generating predictions...\n")
  predictions <- predict_forest(forest, edata, ncores)

  # Step 4: Add predictions to original data
  .dt[, y0_hat := predictions]

  # Report timing if verbose
  if (verbose) {
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "secs")
    cat("impute_forest_dt completed in", round(elapsed, 2), "seconds\n")
  }

  .dt
}

fit_single_tree <- \(.dt, k, min_bucket, max_depth) {
  #' expecting BB weights
  if (!"wts" %in% names(.dt)) .dt[, wts := 1]

  m0 <- rpart(cbind(n / 10^5, y0) ~ .
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

# I think this is not making unnecesary copies
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

forest <- \(.dt, ntrees = 20L
  , k = 10,  min_bucket = 7, max_depth = 30
  , mtry = NULL
  , ncores = 1
  , exvar = FALSE
  , wtype = "exp"
  , shuffle_time = FALSE
) {

  list_of_trees <- mclapply(seq_len(ntrees), function(x) {
    if (x %% 10 == 0) message("Processing tree ", x, " of ", ntrees)

    # Copy the data.table to avoid modifying the original data
    dt_tree <- copy(.dt)

    # Step 1: Shuffle start years
    dt_tree <- shuffle_start(dt_tree)

    # Step 2: Expand features
    dt_tree <- expandx(dt_tree)

    # Step 3: Filter where year < start_year (data.table way)
    dt_tree <- dt_tree[year < start_year]

    # Step 4: Apply bootstrap weights
    if (exvar) {
      dt_tree <- cluster_boot(dt_tree
        , wtype = wtype, shuffle_time = shuffle_time
      )
    }


    # Step 5: Remove response columns (much faster in data.table)
    keep_cols <- setdiff(names(dt_tree), c("y", "y25"))
    dt_tree <- dt_tree[, ..keep_cols]

    # Step 6: Select variables randomly
    if (!is.null(mtry)) {
      dt_tree <- rselect(dt_tree, mtry)
    }

    # Step 7: Fit tree
    fit_single_tree(dt_tree, k, min_bucket, max_depth)

  }, mc.cores = ncores)

  return(list_of_trees)
}

predict_forest_old <- \(forest, .data, ncores) {

  if (!"wts" %in% names(.data))  .data[, wts := 1]
  #.newdata <- .data %>%  expandx #this is redundant
  # predict tree return the estimated rate (pi)
  pi_hat <- mclapply(forest, \(m) predict(m, .data), mc.cores = ncores) %>%
    do.call("cbind", .) %>%
    rowMeans()
  pi_hat * .data$n / 10^5
  # "vector" for Poisson trees it is the estimated response rate
}

predict_forest <- function(forest, .dt, ncores = 1) {
  # Add wts column if not present
  if (!"wts" %in% names(dt)) .dt[, wts := 1]

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
  pi_hat * dt$n / 10^5
}

#cluster_boot(dt1, "mammen") %>% head
#' generate fractional weights
#' which are similar to resampling units (ids)
shuffle_cl <- function(dt, cl = "i", wtype = "mammen") {
  # Make a copy to avoid modifying the original
  result <- copy(dt)

  # Add wts column if it doesn't exist
  if (!"wts" %in% names(result)) {
    result[, wts := 1]
  }

  # Get weight generator from fwb package
  gen_w <- fwb:::make_gen_weights(wtype)

  # Get unique cluster values
  unique_clusters <- unique(result[[cl]])
  n_cl <- length(unique_clusters)

  # Generate new weights
  new_weights <- as.vector(gen_w(n_cl, 1))

  # Create mapping table for fast joining
  wdata <- data.table(
    cluster_val = unique_clusters,
    new_wts = new_weights
  )
  setnames(wdata, "cluster_val", cl)

  # Join weights to the main table (using fast data.table join)
  setkeyv(wdata, cl)
  setkeyv(result, cl)
  result <- result[wdata]

  # Update weights
  k <- nrow(result)
  result[, wts := wts * new_wts]
  result[, new_wts := NULL]

  # Normalize weights
  total_weight <- sum(result$wts)
  result[, wts := wts / total_weight * k]

  # Return the result (will be visible when explicitly printed)
  result
}

cluster_boot <- function(dt
  , cl = "i", wtype = "mammen", shuffle_time = TRUE
) {
  # Make a copy of the input data table
  result <- copy(dt)

  # Apply time shuffling if requested
  if (shuffle_time) {
    result <- shuffle_cl(result, cl = "year", wtype = wtype)
  }

  # Apply cluster shuffling
  result <- shuffle_cl(result, cl = cl, wtype = wtype)

  # Return the result
  return(result)
}

inference_forest <- function(.dt, ntrees = 20,
                                k = 10, min_bucket = 7, max_depth = 30,
                                mtry = NULL, nboot = 20, wtype = "exp",
                                ncores = 1, shuffle_time = FALSE, exvar = FALSE,
                                verbose = TRUE) {

  # Start timer
  t0 <- Sys.time()
  if (verbose) cat("Starting inference with", nboot, "bootstrap iterations\n")

  # Ensure input is a data.table and create a copy to avoid modifying original
  if (!is.data.table(.dt)) {
    result_dt <- as.data.table(.dt)
  } else {
    result_dt <- copy(.dt)
  }

  # Optimize core allocation
  if (ncores > 1) {
    # Use multiple cores for bootstrap iterations, single core for tree building
    ncores_boot <- ncores
    ncores_tree <- 1
  } else {
    # Use single core for bootstrap, multiple cores for tree building
    ncores_boot <- 1
    ncores_tree <- min(20, parallel::detectCores() - 1)
  }

  if (verbose) {
    cat("Using", ncores_boot, "cores for bootstrap iterations and", 
        ncores_tree, "cores for tree building\n")
  }

  # Process bootstrap iterations in parallel
  bootstrap_results <- mclapply(seq_len(nboot), function(i) {
    if (verbose && (i %% 5 == 0 || i == 1 || i == nboot)) {
      cat("Processing bootstrap iteration", i, "of", nboot, "\n")
    }

    # Generate bootstrap sample
    boot_dt <- cluster_boot(result_dt, wtype = wtype
      , shuffle_time = shuffle_time
    )

    # Run forest imputation
    imputed_dt <- impute_forest(boot_dt, ntrees = ntrees,
                                k = k, min_bucket = min_bucket,
                                max_depth = max_depth, mtry = mtry,
                                ncores = ncores_tree, exvar = exvar,
                                verbose = FALSE)

    # Return just the predictions
    imputed_dt$y0_hat
  }, mc.cores = ncores_boot)

  # Combine results into a matrix
  post_matrix <- do.call("cbind", bootstrap_results)

  # Add posterior distribution to result
  # Fast implementation for storing bootstrap results as individual columns
  for (i in 1:nboot) {
    set(result_dt, j = paste0("y0_post_", i), value = post_matrix[, i])
  }

  # Report timing
  t1 <- Sys.time()
  elapsed <- difftime(t1, t0, units = "secs")
  if (verbose) {
    cat("Inference completed in", round(elapsed, 2), "seconds\n")
  } else {
    print(elapsed)
  }

  result_dt
}
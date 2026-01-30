library(rpart)
#' Fit a single Poisson regression tree
#'
#' Fits one Poisson CART model using \code{rpart} to predict \code{y0}
#' (an untreated / admissible outcome proxy) from the available predictors.
#' The model uses \code{method="poisson"} and \code{parms=list(shrink=k)}.
#'
#' @param .dt A data.table containing the response and predictors.
#'   Must include columns \code{y0} and \code{n}. Any additional columns are
#'   treated as predictors.
#' @param k Nonnegative shrinkage parameter passed to \code{rpart} via
#'   \code{parms=list(shrink=k)}. Smaller values imply stronger shrinkage.
#' @param min_bucket Minimum number of observations in a terminal node
#'   (passed to \code{rpart.control(minbucket=...)}).
#' @param max_depth Maximum depth of the fitted tree
#'   (passed to \code{rpart.control(maxdepth=...)}).
#'
#' @details
#' The model formula is \code{cbind(n/1e5, y0) ~ .}. This uses \code{n/1e5}
#' as the exposure component for \code{rpart}'s Poisson method, and returns
#' predictions on the \emph{rate} scale by default (see \code{predict.rpart}).
#'
#' If column \code{wts} is not present, it is created with value 1 for all rows.
#' This is intended to support Bayesian bootstrap / cluster bootstrap weights.
#'
#' @return A fitted \code{rpart} object.
#' @seealso \code{\link[rpart]{rpart}}, \code{\link[rpart]{predict.rpart}}
#' @keywords internal
fit_single_tree <- \(.dt
  , k = 10,  min_bucket = 7, max_depth = 30
  #, k = 10,  min_bucket = 7, max_depth = 30
) {
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
#dt |> fit_single_tree(k = 10, min_bucket = 7, max_depth = 30)

#' Fit a Poisson regression-tree forest under simulated assignments
#'
#' Builds a list of Poisson CART models (one per tree), where each tree is trained
#' on a bootstrapped / resampled version of the data and (optionally) under a
#' simulated treatment assignment. Each tree has its own:
#' \itemize{
#'   \item simulated assignment \code{A_sim} (if \code{shuffle_assignment=TRUE})
#'   \item assignment-dependent features via \code{expand_y0(A_sim = ...)}
#'   \item pre-treatment mask via \code{mask_pre(.dt, A_sim)}
#'   \item cluster bootstrap weights via \code{cluster_boot()} (if \code{shuffle_sample=TRUE})
#'   \item random feature subsample via \code{sample_feature()}
#' }
#'
#' @param .dt A data.table with at least \code{id}, \code{time}, \code{start_year},
#'   the observed outcome \code{y}, and exposure \code{n}, plus any additional
#'   baseline predictors. If assignment simulation depends on hazard-related
#'   variables, those must also be present.
#' @param ntrees Number of trees to fit.
#' @param k,min_bucket,max_depth Tree hyperparameters passed to \code{fit_single_tree()}.
#' @param mtry Number of predictors to sample per tree when bagging. If \code{NULL},
#'   \code{sample_feature()} uses \code{floor(sqrt(p))} where \code{p} is the number
#'   of candidate predictors.
#' @param ncores Number of cores (currently used only for downstream predict; the
#'   tree fitting loop uses \code{lapply} in the current implementation).
#' @param wtype Weight type passed to \code{cluster_boot()} (e.g., \code{"exp"}).
#' @param shuffle_assignment If TRUE, simulates a new assignment \code{A_sim} per tree
#'   via \code{simulate_assignments()} and computes a corresponding pre-treatment mask.
#' @param shuffle_sample If TRUE, draws bootstrap weights per tree via \code{cluster_boot()}.
#' @param shuffle_time,shuffle_lag Parameters forwarded to \code{cluster_boot()} controlling
#'   optional within-cluster time shuffling.
#' @param verbose If TRUE, prints progress messages.
#'
#' @details
#' Per tree, the workflow is:
#' \enumerate{
#'   \item (optional) simulate \code{A_sim} and compute \code{mask_tree = mask_pre(.dt, A_sim)}
#'   \item compute assignment-dependent features: \code{expand_y0(A_sim)} and assign into \code{.dt}
#'   \item restrict to pre-treatment rows for that tree: \code{[mask_tree]}
#'   \item (optional) compute cluster bootstrap weights and store as \code{wts}
#'   \item select a random subset of predictors and fit one \code{rpart} Poisson tree
#' }
#'
#' The returned object is a list with class \code{"dforest"}.
#'
#' @return A list of fitted \code{rpart} models with class \code{"dforest"}.
#' @keywords internal
fit_forest <- \(.dt
  , predictors = NULL
  , ntrees = 20L
  , mtry = NULL
  , ncores = 1
  , wtype = "exp"
  , shuffle_assignment = TRUE
  , shuffle_sample = TRUE
  , shuffle_time = FALSE
  , shuffle_lag = 2
  , verbose = TRUE
  , year_col = "year"
  , start_year_col = "start_year"
  , at_risk_col = "n"
  , id_col = "i"
  , ...
) {

  list_of_trees <- lapply(seq_len(ntrees), \(i) {
    if (i %% 10 == 0) message("Processing tree ", i, " of ", ntrees)
    #r <- sample(1:10000, 1)
    #set.seed(i + r)

    # Step 2: simulate new assigment
    if (shuffle_assignment) {
      A_sim <- .dt |> simulate_assignments()
      mask_tree <- mask_pre(.dt, A_sim)

      # Step 3: compute features that depend on assignment
      art <- .dt |> make_y0_context(A_sim = A_sim) |> materialize()
      tree_dt <- .dt[, names(art) := art][mask_tree]
    } else {
      tree_dt <- copy(.dt)
    }



    # Step 4: bootstrap weights
    if (shuffle_sample) {
      wts <- tree_dt |> cluster_boot(wtype = wtype, shuffle_time = shuffle_time)
      tree_dt[, wts := wts]
    }

    # Step 4: Select variables randomly (bagging)
    keep_cols <- c(year_col, at_risk_col, start_year_col, "y0", "wts")
    keep_cols <- c(keep_cols, get_y0_features(.dt))  |> unique()
    keep_cols <- c(keep_cols, predictors) |> unique()
    #keep_cols <- setdiff(keep_cols, exclude)
    if (i == 1) cat("predictors are", setdiff(keep_cols, c("y0", "wts")), "\n")
    keep_cols <- keep_cols |> sample_feature(mtry = mtry)

    # Step 5: Fit tree
    tree_dt[, ..keep_cols] |>
      fit_single_tree(...)

  })

  class(list_of_trees) <- "dforest"
  list_of_trees
}

#dt |> fit_forest(ntrees = 2)

#' Single forest imputation of untreated outcomes
#'
#' Fits one forest via \code{fit_forest()} and predicts counterfactual/untreated
#' outcomes for all rows in \code{.dt} using \code{predict.dforest()}.
#'
#' @param .dt A data.table.
#' @param ntrees Number of trees.
#' @param ... Additional arguments forwarded to \code{fit_forest()}.
#' @param verbose If TRUE, prints status messages.
#'
#' @return A numeric vector of predicted counts on the original scale.
single_impute_forest <- function(.dt
  , ntrees = 20
  , ...
  , verbose = TRUE
) {

  # Step 1: Build forest
  if (verbose) cat("Building forest with", ntrees, "trees...\n")
  .forest <- .dt |> fit_forest(ntrees = ntrees, ...)

  # Step 2: Predict
  .forest |> predict(.dt, ncores = 1)

}

#dt |> single_impute_forest(ntrees = 2)

#' Multiple forest imputations via repeated forest fits
#'
#' Produces multiple imputations by repeatedly fitting a forest and predicting
#' untreated outcomes, returning an \code{nrow(.dt) x nboot} matrix of draws.
#'
#' This function first computes and attaches the hazard/propensity score
#' \code{cpr} using \code{get_hazard(engine="xgboost")}. It then computes
#' assignment-independent features via \code{expand_y0(A_sim=NULL)} once and
#' reuses them across bootstrap replications. For each bootstrap replication,
#' it fits a forest and predicts outcomes.
#'
#' @param .dt A data.table.
#' @param ntrees Number of trees per forest.
#' @param ... Arguments forwarded to \code{single_impute_forest()} / \code{fit_forest()}.
#' @param verbose If TRUE, prints status and timing.
#' @param nboot Number of forest replications (imputation draws).
#' @param forest_cores Number of cores for parallelizing over \code{nboot}
#'   using \code{parallel::mclapply()}.
#'
#' @return A numeric matrix with \code{nrow(.dt)} rows and \code{nboot} columns.
multiple_impute_forest <- function(
  .dt
  , ntrees = 20
  , ...
  , verbose = TRUE
  , nboot = 20
  , forest_cores = 1
) {

  # Start timer for performance tracking
  start_time <- Sys.time()
  if (verbose) cat("Starting impute_forest_dt with", ntrees, "trees\n")

  forest_dt <- copy(.dt)
  forest_dt |> add_y0_history()

  # Step 3: compute features that depend on assignment
  art <- forest_dt |> make_y0_context() |> materialize()
  forest_dt[, names(art) := art]

  # Step 1: Expand features and add propensity scores
  if (verbose) cat("Expanding features and calculating propensity scores...\n")
  cpr <- get_hazard(forest_dt, engine = "xgboost") #add pro includes expandx
  forest_dt[, cpr := cpr]


  pred <- parallel::mclapply(seq_len(nboot), \(b) {
    forest_dt |> single_impute_forest(ntrees = ntrees, ...)
  }, mc.cores = forest_cores)
  # Combine results
  y0_post_mat <- do.call(cbind, pred)

  # Report timing if verbose
  if (verbose) {
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "secs")
    cat("impute_forest_dt completed in", round(elapsed, 2), "seconds\n")
  }
  y0_post_mat
}

#y0_post <- dt |> multiple_impute_forest(ntrees = 2, nboot = 3, forest_cores = 2)


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
impute_forest <- function(.dt, ..., verbose = TRUE) {
  #' Clears any prior y0_post_* and y0_* summary columns
  #' o avoid stale values, then refits and imputes.
  drop_by_pattern(.dt, "^y0_post_")
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")

  y0_post <-  multiple_impute_forest(.dt, ..., verbose = verbose)
  inference_forest(.dt, y0_post)

  .dt
}


#' Summarize forest posterior draws
#'
#' Computes point predictions from a matrix of posterior draws and attaches them
#' to \code{.dt}. Currently computes only the posterior mean \code{y0_hat} and
#' stores each draw as \code{y0_post_1}, \code{y0_post_2}, ... .
#'
#' @param .dt A data.table modified by reference.
#' @param y0_post A numeric matrix of posterior draws with \code{nrow(.dt)} rows.
#' @return \code{.dt} modified by reference.
inference_forest <- function(.dt, y0_post) {
  nboot <- ncol(y0_post)
  # point prediction on original rows
  y0_hat <- apply(y0_post, 1, mean)
  .dt[, y0_hat := y0_hat]
  .dt[, (paste0("y0_post_", seq_len(nboot))) := data.frame(y0_post)]
  .dt
}

#dt |> impute_forest(ntrees = 2, nboot = 2, forest_cores = 1)



#' Predict from a fitted Poisson regression-tree forest
#'
#' Predicts using each \code{rpart} tree in \code{forest} and averages the
#' resulting predictions across trees. Returns predictions on the count scale.
#'
#' @param forest A \code{"dforest"} object (list of fitted \code{rpart} trees).
#' @param newdata A data.table containing \code{n} and the predictor columns
#'   required by the fitted trees.
#' @param ncores Number of cores for parallel prediction across trees.
#'
#' @details
#' \code{predict.rpart(method="poisson")} returns predictions on the rate scale.
#' This function averages those predicted rates across trees, then converts to
#' counts via \code{pi_hat * n / 1e5}.
#'
#' If \code{wts} is not present in \code{newdata}, it is added (set to 1) for
#' compatibility with training-time conventions.
#'
#' @return A numeric vector of predicted counts.
predict.dforest <- function(forest, newdata, ncores = 1) {
  # Add wts column if not present
  if (!"wts" %in% names(newdata)) newdata[, wts := 1]

  # Use mclapply for parallel prediction
  tree_predictions <- mclapply(forest, function(m) {
    predictions <- predict(m, newdata)
    predictions
  }, mc.cores = ncores)

  # Convert list of predictions to matrix for efficiency
  predictions_matrix <- do.call(cbind, tree_predictions)

  # Calculate row means (average predictions across trees)
  pi_hat <- rowMeans(predictions_matrix)

  # Apply scaling formula - note this returns a vector, not a data.table
  pi_hat * newdata$n / 10^5
}


#' Sample a subset of predictors for bagging
#'
#' Given a set of candidate column names, samples \code{mtry} predictors to keep
#' while always retaining essential columns (exposure, response, and weights).
#'
#' @param cols Character vector of column names.
#' @param mtry Number of non-mandatory predictors to sample. Defaults to
#'   \code{floor(sqrt(length(cols)))}.
#' @param always Columns that are always retained if present.
#'
#' @return Character vector of selected column names.
sample_feature <- function(cols
  , mtry = NULL
  , always = c("n", "y0", "wts")
) {
  if (is.null(mtry)) {
    mtry <- floor(sqrt(length(cols)))
  }
  always <- intersect(cols, always)
  cols <- setdiff(cols, always)
  mtry   <- min(mtry, length(cols))
  keep <- if (mtry < length(cols)) sample(cols, mtry) else cols
  c(always, keep)
}

#sample_feature(c("a", "b", "c", "d", "e"), 3)

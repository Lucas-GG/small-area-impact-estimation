
library(rpart)


impute_forest <- \(.data, ntrees = 20
  , k = 10,  min_bucket = 7, max_depth = 30
  , mtry = Inf
  , ncores = 1
  , exvar = TRUE
  , wtype = "exp"
  , shuffle_time = FALSE
) {
  #ideally we woudl recompute the propensity for eahc tree
  # but that is too expensive
  .edata <- .data  %>% expandx %>% add_prop
  .forest <- forest(.edata
    , ntrees = ntrees
    , k = k
    , min_bucket = min_bucket
    , max_depth = max_depth
    , mtry = mtry
    , ncores = ncores
    , exvar = exvar
    , wtype = "exp"
    , shuffle_time = FALSE
  )

  .data$y0_hat <-  predict_forest(.forest, .edata, ncores)
  .data
}

fit_single_tree_old <- \(.data, k, min_bucket, max_depth) {
  #' expecting BB weights
  if (is.null(.data$wts)) .data$wts <- rep(1, nrow(.data))
  #' making sure there is no outcome info
  .data$y <- NULL
  #' y25 is not treated as exogenous
  .data$y25 <- NULL

  m0 <- rpart(cbind(n / 10^5, y0) ~ ., data = .data
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

fit_single_tree <- function(.data, k, min_bucket, max_depth) {
  fit_single_tree_cpp(.data, k, min_bucket, max_depth)
}

cppFunction('
Rcpp::List fit_single_tree_cpp(Rcpp::DataFrame data, double k, int min_bucket, int max_depth) {
  Rcpp::NumericVector n = data["n"];
  Rcpp::NumericVector y0 = data["y0"];
  Rcpp::NumericVector wts = data.containsElementNamed("wts") ? data["wts"] : Rcpp::NumericVector(data.nrows(), 1.0);

  // Remove unnecessary columns
  //data.erase("y");
  //data.erase("y25");

  // Fit the rpart model
  Rcpp::Environment rpart_env = Rcpp::Environment::namespace_env("rpart");
  Rcpp::Function rpart = rpart_env["rpart"];
  Rcpp::Function rpart_control = rpart_env["rpart.control"];

  Rcpp::List control = rpart_control(Rcpp::Named("cp") = -INFINITY, Rcpp::Named("minbucket") = min_bucket, Rcpp::Named("maxdepth") = max_depth);
  Rcpp::List model = rpart(Rcpp::Named("formula") = Rcpp::Formula("cbind(n / 10^5, y0) ~ ."),
                           Rcpp::Named("data") = data,
                           Rcpp::Named("method") = "poisson",
                           Rcpp::Named("parms") = Rcpp::List::create(Rcpp::Named("shrink") = k),
                           Rcpp::Named("weights") = wts,
                           Rcpp::Named("control") = control);

  return model;
}
')

#'select subset of features for forest
rselect <- \(.data, mtry = NULL) {
  y <- .data %>% select(any_of(c("n", "y0", "wts", "y0_hat")))
  x <- .data %>% select(- any_of(c("n", "y0", "wts", "y0_hat")))
  h <- floor(sqrt(ncol(x)))
  if (is.null(mtry)) mtry <- h
  if (mtry == Inf) {
    return(.data)
  } else {
    x <- x %>% sample(mtry)
    .data <- data.frame(y, x)
    return(.data)
  }
}

#' this forest is base on shuflyng the start date
forest <- \(.data, ntrees = 20L
  , k = 10,  min_bucket = 7, max_depth = 30
  , mtry = Inf
  , ncores = 1
  , exvar = TRUE
  , wtype = "exp"
  , shuffle_time = FALSE
) {
  if (!exvar) cluster_boot <- \(.data, ...) .data

  list_of_trees <- mclapply(seq_len(ntrees), \(x) {
    print(x)
    .data %>%
      shuffle_start %>% # this set to missing y0[yr>yr1]
      expandx %>% #this recomputed propesnity every time
      filter(year < start_year) %>%
      cluster_boot(wtype = wtype, shuffle_time = shuffle_time) %>%
      select(-y, -y25) %>%
      rselect(mtry) %>%
      fit_single_tree(k, min_bucket, max_depth)
  }, mc.cores = ncores)
  return(list_of_trees)
}

predict_forest <- \(forest, .data, ncores) {
  if (is.null(.data$wts)) .data$wts <- rep(1, nrow(.data))
  #.newdata <- .data %>%  expandx #this is redundant
  # predict tree return the estimated rate (pi)
  pi_hat <- mclapply(forest, \(m) predict(m, .data), mc.cores = ncores) %>%
    do.call("cbind", .) %>%
    rowMeans()
  pi_hat * .data$n / 10^5
  # "vector" for Poisson trees it is the estimated response rate
}



#cluster_boot(dt1, "mammen") %>% head
#' generate fractional weights
#' which are similar to resampling units (ids)
shuffle_cl <- \(.data, wtype, cl = "i") {

  if (is.null(.data$wts)) .data$wts <- rep(1, nrow(.data))

  gen_w <- fwb:::make_gen_weights(wtype)
  .data$cl <- .data[[cl]]
  n_cl <- n_distinct(.data$cl)

  .wdata <- data.frame(
    cl = unique(.data$cl)
    , new_wts = as.vector(gen_w(n_cl, 1))
  )
  .data <- left_join(.data, .wdata, by = "cl") %>%
    mutate(wts = wts * new_wts, new_wts = NULL)
  k <- nrow(.data)
  .data <- .data %>%
    mutate(
      wts = wts / sum(wts) * k
      , cl = NULL
    )
  .data
}

cluster_boot <- \(.data, cl = "i", wtype = "mammen"
  , shuffle_time = TRUE
) {

  if (shuffle_time) {
    .data <- .data %>% shuffle_cl(wtype, cl = "year")
  }
  .data <- .data %>% shuffle_cl(wtype, cl = "i")

  return(.data)
}

inference_forest <- \(.data, ntrees = 20
  , k = 10,  min_bucket = 7, max_depth = 30
  , mtry = Inf
  , nboot = 20, wtype = "exp"
  , ncores = 1
  , shuffle_time = FALSE
  , exvar = TRUE
) {
  t0 <- Sys.time()
  .data$y0_post <- mclapply(seq_len(nboot), \(i) {
    .data %>%
      cluster_boot(wtype = wtype, shuffle_time = shuffle_time) %>%
      impute_forest(ntrees
        , k, min_bucket, max_depth, mtry, ncores, exvar
      ) %>%
      pull(y0_hat)
  }, mc.cores = ncores) %>%
    do.call("cbind", .)
  .data <- .data %>% post_summary

  t1 <- Sys.time()
  print(t1 - t0)
  return(.data)
}

mrpois <- \(m) rpois(prod(dim(m)), m) %>% matrix(ncol = ncol(m))


post_summary <- \(.data, include_mean = FALSE) {
  if (is.null(.data$y0_post)) stop("No posterior samples")

  mp <- .data$y0_post %>% mrpois
  if (include_mean) .data$y0_hat <- rowMeans(mp)
  .data$y0_sd <- dapply(mp, fsd, MARGIN = 1)
  .data$y0_pi <- mp == replicate(ncol(mp), .data$y) %>% rowMeans
  .data[, c("y0_lower", "y0_upper")] <-  
    dapply(mp, fquantile, MARGIN = 1, c(0.025, 0.975))
  .data
}


post_summary_old <- \(.data, include_mean = FALSE) {
  if (is.null(.data$y0_post)) stop("No posterior samples")
  if (include_mean) .data$y0_hat <- apply(.data$y0_post, 1, mean)

  .data$y0_sd <- dapply(.data$y0_post, sd)
  .data$eta_mean <- apply(log(.data$y0_post), 1, mean)
  .data$eta_sd <- apply(log(.data$y0_post), 1, sd)
  .data$y0_pi <- with(.data, pi_poilog(y, eta_mean, eta_sd))
  .data[, c("y0_lower", "y0_upper")] <- with(.data, ci_poilog(eta_mean, eta_sd))
  .data
}
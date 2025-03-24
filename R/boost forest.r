#------------------------------------------------------------------------------
#' Thus, in every step of the Poisson boosting machine
#' we receive updated working weights (w(m) i )i=1,...,n,
#' playing the role of the volumes in the Poisson model, and
#' replacing the role of the working responses in the generic GBM.
#' the logarithms of the working responses play the role of fixed offsets
#' that are enhanced by a next regression model/boosting step.
boost_single <- \(.idata
  , k, min_bucket, max_depth
) {

  tree <- rpart(cbind(y0_hat, y0) ~ .
    , data = .idata
    , method = "poisson"
    , parms = list(shrink = k)
    , control = rpart.control(cp = -Inf,
      , minbucket = min_bucket
      , maxdepth = max_depth
    )
  )
  tree
}
impute_boost_forest <- \(.idata, ntrees = 20
  , k = 1,  min_bucket = 7, max_depth = 30
  , beta = 1
) {
  .idata$y0_hat <-
    boost_forest(.idata, ntrees, k, min_bucket, max_depth) %>%
    predict_boost_forest(.idata) * .idata$y0_hat * beta +
    .idata$y0_hat * (1 - beta)
  .idata
}
predict_boost_forest <- \(bforest, .idata) {
  lapply(bforest, \(m) {
    .idata %>%
      add_avg %>%
      predict(m, .)
  }) %>%
    do.call("cbind", .) %>%
    apply(1, mean)
}
boost_forest <- \(.idata, ntrees = 2L
  , k = 1,  min_bucket = 7, max_depth = 30
  , mtry = NULL
) {
  lapply(seq_len(ntrees), \(x) {
    .idata %>%
      shuffle_start %>% # this set to missing y0[yr>yr1]
      filter(year < start_year) %>%
      add_avg %>%
      select(-y) %>%
      rselect(mtry) %>%
      boost_single(k, min_bucket, max_depth)
  })
}
#dt1 %>% boost_forest(2) %>% predict_boost_forest(dt1) %>% summary
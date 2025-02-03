#shuffle_start_mar <- \(.gdata) {
#  p0 <- cumsum(prop.table(table(.gdata$start_year)))
#  p0 <- p0[- length(p0)]
#  z0 <- qlogis(p0)
#  q <- length(z0)
#  n <- nrow(.gdata)
#  mz0 <- matrix(z0, n, q, byrow = TRUE)
#  eta <- drop(scale(ntile(.gdata$my, 5)))
#  cumpr <- plogis(mz0 + eta)
#  p <- t(apply(cumpr, 1, \(x) diff(c(0, x, 1))))
#  start_year <- apply(p, 1, \(pi) rmultinom(1, 1, pi)) %>%
#    apply(., 2, \(v) {
#      levels(ordered(.gdata$start_year))[which(v == 1)] %>%
#        as.numeric
#    })
#  .gdata$start_year <- start_year
#  return(.gdata)
#}

shuffle_start_mar <- \(.gdata) {
  v0  <- table(.gdata$start_year)
  v0 <- v0[-length(v0)]
  s <- NULL
  .gdata$size <- rank(.gdata$my)
  s <<- NULL
  sl <- lapply(seq_len(length(v0)), \(j) {
    a <- with(.gdata[!.gdata$i %in% s, ], sample(i, v0[j], prob = size))
    s <<- c(s, a)
    data.frame(i = a, start_year = as.numeric(names(v0))[j])
  }) %>% bind_rows
  .gdata$start_year <- NULL
  .gdata <- left_join(.gdata, sl, by = "i") %>%
    mutate(start_year = ifelse(is.na(start_year), Inf, start_year))
  return(.gdata)
}

shuffle_start <- \(.data, type = "mcar", a = 3) {
  .gdata <- .data %>%
    mutate(
      tots = ifelse(year < min(yr1) - a, Y + y25, NA)
      , totn = ifelse(year < min(yr1) - a, n + n25, NA)
    ) %>%
    group_by(i, start_year, st) %>%
    reframe(
      my = sum(tots, na.rm = TRUE) / sum(totn, na.rm = TRUE) * 10^5
      , start_year = first(start_year) - a
    )
  if (type == "mar") {
    .gdata <-  .gdata %>%
      shuffle_start_mar
  } else if (type == "mcar") {
    .gdata <- .gdata %>%
      mutate(start_year = sample(start_year))
  } else if (type == "mnar") {
    .gdata <- .gdata %>%
      group_by(st) %>%
      mutate(start_year = sample(start_year)) %>%
      ungroup
  }
  .data <- left_join(
    select(.data, - start_year)
    , select(.gdata, c(i, start_year))
    , by = "i"
  ) %>%
    arrange(i, year) %>%
    mutate(
      y0 = ifelse(year < start_year, y0, NA)
      , y25 = ifelse(year < start_year, y25, NA)
    )
  return(.data)
}

#df %>%
#  shuffle_start("MAR") %>%
#  mutate(
#    y = ifelse(year <= 2002, Y, NA)
#    , n = ifelse(year <= 2002, n, NA)
#  ) %>%
#  group_by(i, start_year, st) %>%
#  reframe(
#    my = sum(y, na.rm = TRUE) / sum(n, na.rm = TRUE) * 10^5
#  )  %>%
#  reframe(cor(my, start_year, method = "spearman"))

# lm(rank(start_year) ~ st, data = .) %>%
#  summary %>%
#  .[["adj.r.squared"]]

#  reframe(cor(my, start_year, method = "spearman"))



add_avg <- \(.data) {
  .data <- .data %>%
    #    group_by(i) %>%
    #    mutate(
    #      m_slp = ifelse(sum(!is.na(y0)) >1, coef(lm(y0 ~ year))[2], 0)
    #    ) %>%
    #    ungroup %>%
    mutate(
      mi_y0 = ave(y0, i, FUN = \(x) mean(x, na.rm = TRUE))
      , mt_y0 = ave(y0, year, FUN = \(x) mean(x, na.rm = TRUE))
      , mst_y0 = ave(y0, st, FUN = \(x) mean(x, na.rm = TRUE))
      , mc_y0 = ave(y0, start_year, FUN = \(x) mean(x, na.rm = TRUE)) #cohort
      , mi_y25 = ave(y25, i, FUN = \(x) mean(x, na.rm = TRUE))
      , mt_y25 = ave(y25, year, FUN = \(x) mean(x, na.rm = TRUE))
      , mst_y25 = ave(y25, st, FUN = \(x) mean(x, na.rm = TRUE))
      , mc_y25 = ave(y25, start_year, FUN = \(x) mean(x, na.rm = TRUE)) #cohort
      #, mtr_y0 = ave(y0, paste0(year, region), FUN = \(x) mean(x, na.rm = TRUE))
    )
  return(.data)
}

#dt0 %>% shuffle_start %>% add_avg %>% summary
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

fit_single_tree <- \(.data
  , k, min_bucket, max_depth
) {
  if (is.null(.data$wts)) .data$wts <- rep(1, nrow(.data))
  .data$Y <- NULL # this was eliminated before, but just in case
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

#dt1 <- dt0 %>% shuffle_start %>% add_avg
#yhat <- dt1 %>% fit_single_tree(.3, 50, 30) %>% predict(dt1)
#plot(dt0$Y, yhat)
#table(dt0$Y)

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
#dt1 %>% select(-Y) %>% boost_single(.3, 500, 30)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#select subset of features for forest
rselect <- \(.data, mtry = NULL) {
  if (is.null(mtry)) mtry <- floor(sqrt(ncol(.data) - 2))

  y <- .data %>% select(c(n, y0))
  x <- .data %>% select(-c(n, y0)) %>% sample(mtry)
  cbind(y, x)
}

irselect <- \(.idata, mtry = NULL) {
  if (is.null(mtry)) mtry <- floor(sqrt(ncol(.idata) - 2))

  y <- .idata %>% select(c(y0_hat, y0))
  x <- .idata %>% select(-c(y0_hat, y0)) %>% sample(mtry)
  cbind(y, x)
}

forest <- \(.data, ntrees = 20L
  , k = 1,  min_bucket = 7, max_depth = 30
  , mtry = NULL
) {
  lapply(seq_len(ntrees), \(x) {
    .data %>%
      shuffle_start %>% # this set to missing y0[yr>yr1]
      add_avg %>%
      filter(year < start_year) %>%
      select(-Y) %>%
      rselect(mtry) %>%
      fit_single_tree(k, min_bucket, max_depth)
  })#, mc.cores = 20)
}

#dt0 %>% shuffle_start %>% add_avg %>% fit_single_tree(1,1,30)
#dt1 <- dt0 %>% shuffle_start %>% add_avg
#dt1 %>% forest %>% predict_forest(dt1)

boost_forest <- \(.idata, ntrees = 2L
  , k = 1,  min_bucket = 7, max_depth = 30
  , mtry = NULL
) {
  lapply(seq_len(ntrees), \(x) {
    .idata %>%
      shuffle_start %>% # this set to missing y0[yr>yr1]
      filter(year < start_year) %>%
      add_avg %>%
      select(-Y) %>%
      irselect(mtry) %>%
      boost_single(k, min_bucket, max_depth)
  })
}
#dt1 %>% boost_forest(2) %>% predict_boost_forest(dt1) %>% summary


#------------------------------------------------------------------------------

#predict_forest(forest(2), dt0 %>% shuffle_start)
predict_forest <- \(forest, .data) {
  if (is.null(.data$wts)) .data$wts <- rep(1, nrow(.data))
  .newdata <- .data %>%  add_avg
  # predict tree return the estimated rate (pi)
  pi_hat <- lapply(forest, \(m) predict(m, .newdata)) %>%
    do.call("cbind", .) %>%
    apply(1, mean)
  pi_hat * (.data$n / 10^5)
  # "vector" for Poisson trees it is the estimated response rate
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

impute_forest <- \(.data, ntrees = 20
  , k = 1,  min_bucket = 7, max_depth = 30
) {
  .data$y0_hat <-
    forest(.data, ntrees, k, min_bucket, max_depth) %>%
    predict_forest(.data)
  .data
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

inferece_forest <- \(.idata, ntrees = 20
  , k = 1,  min_bucket = 7, max_depth = 30
  , nboot = 20
) {
  t0 <- Sys.time()
  .idata$y0_post <- lapply(seq_len(nboot), \(i) {
    .wdata <- left_join(
      .idata
      , data.frame(i = unique(.idata$i), wts = bb(n_distinct(.idata$i)))
      , by = "i"
    ) %>%
      mutate(wts = wts / sum(wts) * nrow(.idata))
    .wdata %>%
      forest(ntrees, k, min_bucket, max_depth) %>%
      predict_forest(.wdata)
  }) %>%
    do.call("cbind", .)

  #eta <- log(cbind(.idata$y0_post, y0_post))
  #r <- mclapply(2:1400, \(j) apply(eta[, 1:j], 1, sd), mc.cores = 20) %>% do.call("cbind", .)
  #r <- r - apply(eta, 1, sd)
  #d <- apply(r, 1, diff) %>% t
  #rd <- abs(d) / r[, -1399]
  #plot(c(0, 1500), c(0, .01), type = "n", xlab = "", ylab = "")
  #abline(h = c(.01, .001))
  #abline(v = c(200, 400))
  #points(20:1398, apply(rd, 2, mean) [20:1398], type = "l")

  .idata$y0_sd <- apply(.idata$y0_post, 1, sd)
  .idata$eta_mean <- log(.idata$y0_hat)
  .idata$eta_sd <- apply(log(.idata$y0_post), 1, sd)
  .idata$y0_pi <- with(.idata, pi_poilog(Y, eta_mean, eta_sd))
  .idata[, c("y0_lower", "y0_upper")] <-
    with(.idata, ci_poilog(eta_mean, eta_sd))
  t1 <- Sys.time()
  print(t1 - t0)
  return(.idata)
}


#n <- 1000
#eta <- rnorm(n)
#x <- rpois(n, exp(eta))
#x %>% quantile(c(.025, .975))

#ci_poilog(0, 1)

#test <-  dt0 %>%
#  shuffle_start %>%
#  add_avg %>%
#  impute_forest %>%

#filter(placebo_dt, replication == 1) %>% summary
#filter(placebo_dt, replication == 2) %>% summary

#filter(placebo_dt, replication == 1) %>%
#  inferece_forest %>%
#  filter(is.na(y0)) %>%
#  scoring

#filter(placebo_dt, replication == 2) %>%
#  inferece_forest %>%
#  filter(is.na(y0)) %>%
#  scoring



#test %>%
#  filter(is.na(y0)) %>%
#  scoring

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
assess_impute <- \(.idata) {
  .idata %>%
    filter(year >= start_year) %>%
    mutate(
      , d = Y - y0_hat
      , ae = ae(Y, y0_hat)
      , se = se(Y, y0_hat)
      , d2 = d2(Y, y0_hat)
      #, lns = log(y0_pi)
      #, ds = ds(Y, y0_hat, y0_sd)
      #, cvg = cvg(Y, y0_lower, y0_upper)
      #, intv = intv(Y, y0_lower, y0_upper)
      #, sigpred = y0_hat + y0_sd^2
    ) %>%
    reframe(
      across(c(d, ae, se, d2), mean)
      , se = sqrt(se)
      #, sigpred = sqrt(sigpred)
      , pcor = cor(Y, y0_hat)
      , scor = cor(Y, y0_hat, method = "spearman")
    )
}

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#Mixed effect (ME) model
#a secular trend, with random intercept and slopes. 

#library(lme4)
#impute_me <- \(.data) {
  #.data <- dt0 %>% shuffle_start
#  .data$t <- scale(.data$year)
#  fm <- formula(y0 ~ offset(log(n / 10^5)) + t + (t | i))
#  glmer(fm, data = .data, family = "poisson") %>%
#    predict(newdata = .data,  type = "response")  ->
#    .data$y0_hat
#  .data
#}


#ETWFE
#This can be improved adding a predicitve score
#library(fixest)
#impute_etwfe <- \(.data) {
#  .data <- dt0 %>% shuffle_start
#  fm <- formula(y0 ~ offset(log(n / 10^5)) +
#      i(start_year, i.year, ref = Inf) | start_year + year
#  )
#  suppressMessages(fepois(fm, .data, vcov = ~ i)) %>%
#    predict(newdata = .data,  type = "response") ->
#    .data$y0_hat
#  .data
#}
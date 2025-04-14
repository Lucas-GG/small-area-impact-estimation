
mrpois <- \(.m) {
  .m <- as.matrix(.m)
  rpois(prod(dim(.m)), .m) %>% matrix(ncol = ncol(.m))
}

bounds_sim <- \(.m){
  apply(mrpois(.m), 1, quantile, c(.025, .975)) %>% t
}

eval <- \(.dt
  , pi_type = "poilog"
  , bounds_type = "sim"
  , include_mean = FALSE
) {
  results <- copy(.dt)
  if (is.null(results$y0_post)) { #y0_post is no a nested data.frame
    y0_post_cols <- grep("y0_post", names(results), value = TRUE)
    y0_post <- results[, ..y0_post_cols] %>% as.matrix
  } else {
    y0_post <- results$y0_post  %>% as.matrix
  }
  dim(y0_post)

  # Calculate all values first
  y0_hat_val <- if (include_mean) rowMeans(y0_post) else NULL
  y0_sd_val <- fsd(t(y0_post))
  length(y0_hat_val)
  length(y0_sd_val)

  # Calculate PI values
  if (pi_type == "interval") {
    y0_pi <- rowMeans(replicate(ncol(y0_post), results$y) == y0_post)
  }

  if (pi_type == "poilog") {
    # Poisson Lognormal
    log_y0_post  <- log(y0_post)
    eta_mean <- rowMeans(log_y0_post)
    eta_sd   <- fsd(t(log_y0_post))
    y0_pi <- pi_poilog(results$y, eta_mean, eta_sd)
  }

  if (pi_type == "kernel") {
    #univariate local polynomial kernel density estimators
    F_y0 <- apply(y0_post, 1, kde1d:::kde1d, type = "d")
    y0_pi <- sapply(seq_len(length(F_y0))
      , \(i) kde1d:::dkde1d(results$y[i], F_y0[[i]])
    )
    y0_pi[y0_pi == 0] <- min(y0_pi[y0_pi > 0])
  }

  # Calculate bounds
  if (bounds_type == "sim") {
    bounds <- bounds_sim(y0_post)
  }

  if (bounds_type == "poilog") {
    log_y0_post  <- log(y0_post)
    eta_mean     <- rowMeans(log_y0_post)
    eta_sd       <- fsd(t(log_y0_post))
    bounds       <- with(.dt, ci_poilog(eta_mean, eta_sd))
  }

  if (pi_type == "interval") {
    bounds <- apply(y0_post, 1, quantile, c(.025, .975)) %>% t
  }

  if (include_mean) results[, y0_hat := y0_hat_val]
  results[, `:=`(
    y0_sd = y0_sd_val
    , y0_pi = y0_pi
    , y0_lower = bounds[, 1]
    , y0_upper = bounds[, 2]
  )]

  results
}



#change in number
ae <- \(y, y0) abs(y - y0)
se <- \(y, y0) (y - y0)^2

#absolute risk change
sae <- \(y, y0, n) abs((y - y0) / n * 10^5)
sse <- \(y, y0, n) ((y - y0) /  n * 10^5)^2

#relative risk change
rae <- \(y, y0) abs((y / y0) - 1)
rse <- \(y, y0) ((y / y0) - 1)^2

#coverage
cvg <- \(y, l, u) as.numeric(l <= y & y <= u)

#interval score 
intv <- \(y, l, u, α = .05) {
  (u - l) + 2 / α * (l - y) * (y < l) + 2 / α * (y - u) * (y > u)
}


#Poisson deviance
d2 <- \(y, y0) {
  ifelse(y > 0
    , 2 * (y * log(y / y0) + y0 - y)
    , 2 * y0
  )
}

#Dawid–Sebastiani score
ds <- \(y, y0, y0_sd) {
  sigpred <- sqrt(y0 + y0_sd^2)
  ((y - y0) / sigpred) ^ 2 + 2 * log(sigpred)
}

score <- \(m) {
  m %>%
    mutate(y = y
      , d = y - y0_hat
      , ae = ae(y, y0_hat)
      , se = se(y, y0_hat)
      , r = (y / y0_hat) - 1
      , rae = rae(y, y0_hat)
      , rse = rse(y, y0_hat)
      , dr = (y - y0_hat) / n * 10^5
      , sae = sae(y, y0_hat, n)
      , sse = sse(y, y0_hat, n)
      , d2 = d2(y, y0_hat)
      , lns = - log(y0_pi)
      , ds = ds(y, y0_hat, y0_sd)
      , cvg = cvg(y, y0_lower, y0_upper)
      , intv = intv(y, y0_lower, y0_upper)
      , sigpred = y0_hat + y0_sd^2
    )
}

summ_score <- \(.data) {
  .data %>%
    reframe(
      across(c(d:sigpred)
        , \(x) mean(x, na.rm = TRUE)
      )
      , pcor = cor(y, y0_hat)
      , scor = cor(y, y0_hat, method = "spearman")
    ) %>%
    mutate(
      se = sqrt(se)
      , rse = sqrt(rse)
      , sse = sqrt(sse)
      , sigpred = sqrt(sigpred)
    )
}

summ_score_dt <- function(.data) {
  # Convert to data.table if not already
  dt <- as.data.table(.data)

  # Get group columns from the data
  group_cols <- key(dt)
  if (is.null(group_cols)) {
    # If no grouping is set via keys, pass through ungrouped
    result <- dt[, lapply(.SD, mean, na.rm = TRUE), 
                .SDcols = c("d", "ae", "se", "r", "rae", "rse", "dr", "sae", "sse", 
                            "d2", "lns", "ds", "cvg", "intv", "sigpred")]

    result[, `:=`(
      pcor = cor(dt$y, dt$y0_hat),
      scor = cor(dt$y, dt$y0_hat, method = "spearman")
    )]
  } else {
    # With grouping
    result <- dt[, c(lapply(.SD, mean, na.rm = TRUE),
                    list(pcor = cor(y, y0_hat),
                         scor = cor(y, y0_hat, method = "spearman"))), 
                by = group_cols,
                .SDcols = c("d", "ae", "se", "r", "rae", "rse", "dr", "sae", "sse", 
                            "d2", "lns", "ds", "cvg", "intv", "sigpred")]
  }

  # Take square roots
  result[, `:=`(
    se = sqrt(se),
    rse = sqrt(rse),
    sse = sqrt(sse),
    sigpred = sqrt(sigpred)
  )]

  result
}

#aggregate
agg <- \(.idata, ...) {
  .data <- .idata %>%
    group_by(...) %>%
    reframe(across(c(n, y, y0_hat), sum))

  .data$y0_post <- .idata$y0_post %>%
    split(select(.idata, ...), drop = TRUE, lex.order = TRUE) %>%
    lapply(matrix, ncol = ncol(.idata$y0_post)) %>%
    lapply(\(m) apply(m, 2, sum)) %>%
    do.call(rbind, .)
  .data
}



score0 <- \(m) {
  m %>%
    mutate(y = y
      , d = y - y0_hat
      , ae = ae(y, y0_hat)
      , se = se(y, y0_hat)
      , r = (y / y0_hat) - 1
      , rae = rae(y, y0_hat)
      , rse = rse(y, y0_hat)
      , dr = (y - y0_hat) / n * 10^5
      , sae = sae(y, y0_hat, n)
      , sse = sse(y, y0_hat, n)
      , d2 = d2(y, y0_hat)
      #, lns = log(y0_pi)
      #, ds = ds(y, y0_hat, y0_sd)
      #, cvg = cvg(y, y0_lower, y0_upper)
      #, intv = intv(y, y0_lower, y0_upper)
      #, sigpred = y0_hat + y0_sd^2
    ) %>%
    reframe(
      across(c(d:d2)
        , \(x) mean(x, na.rm = TRUE)
      )
      , pcor = cor(y, y0_hat)
      , scor = cor(y, y0_hat, method = "spearman")
    ) %>%
    mutate(
      se = sqrt(se)
      , rse = sqrt(rse)
      , sse = sqrt(sse) # this is the root mean square error (it include the bias)
      #, sigpred = sqrt(sigpred)
    )
}
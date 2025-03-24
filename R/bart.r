

library(dbarts)
impute_bart <- \(.data) {
  t0 <- Sys.time()
  train <- filter(.data, year < first_year)
  test <- .data
  m <- bart2(I(y0 / n * 10^5) ~ . - y
    , data = train, test = test
    , n.samples = 2000L, n.burn = 2000L
    , n.chains = 16L
    , n.trees = 200
    , seed = 0203
  )
  n <- .data$n

  .data$y0_hat <- m$yhat.test.mean * n / 10^5

  .data$y0_post <- apply(m$yhat.test, 3, c) %>% t
  .data$y0_post <- .data$y0_post * n / 10^5
  .data$y0_sd <- apply(.data$y0_post, 1, sd)
  #ecdf(.data$y0_post[1, ]) %>% plot
  #hist(.data$y0_post[1, ])
  .data$eta_mean <- log(.data$y0_hat)
  .data$eta_sd <- apply(log(.data$y0_post), 1, sd)
  .data$y0_pi <- with(.data, pi_poilog(y, eta_mean, eta_sd))
  .data[, c("y0_lower", "y0_upper")] <-
    with(.data, ci_poilog(eta_mean, eta_sd))
  t1 <- Sys.time()
  print(t1 - t0)

  return(.data)
}

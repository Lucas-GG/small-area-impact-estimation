
mrpois <- \(m) rpois(prod(dim(m)), m) %>% matrix(ncol = ncol(m))


eval_sim <- \(.data, include_mean = FALSE) {
  if (is.null(.data$y0_post)) stop("No posterior samples")

  mp <- .data$y0_post %>% mrpois
  if (include_mean) .data$y0_hat <- rowMeans(mp)
  .data$y0_sd <- dapply(mp, fsd, MARGIN = 1)
  .data$y0_pi <- mp == replicate(ncol(mp), .data$y) %>% rowMeans
  .data[, c("y0_lower", "y0_upper")] <-  
    dapply(mp, fquantile, MARGIN = 1, c(0.025, 0.975))
  .data
}


eval_poilog <- \(.data, include_mean = FALSE) {
  if (is.null(.data$y0_post)) stop("No posterior samples")
  if (include_mean) .data$y0_hat <- apply(.data$y0_post, 1, mean)

  .data$y0_sd <- dapply(.data$y0_post, sd)
  .data$eta_mean <- apply(log(.data$y0_post), 1, mean)
  .data$eta_sd <- apply(log(.data$y0_post), 1, sd)
  .data$y0_pi <- with(.data, pi_poilog(y, eta_mean, eta_sd))
  .data[, c("y0_lower", "y0_upper")] <- with(.data, ci_poilog(eta_mean, eta_sd))
  .data
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
      , lns = log(y0_pi)
      , ds = ds(y, y0_hat, y0_sd)
      , cvg = cvg(y, y0_lower, y0_upper)
      , intv = intv(y, y0_lower, y0_upper)
      , sigpred = y0_hat + y0_sd^2
    )
}

summ_score <- \(.data) {
  .data%>%
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
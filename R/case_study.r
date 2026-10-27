# Case study: aggregate effects, heterogeneity test p-values against the placebo null, and the
# effect-modification model with its placebo calibration. Requires performance.r, injection.r.

#' Cumulative county counts through event year s, from the post-adoption rows (s >= 0) of a
#' fitted panel: observed y, latent predictive draws (un-noised) and population n. Returns one
#' row per county and event year with the draws as a matrix column set.
case_cumulative <- function(fits, s_max = 4L) {
  draws <- grep("^y0_post_", names(fits), value = TRUE)
  f <- fits[is.finite(start_year) & s >= 0 & s <= s_max & is.finite(y0_hat)]
  setorderv(f, c("spec", "model", "nc", "i", "s"))
  f[, c("y", "n", "y0_hat", draws) := lapply(.SD, cumsum), by = .(spec, model, nc, i), .SDcols = c("y", "n", "y0_hat", draws)]
  f
}

#' Aggregate effect by event year: mean over counties of the cumulative observed minus predicted
#' count, and the rate ratio sum(O) / sum(mu), with 95% intervals from Poisson-noised predictive
#' draws (as in the placebo scoring). `balanced` keeps counties observed at every event year 0..s_max.
case_effects <- function(cum, seed = 20261022L, balanced = FALSE, s_max = 4L) {
  draws <- grep("^y0_post_", names(cum), value = TRUE)
  d <- copy(cum)
  if (balanced) d <- d[d[, .I[max(s) == s_max], by = .(spec, model, nc, i)]$V1]
  # Poisson noise with a seed per specification, model and county group, so a model's interval
  # depends only on its own fit (one stream over all models made it shift with the others').
  key <- d[, paste(spec, model, nc)]
  for (g in unique(key)) {
    rows <- which(key == g)
    set.seed(seed + sum(utf8ToInt(g)))
    d[rows, (draws) := lapply(.SD, function(lam) rpois(.N, lam)), .SDcols = draws]
  }
  d[, {
    pred <- colMeans(as.matrix(.SD)); obs <- mean(y)
    eff <- obs - pred; rr <- sum(y) / colSums(as.matrix(.SD))
    .(counties = .N, observed = obs, predicted = mean(y0_hat),
      effect = obs - mean(y0_hat), effect_lo = quantile(eff, 0.025), effect_hi = quantile(eff, 0.975),
      rr = sum(y) / sum(y0_hat), rr_lo = quantile(rr, 0.025), rr_hi = quantile(rr, 0.975))
  }, by = .(spec, model, nc, s), .SDcols = draws]
}

#' p-value of an observed statistic against its placebo null: the share of null replications at
#' least as large, with the +1 correction, (1 + #{null >= obs}) / (1 + R).
placebo_p <- function(obs, null) (1 + sum(null >= obs)) / (1 + length(null))

#' Heterogeneity statistics on random subsets of `size` counties (matched to the placebo
#' evaluation sets), B times with a fixed seed; returns the statistics per subset.
heterogeneity_subsets <- function(h, size = 45L, B = 200L, seed = 20261023L, O = "y") {
  set.seed(seed)
  rbindlist(lapply(seq_len(B), function(b) {
    k <- sample(nrow(h), min(size, nrow(h)))
    heterogeneity_stats(h[k][, r := b], O = O, by = "r")
  }))
}

#' Effect-modification model over predictive draws. For each draw d of the counterfactual expected
#' counts mu^(d): quasi-Poisson regression O ~ x [+ z], offset log mu^(d); the dispersion is floored
#' at 1 for the SEs. Uncertainty: n_imp normal draws N(beta_d, SE_d^2) per d, pooled; Rubin's rules
#' as a check. `x` and `z` are numeric vectors over counties; returns a one-row summary.
em_fit <- function(O, D, x, z = NULL, n_imp = 5L, seed = 1L) {
  X <- cbind(1, x, z)
  per <- t(vapply(seq_len(ncol(D)), function(d) {
    fit <- stats::glm.fit(X, O, offset = log(D[, d]), family = stats::poisson())
    mu <- fit$fitted.values
    phi <- sum((O - mu)^2 / mu) / (length(O) - ncol(X))
    se <- sqrt(diag(chol2inv(qr.R(fit$qr))))[2]
    c(beta = unname(fit$coefficients[2]), se = se, phi = phi, converged = fit$converged)
  }, numeric(4)))
  se_floor <- per[, "se"] * sqrt(pmax(per[, "phi"], 1))
  set.seed(seed)
  tilde <- rnorm(length(se_floor) * n_imp, rep(per[, "beta"], each = n_imp), rep(se_floor, each = n_imp))
  Dn <- nrow(per); within <- mean(se_floor^2); between <- stats::var(per[, "beta"])
  total <- within + (1 + 1 / Dn) * between
  data.table(beta = mean(tilde), lo = unname(stats::quantile(tilde, 0.025)), hi = unname(stats::quantile(tilde, 0.975)),
    beta_hat_mean = mean(per[, "beta"]), rubin_se = sqrt(total), rubin_lo = mean(per[, "beta"]) - 1.96 * sqrt(total),
    rubin_hi = mean(per[, "beta"]) + 1.96 * sqrt(total), share_between = (1 + 1 / Dn) * between / total,
    phi_median = stats::median(per[, "phi"]), phi_min = min(per[, "phi"]), phi_max = max(per[, "phi"]),
    share_phi_below_1 = mean(per[, "phi"] < 1), nonconverged = sum(per[, "converged"] == 0), counties = length(O))
}

#' Long-training exposure transformations: percentile rank among the analysed counties scaled to
#' [0, 1] (ties averaged; primary), and the raw share in units of 10 percentage points (secondary).
em_rank <- function(share) (data.table::frank(share, ties.method = "average") - 1) / (length(share) - 1)
em_per10 <- function(share) share * 10

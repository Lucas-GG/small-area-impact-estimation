# Placebo performance: per-replication aggregate scores and their summaries.
# Requires scoring.r (score_avg) and data.table.

#' Poisson noise on the predictive draws (latent means), by reference. Rows are put in a fixed
#' order and each scenario x county group x model gets its own seed (a stable hash of its label),
#' so a model's noise, and so its intervals and scores, do not depend on which other models are
#' in the table. Used by aggregate_placebo() and injection_base().
add_poisson_noise <- function(fits, draws, seed = 0203) {
  setorderv(fits, c("dgp", "nc", "model", "r", "i", "s"))
  key <- paste(fits$dgp, fits$nc, fits$model)
  for (g in unique(key)) {
    rows <- which(key == g)
    set.seed(seed + stable_hash(g))
    fits[rows, (draws) := lapply(.SD, function(lam) rpois(.N, lam)), .SDcols = draws]
  }
  invisible(fits)
}

#' Stable integer hash of a string (for per-group seeds).
stable_hash <- function(s) {
  h <- 0
  for (x in utf8ToInt(s)) h <- (h * 31 + x) %% 1000003
  as.integer(h)
}

#' Score a prepared placebo sample (03a). Returns two per-replication ledgers:
#'   aggregate: the mean across counties of the cumulative observed count, prediction
#'     and predictive draws, scored (one row per county group, scenario, model, event
#'     year and replication);
#'   county: county-level effect metrics (see county_effect_metrics).
#' Draws are latent means; Poisson noise turns them into predictive counts.
aggregate_placebo <- function(fits, seed = 0203, tau = 0) {
  draws <- grep("^y0_post_", names(fits), value = TRUE)
  add_poisson_noise(fits, draws, seed)
  fits[, c("n", "y", "y0_hat", draws) := lapply(.SD, cumsum),
    by = .(dgp, nc, model, r, i), .SDcols = c("n", "y", "y0_hat", draws)]
  county <- county_effect_metrics(fits, draws, tau = tau)
  agg <- fits[, c(lapply(.SD, mean), list(
    counties = .N,
    # Poisson noise floor: variance of the county-mean count if each expected count
    # were known exactly, estimated with the observed counts (E[y] = mu).
    noise_var = sum(y) / .N^2,
    pcor = cor(y, y0_hat, use = "complete.obs"),
    scor = cor(y, y0_hat, use = "complete.obs", method = "spearman")
  )), by = .(nc, dgp, model, s, r), .SDcols = c("y", "y0_hat", draws)]
  score_avg(agg)
  list(aggregate = agg[, setdiff(names(agg), draws), with = FALSE], county = county)
}

#' County-level effect metrics within each replication, for cumulative effects
#' through event year s. The estimated effect is observed minus predicted
#' untreated count; its 95% interval is observed minus the 97.5% and 2.5%
#' quantiles of the predictive draws. `tau` is the true county effect: 0 in the
#' placebo, the injected effect in simulations (a scalar or a column name).
#' Counties are weighted equally.
county_effect_metrics <- function(fits, draws, tau = 0, level = 0.95) {
  a <- (1 - level) / 2
  dm <- as.matrix(fits[, ..draws])
  q <- matrixStats::rowQuantiles(dm, probs = c(a, 1 - a))
  truth <- if (is.character(tau)) fits[[tau]] else rep(tau, nrow(fits))
  # CRPS of the predictive distribution of the effect, y - draws, at the true effect;
  # equal to the CRPS of the draws at y - truth (translation invariance).
  crps <- crps_rows(dm, fits$y - truth)
  rm(dm)
  z <- fits[, .(nc, dgp, model, s, r, i, y_obs = y, est = y - y0_hat, lower = y - q[, 2], upper = y - q[, 1])]
  z[, `:=`(truth = truth, crps = crps)]
  z[, .(
    counties = .N,
    missing = sum(!is.finite(est)),
    # 1. county-effect RMSE (root-PEHE)
    root_pehe = sqrt(mean((est - truth)^2, na.rm = TRUE)),
    # 2. centered heterogeneity RMSE: between-county differences net of aggregate error
    centered_rmse = sqrt(mean(((est - mean(est, na.rm = TRUE)) - (truth - mean(truth)))^2, na.rm = TRUE)),
    # 3. county-level interval coverage and 4. mean interval width
    county_coverage = mean(lower <= truth & truth <= upper, na.rm = TRUE),
    county_width = mean(upper - lower, na.rm = TRUE),
    # 5. county CRPS: accuracy, calibration and sharpness of the county predictive distribution
    county_crps = mean(crps, na.rm = TRUE),
    # Noise floor: county-effect RMSE if each expected count were known exactly
    # (error y - mu, variance mu), estimated with the observed counts (E[y] = mu).
    noise_floor = sqrt(mean(y_obs))
  ), by = .(nc, dgp, model, s, r)]
}

#' Sample CRPS for each row of a draw matrix against its observation, in chunks:
#' mean |x - obs| minus half the mean pairwise distance, via within-row ranks
#' (identical to crps_sample() row by row).
crps_rows <- function(m, obs, chunk = 2e5) {
  S <- ncol(m)
  out <- numeric(nrow(m))
  for (start in seq(1, nrow(m), by = chunk)) {
    k <- start:min(nrow(m), start + chunk - 1)
    mk <- m[k, , drop = FALSE]
    R <- matrixStats::rowRanks(mk, ties.method = "average")
    out[k] <- rowMeans(abs(mk - obs[k])) - rowSums((2 * R - S - 1) * mk) / S^2
  }
  out
}

#' Summary of the county-level metrics across replications (means with MC SEs).
county_summary <- function(county, by = c("specification", "nc", "dgp", "model", "s")) {
  metrics <- c("root_pehe", "centered_rmse", "county_coverage", "county_width", "county_crps", "noise_floor")
  county[, rbindlist(c(lapply(metrics, function(m) c(list(metric = m), mc_mean(get(m)))),
    list(list(metric = "counties", estimate = mean(counties), mc_se = NA_real_, n = .N),
         list(metric = "missing", estimate = sum(missing), mc_se = NA_real_, n = .N)))), by = by]
}

# Simulation mean with Monte Carlo SE.
mc_mean <- function(x) {
  x <- x[is.finite(x)]
  list(estimate = mean(x), mc_se = if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_, n = length(x))
}

# RMSE from unrounded squared errors, with a delta-method Monte Carlo SE.
mc_rmse <- function(sq) {
  sq <- sq[is.finite(sq)]
  r <- sqrt(mean(sq))
  list(estimate = r, mc_se = if (r > 0 && length(sq) > 1) sd(sq) / (2 * r * sqrt(length(sq))) else NA_real_, n = length(sq))
}

#' Summary of the per-replication scores, one row per metric with its MC SE.
performance_summary <- function(ledger, by = c("specification", "nc", "dgp", "model", "s")) {
  ledger[, {
    rows <- list(
      RMSE = mc_rmse(sq_diff), bias = mc_mean(diff), MAE = mc_mean(abs_diff),
      relative_RMSE = mc_rmse(sq_rel_diff), coverage = mc_mean(cvg), interval_score = mc_mean(intv),
      CRPS = mc_mean(crps), deviance = mc_mean(deviance), noise_floor_RMSE = mc_rmse(noise_var))
    rbindlist(lapply(names(rows), function(m) c(list(metric = m), rows[[m]])))
  }, by = by]
}

#' Paired percent difference of a model's metric from a reference, with a
#' delta-method MC SE from the paired replications. RMSE-type metrics compare
#' sqrt(mean squared error); others compare means.
paired_change <- function(x, b, rmse = FALSE) {
  ok <- is.finite(x) & is.finite(b); x <- x[ok]; b <- b[ok]
  power <- if (rmse) .5 else 1
  ratio <- (mean(x) / mean(b))^power
  influence <- if (rmse) .5 * ratio * (x / mean(x) - b / mean(b)) else (x - mean(x)) / mean(b) - mean(x) * (b - mean(b)) / mean(b)^2
  list(percent_change = 100 * (ratio - 1), mc_se = 100 * sd(influence) / sqrt(length(x)), n = length(x))
}

change_metrics <- c(RMSE = "sq_diff", MAE = "abs_diff", interval_score = "intv", CRPS = "crps")

#' Percent change of each model relative to a reference model, paired by replication.
relative_to_model <- function(ledger, reference = "TWFE", by = c("specification", "nc", "dgp", "s")) {
  ref <- ledger[model == reference, c(by, "r", unname(change_metrics)), with = FALSE]
  setnames(ref, unname(change_metrics), paste0("ref_", change_metrics))
  z <- merge(ledger, ref, by = c(by, "r"))
  z[, rbindlist(lapply(names(change_metrics), function(m) {
    col <- change_metrics[[m]]
    c(list(metric = m), paired_change(get(col), get(paste0("ref_", col)), rmse = m == "RMSE"))
  })), by = c(by, "model")]
}

#' Percent change from one specification to another, per model, paired by replication
#' (e.g. adding covariates, or varying covariate slopes instead of common ones).
specification_change <- function(ledger, from = "Without covariates", to = "With covariates",
  by = c("nc", "dgp", "model", "s")) {
  w <- ledger[specification == from, c(by, "r", unname(change_metrics)), with = FALSE]
  setnames(w, unname(change_metrics), paste0("ref_", change_metrics))
  z <- merge(ledger[specification == to], w, by = c(by, "r"))
  z[, rbindlist(lapply(names(change_metrics), function(m) {
    col <- change_metrics[[m]]
    c(list(metric = m), paired_change(get(col), get(paste0("ref_", col)), rmse = m == "RMSE"))
  })), by = by]
}

#' Worst-case regret: each model's RMSE over the best model's RMSE in the same
#' scenario, then the worst scenario. Computed within each specification.
worst_case_regret <- function(perf, by = c("specification", "nc", "s")) {
  r <- perf[metric == "RMSE", .(dgp, model, rmse = estimate), by = by]
  r[, regret := rmse / min(rmse) - 1, by = c(by, "dgp")]
  r[, .(worst_regret = max(regret), worst_scenario = dgp[which.max(regret)]), by = c(by, "model")]
}

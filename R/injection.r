# Injected-effect extension of the placebo study.
# Models were fit only on cells untreated under the simulated schedule, so the saved
# predictive draws stay valid; only the "observed" outcomes of pseudo-treated cells change.
# Requires scoring.r (score_avg), performance.r and aux.r (set_y0).

#' County-year expected counts for the increase rule and the size measure: a Poisson
#' GLMM with year effects and a county random intercept, fitted on all real untreated
#' county-years; predictions (with the county effect) for every county-year of the panel.
injection_mu <- function(panel) {
  dt <- data.table::copy(data.table::as.data.table(panel))
  set_y0(dt)
  fit <- lme4::glmer(y0 ~ factor(year) + (1 | i) + offset(log(n)), data = dt[!is.na(y0)], family = poisson)
  dt[, mu_hat := stats::predict(fit, newdata = dt, type = "response", allow.new.levels = FALSE)]
  stopifnot(all(is.finite(dt$mu_hat) & dt$mu_hat > 0))
  size <- dt[, .(size = mean(mu_hat), nc = ifelse(urb[1] != "6.NonCore", "Core", "NonCore")), by = i]
  # Quartile of mean annual expected count, within county group; fixed across simulations.
  size[, size_quartile := cut(size, quantile(size, 0:4 / 4), include.lowest = TRUE, labels = FALSE), by = nc]
  list(mu = dt[, .(i, year, mu_hat)], size = size, fit = fit)
}

#' Everything that does not depend on the injected outcomes, computed once per
#' specification exactly as aggregate_placebo() does (same seed, row order, noise).
injection_base <- function(fits, seed = 0203) {
  draws <- grep("^y0_post_", names(fits), value = TRUE)
  add_poisson_noise(fits, draws, seed)
  fits[, y_annual := y]
  fits[, c("n", "y", "y0_hat", draws) := lapply(.SD, cumsum),
    by = .(dgp, nc, model, r, i), .SDcols = c("n", "y", "y0_hat", draws)]
  q <- matrixStats::rowQuantiles(as.matrix(fits[, ..draws]), probs = c(0.025, 0.975))
  county <- fits[, .(nc, dgp, model, s, r, i, y0_hat, q_lo = q[, 1], q_hi = q[, 2])]
  agg <- fits[, c(lapply(.SD, mean), list(counties = .N)), by = .(nc, dgp, model, s, r),
    .SDcols = c("y0_hat", draws)]
  # Observed annual outcomes are the same for every model: one row per cell.
  cells <- unique(fits[, .(dgp, nc, r, i, s, year, y = y_annual)])
  stopifnot(!anyDuplicated(cells[, .(dgp, r, i, s)]))
  list(county = county, agg = agg, cells = cells, draws = draws)
}

#' Scenario grid of the injected-effect study (Appendix C.3). RR 0.5, 0.6 and 1.5 are detectability
#' benchmarks, not plausible program effects.
injection_scenarios <- function() {
  g <- data.table::rbindlist(list(
    data.table::data.table(scenario = "E0", type = "uniform", rr = 1, sigma = NA_real_),
    data.table::data.table(scenario = "E1", type = "uniform", rr = c(0.9, 0.8, 0.7, 0.6, 0.5), sigma = NA_real_),
    data.table::data.table(scenario = "E2", type = "uniform", rr = c(1.1, 1.2, 1.3, 1.5), sigma = NA_real_),
    data.table::data.table(scenario = "E3", type = "lognormal", rr = 0.8, sigma = c(0.25, 0.5)),
    data.table::data.table(scenario = "E4", type = "quartile", rr = 0.75, sigma = NA_real_),
    data.table::data.table(scenario = "E5", type = "dynamic", rr = 0.7, sigma = NA_real_)
  ))
  g[, label := data.table::fcase(
    scenario == "E0", "Null",
    scenario %in% c("E1", "E2"), paste0("Uniform RR ", rr, ifelse(rr %in% c(0.5, 0.6, 1.5), " (benchmark)", "")),
    scenario == "E3", paste("Heterogeneous, sigma", sigma),
    scenario == "E4", "Size quartiles, RR 0.6-0.9",
    scenario == "E5", "Dynamic, RR 1.0 to 0.7")]
  g[, id := paste0(scenario, "_", seq_len(.N)), by = scenario][]
}

#' Injected outcomes for one scenario, the same for every model (common random numbers).
#' Reductions: y* ~ Binomial(y, RR) (exact for Poisson outcomes). Increases:
#' y* = y + Poisson((RR - 1) mu_hat), independent of y (exact if mu_hat = mu).
#' `target` is the county-level RR target: RR_i, or for E5 the injected RRs weighted
#' by the county's untreated counts, sum(RR_t y_t) / sum(y_t), cumulated over event years.
inject_outcomes <- function(cells, sc, mu, size, seed) {
  x <- copy(cells)
  setorderv(x, c("dgp", "r", "i", "s"))
  set.seed(seed)
  county_rr <- unique(x[, .(dgp, r, i)])
  if (sc$type %in% c("uniform", "dynamic")) county_rr[, rr_i := sc$rr]
  if (sc$type == "lognormal") county_rr[, rr_i := exp(rnorm(.N, log(sc$rr), sc$sigma))]
  if (sc$type == "quartile") county_rr[size, rr_i := c(0.6, 0.7, 0.8, 0.9)[i.size_quartile], on = "i"]
  stopifnot(!anyNA(county_rr$rr_i))
  x[county_rr, rr_i := i.rr_i, on = .(dgp, r, i)]
  x[, rr := if (sc$type == "dynamic") 1 - (1 - sc$rr) * s / 4 else rr_i]
  x[mu, mu_hat := i.mu_hat, on = .(i, year)]
  stopifnot(!anyNA(x$mu_hat))
  x[, rule := fifelse(rr < 1, "thinning", fifelse(rr > 1, "increase", "none"))]
  x[, y_star := y]
  x[rr < 1, y_star := rbinom(.N, y, rr)]
  x[rr > 1, y_star := y + rpois(.N, (rr - 1) * mu_hat)]
  x[, `:=`(y_cum = cumsum(y), y_star_cum = cumsum(y_star),
    target = if (sc$type == "dynamic") cumsum(rr * y) / cumsum(y) else rr_i), by = .(dgp, r, i)]
  x
}

#' Score one scenario with the existing estimator: the per-replication aggregate
#' ledger (as in 03b, with y*) and the county table with its targets.
score_injection <- function(base, inj) {
  county <- copy(base$county)
  county[inj, `:=`(y = i.y_star_cum, y_obs = i.y_cum, target = i.target), on = .(dgp, r, i, s)]
  stopifnot(!anyNA(county$y))
  y_agg <- county[, .(y = mean(y), noise_var = sum(y) / .N^2, y_obs_mean = mean(y_obs),
    pcor = cor(y, y0_hat, use = "complete.obs"), scor = cor(y, y0_hat, use = "complete.obs", method = "spearman")),
    by = .(nc, dgp, model, s, r)]
  agg <- merge(base$agg, y_agg, by = c("nc", "dgp", "model", "s", "r"))
  score_avg(agg)
  list(aggregate = agg[, setdiff(names(agg), base$draws), with = FALSE], county = county)
}

#' Per-replication aggregate metrics. The interval for the effect y* - Y_hat excludes
#' zero when y* falls below (reduction) or above (increase) the predictive interval.
#' Aggregate RR = sum(y*)/sum(Y_hat), interval y*/upper .. y*/lower; its target is the
#' injected RR (uniform) or the realized ratio sum(y*)/sum(y) (E3-E5), which is also
#' reported for every scenario.
injection_aggregate_metrics <- function(agg, sc) {
  agg[, .(nc, dgp, model, s, r,
    effect = y - y0_hat,
    # the complement of cvg(): same closed interval and tolerance, so E0 fp = 1 - coverage
    below = as.numeric(y < y0_lower - 1e-9 * pmax(1, abs(y))), above = as.numeric(y > y0_upper + 1e-9 * pmax(1, abs(y))),
    rr_hat = y / y0_hat,
    rr_realized = y / y_obs_mean,
    rr_target = if (sc$type == "uniform") sc$rr else y / y_obs_mean,
    rr_lower = y / y0_upper, rr_upper = ifelse(y0_lower > 0, y / y0_lower, Inf))]
}

#' County-level rate-ratio metrics per replication (all counties and by size quartile).
#' Counties with no observed events in the window are excluded whatever their RR_i or
#' rule; draws summing to zero give an unbounded upper limit, counted separately.
#' The estimator treats y* as fixed, so county intervals omit outcome Poisson noise.
injection_county_metrics <- function(county, size) {
  x <- county[y_obs > 0]
  x[, `:=`(est = y / y0_hat, lo = y / q_hi, hi = ifelse(q_lo > 0, y / q_lo, Inf))]
  x[size, size_quartile := paste0("Q", i.size_quartile), on = "i"]
  summarise <- function(d, by) d[, .(counties = .N, unbounded = sum(!is.finite(hi)),
    rr_bias = mean(est - target), rr_rmse = sqrt(mean((est - target)^2)),
    rr_coverage = mean(lo <= target & target <= hi),
    iqr_est = IQR(est), iqr_target = IQR(target),
    cor_est_target = if (isTRUE(sd(target) > 0)) cor(est, target, method = "spearman") else NA_real_), by = by]
  excluded <- county[y_obs == 0, .(zero_excluded = .N), by = .(nc, dgp, model, s, r)]
  all <- merge(summarise(x, c("nc", "dgp", "model", "s", "r")), excluded, by = c("nc", "dgp", "model", "s", "r"), all.x = TRUE)
  all[is.na(zero_excluded), zero_excluded := 0L][, size_quartile := "All"]
  rbind(all, summarise(x, c("nc", "dgp", "model", "s", "r", "size_quartile"))[, zero_excluded := NA_integer_], fill = TRUE)
}

#' Means across replications with MC SEs for every numeric metric column.
injection_summary <- function(d, by) {
  metrics <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], c(by, "r"))
  d[, rbindlist(lapply(metrics, function(m) {
    v <- get(m); v <- v[is.finite(v)]
    list(metric = m, estimate = if (length(v)) mean(v) else NA_real_,
      mc_se = if (length(v) > 1) sd(v) / sqrt(length(v)) else NA_real_, n = length(v))
  })), by = by]
}

#' Minimum detectable effect: the RR at which power reaches `target`, interpolated
#' linearly on the grid moving away from RR = 1, or NA if never reached.
#' `rr` and `power` cover one direction, including RR = 1.
mde <- function(rr, power, target = 0.8) {
  o <- order(abs(rr - 1)); rr <- rr[o]; power <- power[o]
  k <- which(power >= target)[1]
  if (is.na(k)) return(NA_real_)
  if (k == 1) return(rr[1])
  rr[k - 1] + (target - power[k - 1]) / (power[k] - power[k - 1]) * (rr[k] - rr[k - 1])
}

#' MDE table: primary (correct-sign interval test) and size-adjusted (reject when the
#' estimated effect falls outside the 2.5-97.5% range of the model's E0 estimates),
#' separately for reductions and increases, with the E0 false-positive rate.
mde_table <- function(metrics, grid, by = c("specification", "nc", "dgp", "model", "s"), fp_limit = 0.10) {
  u <- merge(metrics, grid[type == "uniform", .(id, rr_grid = rr)], by = "id")
  null <- u[rr_grid == 1, .(fp_rate = mean(below + above), q_lo = quantile(effect, 0.025), q_hi = quantile(effect, 0.975)), by = by]
  u <- merge(u, null, by = by)
  pw <- u[, .(primary_reduction = mean(below), primary_increase = mean(above),
    adjusted_reduction = mean(effect < q_lo), adjusted_increase = mean(effect > q_hi)), by = c(by, "rr_grid")]
  out <- pw[, .(
    mde_reduction = mde(rr_grid[rr_grid <= 1], primary_reduction[rr_grid <= 1]),
    mde_increase = mde(rr_grid[rr_grid >= 1], primary_increase[rr_grid >= 1]),
    mde_reduction_adjusted = mde(rr_grid[rr_grid <= 1], adjusted_reduction[rr_grid <= 1]),
    mde_increase_adjusted = mde(rr_grid[rr_grid >= 1], adjusted_increase[rr_grid >= 1])), by = by]
  out <- merge(out, null[, c(by, "fp_rate"), with = FALSE], by = by)
  out[, interpretable := fp_rate <= fp_limit]
  list(mde = out, power = pw)
}

#' County inputs to the heterogeneity tests, cumulative through event year `s_max`, from
#' the un-noised predictive draws (latent means): mu = predicted untreated count (posterior or
#' bootstrap mean), var_mu = variance of the cumulative latent draws (posterior variance for
#' the INLA models; bootstrap / ensemble variance for TWFE, ETWFE and GB-Ens), y = observed
#' count. One row per county observed at s_max (rows are contiguous event years 0..s).
heterogeneity_inputs <- function(fits, s_max = 4L) {
  draws <- grep("^y0_post_", names(fits), value = TRUE)
  f <- fits[s <= s_max]
  f <- f[f[, .I[max(s) == s_max], by = .(dgp, nc, model, r, i)]$V1]
  cum <- f[, lapply(.SD, sum), by = .(dgp, nc, model, r, i), .SDcols = c("y", "y0_hat", draws)]
  D <- as.matrix(cum[, ..draws])
  out <- cum[, .(dgp, nc, model, r, i, y, mu = y0_hat)]
  out[, var_mu := matrixStats::rowVars(D)]
  out[is.finite(mu) & mu > 0]
}

#' Between-county heterogeneity statistics per replication (cumulative counts O_i, predicted
#' mu_i, var_mu_i). r = sum(O) / sum(mu) is the common ratio. Primary: Dean-type score test
#' for extra-Poisson variation of the form Var(O) = mu + tau mu^2,
#'   T = sum[(O - r mu)^2 - r mu - r^2 var_mu] / sqrt(2 sum r^2 mu^2),
#' which weights counties by mu^2. Secondary: Pearson Q* = Q / (n - 1),
#'   Q = sum (O - r mu)^2 / (r mu + r^2 var_mu).
#' Diagnostics: counties with mu < 1, smallest denominator of Q, and the share of the T
#' numerator and of the squared residuals from the top 10% of counties by mu.
heterogeneity_stats <- function(h, O = "y", by = c("specification", "nc", "dgp", "model", "r")) {
  by <- intersect(by, names(h))
  dean <- function(o, mu, v) {
    rr <- sum(o) / sum(mu); num <- (o - rr * mu)^2 - rr * mu - rr^2 * v
    list(T = sum(num) / sqrt(2 * sum(rr^2 * mu^2)), num = num)
  }
  outcome_col <- O   # a column named like the argument would shadow it inside h[...]
  h[, {
    o <- get(outcome_col); rr <- sum(o) / sum(mu); res <- o - rr * mu
    d <- dean(o, mu, var_mu); num <- d$num
    # T(-1): drop the county with the largest numerator term, recompute r and T on the rest.
    k <- which.max(num)
    top <- mu >= stats::quantile(mu, 0.9)
    .(n = .N, n_mu_lt1 = sum(mu < 1), ratio = rr,
      T = d$T,
      T_minus1 = dean(o[-k], mu[-k], var_mu[-k])$T,
      Qstar = sum(res^2 / (rr * mu + rr^2 * var_mu)) / (.N - 1),
      min_denominator = min(rr * mu + rr^2 * var_mu),
      top10_share_T_numerator = sum(num[top]) / sum(num),
      top10_share_sq_residuals = sum(res[top]^2) / sum(res^2))
  }, by = by]
}

#' Placebo calibration: replications split in half with a fixed seed; the 95th percentile of
#' the null statistic in one half is the critical value, its rejection rate in the other half
#' the empirical size; both directions, averaged. Returns the two critical values per group.
heterogeneity_calibration <- function(null, stat = "T", seed = 2032027L, level = 0.95,
  by = c("specification", "nc", "dgp", "model")) {
  by <- intersect(by, names(null))
  reps <- sort(unique(null$r))
  set.seed(seed)
  half_a <- sample(reps, length(reps) %/% 2)
  null[, {
    x <- get(stat); a <- r %in% half_a
    crit_a <- stats::quantile(x[a], level, names = FALSE); crit_b <- stats::quantile(x[!a], level, names = FALSE)
    size_ab <- mean(x[!a] > crit_a); size_ba <- mean(x[a] > crit_b); p <- (size_ab + size_ba) / 2
    .(statistic = stat, crit_a = crit_a, crit_b = crit_b, size_a_to_b = size_ab, size_b_to_a = size_ba,
      size = p, size_se = sqrt(p * (1 - p) / .N), n_reps = .N,
      null_median = stats::median(x), null_p95 = stats::quantile(x, 0.95, names = FALSE))
  }, by = by]
}

#' Rejection rates of an injected-scenario statistic at the placebo critical values (both
#' split directions, averaged): the power of the calibrated test. Unpaired: every injected
#' replication is compared with the E0 critical values.
heterogeneity_power <- function(stats, cal, by = c("specification", "nc", "dgp", "model")) {
  by <- intersect(by, names(stats))
  m <- merge(stats, cal[, c(by, "statistic", "crit_a", "crit_b"), with = FALSE], by = by, allow.cartesian = TRUE)
  m[, x := fifelse(statistic == "T", T, fifelse(statistic == "T_minus1", T_minus1, Qstar))]
  m[, .(power = (mean(x > crit_a) + mean(x > crit_b)) / 2, n_reps = .N), by = c(by, "statistic")][,
    power_se := sqrt(power * (1 - power) / n_reps)][]
}

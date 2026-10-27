# Prior settings and backdated fits for supplement S7 (analysis/S7_priors/; rules in s4_decision.yml).
# Requires the core helpers in R/ (set_y0, fit_blind, impute_inla, fit_nested,
# aggregate_placebo, crps_rows, crps_sample) and data.table.

#' The 10 settings of the nested RW-IID model: the 3 x 3 grid of PC-prior upper
#' bounds (U_dev on the deviation walks, U_cty on the county effects) and the
#' structural-only point (current priors with the nested structure).
s4_settings <- function() {
  grid <- CJ(U_dev = c(1, 0.2, 0.05), U_cty = c(1, 0.3, 0.1), sorted = FALSE)
  grid[, structural_only := FALSE]
  z <- rbind(grid, data.table(U_dev = NA_real_, U_cty = NA_real_, structural_only = TRUE))
  z[, setting := fifelse(structural_only, "structural", sprintf("dev%s_cty%s", U_dev, U_cty))]
  # Tie order: tighter prior first (smaller U_dev, then smaller U_cty); structural-only loosest.
  z[, tie_order := frank(list(structural_only, U_dev, U_cty), ties.method = "first", na.last = TRUE)]
  z[]
}

#' Backdate one dataset: mask y0 for years mask_year - h ... mask_year - 1 in every
#' county with a finite mask_year, keeping the cohort label unchanged. Marks the
#' held-out cells (observed y0 before masking) and the scored counties (at least 2
#' observed years before the window and at least 1 held-out observed year). Unscored
#' counties stay in the fit. Returns a modified copy.
backdate_mask <- function(dt, h = 3L) {
  z <- copy(dt)
  window <- is.finite(z$mask_year) & z$year >= z$mask_year - h & z$year < z$mask_year
  z[, heldout := window & !is.na(y0)]
  z[window, y0 := NA]
  z[, scored := sum(year < mask_year - h & !is.na(y0)) >= 2 & any(heldout), by = i]
  z
}

#' Backdated scores of one fit: Poisson noise on the latent draws of the held-out
#' cells, summed over each scored county's held-out years. Aggregate: CRPS of the
#' noncore mean cumulative count; county: mean county CRPS. Counties with a missing
#' prediction (zero-event cohort) are left out and counted.
backdated_scores <- function(z, draws, seed) {
  ho <- z[heldout & scored]
  setorder(ho, i, year)
  D <- as.matrix(ho[, ..draws])
  set.seed(seed)
  D[] <- rpois(length(D), D)
  Dc <- rowsum(D, ho$i, reorder = FALSE)
  obs <- rowsum(ho$y, ho$i, reorder = FALSE)[, 1]
  noncore <- ho[, .(nc = urb[1] == "6.NonCore"), by = i]$nc
  ok <- rowSums(!is.finite(Dc)) == 0
  list(bd_agg_crps = crps_sample(colMeans(Dc[ok & noncore, , drop = FALSE]), mean(obs[ok & noncore])),
    bd_county_crps = mean(crps_rows(Dc[ok, , drop = FALSE], obs[ok])),
    bd_counties = sum(ok), bd_noncore = sum(ok & noncore), bd_missing = sum(!ok))
}

#' One placebo dataset (scenario k, replication r): main and backdated fits of each
#' setting, plus the production RW-IID predictions (`old`, post rows with y0_hat and
#' draws) as a reference. Post-period rows are event years 0-4 before actual
#' adoption, kept where every model predicts. Aggregate and county ledgers come
#' from aggregate_placebo() per half, with the same noise seed for every model
#' (common random numbers). Returns the ledgers by half, the backdated scores and
#' fit diagnostics.
s4_run_one <- function(panel, a, k, r, dgp, settings, split, old, covariates = NULL,
                       h = 3L, nboot = 200L, seed = 203L + (k - 1L) * 10000L + r) {
  draws <- paste0("y0_post_", seq_len(nboot))
  dt <- copy(panel)
  set_y0(dt, A_sim = a)
  setorder(dt, i, year)
  post <- dt$year >= dt$start_year_sim & dt$year <= dt$start_year_sim + 4 & dt$year < dt$start_year
  keys <- dt[post, .(i, year, urb, y, n, s = year - start_year_sim)]
  bd <- backdate_mask(dt, h)
  preds <- list()
  backdated <- list()
  diagnostics <- list()
  for (j in seq_len(nrow(settings))) {
    st <- settings[j]
    env <- new.env()
    fit_fun <- function(x) {
      m <- fit_nested(x, U_dev = st$U_dev, U_cty = st$U_cty, structural_only = st$structural_only,
        covariates = covariates)
      env$hyper <- m$summary.hyperpar[, "mean"]
      m
    }
    for (version in c("main", "backdated")) {
      x <- if (version == "main") dt else bd
      t0 <- Sys.time()
      set.seed(seed + j)
      res <- tryCatch({fit_blind(x, impute_inla, fit_fun = fit_fun, nboot = nboot); NULL},
        error = function(e) conditionMessage(e))
      diagnostics[[length(diagnostics) + 1L]] <- data.table(setting = st$setting, version = version,
        seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), error = if (is.null(res)) NA_character_ else res,
        prec_dev = if (is.null(res)) env$hyper[1] else NA_real_, prec_cty = if (is.null(res)) env$hyper[2] else NA_real_)
      if (!is.null(res)) next
      if (version == "main") preds[[st$setting]] <- x[post, c("y0_hat", draws), with = FALSE]
      else backdated[[st$setting]] <- c(list(setting = st$setting), backdated_scores(x, draws, seed))
    }
  }
  # Production RW-IID (no year effect), matched on county and year.
  if (!is.null(old)) preds[["old_rw_iid"]] <- old[keys, on = .(i, year), c("y0_hat", draws), with = FALSE]
  # Common sample: rows with a prediction from every model that fitted.
  ok <- Reduce(`&`, lapply(preds, function(p) is.finite(p$y0_hat)))
  half <- split[keys, on = .(id = i), half]
  stopifnot(!anyNA(half[ok]))
  ledgers <- lapply(c(A = "A", B = "B"), function(hh) {
    rows <- ok & half == hh
    out <- lapply(names(preds), function(model) {
      f <- cbind(keys[rows], preds[[model]][rows])
      f[, `:=`(model = model, dgp = dgp, r = r, nc = ordered(urb != "6.NonCore", c(FALSE, TRUE), c("NonCore", "Core")))]
      aggregate_placebo(f, seed = seed)
    })
    list(aggregate = rbindlist(lapply(out, `[[`, "aggregate")), county = rbindlist(lapply(out, `[[`, "county")))
  })
  backdated <- rbindlist(backdated)
  backdated[, `:=`(dgp = dgp, r = r)]
  diagnostics <- rbindlist(diagnostics)
  diagnostics[, `:=`(dgp = dgp, r = r, post_rows = sum(post), common_rows = sum(ok))]
  list(A = ledgers$A, B = ledgers$B, backdated = backdated, diagnostics = diagnostics)
}

#' Apply the selection rules to the backdated scores of each replication. R1: lowest
#' backdated aggregate CRPS; R2: lowest backdated county CRPS; R3: lowest backdated
#' county CRPS among settings within 5% of the best aggregate CRPS. Ties go to the
#' tighter prior (tie_order).
s4_select <- function(backdated, settings = s4_settings()) {
  z <- merge(backdated[is.finite(bd_agg_crps) & is.finite(bd_county_crps)],
    settings[, .(setting, tie_order)], by = "setting")
  z[, {
    near <- bd_agg_crps <= 1.05 * min(bd_agg_crps)
    .(rule = c("R1", "R2", "R3"), selected = c(
      setting[order(bd_agg_crps, tie_order)[1]],
      setting[order(bd_county_crps, tie_order)[1]],
      setting[near][order(bd_county_crps[near], tie_order[near])[1]]))
  }, by = .(specification, dgp, r)]
}

#' Add each rule's ledger (its selected setting's rows in each replication) to a ledger.
s4_rule_ledger <- function(L, selections) {
  z <- merge(selections, L, by.x = c("specification", "dgp", "r", "selected"),
    by.y = c("specification", "dgp", "r", "model"))
  z[, model := rule][, rule := NULL]
  rbindlist(list(L, z), use.names = TRUE, fill = TRUE)
}

#' Ratio of a candidate's metric to a reference's, pooled over scenarios (mean over
#' scenarios of the per-scenario mean over replications; for RMSE the root of the
#' pooled MSE ratio), with a delta-method SE paired by replication.
s4_pooled_ratio <- function(x, b, dgp, rmse = FALSE) {
  mx <- mean(tapply(x, dgp, mean)); mb <- mean(tapply(b, dgp, mean))
  ratio <- if (rmse) sqrt(mx / mb) else mx / mb
  infl <- if (rmse) ratio / 2 * (x / mx - b / mb) else (x - ratio * b) / mb
  se <- sqrt(sum(tapply(infl, dgp, function(v) var(v) / length(v)))) / length(unique(dgp))
  list(ratio = ratio, se = se, upper = ratio + 1.96 * se, n = length(x))
}

#' The backdating procedure in one dataset, as an analyst would run it: backdated fits
#' of every setting, the selection rule on their backdated scores, then the main fit of
#' the selected setting only. Predictions are attached to dt by reference (fit_blind);
#' returns the backdated scores with the selected setting flagged.
s4_procedure <- function(dt, settings = s4_settings(), rule = "R3", covariates = NULL,
                         h = 3L, nboot = 200L, seed = 1L) {
  draws <- paste0("y0_post_", seq_len(nboot))
  nested <- function(st) function(x) fit_nested(x, U_dev = st$U_dev, U_cty = st$U_cty,
    structural_only = st$structural_only, covariates = covariates)
  bd <- backdate_mask(dt, h)
  scores <- rbindlist(lapply(seq_len(nrow(settings)), function(j) {
    set.seed(seed + j)
    ok <- tryCatch({fit_blind(bd, impute_inla, fit_fun = nested(settings[j]), nboot = nboot); TRUE},
      error = function(e) FALSE)
    if (ok) c(list(setting = settings$setting[j]), backdated_scores(bd, draws, seed))
  }))
  scores[, `:=`(specification = "", dgp = "", r = 0L)]
  chosen <- rule
  pick <- s4_select(scores, settings)[rule == chosen, selected]
  set.seed(seed)
  fit_blind(dt, impute_inla, fit_fun = nested(settings[setting == pick]), nboot = nboot)
  scores[, `:=`(specification = NULL, dgp = NULL, r = NULL, selected = setting == pick)][]
}

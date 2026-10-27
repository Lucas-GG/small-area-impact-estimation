# Functions for the separate backdated diagnostic; no main fits are run here.
find_backdate_root <- function(start = getwd()) {
  path <- normalizePath(start)
  repeat {
    if (file.exists(file.path(path, "R", "fe.r")) &&
        file.exists(file.path(path, "data", "simulation", "simulated_assignments_with_meta.rds"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Cannot find repository from ", start)
    path <- parent
  }
}

backdate_assignment <- function(a, h) {
  stopifnot(h >= 1L, h == as.integer(h), !is.null(names(a)), !anyNA(a))
  a[is.finite(a)] <- a[is.finite(a)] - h
  a
}

backdate_panel <- function(panel, a, h) {
  # Start from raw columns, never predictions or previously generated features.
  z <- copy(panel)
  set_y0(z, A_sim = backdate_assignment(a, h))
  masked <- is.na(z$y0)
  stopifnot(all(is.na(z$y0[z$year >= z$start_year_sim])))
  # The original y is kept outside model inputs for truth summaries only.
  z
}

backdate_eligibility <- function(panel, a, h, sim, scenario, seed) {
  units <- data.table(id = as.integer(names(a)), x = unname(a))[is.finite(x)]
  z <- merge(units, panel[, .(id = i, year, y, start_year)], by = "id",
             allow.cartesian = TRUE)
  out <- z[, .(
    n_pre = sum(year < x - h & year < start_year & is.finite(y)),
    positive_pre = sum(y[year < x - h & year < start_year], na.rm = TRUE) > 0,
    n_ho = sum(year >= x - h & year < x & year < start_year & is.finite(y)),
    n_post = sum(year >= x & year <= x + 4 & year < start_year & is.finite(y)),
    y_ho = if (sum(year >= x - h & year < x & year < start_year & is.finite(y)) == h)
      as.numeric(sum(y[year >= x - h & year < x])) else NA_real_
  ), by = .(id, x)]
  out[, `:=`(sim = sim, scenario = scenario, seed = seed, backdate_years = h,
              fake_start = x - h)]
  out
}

with_backdate_seed <- function(seed, expr) {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  force(expr)
}

check_backdate_features <- function(z, seed = 203L, predictors = NULL) {
  # Perturb unavailable truth: no derived feature, label or early-stop mask may change.
  feature_inputs <- function(dt) {
    build_y0_features(dt)
    dt[, c("y0", get_y0_features(dt)), with = FALSE]
  }
  changed <- copy(z)
  changed[is.na(y0), y := y + 100000L]
  left <- copy(z); right <- copy(changed)
  stopifnot(isTRUE(all.equal(feature_inputs(left), feature_inputs(right))))
  # Reuse the usual hazard and bootstrap feature/validation construction.
  left[, cpr := get_hazard(left, engine = "xgboost", predictors = predictors)]
  right[, cpr := get_hazard(right, engine = "xgboost", predictors = predictors)]
  stopifnot(isTRUE(all.equal(left$cpr, right$cpr)))
  p1 <- with_backdate_seed(seed, pre_xgb(left, predictors = predictors))
  p2 <- with_backdate_seed(seed, pre_xgb(right, predictors = predictors))
  stopifnot(isTRUE(all.equal(p1, p2)), any(p1$mask & !is.na(p1$data$y0)),
            any(!p1$mask & !is.na(p1$data$y0)))
  invisible(TRUE)
}

# Covariate specifications match 2c exactly; the default preserves the no-covariate branch.
backdate_covariate_panel <- function(original) {
  z <- copy(as.data.table(original))
  setnames(z, tolower(names(z)))
  z[, intersect(c("y0", "y0_25"), names(z)) := NULL]
  z[, nc := ordered(urb != "6.NonCore", levels = c(FALSE, TRUE), labels = c("NonCore", "Core"))]
  predictors <- c("uemp", "pov", "mhinc", "aian")
  stopifnot(all(predictors %in% names(z)))
  z[, (predictors) := lapply(.SD, scale), .SDcols = predictors]
  z
}

backdate_specs <- function(with_covariates = FALSE) {
  if (with_covariates) return(list(
    TWFE = list(fun = impute_fe, args = list(fm = y0 ~ uemp + pov + mhinc + aian | i + year)),
    ETWFE = list(fun = impute_fe, args = list(fm = y0 ~ uemp + pov + mhinc + aian | cell + year)),
    # Nested RW-IID with the default priors, as in the main results (S4 notes).
    `RW-IID` = list(fun = impute_inla, args = list(fit_fun = function(x) fit_nested(x, structural_only = TRUE,
      covariates = c("uemp", "pov", "mhinc", "aian")))),
    `GB-Ens` = list(fun = impute_xgb, args = list(predictors = c("uemp", "pov", "mhinc", "aian")))))
  list(
    TWFE = list(fun = impute_fe, args = list(fm = y0 ~ 1 | i + year)),
    ETWFE = list(fun = impute_fe, args = list(fm = y0 ~ 1 | cell + year)),
    # Nested RW-IID with default priors (S4 decision without covariates).
    `RW-IID` = list(fun = impute_inla, args = list(fit_fun = function(x) fit_nested(x, structural_only = TRUE))),
    `GB-Ens` = list(fun = impute_xgb, args = list()))
}

backdate_predict <- function(z, spec, nboot) {
  # FE estimation implicitly drops NA response rows; predictions retain all rows.
  do.call(fit_blind, c(list(dt = z, .impute = spec$fun, nboot = nboot), spec$args))
  z
}

summarize_backdate <- function(fitted, eligibility, model, nboot, ppp_seed) {
  draw_cols <- paste0("y0_post_", seq_len(nboot))
  stopifnot(all(c("y0_hat", draw_cols) %in% names(fitted)))
  base <- copy(eligibility)
  value <- lapply(seq_len(nrow(base)), function(j) {
    e <- base[j]
    rows <- fitted[i == e$id & year >= e$fake_start & year < e$x]
    result <- list(mu_ho_mean = NA_real_, mu_ho_var = NA_real_, ppp_ho = NA_real_, ppp_ho_upper = NA_real_,
                   ppp_ho_lower = NA_real_, tail_convention = "mid-p",
                   n_valid_draws = 0L, heldout_status = "incomplete_heldout")
    if (e$n_ho != e$backdate_years || nrow(rows) != e$backdate_years) return(result)
    if (e$n_pre < 2L) {result$heldout_status <- "insufficient_history"; return(result)}
    if (any(!is.finite(rows$y0_hat) | rows$y0_hat < 0)) {
      result$heldout_status <- "nonestimable_point_prediction"; return(result)
    }
    result$mu_ho_mean <- sum(rows$y0_hat)
    draws <- as.matrix(rows[, ..draw_cols])
    valid <- colSums(!is.finite(draws) | draws < 0) == 0L
    result$n_valid_draws <- as.integer(sum(valid))
    # Do not silently condition diagnostics on only the successful bootstrap draws.
    if (!all(valid) || nboot < 2L) {
      result$heldout_status <- "incomplete_draws"; return(result)
    }
    lambda <- colSums(draws)
    result$mu_ho_var <- var(lambda)
    predictive_counts <- with_backdate_seed(ppp_seed + j,
      rpois(length(lambda), lambda))
    result$ppp_ho_upper <- mean(predictive_counts > e$y_ho) + 0.5 * mean(predictive_counts == e$y_ho)
    result$ppp_ho_lower <- mean(predictive_counts < e$y_ho) + 0.5 * mean(predictive_counts == e$y_ho)
    result$ppp_ho <- result$ppp_ho_upper # Preserve the original upper-tail field.
    result$heldout_status <- "ok"
    result
  })
  cbind(base, rbindlist(value))[, model := model][]
}

run_backdate_one <- function(panel, a, h, sim, scenario, seed, model, spec, nboot,
                             audit_features = TRUE) {
  e <- backdate_eligibility(panel, a, h, sim, scenario, seed)
  z <- backdate_panel(panel, a, h)
  warnings <- character()
  result <- tryCatch(withCallingHandlers(with_backdate_seed(seed, {
    if (model == "GB-Ens" && audit_features) check_backdate_features(z, seed, predictors = spec$args$predictors)
    # Feature auditing must not consume the random stream of the actual fit.
    set.seed(seed)
    fitted <- backdate_predict(z, spec, nboot)
    summarize_backdate(fitted, e, model, nboot, seed + 1000000L)
  }), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(err) {
    e[, `:=`(model = model, mu_ho_mean = NA_real_, mu_ho_var = NA_real_,
              ppp_ho = NA_real_, ppp_ho_upper = NA_real_, ppp_ho_lower = NA_real_,
              tail_convention = "mid-p", n_valid_draws = 0L, heldout_status = "fit_failed",
              fit_error = conditionMessage(err))]
    e
  })
  if (!"fit_error" %in% names(result)) result[, fit_error := NA_character_]
  result[, fit_warnings := paste(unique(warnings), collapse = " | ")]
  result
}

read_main_summary <- function(path, panel, assignments, scenario_index, model) {
  saved <- readRDS(path)
  cols <- c("r", "i", "year", "start_year_sim", "y", "n", "y0_hat")
  stopifnot(all(cols %in% names(saved)))
  main <- as.data.table(saved)[, ..cols]
  rm(saved); gc()
  stopifnot(!anyDuplicated(main[, .(r, i, year)]),
            setequal(unique(main$r), seq_len(dim(assignments)[2])))
  # Check assignment dates and observed input values for every saved prediction row.
  ii <- match(as.character(main$i), dimnames(assignments)[[1]])
  expected_x <- assignments[cbind(ii, main$r, rep(scenario_index, nrow(main)))]
  stopifnot(!anyNA(ii), isTRUE(all.equal(main$start_year_sim, expected_x,
                                      check.attributes = FALSE)))
  source_row <- panel[main, on = .(i, year), .(y, n)]
  stopifnot(isTRUE(all.equal(main$y, source_row$y, check.attributes = FALSE)),
            isTRUE(all.equal(main$n, source_row$n, check.attributes = FALSE)))
  main[, x := start_year_sim]
  main <- main[is.finite(x)]
  main[, {
    post <- year >= x[1L] & year <= x[1L] + 4
    pre <- year < x[1L]
    pre_valid <- pre & is.finite(y) & is.finite(y0_hat)
    pearson <- pre_valid & y0_hat > 0
    n_post <- sum(post & is.finite(y))
    n_pred <- sum(post & is.finite(y0_hat) & y0_hat >= 0)
    list(y_post = if (n_post == 5L) as.numeric(sum(y[post])) else NA_real_,
         mu_post_mean = if (n_post == 5L && n_pred == 5L) sum(y0_hat[post]) else NA_real_,
         main_n_post = as.integer(n_post), main_n_post_pred = as.integer(n_pred),
         pre_rmse = if (any(pre_valid)) sqrt(mean((y[pre_valid] - y0_hat[pre_valid])^2)) else NA_real_,
         pre_rms_pearson = if (any(pearson))
           sqrt(mean((y[pearson] - y0_hat[pearson])^2 / y0_hat[pearson])) else NA_real_,
         pre_mu_mean = if (any(pre_valid)) mean(y0_hat[pre_valid]) else NA_real_,
         pre_n = sum(pre_valid), pre_pearson_n = sum(pearson))
  }, by = .(sim = r, id = i, x)][, model := model][]
}

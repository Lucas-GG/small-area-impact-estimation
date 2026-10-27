
library(INLA)



#' Fit an INLA model and attach posterior predictive draws and summaries to the data.table.
#' Clears any prior y0_post_* and y0_* summary columns to avoid stale values, then refits and imputes.
impute_inla <- \(.dt, fit_fun = fit_ar, nboot = 200) {
  # drop posterior draws
  drop_by_pattern(.dt, "^y0_post_")
  # drop derived summaries
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")

  # A cohort with no observed events has no finite effect (fixest drops such
  # groups in ETWFE). Leave it out of the fit and report its rows as NA.
  events <- .dt[, .(events = sum(y0, na.rm = TRUE)), by = cohort]
  dead <- .dt$cohort %in% events[events == 0, cohort]
  y0_saved <- .dt$y0
  if (any(dead)) .dt[dead, y0 := NA]
  # Cohort x group cells for fit_nested(), from the training rows (a cohort that is dead
  # overall stays out of the fit, as before).
  if ("ncf" %in% names(.dt)) {
    cell <- cohort_cell(.dt$cohort, .dt$ncf, .dt$y0)
    .dt[, cell := as.vector(cell)]
    data.table::setattr(.dt, "cell_merges", attr(cell, "merged"))
  }

  m <-  fit_fun(.dt)
  inference_inla(.dt, m, nboot)

  if (any(dead)) {
    .dt[, y0 := y0_saved]
    pred_cols <- grep("^y0_(post_\\d+|hat)$", names(.dt), value = TRUE)
    .dt[dead, (pred_cols) := NA_real_]
  }
  .dt
}

#' Generate posterior predictive draws for the latent predictor and attach them to the data.table.
#' Adds y0_post_1:nboot (exp(Predictor)) and a point summary y0_hat from the fitted INLA object.
inference_inla <- \(.dt, .m, nboot = 200) {
  # Seed INLA's sampler from R's RNG so set.seed() makes the draws reproducible
  # (the default seed = 0L draws a random seed).
  y0_post <- inla.posterior.sample(nboot, .m, seed = as.integer(runif(1) * .Machine$integer.max))
  y0_post <- inla.posterior.sample.eval("Predictor", y0_post) |>
    remove_names() |>
    exp()

  y0_post_cols <- paste0("y0_post_", 1:nboot)
  .dt[, (y0_post_cols) := as.data.frame(y0_post)]

  y0_hat   <- .m$summary.fitted.values$mean
  .dt[, y0_hat := y0_hat]

  .dt
}

#' Default INLA compute options used across model fits (fit criteria + posterior sampling support).
#' Enables DIC/CPO/WAIC and predictor configuration needed for posterior sampling of the linear predictor.
ctrl_comp <- list(
  dic = TRUE, cpo = TRUE, waic = TRUE
  , config = TRUE, return.marginals.predictor = TRUE
)

#' Pattern-mixture Poisson model with cohort fixed effects and cohort-specific RW1 time trends.
#' Uses replicate=cid RW1 to allow each start-year cohort its own smooth deviation over calendar time.
fit_pm <- \(.dt) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}


fit_pm_c <- \(.dt) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  .dt$c4 <- .dt$c3 <- .dt$c2 <- .dt$c1 <- .dt$cid
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$t4 <- .dt$t3 <- .dt$t2 <- .dt$t1 <- .dt$time

  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + uemp + pov + mhinc + aian
    + f(c1, uemp, model = "iid", constr = TRUE)
    + f(c2, pov, model = "iid", constr = TRUE)
    + f(c3, mhinc, model = "iid", constr = TRUE)
    + f(c4, aian, model = "iid", constr = TRUE)
    + f(t1, uemp, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t2, pov, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t3, mhinc, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t4, aian, model = "rw1", constr = TRUE, scale.model = TRUE)
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}


#' Pattern-mixture model + unit random intercept (iid) for extra between-county heterogeneity.
#' Adds f(i, iid) on top of cohort FE and cohort-specific RW1 trends.
fit_me <- \(.dt) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()

  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + f(i, model = "iid")
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}


#' Pattern-mixture model + unit random intercept (iid) for extra between-county heterogeneity.
#' Adds f(i, iid) on top of cohort FE and cohort-specific RW1 trends.
fit_me_c <- \(.dt) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  .dt$c4 <- .dt$c3 <- .dt$c2 <- .dt$c1 <- .dt$cid
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$t4 <- .dt$t3 <- .dt$t2 <- .dt$t1 <- .dt$time

  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + f(i, model = "iid")
    + uemp + pov + mhinc + aian
    + f(c1, uemp, model = "iid", constr = TRUE)
    + f(c2, pov, model = "iid", constr = TRUE)
    + f(c3, mhinc, model = "iid", constr = TRUE)
    + f(c4, aian, model = "iid", constr = TRUE)
    + f(t1, uemp, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t2, pov, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t3, mhinc, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t4, aian, model = "rw1", constr = TRUE, scale.model = TRUE)
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}

#' Pattern-mixture model + spatial BYM2 random effect for county-level structured/unstructured variation.
#' Adds f(i, bym2, graph=W) to borrow strength across neighboring areas.
fit_sp <- \(.dt, graph = W) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()

  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + f(i, model = "bym2", graph = graph)
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}

#' Spatial pattern-mixture model with additional observation-level iid noise.
#' Adds f(id, iid) to absorb extra-Poisson variability not captured by spatial/cohort-time structure.
fit_spo <- \(.dt, graph = W) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$id <- seq_len(nrow(.dt))

  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + f(i, model = "bym2", graph = graph)
    + f(id, model = "iid")
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}


#' Reduced pattern-mixture model with unit-specific RW1 trends (no spatial term).
#' Uses replicate=i RW1 to let each county have its own smooth temporal trajectory (random slope-like).
fit_ar <- \(.dt) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$time2   <- .dt$time
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time2, model = "rw1"
      , replicate = i
      , constr = TRUE
      , scale.model = TRUE
    )
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}


#' Combined model: cohort-specific RW1 + spatial BYM2 + unit-specific RW1 trends.
#' Captures cohort-level dynamics, spatial dependence, and idiosyncratic county time trends simultaneously.
fit_iii <- \(.dt, graph = W) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$time2   <- .dt$time
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + f(i, model = "bym2", graph = graph)
    + f(time2, model = "rw1"
      , replicate = i
      , constr = TRUE, scale.model = TRUE
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}

#' Spatiotemporal model with spatial BYM2 evolving over time (AR1 across time groups).
#' Lets spatial risk vary by year via group=time with AR1 dependence, plus cohort-time RW1 structure.
fit_sptiv <- \(.dt, graph = W) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + f(i, model = "bym2", graph = graph
      , group = time, control.group = list(model = "ar1")
    )
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}



#' Spatiotemporal model with spatial BYM2 evolving over time (AR1 across time groups).
#' Lets spatial risk vary by year via group=time with AR1 dependence, plus cohort-time RW1 structure.
fit_sptiv_c <- \(.dt, graph = W) {
  pc_c <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.1)))
  pc_t <- list(prec = list(prior = "pc.prec", param = c(0.03, 0.1)))
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  .dt$c4 <- .dt$c3 <- .dt$c2 <- .dt$c1 <- .dt$cid
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$t4 <- .dt$t3 <- .dt$t2 <- .dt$t1 <- .dt$time
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + f(i, model = "bym2", graph = graph
      , group = time, control.group = list(model = "ar1")
    )
    + uemp + pov + mhinc + aian
    + f(c1, uemp, model = "iid", constr = TRUE)
    + f(c2, pov, model = "iid", constr = TRUE)
    + f(c3, mhinc, model = "iid", constr = TRUE)
    + f(c4, aian, model = "iid", constr = TRUE)
    + f(t1, uemp, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t2, pov, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t3, mhinc, model = "rw1", constr = TRUE, scale.model = TRUE)
    + f(t4, aian, model = "rw1", constr = TRUE, scale.model = TRUE)
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}





#' Random-slope-style extension: spatial random intercept + spatial random slope over centered time.
#' Adds a second BYM2 term on (ii, time2) to allow spatially-structured differential trends.
fit_rs <- \(.dt, graph = W) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$cohort |> factor()
  # cid    : numeric cohort index for RW replication
  .dt$cid <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$ii <- .dt$i
  .dt$time2 <- as.integer(as.factor(.dt$time)) |> scale(scale = FALSE)
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cid
      , constr = TRUE       # impose sum-to-zero per RW
      , scale.model = TRUE   # optional, rescales for stability
      , hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    )
    + f(i, model = "bym2", graph = graph)
    + f(ii, time2, model = "bym2", graph = graph)
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = ctrl_comp
  )
  model
}
#' Nested hierarchical model: cohort effects, common year effects,
#' cohort deviation random walks for adopting cohorts only (never-treated counties carry
#' the common path), and county effects. Every prior is explicit. With tight deviation
#' and county priors the model approaches ETWFE; with a loose county prior, roughly TWFE.
#' `structural_only` uses the current priors (pc.prec(1, 0.1) walks, loggamma(1, 5e-5)
#' county effects). `covariates` enter as common fixed coefficients.
#' Cohort x group: with a `cell` column (impute_inla), the fixed effects
#' are per cohort x group cell (never-treated cells included) and, with walks = "cell", there is
#' one deviation walk per adopting cell; walks = "cohort" keeps one walk per adopting cohort
#' (the Table S9 intermediate). All walks share one precision; year effects are common.
fit_nested <- function(.dt, U_dev = 1, U_cty = 1, structural_only = FALSE, iid = TRUE,
  covariates = NULL, spatial = FALSE, graph = if (exists("W")) W else NULL, walks = c("cell", "cohort")) {
  walks <- match.arg(walks)
  cell <- if ("cell" %in% names(.dt)) .dt$cell else as.character(.dt$cohort)
  .dt$fcohort <- factor(cell)
  .dt$fyear <- factor(.dt$year)
  adopting <- is.finite(.dt$cohort)
  .dt$time_dev <- ifelse(adopting, as.numeric(factor(.dt$year)), NA)
  .dt$cid_dev <- 1
  walk_id <- if (walks == "cell") cell else as.character(.dt$cohort)
  .dt$cid_dev[adopting] <- as.numeric(factor(walk_id[adopting]))
  hyper_dev <- list(prec = list(prior = "pc.prec", param = if (structural_only) c(1, 0.1) else c(U_dev, 0.01)))
  hyper_cty <- list(prec = if (structural_only) list(prior = "loggamma", param = c(1, 5e-05)) else
    list(prior = "pc.prec", param = c(U_cty, 0.01)))
  hyper_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), phi = list(prior = "pc", param = c(0.5, 0.5)))
  hyper_ar1 <- list(rho = list(prior = "normal", param = c(0, 0.15)))
  terms <- c("offset(log(n)) - 1", "fcohort", "fyear", covariates,
    'f(time_dev, model = "rw1", replicate = cid_dev, constr = TRUE, scale.model = TRUE, hyper = hyper_dev)')
  if (iid) terms <- c(terms, 'f(i, model = "iid", hyper = hyper_cty)')
  if (spatial) {
    .dt$i_sp <- .dt$i
    .dt$time_sp <- as.numeric(factor(.dt$year))
    terms <- c(terms, 'f(i_sp, model = "bym2", graph = graph, hyper = hyper_bym2, group = time_sp, control.group = list(model = "ar1", hyper = hyper_ar1))')
  }
  fm <- stats::as.formula(paste("y0 ~", paste(terms, collapse = " + ")), env = environment())
  inla(fm, family = "poisson", data = .dt,
    control.predictor = list(link = 1), control.compute = ctrl_comp,
    control.fixed = list(mean = 0, prec = 0.001))
}

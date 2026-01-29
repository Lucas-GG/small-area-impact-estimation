
library(INLA)



#' Fit an INLA model and attach posterior predictive draws and summaries to the data.table.
#' Clears any prior y0_post_* and y0_* summary columns to avoid stale values, then refits and imputes.
impute_inla <- \(.dt, fit_fun = fit_ar, nb = 200) {
  # drop posterior draws
  drop_by_pattern(.dt, "^y0_post_")
  # drop derived summaries
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")

  m <-  fit_fun(.dt)
  inference_inla(.dt, m, nb)
}

#' Generate posterior predictive draws for the latent predictor and attach them to the data.table.
#' Adds y0_post_1:nb (exp(Predictor)) and a point summary y0_hat from the fitted INLA object.
inference_inla <- \(.dt, .m, nb = 200) {
  y0_post <- inla.posterior.sample(nb, .m)
  y0_post <- inla.posterior.sample.eval("Predictor", y0_post) |>
    remove_names() |>
    exp()

  y0_post_cols <- paste0("y0_post_", 1:nb)
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
#' Uses replicate=cohort RW1 to allow each start-year cohort its own smooth deviation over calendar time.
fit_pm <- \(.dt) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$start_year |> factor()
  # cohort : numeric cohort index for RW replication
  .dt$cohort <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cohort
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

#' Pattern-mixture model + unit random intercept (iid) for extra between-county heterogeneity.
#' Adds f(i, iid) on top of cohort FE and cohort-specific RW1 trends.
fit_me <- \(.dt) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$start_year |> factor()
  # cohort : numeric cohort index for RW replication
  .dt$cohort <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()

  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cohort
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

#' Pattern-mixture model + spatial BYM2 random effect for county-level structured/unstructured variation.
#' Adds f(i, bym2, graph=W) to borrow strength across neighboring areas.
fit_sp <- \(.dt, graph = W) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$start_year |> factor()
  # cohort : numeric cohort index for RW replication
  .dt$cohort <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()

  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cohort
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
  .dt$fcohort <- .dt$start_year |> factor()
  # cohort : numeric cohort index for RW replication
  .dt$cohort <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$id <- seq_len(nrow(.dt))

  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cohort
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
  .dt$fcohort <- .dt$start_year |> factor()
  # cohort : numeric cohort index for RW replication
  .dt$cohort <- .dt$fcohort |> factor() |> as.numeric()
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
  .dt$fcohort <- .dt$start_year |> factor()
  # cohort : numeric cohort index for RW replication
  .dt$cohort <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$time2   <- .dt$time
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cohort
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
  .dt$fcohort <- .dt$start_year |> factor()
  # cohort : numeric cohort index for RW replication
  .dt$cohort <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cohort
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

#' Random-slope-style extension: spatial random intercept + spatial random slope over centered time.
#' Adds a second BYM2 term on (ii, time2) to allow spatially-structured differential trends.
fit_rs <- \(.dt, graph = W) {
  # fcohort: cohort fixed effects
  .dt$fcohort <- .dt$start_year |> factor()
  # cohort : numeric cohort index for RW replication
  .dt$cohort <- .dt$fcohort |> factor() |> as.numeric()
  #time   : numeric time index for RW
  .dt$time   <- .dt$year |> factor() |> as.numeric()
  .dt$ii <- .dt$i
  .dt$time2 <- as.integer(as.factor(.dt$time)) |> scale(scale = FALSE)
  model <- inla(y0 ~ offset(log(n)) - 1
    + fcohort
    + f(time, model = "rw1"
      , replicate = cohort
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
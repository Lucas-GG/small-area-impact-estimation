
library(INLA)

fit_ime <- \(.dt) {
  model <- inla(y0 ~ offset(log(n)) - 1
    + f(i, model = "iid")
    + f(time, model = "rw1")
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = list(config = TRUE)
  )
  model
}

fit_ar <- \(.dt) {
  model <- inla(y0 ~ offset(log(n)) - 1
    + f(time, model = "rw1")
    + f(i, model = "iid", group = time, control.group = list(model = "ar1"))
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = list(config = TRUE)
  )
  model
}

fit_spt <- \(.dt) {
  model <- inla(y0 ~ offset(log(n)) - 1
    + f(time, model = "rw1")
    + f(i, model = "bym2", graph = graph
      , group = time, control.group = list(model = "ar1")
    )
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = list(config = TRUE)
  )
  model
}

impute_inla <- \(.dt, fit_fun = fit_ar, nb = 200) {
  result <- copy(.dt)
  result[, time := as.integer(as.factor(year))]
  m <-  fit_fun(result)
  result[, time := NULL]
  inference_inla(result, m, nb)
}

inference_inla <- \(.dt, .m, nb = 200) {
  y0_post <- inla.posterior.sample(nb, .m) %>%
    inla.posterior.sample.eval("Predictor", .) %>%
    remove_names %>%
    exp

  y0_post_cols <- paste0("y0_post_", 1:nb)
  .dt[, (y0_post_cols) := as.data.frame(y0_post)]

  y0_hat   <- .m$summary.fitted.values$mean

  .dt[, y0_hat := y0_hat]
  .dt
}
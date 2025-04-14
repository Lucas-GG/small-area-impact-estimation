
library(INLA)

fit_imex <- \(.dt) {
  model <- inla(y0 ~ offset(log(x)) - 1
    + f(i, model = "iid")
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = list(config = TRUE)
  )
  model
}

fit_arx <- \(.dt) {
  model <- inla(y0 ~ offset(log(x)) - 1
    + f(i, model = "iid", group = time, control.group = list(model = "ar1"))
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = list(config = TRUE)
  )
  model
}

fit_sptx <- \(.dt) {
  model <- inla(y0 ~ offset(log(x)) - 1
    + f(i, model = "besagproper2", graph = graph
      , group = time, control.group = list(model = "ar1")
    )
    , family = "poisson"
    , data = .dt
    , control.predictor = list(link = 1)
    , control.compute = list(config = TRUE)
  )
  model
}

#dt %>% mutate(x = y0_post_1, time = as.numeric(factor(year))) %>% fit_spt %>% summary
#dt %>% combx_inla(nb = 20, ncores = 1, fit_fun = fit_spt)%>% summary


#combine inla ran on diffrent draws of x
combx_inla <- function(.dt, nb = 20, ncores = 1, fit_fun = fit_arx) {
  result <- copy(.dt)
  # Create a temporary copy for modifications
  #y0_post_cols <- grep("y0_post", names(result), value = TRUE)
  y0_post_cols <- paste0("y0_post_", round(seq(1, 200, length.out  = nb)))

  # Process just the first two posterior columns as in your example
  .dt0 <- .dt[, .(i, year, n, y0)]
  .dt0[, time := as.numeric(factor(year))]

  m_list <- mclapply(y0_post_cols, \(b) {
    .dt0[, x := .dt[[b]]]
    fit_fun(.dt0)
  }, mc.cores = ncores)

  m <- inla.merge(m_list)
  m_list <- NULL
  m
}
#.m <- m
#.dt <- dt

inferece_inla_old <- \(.dt, .m, nb = 200) {

  y0_post <- inla.posterior.sample(nb, .m) %>%
    inla.posterior.sample.eval("Predictor", .) %>%
    remove_names %>%
    exp
  y0_post_cols <- paste0("y0_post_", 1:nb)
  .dt[, (y0_post_cols) := as.data.frame(y0_post)]


#  if (length(.m$marginals.linear.predictor) != 0) {
#    y0_hat   <- sapply(.m$marginals.linear.predictor
#      , \(m) inla.emarginal(exp, m)
#    )
#    y0_2   <- sapply(
#      .m$marginals.linear.predictor, \(m) inla.emarginal(\(x) exp(x)^2L, m)
#    )
#    y0_sd    <- sqrt(y0_2 - y0_hat^2L)
#  } else 
  if (is.null(.m$summary.fitted.values)) {
    y0_hat   <- rowMeans(y0_post)
    y0_sd    <- fsd(t(y0_post))
  } else {
    y0_hat   <- .m$summary.fitted.values$mean
    y0_sd    <- .m$summary.fitted.values$sd
  }

  eta_mean <- .m$summary.linear.predictor$mean
  eta_sd   <- .m$summary.linear.predictor$sd
  y0_pi    <- with(.dt, pi_poilog(y, eta_mean, eta_sd))
  ci       <- with(.dt, ci_poilog(eta_mean, eta_sd))

  .dt[, ':='(
    y0_hat = y0_hat
    , y0_sd  = y0_sd
    , y0_pi   = y0_pi
    , y0_lower = ci[, 1]
    , y0_upper = ci[, 2]
  )]
  .dt

}


ar_inference_old <- function(.dt, mcores = 1) {
  result <- copy(.dt)
  # Create a temporary copy for modifications
  y0_post_cols <- grep("y0_post", names(result), value = TRUE)

  # Process just the first two posterior columns as in your example
  .dt0 <- .dt[, .(i, year, n, y0)]
  .dt0[, time := as.numeric(factor(year))]

  y0_post <- mclapply(y0_post_cols, \(b) {
    .dt0[, x := .dt[[b]]]
    # Apply cluster_boot
    #.dt0 <- cluster_boot(.dt0)

    # Get predictions
    predicted_values <- impute_ar(.dt0)
    predicted_values
  }, mc.cores = mcores)

  lapply(seq_along(y0_post), \(i) {
    set(result, j = y0_post_cols[i], value = y0_post[[i]])
  })

  result
}
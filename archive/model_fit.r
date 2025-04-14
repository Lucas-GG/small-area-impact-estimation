#***------------------------------------------------------------***
#***fit functions***
#***------------------------------------------------------------***

library(fixest)


#.data <- NULL
#.data <- filter(placebo_dt, replication == 1)
#.data %>% impute_etwfe
#.data %>% fit_fe(3) %>% score
boot_fe <- \(.data, impute, nboot = 200) {
  lapply(1:nboot, \(i) {
    .data$wts <- NULL
    .wdata <- left_join(
      .data
      , data.frame(i = unique(.data$i), wts = bb(n_distinct(.data$i)))
      , by = "i"
    )
    .wdata %>% mutate(wts = .wdata$wts / sum(.wdata$wts) * nrow(.wdata)) %>%
      impute
  }) %>%
    do.call("cbind", .)
}

impute_fe <- \(.data) {
  if (is.null(.data$wts)) .data$wts <- 1
  fepois(fm
    , vcov = ~ i, offset = ~ log(n)
    , .data, weights = .data$wts
  ) %>% predict(.data, type = "response")
}

fit_fe <- \(.data, nboot = 200) {
  t0 <- Sys.time()
  .data$y0_hat <- .data %>% impute_fe
  .data$y0_post <- .data %>% boot_fe(impute_fe, nboot)

  .data$y0_sd <- apply(.data$y0_post, 1, sd)
  .data$eta_mean <- log(.data$y0_hat)
  .data$eta_sd <- apply(log(.data$y0_post), 1, sd)
  .data$y0_pi <- with(.data, pi_poilog(Y, eta_mean, eta_sd))
  .data[, c("y0_lower", "y0_upper")] <-
    with(.data, ci_poilog(eta_mean, eta_sd))
  t1 <- Sys.time()
  print(t1 - t0)
  .data
}

#------------------------------------------------------------------#
#ME ---------------------------------------------------------------#
library(lme4)

#fm <- formula(Y0 ~ offset(log(n)) + x + (x - 1 | id) + (1 | fyr))
fm <- formula(y0 ~ offset(log(n)) + scale(year) + (scale(year) | i))

boot_me <- \(.data, nboot = 200) {
  lapply(1:nboot, \(i) {
    .data$wts <- NULL
    .wdata <- left_join(
      .data
      , data.frame(i = unique(.data$i), wts = bb(n_distinct(.data$i)))
      , by = "i"
    )
    .wdata %>% mutate(wts = .wdata$wts / sum(.wdata$wts) * nrow(.wdata)) %>%
      impute_me
  }) %>%
    do.call("cbind", .)
}

#.data <- filter(placebo_dt, replication == 1)
#.data$wts <- 1
impute_me <- \(.data) {
  if (is.null(.data$wts)) .data$wts <- 1
  glmer(fm
    , family = "poisson"
    , .data, weights = wts
    , control = glmerControl(optimizer = "bobyqa")
  ) %>% predict(.data, type = "response")
}

fit_me <- \(.data, nboot = 200) {
  t0 <- Sys.time()
  .data$y0_hat <- .data %>% impute_me
  .data$y0_post <- .data %>% boot_me(nboot)
  .data$y0_sd <- apply(.data$y0_post, 1, sd)
  .data$eta_mean <- log(.data$y0_hat)
  .data$eta_sd <- apply(log(.data$y0_post), 1, sd)
  .data$y0_pi <- with(.data, pi_poilog(Y, eta_mean, eta_sd))
  .data[, c("y0_lower", "y0_upper")] <-
    with(.data, ci_poilog(eta_mean, eta_sd))
  .data$method <- "me"
  t1 <- Sys.time()
  print(t1 - t0)
  .data
}

#filter(placebo_dt, replication == 1) %>% fit_me(2) %>% score
#filter(placebo_dt, replication == 2) %>% fit_me(2) %>% score
#.data %>% impute_me
#.data %>% boot_me(3)
#.data %>% fit_me(2) %>% filter(is.na(y0)) %>% score
#.data %>% fit_me(3) %>% filter(is.na(y0)) %>% score

#------------------------------------------------------------------#

fit_spwe <- \(.data
  , nboot = 200, compute_we = TRUE
) {
  t0 <- Sys.time()
  if (compute_we)
    graph <- W_est(.data, W)
  else
    graph <- W
  .data$time <- as.numeric(factor(.data$year))
  .data$time2 <- .data$time
  model <- inla(y0 ~ offset(log(n))
    + f(time2, model = "rw1")
    + f(i, model = "besagproper2", graph = graph
      , group = time, control.group = list(model = "ar1")
    )
    , family = "poisson"
    , data = .data
    , control.predictor = list(link = 1)
    , control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE
      , config = TRUE
    )
  )
  #summary(model)
  .data <- append_results(model, .data, nboot)
  .data$method <- "spwe_inla"
  t1 <- Sys.time()
  print(t1 - t0)
  return(.data)
}
#------------------------------------------------------------------#
#------------------------------------------------------------------#
#------------------------------------------------------------------#

fit_spx <- \(.data
  , nboot = 200, compute_we = TRUE
) {
  t0 <- Sys.time()
  if (compute_we)
    graph <- W_est(.data, W)
  else
    graph <- W
  .data$time <- as.numeric(factor(.data$year))
  model <- inla(y0 ~ offset(log(y0_hat))  - 1
    + f(i, model = "besagproper2", graph = graph
      , group = time, control.group = list(model = "ar1")
    )
    , family = "poisson"
    , data = .data
    , control.predictor = list(link = 1)
    , control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE
      , config = TRUE
    )
  )
  .data <- append_results(model, .data, nboot)
  .data$method <- "spx_inla"
  t1 <- Sys.time()
  print(t1 - t0)
  return(.data)
}

#------------------------------------------------------------------#
#------------------------------------------------------------------#
#------------------------------------------------------------------#

fit_me_inla <- \(.data
  , nboot = 200, compute_we = TRUE
) {
  t0 <- Sys.time()
  if (compute_we)
    graph <- W_est(.data, W)
  else
    graph <- W
  model <- inla(y0 ~ offset(log(y0_hat))  - 1
    + f(i, model = "iid")
    , family = "poisson"
    , data = .data
    , control.predictor = list(link = 1)
    , control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE
      , config = TRUE
    )
  )
  .data <- append_results(model, .data, nboot)
  .data$method <- "me_inla"
  t1 <- Sys.time()
  print(t1 - t0)
  return(.data)
}
#------------------------------------------------------------------#
#------------------------------------------------------------------#
#------------------------------------------------------------------#
fit_sx <- \(.data
  , nboot = 200, compute_we = TRUE
) {
  t0 <- Sys.time()
  if (compute_we)
    graph <- W_est(.data, W)
  else
    graph <- W
  .data$time <- as.numeric(factor(.data$year))
  .data$time2 <- .data$time
  model <- inla(y0 ~ offset(log(y0_hat))  - 1
    + f(i, model = "besagproper2", graph = graph)
    , family = "poisson"
    , data = .data
    , control.predictor = list(link = 1)
    , control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE
      , config = TRUE
    )
  )
  #summary(model)
  .data <- append_results(model, .data, nboot)
  .data$method <- "sx_inla"
  t1 <- Sys.time()
  print(t1 - t0)
  return(.data)
}
#------------------------------------------------------------------#
#------------------------------------------------------------------#
#inla ar
#------------------------------------------------------------------#
fit_arx <- \(.data
  , nboot = 200, compute_we = TRUE
) {
  t0 <- Sys.time()
  .data$time <- as.numeric(factor(.data$year))
  model <- inla(y0 ~ offset(log(y0_hat)) - 1
    + f(i, model = "iid", group = time, control.group = list(model = "ar1"))
    , family = "poisson"
    , data = .data
    , control.predictor = list(link = 1)
    , control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE
      , config = TRUE
    )
  )
  #summary(model)
  .data <- append_results(model, .data, nboot)
  .data$method <- "arx_inla"
  t1 <- Sys.time()
  print(t1 - t0)
  return(.data)
}

#------------------------------------------------------------------#
#------------------------------------------------------------------#
#rirw
#------------------------------------------------------------------#

#.data <- filter(placebo_dt, replication == 2)
#.data <- NULL
fit_rirw_inla <- \(.data
  , nboot = 200
) {
  t0 <- Sys.time()
  .data_inla <- .data %>%
    mutate(i1 = i
      , i2 = i + max(.data$i)
      , t1 = as.numeric(factor(year))
    )

  model <- inla(y0 ~ offset(log(n))
    + f(t1, model = "rw1")
    + f(i1, model = "iid")
    , family = "poisson"
    , data = .data_inla
    , control.predictor = list(link = 1)
    , control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE
      , config = TRUE
    )
  )
  summary(model)
  .data <- append_results(model, .data, nboot)
  .data$method <- "me_rirw_inla"
  t1 <- Sys.time()
  print(t1 - t0)
  return(.data)
}
#------------------------------------------------------------------#
#------------------------------------------------------------------#
#mixed effect with inla 
#------------------------------------------------------------------#
fit_me_inla <- \(.data
  , nboot = 200
) {
  t0 <- Sys.time()
  .data_inla <- .data %>%
    mutate(i1 = i
      , i2 = i + max(.data$i)
      , t1 = as.numeric(factor(year))
    )

  model <- inla(y0 ~ offset(log(n))+ t1
    + f(i1, model = "iid2d", n = 2 * max(.data$i), constr = TRUE)
    + f(i2, t1, copy = "i1")
    , family = "poisson"
    , data = .data_inla
    , control.predictor = list(link = 1)
    , control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE
      , config = TRUE
    )
  )
  summary(model)
  .data <- append_results(model, .data, nboot)
  .data$method <- "me_inla"
  t1 <- Sys.time()
  print(t1 - t0)
  return(.data)
}

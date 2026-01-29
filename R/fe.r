library(fixest)

#dt <- dt[, x := dt[["y0_post_1"]]]
#a <- impute_fe(dt, fm_twfe)
#dt[is.na(a), ]
#help(fepois)


fm_twfe <- formula(y0 ~ 1 | i + year)
fm_etwfe <- formula(y0 ~ 1 | start_year + year)
fm_etwfe_c <- formula(y0 ~ unemplyment + state | start_year + year)

#' Fit a Poisson fixed-effects imputer and attach bootstrap-based posterior draws.
#' Clears any prior y0_post_* and y0_* summary columns, then imputes y0_hat and y0_post_1:nb via clustered bootstrap.
impute_fe <- function(.dt, fm = fm_twfe, nb = 200, ncores = 1) {
  drop_by_pattern(.dt, "^y0_post_")
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")

  m <- fit_fe(.dt, fm)
  inference_fe(.dt, m, nb = nb, ncores = ncores)
}

#' Fit a Poisson FE model (fixest) for bootstrap-based imputation.
#' Ensures a default weight column and uses log(n) as the offset (exposure).
fit_fe <- function(.dt, fm) {
  if (!"wts" %in% names(.dt)) .dt[, wts := 1]

  fepois(
    fm,
    vcov   = ~ i,
    offset = ~ log(n),
    data   = .dt,
    weights = .dt$wts
  )
}

#' Attach clustered-bootstrap draws and point predictions from a fitted FE model.
#' Adds y0_hat and y0_post_1:nb, where each draw is a prediction on the original rows.
inference_fe <- function(.dt, m, nb = 20, ncores = 1) {

  # point prediction on original rows
  y0_hat <- predict(m, newdata = .dt, type = "response")

  # bootstrap draws: refit on bootstrap sample, predict on original rows
  y0_post <- parallel::mclapply(seq_len(nb), \(b) {
    dtb <- cluster_boot(.dt)
    mb  <- fit_fe(dtb, m$fml)
    predict(mb, newdata = .dt, type = "response")
  }, mc.cores = ncores)

  .dt[, y0_hat := y0_hat]
  .dt[, (paste0("y0_post_", seq_len(nb))) := y0_post]
  .dt
}

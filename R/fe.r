library(fixest)

#dt <- dt[, x := dt[["y0_post_1"]]]
#a <- impute_fe(dt, fm_twfe)
#dt[is.na(a), ]
#help(fepois)


fm_twfe <- formula(y0 ~ 1 | i + year)
# ETWFE: cohort x group cells (`cell`, see cohort_cell()) and common year effects.
fm_etwfe <- formula(y0 ~ 1 | cell + year)
fm_etwfe_c <- formula(y0 ~ uemp + pov + mhinc + aian | cell + year)

#' Fit a Poisson fixed-effects imputer and attach bootstrap-based posterior draws.
#' Clears any prior y0_post_* and y0_* summary columns, then imputes y0_hat and y0_post_1:nboot via clustered bootstrap.
impute_fe <- function(.dt, fm = fm_twfe, nboot = 200, ncores = 1) {
  drop_by_pattern(.dt, "^y0_post_")
  drop_by_pattern(.dt, "^y0_(hat|sd|pi|lower|upper)$")

  # Cohort x group cells from the training rows; the bootstrap weights are never zero,
  # so the point fit and every refit use the same cells.
  if ("cell" %in% all.vars(fm)) {
    cell <- cohort_cell(.dt$cohort, .dt$ncf, .dt$y0)
    .dt[, cell := as.vector(cell)]
    data.table::setattr(.dt, "cell_merges", attr(cell, "merged"))
  }
  m <- fit_fe(.dt, fm)
  inference_fe(.dt, m, nboot = nboot, ncores = ncores)
}


#' Fit a Poisson FE model (fixest) for bootstrap-based imputation.
#' Ensures a default weight column and uses log(n) as the offset (exposure).
fit_fe <- function(.dt, fm, wts = NULL) {
  if (is.null(wts)) wts <- rep(1L, nrow(.dt))

  fepois(
    fm,
    vcov   = ~ i,
    offset = ~ log(n),
    data   = .dt,
    weights = wts
  )
}

#' Attach clustered-bootstrap draws and point predictions from a fitted FE model.
#' Adds y0_hat and y0_post_1:nboot, where each draw is a prediction on the original rows.
inference_fe <- function(.dt, m, nboot = 20, ncores = 1) {

  # point prediction on original rows
  y0_hat <- predict(m, newdata = .dt, type = "response")

  # bootstrap draws: refit on bootstrap sample, predict on original rows
  y0_post <- parallel::mclapply(seq_len(nboot), \(b) {
    wts <- cluster_boot(.dt)
    mb  <- fit_fe(.dt, formula(m), wts)
    post <- predict(mb, newdata = .dt, type = "response")
    post
  }, mc.cores = ncores)

  .dt[, y0_hat := y0_hat]
  .dt[, (paste0("y0_post_", seq_len(nboot))) := y0_post]
  .dt
}

#w1 <- cluster_boot(dt)
#w2 <- cluster_boot(dt)
#identical(w1, w2)
#sum(w1)
#summary(w2)
#m1 <- fit_fe(dt, fm_twfe, w1)
#m2 <- fit_fe(dt, fm_twfe, w2)
#p1 <- predict(m1, newdata = dt, type = "response")
#p2 <- predict(m2, newdata = dt, type = "response")

#all.equal(p1, p2)
#head(p1)
#head(p2)

#formula(m1)
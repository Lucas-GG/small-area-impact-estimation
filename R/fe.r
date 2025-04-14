library(fixest)

#dt <- dt[, x := dt[["y0_post_1"]]]
#a <- impute_fe(dt, fm_twfe)
#dt[is.na(a), ]

fm_twfe <- formula(y0 ~ 1 | i + year)
fm_etwfe <- formula(y0 ~ 1 | start_year + year)

impute_fe <- \(.dt, fm) {
  if (!"wts" %in% names(.dt)) .dt[, wts := 1]
  fepois(fm
    , vcov = ~ i
    , offset = ~ log(x) #the offset is yhat
    , .dt
    , weights = .dt$wts
  ) %>%
    predict(.dt, type = "response")
}


bbx_inference <- function(.dt
  , ncores = 1
  , impute_fun = impute_fe, fm = fm_twfe
) {
  result <- copy(.dt)
  # Create a temporary copy for modifications
  y0_post_cols <- grep("y0_post", names(result), value = TRUE)

  # Process just the first two posterior columns as in your example
  .dt0 <- .dt[, .(i, year, n, y0, start_year, y0_hat)][, x := y0_hat]

  y0_hat_val <- impute_fun(.dt0, fm)

  y0_post <- mclapply(y0_post_cols, \(b) {

    .dt0[, x := .dt[[b]]]
    # Apply cluster_boot
    .dt0 <- cluster_boot(.dt0)

    # Get predictions
    predicted_values <- impute_fun(.dt0, fm)
    predicted_values
  }, mc.cores = ncores)


  result[, y0_hat := y0_hat_val]
  result[, (y0_post_cols) := y0_post]
  result
}

bbx_inference_nox <- function(.dt
  , mcores = 1
  , impute_fun = impute_fe
  , fm = fm_twfe
) {
  result <- copy(.dt)
  # Create a temporary copy for modifications
  y0_post_cols <- grep("y0_post", names(result), value = TRUE)

  # Process just the first two posterior columns as in your example
  .dt0 <- .dt[, .(i, year, n, y0, start_year)][, x := n]

  y0_hat_val <- impute_fun(.dt0, fm)

  y0_post <- mclapply(y0_post_cols, \(b) {
    # Apply cluster_boot
    .dt0 <- cluster_boot(.dt0)

    # Get predictions
    predicted_values <- impute_fun(.dt0, fm)
    predicted_values
  }, mc.cores = mcores)

  result[, y0_hat := y0_hat_val]
  result[, (y0_post_cols) := y0_post]
  result
}
library(data.table)


#' shuffle_cl()
#' Apply fractional bootstrap weights at a cluster level
#' (e.g., by unit i or by year).
#' Multiplies existing wts (or creates wts := 1) by new cluster weights
#' and normalizes weights to sum to nrow(dt).
#' Returns: a modified copy of dt with updated wts. Does NOT modify input dt.
shuffle_cl <- function(dt, cl = "i", wtype = "exp") {

  gen_w <- fwb:::make_gen_weights(wtype)
  clv <- dt[[cl]]
  ucl <- unique(clv)
  w   <- gen_w(length(ucl), 1) |> as.vector()
  pos <- match(clv, ucl)         # integer positions, length N
  w_row <- w[pos]

  # Update weights
  if ("wts" %in% colnames(dt)) {
    wts <- dt$wts * w_row
  }  else {
    wts <- w_row
  }

  # Return the result (will be visible when explicitly printed)
  wts
}


#' cluster_boot()
#' Convenience wrapper around shuffle_cl():
#' optionally shuffle by year first, then by cluster cl (default i).
#' Returns: a modified copy of dt with updated wts. Does NOT modify input dt.
cluster_boot <- function(dt
  , cl = "i", wtype = "exp", shuffle_time = FALSE
) {
  k <- nrow(dt)
  if (!"wts" %in% colnames(dt)) wts <- rep(1, k) else wts <- dt$wts

  # Apply time shuffling if requested
  if (shuffle_time) wts <- wts * shuffle_cl(dt, cl = "year", wtype = wtype)

  # Apply cluster shuffling
  wts <- wts * shuffle_cl(dt, cl = cl, wtype = wtype)
  # Return the result

  # Normalize weights
  total_weight <- sum(wts)
  wts <- wts / total_weight * k
  wts
}

if (FALSE) {
  # Example usage
  dt <- data.table(
    i = rep(1:5, each = 4),
    year = rep(2001:2004, times = 5),
    y0 = rpois(20, lambda = 5),
    n = rep(10, 20)
  )
#  dt[, wts := 1]

  # Apply cluster bootstrap
  test <- dt |> cluster_boot(cl = "i", wtype = "exp", shuffle_time = TRUE)
  print(test)
  sum(test)
}

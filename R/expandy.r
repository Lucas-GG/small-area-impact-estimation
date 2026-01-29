set_y0 <- function(dt
  , y = "y"
  , id = "i"
  , time = "year"
  , start_year_col = "start_year"
  , A_sim = NULL
) {

  if (is.null(A_sim)) {
    ok_row <- mask_pre(dt, start_year_col = start_year_col)
  } else {
    ok_row <- mask_pre(dt, A_sim, start_year_col = start_year_col)
  }

  # create masked outcome (by reference)
  dt[, y0 := fifelse(ok_row, (y), NA)]
  dt
}

#' Construct within-unit history features from a (masked) outcome series
#'
#' Adds within-unit lagged summaries of an outcome column. Intended to be used
#' with \code{y0}, a masked untreated-outcome proxy, but can be applied to any
#' numeric outcome series with missing values.
#'
#' @param dt A \code{data.table} modified by reference.
#' @param y Name of the outcome column used to build lags and summaries.
#' @param id Name of the unit identifier column.
#' @param time Name of the time index column. Must be sortable within unit.
#' @param Lset Integer vector of window lengths. For each \code{L} in \code{Lset},
#'   the function computes \code{ybar_L}, \code{ydot_L}, and \code{czero_L}.
#'
#' @details
#' Uses \code{collapse::flag()} to compute a matrix of lags for each row within
#' unit and time ordering. Summaries are computed over the first \code{L} columns
#' of the lag matrix, corresponding to the most recent \code{L} lags.
#'
#' Missing values in the lag matrix are handled with \code{na.rm=TRUE} for means
#' and are excluded from the slope calculation via \code{fast_slope_from_lags()}.
#'
#' @return \code{dt}, invisibly, with history feature columns added by reference.
#' @seealso \code{\link[collapse]{flag}}, \code{\link{fast_slope_from_lags}}
add_y0_history <- function(dt, id = "i", time = "year", conf = getOption("xpn.conf")) {
  Lset <- conf$Lset
  maxL <- max(Lset)
  Ylags_all <- collapse::flag(dt$y0, 1:maxL, dt[[id]], dt[[time]])  # N x maxL

  dt[, h_y_l1 := Ylags_all[, 1]]

  for (L in Lset) {
    Ylags <- Ylags_all[, 1:L, drop = FALSE]
    dt[, (paste0("ybar_", L))  := rowMeans(Ylags, na.rm = TRUE)]
    dt[, (paste0("ydot_", L))  := fast_slope_from_lags(Ylags, L)]
    dt[, (paste0("czero_", L)) := rowMeans(Ylags == 0, na.rm = TRUE)]
  }
  dt
}

#' Compute a local linear trend (slope) from lagged outcomes
#'
#' Computes a per-row slope of a local linear regression of lagged outcomes on
#' equally spaced lag indices, using closed-form OLS calculations for speed.
#'
#' @param Ylags Numeric matrix of lagged outcomes with \code{nrow} equal to the
#'   number of observations and \code{ncol} equal to the maximum lag length.
#'   Columns should correspond to consecutive lags (most recent first).
#' @param L Window length. Only the first \code{L} columns of \code{Ylags} are used.
#' @param min_obs Minimum number of non-missing lag values required to compute
#'   a slope. Rows with fewer observed lags return \code{NA}.
#'
#' @details
#' The slope is computed using centered weights \code{w = 1:L - (L+1)/2}.
#' Missing outcomes are excluded by zeroing their contribution in the numerator
#' and denominator using an observation mask.
#'
#' @return Numeric vector of length \code{nrow(Ylags)} with per-row slopes.
#' @keywords internal
fast_slope_from_lags <- function(Ylags, L, min_obs = 2) {
  w <- 1:L - (L + 1) / 2
  ok <- !is.na(Ylags)

  num <- rowSums(Ylags * rep(w, each = nrow(Ylags)), na.rm = TRUE)
  den <- rowSums(ok * rep(w^2, each = nrow(Ylags)), na.rm = TRUE)

  out <- num / den
  out[rowSums(ok) < min_obs | den == 0] <- NA_real_
  out
}

if (FALSE) {
  library(data.table)
  library(Matrix)
  library(ggplot2)

  set.seed(1)
  N <- 10
  Tt <- 10
  years <- 2001:2010

  dt <- data.table(
    i = rep(1:N, each = Tt),
    year = rep(years, times = N)
  )

  # unit-specific start year (Inf = never treated)
  start_per_id <- rep(2006:2010, length.out = N)
  start_per_id[start_per_id <= 2008] <- Inf
  dt[, start_year := rep(start_per_id, each = Tt)]

  # covariates
  dt[, x_unit := rnorm(N)[i]]
  dt[, x_time := 0.2 * (year - min(years)) + rnorm(.N, sd = 0.5)]
  dt[, x_bin  := rbinom(.N, 1, plogis(x_unit + 0.1 * x_time))]

  # ring adjacency G (i matches row/col)
  G <- Matrix(0, nrow = N, ncol = N, sparse = TRUE)
  for (k in 1:N) {
    G[k, ifelse(k == 1, N, k - 1)] <- 1
    G[k, ifelse(k == N, 1, k + 1)] <- 1
  }
  diag(G) <- 0

  # ---- simulate an "untreated outcome" y with mild spatial correlation ----
  # baseline linear predictor
  alpha_i <- rnorm(N, sd = 0.4)                 # unit effect
  beta_t  <- 0.05 * (dt$year - min(years))      # common time trend

  lp0 <- 0.3 + alpha_i[dt$i] + beta_t + 0.4 * dt$x_unit + 0.2 * dt$x_time + 0.2 * dt$x_bin

  # add a simple neighbor component: mean neighbor x_unit (time-invariant, easy)
  neigh_xunit <- as.numeric(G %*% alpha_i) / pmax(as.numeric(rowSums(G)), 1)
  lp0 <- lp0 + 0.3 * neigh_xunit[dt$i]

  # Poisson outcome
  dt[, y := rpois(.N, lambda = exp(lp0))]
  dt[, y25 := rpois(.N, lambda = exp(lp0 * 5))]

  # ---- set y missing after treatment starts (observed Y(0) only when untreated) ----
  conf <- list(Lset = choose_lset(dt))
  options(xpn.conf = conf)
  dt |> expand_y()
  colnames(dt)

}





#' Remove previously created within-history feature columns
#'
#' Deletes columns created by \code{add_within_history()} (and related helpers)
#' from \code{dt}. This prevents stale history features from persisting across
#' repeated feature-expansion calls.
#'
#' Columns removed match these patterns:
#' \itemize{
#'   \item \code{^h_y_l\\d+$} (lag features)
#'   \item \code{^ybar_\\d+$} (rolling means)
#'   \item \code{^ydot_\\d+$} (rolling slopes)
#'   \item \code{^czero_\\d+$} (zero shares)
#' }
#'
#' @param dt A \code{data.table} modified by reference.
#' @return \code{dt}, invisibly.
#' @keywords internal
reset_y0_features <- function(dt) {
  pats <- c(
    "^h_y_l\\d+$",     # lagged outcomes
    "^ybar_\\d+$",      # within-unit means
    "^ydot_\\d+$",      # within-unit slopes
    "^czero_\\d+$"     # zero shares
  )

  drop <- unique(unlist(lapply(
    pats, \(p) grep(p, names(dt), value = TRUE)
  )))

  if (length(drop)) data.table::set(dt, j = drop, value = NULL)
  dt
}

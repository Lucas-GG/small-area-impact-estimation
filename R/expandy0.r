# Expand features for nonparametric hazard and outcome models.
# Constructs within-unit history summaries and augments them with
# spatial, cohort-level, and contemporaneous control averages.
# Feature expansion following Appendix 1: history-based, spatial,
# cohort-level, and contemporaneous summaries of Y_it(0).
make_y0_context <- function(dt
  , A_sim = NULL
  , id = "i"
  , time = "year"
  , add_cohort = TRUE
  , cohort_col = "start_year"
  , conf = getOption("xpn.conf")
) {
  Lset <- conf$Lset
  G <- conf$graph

  if (is.null(A_sim)) {
    ok_row <- !is.na(dt$y0)
  } else {
    ok_row <- mask_pre(dt, A_sim)
  }

  hist_vars <- lapply(Lset, \(L) c(paste0("ybar_", L), paste0("ydot_", L))) |>
    unlist()

  # contemporaneous control mean
  ybar_ctrl <- dt |> contemp_control_mean(time = time, ok_row = ok_row)

  # cohort means by (A_row, year)
  coh <- NULL
  if (add_cohort) {
    coh <- dt |>
      cohort_means(cohort = dt[[cohort_col]], time = time
        , vars = hist_vars, ok_row = ok_row
      )
  }
  # spatial neighbors restricted to controls
  nbr <- NULL
  if (!is.null(G)) {
    nbr <- dt |>
      spatial_neighbors_means(G = G
        , time = time, ok_row = ok_row, vars = hist_vars
      )
  }

  list(
    ok_row = ok_row,
    cohort_row = dt[[cohort_col]],
    ybar_ctrl = ybar_ctrl,
    coh = coh,   # list of vectors (or NULL)
    nbr = nbr    # matrix (or NULL)
  )
}



# Compute cohort-level averages among untreated units.
# Aggregates within-history features by assignment cohort and time,
# using only units with observed untreated outcomes.
cohort_means <- function(dt
  , cohort
  , time = "year", vars, ok_row = NULL
) {

  g <- collapse::GRP(list(cohort, dt[[time]]))

  out <- vector("list", length(vars))
  names(out) <- paste0(vars, "_coh")

  for (v in vars) {
    x <- dt[[v]]
    x[!ok_row] <- NA_real_
    out[[paste0(v, "_coh")]] <-
      collapse::fbetween(x, g, na.rm = TRUE, fill = TRUE)
  }
  out
}


# Compute cohort-level averages among untreated units.
# Aggregates within-history features by assignment cohort and time,
# using only units with observed untreated outcomes.
contemp_control_mean <- function(dt, time = "year", ok_row) {
  x <- dt$y0
  x[!ok_row] <- NA_real_
  collapse::fbetween(x, dt[[time]], na.rm = TRUE, fill = TRUE)
}



# Compute spatially weighted neighbor averages of history features.
# Uses row-standardized contiguity weights restricted to untreated
# neighbors, recomputed separately for each time period.
spatial_neighbors_means <- function(dt, G, time = "year", ok_row, vars) {
  stopifnot(inherits(G, "Matrix"))
  N <- nrow(dt)
  p <- length(vars)

  outM <- matrix(NA_real_, N, p)
  colnames(outM) <- paste0(vars, "_nbr")

  idx_by_year <- split(seq_len(N), dt[[time]])

  for (idx in idx_by_year) {
    o <- as.numeric(ok_row[idx])
    d <- as.numeric(G %*% o)

    for (j in seq_along(vars)) {
      vv  <- dt[[vars[j]]][idx]
      num <- as.numeric(G %*% (o * vv))
      outM[idx, j] <- ifelse(d > 0, num / d, 0)
    }
  }

  outM
}
#spatial_neighbors_means(dt, G, vars = "ybar_3", ok_row = mask_pre(dt))


reset_y0_features <- function(dt) {
  pats <- c(
    "_nbr$",            # spatial neighbors
    "_coh$",            # cohort averages
    "^ybar_ctrl$"       # contemporaneous control mean
  )

  drop <- unique(unlist(lapply(
    pats, \(p) grep(p, names(dt), value = TRUE)
  )))

  if (length(drop)) data.table::set(dt, j = drop, value = NULL)
  dt
}

materialize <- function(
  art  # output of expand_y0()
  , idx = NULL # optional row subset
  , include = c("ybar_ctrl", "coh", "nbr")
  , features = NULL
) {   # optional subset of feature names

  if (is.null(idx)) idx <- seq_along(art$ybar_ctrl)
  include <- match.arg(include, several.ok = TRUE)

  out <- data.table::data.table(row = idx)

  # helper to check/select names
  want <- function(nms) {
    if (is.null(features)) nms else intersect(nms, features)
  }

  # ybar_ctrl
  if ("ybar_ctrl" %in% include && !is.null(art$ybar_ctrl)) {
    if (is.null(features) || "ybar_ctrl" %in% features) {
      out[, ybar_ctrl := art$ybar_ctrl[idx]]
    }
  }

  # cohort means: list of vectors
  if ("coh" %in% include && !is.null(art$coh)) {
    nms <- want(names(art$coh))
    for (nm in nms) {
      out[[nm]] <- art$coh[[nm]][idx]
    }
  }

  # neighbor means: matrix
  if ("nbr" %in% include && !is.null(art$nbr)) {
    nms <- want(colnames(art$nbr))
    if (length(nms)) {
      j <- match(nms, colnames(art$nbr))
      # bind selected columns
      tmp <- art$nbr[idx, j, drop = FALSE]
      for (k in seq_along(nms)) out[[nms[k]]] <- tmp[, k]
    }
  }

  out[, row := NULL]
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
  start_per_id[start_per_id > 2008] <- Inf
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
  conf <- list(Lset = choose_lset(dt), graph = G)
  options(xpn.conf = conf)
  dt |> expand_y()
  art <- dt |> expand_y0() |> materialize()
  dt[, names(art) := art] |> filter(i == 1)



  A_sim <- simulate_assignments(dt)
  art <- dt |> expand_y0(A_sim = A_sim, G = G)
  art |> materialize() |> head()
  art <- dt |> expand_y0(A_sim = A_sim, G = G, add_cohort = FALSE)
  art |> materialize() |> head()
  str(art)

}
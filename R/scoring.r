library(data.table)

# -------------------------------------------------------------------
# Poisson draw helpers
# -------------------------------------------------------------------

#' mrpois()
#' Draw independent Poisson variates with mean matrix .m; returns matrix same dim.
mrpois <- function(.m) {
  m <- as.matrix(.m)
  matrix(stats::rpois(length(m), lambda = as.numeric(m)),
         nrow = nrow(m), ncol = ncol(m))
}

#' bounds_sim()
#' Simulation-based 95% bounds per row using one Poisson draw per cell.
#' (Matches your original behavior; not "many simulations".)
bounds_sim <- function(.m, probs = c(.025, .975)) {
  out <- apply(mrpois(.m), 1, stats::quantile, probs = probs, na.rm = TRUE)
  t(out)
}

# -------------------------------------------------------------------
# Posterior evaluation (data.table)
# -------------------------------------------------------------------

#' eval_post_dt()
#'
#' Purpose:
#'   From observed y and posterior draws y0_post, compute rowwise:
#'     - y0_hat (optional mean)
#'     - y0_sd  (posterior SD of y0_post)
#'     - y0_pi  (predictive mass/density at observed y)
#'     - y0_lower/y0_upper (interval bounds)
#'
#' Data contract (.dt must be long or wide; one row per time-unit):
#'   Required:
#'     - y : observed count (numeric/integer)
#'   Posterior draws (one of):
#'     (A) matrix-like column `y0_post` coercible via as.matrix()
#'     (B) multiple columns named y0_post* (e.g., y0_post_1, y0_post_2, ...)
#'
#' Notes:
#'   - By default returns a copy. Set copy_dt=FALSE to modify by reference.
#'   - Requires fsd(), pi_poilog(), ci_poilog() in your environment for "poilog".
eval_post <- function(.dt,
                      pi_type = c("poilog", "interval", "kernel"),
                      bounds_type = c("sim", "poilog", "interval"),
                      include_mean = FALSE,
                      copy_dt = FALSE) {

  pi_type     <- match.arg(pi_type)
  bounds_type <- match.arg(bounds_type)

  dt <- if (copy_dt) copy(.dt) else .dt
  stopifnot(is.data.table(dt))
  if (!("y" %in% names(dt))) stop("eval_post(): dt must contain column 'y'.")

  # --- Extract posterior draws matrix ---
  y0_post_mat <- NULL

  if ("y0_post" %in% names(dt) && !is.null(dt[["y0_post"]])) {
    y0_post_mat <- as.matrix(dt[["y0_post"]])
  } else {
    y0_post_cols <- grep("^y0_post", names(dt), value = TRUE)
    if (length(y0_post_cols) == 0) {
      stop("eval_post(): expected column 'y0_post' or columns starting with 'y0_post'.")
    }
    y0_post_mat <- as.matrix(dt[, ..y0_post_cols])
  }

  if (nrow(y0_post_mat) != nrow(dt)) {
    stop("eval_post_dt(): y0_post matrix must have nrow equal to nrow(dt).")
  }

  y_vec <- dt[["y"]]

  # Optional mean
  if (include_mean) {
    dt[, y0_hat := rowMeans(y0_post_mat)]
  }

  # Posterior SD per row
  y0_sd_val <- fsd(t(y0_post_mat))

  # --- y0_pi ---
  y0_pi <- switch(
    pi_type,

    "interval" = {
      # fraction of draws equal to y
      rowMeans(replicate(ncol(y0_post_mat), y_vec) == y0_post_mat)
    },

    "poilog" = {
      log_y0 <- log(y0_post_mat)
      eta_mean <- rowMeans(log_y0)
      eta_sd   <- fsd(t(log_y0))
      pi_poilog(y_vec, eta_mean, eta_sd)
    },

    "kernel" = {
      # KDE per row using internal kde1d objects
      F_y0 <- apply(y0_post_mat, 1, kde1d:::kde1d, type = "d")
      out  <- vapply(seq_along(F_y0),
                     function(i) kde1d:::dkde1d(y_vec[i], F_y0[[i]]),
                     numeric(1))
      if (any(out == 0, na.rm = TRUE)) {
        out[out == 0] <- min(out[out > 0], na.rm = TRUE)
      }
      out
    }
  )

  # --- bounds ---
  bounds <- switch(
    bounds_type,

    "sim" = bounds_sim(y0_post_mat),

    "poilog" = {
      log_y0 <- log(y0_post_mat)
      eta_mean <- rowMeans(log_y0)
      eta_sd   <- fsd(t(log_y0))
      ci_poilog(eta_mean, eta_sd)  # should return n x 2
    },

    "interval" = {
      t(apply(y0_post_mat, 1, stats::quantile, probs = c(.025, .975), na.rm = TRUE))
    }
  )

  dt[, `:=`(
    y0_sd    = y0_sd_val,
    y0_pi    = y0_pi,
    y0_lower = bounds[, 1],
    y0_upper = bounds[, 2]
  )]

  dt
}

# -------------------------------------------------------------------
# Scoring rules (vectorized helpers)
# -------------------------------------------------------------------

ae  <- function(y, y0) abs(y - y0)
se  <- function(y, y0) (y - y0)^2

sae <- function(y, y0, n) abs((y - y0) / n * 1e5)
sse <- function(y, y0, n) ((y - y0) / n * 1e5)^2

rae <- function(y, y0) abs((y / y0) - 1)
rse <- function(y, y0) ((y / y0) - 1)^2

cvg <- function(y, l, u) as.numeric(l <= y & y <= u)

intv <- function(y, l, u, alpha = .05) {
  (u - l) + 2 / alpha * (l - y) * (y < l) + 2 / alpha * (y - u) * (y > u)
}

d2 <- function(y, y0) {
  ifelse(y > 0, 2 * (y * log(y / y0) + y0 - y), 2 * y0)
}

ds <- function(y, y0, y0_sd) {
  sigpred <- sqrt(y0 + y0_sd^2)
  ((y - y0) / sigpred)^2 + 2 * log(sigpred)
}

# -------------------------------------------------------------------
# Add rowwise scores to a data.table
# -------------------------------------------------------------------

#' score_dt()
#'
#' Data contract:
#'   dt must contain y, n, y0_hat.
#'   For full scoring also needs y0_pi, y0_sd, y0_lower, y0_upper.
#'
#' Returns:
#'   dt with additional score columns (modifies in-place unless copy_dt=TRUE).
score <- function(dt, idx, copy_dt = FALSE) {
  if (copy_dt) dt <- copy(dt)

  # apply subset if provided
  if (!missing(idx)) {
    dt <- dt[eval(substitute(idx))]
  }

  # ... exist
  req <- c("y", "n", "y0_hat")
  miss <- setdiff(req, names(dt))
  if (length(miss)) {
    stop("score_dt(): missing columns: ", paste(miss, collapse = ", "))
  }

  # Optional columns (used if present)
  has_pi <- "y0_pi" %in% names(dt)
  has_ci <- all(c("y0_lower", "y0_upper") %in% names(dt))
  has_sd <- "y0_sd" %in% names(dt)

  dt[, `:=`(
    d   = y - y0_hat,
    ae  = ae(y, y0_hat),
    se  = se(y, y0_hat),
    r   = (y / y0_hat) - 1,
    rae = rae(y, y0_hat),
    rse = rse(y, y0_hat),
    dr  = (y - y0_hat) / n * 1e5,
    sae = sae(y, y0_hat, n),
    sse = sse(y, y0_hat, n),
    d2  = d2(y, y0_hat)
  )]

  if (has_pi) {
    dt[, lns := -log(y0_pi)]
  }

  if (has_sd) {
    dt[, ds := ds(y, y0_hat, y0_sd)]
    dt[, sigpred := y0_hat + y0_sd^2]
  }

  if (has_ci) {
    dt[, cvg := cvg(y, y0_lower, y0_upper)]
    dt[, intv := intv(y, y0_lower, y0_upper)]
  }

  dt
}

# -------------------------------------------------------------------
# Summaries (data.table)
# -------------------------------------------------------------------

#' summ_score_dt()
#'
#' If dt is keyed, summarizes by key(dt). Otherwise returns overall means.
#' Returns RMSE-style versions for se/rse/sse, and sqrt(sigpred).
summ_score <- function(dt) {
  x <- as.data.table(dt)

  # choose columns that exist
  score_cols <- intersect(
    c("d","ae","se","r","rae","rse","dr","sae","sse","d2","lns","ds","cvg","intv","sigpred"),
    names(x)
  )

  group_cols <- key(x)
  if (is.null(group_cols) || length(group_cols) == 0) {
    out <- x[, lapply(.SD, mean, na.rm = TRUE), .SDcols = score_cols]
    out[, `:=`(
      pcor = stats::cor(x$y, x$y0_hat, use = "pairwise"),
      scor = stats::cor(x$y, x$y0_hat, method = "spearman", use = "pairwise")
    )]
  } else {
    out <- x[, c(lapply(.SD, mean, na.rm = TRUE),
                 list(pcor = stats::cor(y, y0_hat, use = "pairwise"),
                      scor = stats::cor(y, y0_hat, method = "spearman", use = "pairwise"))),
             by = group_cols,
             .SDcols = score_cols]
  }

  # convert MSE-like to RMSE-like where present
  if ("se"  %in% names(out)) out[, se := sqrt(se)]
  if ("rse" %in% names(out)) out[, rse := sqrt(rse)]
  if ("sse" %in% names(out)) out[, sse := sqrt(sse)]
  if ("sigpred" %in% names(out)) out[, sigpred := sqrt(sigpred)]

  out
}

# -------------------------------------------------------------------
# Aggregation of y0_post by groups (data.table)
# -------------------------------------------------------------------

#' agg_dt()
#'
#' Aggregate n, y, y0_hat by groups and also aggregate posterior draws y0_post
#' by summing draws within groups (columnwise).
#'
#' Data contract:
#'   - .dt has columns: n, y, y0_hat
#'   - posterior draws stored either in column y0_post (matrix-like per row)
#'     OR as y0_post* columns (wide draws).
agg_dt <- function(.dt, ..., copy_dt = TRUE) {
  x <- if (copy_dt) copy(as.data.table(.dt)) else as.data.table(.dt)

  grp <- as.character(substitute(list(...)))[-1L]
  if (length(grp) == 0) stop("agg_dt(): provide grouping columns via ...")

  # Extract y0_post matrix (same logic as eval_post_dt)
  if ("y0_post" %in% names(x) && !is.null(x[["y0_post"]])) {
    y0_post_mat <- as.matrix(x[["y0_post"]])
  } else {
    y0_post_cols <- grep("^y0_post", names(x), value = TRUE)
    if (length(y0_post_cols) == 0) stop("agg_dt(): no y0_post found.")
    y0_post_mat <- as.matrix(x[, ..y0_post_cols])
  }

  # Sum n,y,y0_hat by group
  out <- x[, .(
    n = sum(n, na.rm = TRUE),
    y = sum(y, na.rm = TRUE),
    y0_hat = sum(y0_hat, na.rm = TRUE)
  ), by = grp]

  # Aggregate posterior draws by group: sum draws across rows in group
  # Strategy: split row indices by group, then colSums on each block.
  gidx <- x[, .I, by = grp]
  idx_list <- split(gidx$I, gidx[, do.call(paste, c(.SD, sep = "\r")) , .SDcols = grp])

  y0_post_agg <- lapply(idx_list, function(ii) {
    colSums(y0_post_mat[ii, , drop = FALSE], na.rm = TRUE)
  })
  y0_post_agg <- do.call(rbind, y0_post_agg)

  # Attach draws as a matrix column to out (keeps it compact)
  out[, y0_post := y0_post_agg]

  out[]
}




.pretty_labels <- c(
  d      = "Mean error",
  ae     = "Mean absolute error",
  se     = "Mean squared error",
  r      = "Mean relative error",
  rae    = "Mean absolute relative error",
  rse    = "Mean squared relative error",
  dr     = "Mean rate difference (per 100k)",
  sae    = "Mean abs. rate error (per 100k)",
  sse    = "Mean sq. rate error (per 100k)",
  d2     = "Poisson deviance",
  pcor   = "Pearson correlation",
  scor   = "Spearman correlation"
)

pretty_fit <- function(x, digits = 3, labels = NULL) {
  dt <- data.table::as.data.table(x)

  # metric columns are numeric columns (you can tweak this rule if needed)
  metric_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  group_cols  <- setdiff(names(dt), metric_cols)

  out <- data.table::melt(
    dt,
    id.vars = group_cols,
    measure.vars = metric_cols,
    variable.name = "metric",
    value.name = "value"
  )

  out[, value := round(value, digits)]

  if (!is.null(labels)) {
    out[, label := labels[as.character(metric)]]
    out[is.na(label), label := as.character(metric)]
    # Put label first for readability
    setcolorder(out, c(group_cols, "label", "metric", "value"))
  } else {
    setcolorder(out, c(group_cols, "metric", "value"))
  }

  out[]
}

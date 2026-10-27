library(data.table)

#' Simulate "sticky" (mostly-backward) reassigned adoption years per unit.
#'
#' Purpose:
#'   Given a panel dt with unit id i, calendar year, and observed start_year,
#'   generate a simulated start year A_sim for each unit.
#'
#' Typical use:
#'   - placebo / permutation / resampling in staggered adoption settings
#'   - creating alternative (plausible) assignment vectors for training / validation
#'
simulate_sticky_assignments <- function(
  dt,
  shift = 1,
  shift_type = c("poisson", "deterministic"),
  clamp_years = TRUE, # clamp A_sim to [min(year), max(year)]
  snap_to_available = TRUE, # snap to closest observed panel year (handles gaps)
  haz_col = "cpr",
  use_hazard = TRUE,           # multiply sticky weights by hazard weights
  max_shift = NULL,            # backward window size (years)
  hazard_power = 1,            # temper hazard effect (0.5, 1, 2, ...)
  eps_haz = 1e-8,              # avoid zero weights if hazard is zero everywhere
  adopt_col = "mask_year"      # centre: first masked year, min(actual, simulated)
) {
  shift_type <- match.arg(shift_type)
  if (!inherits(dt, "data.table")) dt <- as.data.table(dt)

  # --- Block 1: panel year support (used for clamping/snapping)
  years_avail <- sort(unique(dt$year))
  miny <- min(years_avail)
  maxy <- max(years_avail)

  snap_year <- function(y) {
    if (!snap_to_available) return(y)
    years_avail[which.min(abs(years_avail - y))]
  }

  clamp_year <- function(y) {
    if (!clamp_years) return(y)
    min(max(y, miny), maxy)
  }

  # Default backward window: at least 6, or about 2*shift
  if (is.null(max_shift)) max_shift <- max(6L, as.integer(2 * shift))

  # --- Block 2: unit-level observed starts
  # One row per unit i with its observed adoption year A_obs
  unit_dt <- unique(dt[, .(i, A_obs = get(adopt_col))])

  # --- Block 3: optional hazard lookup table h(i, year)
  # If use_hazard and haz_col exists, we build a keyed lookup haz_dt.
  has_haz <- use_hazard && (haz_col %in% names(dt))
  if (has_haz) {
    haz_dt <- unique(dt[, .(i, year, h = get(haz_col))])
    setkey(haz_dt, i, year)
  }

  # --- Block 4: candidate construction (ONLY backward from A_obs)
  make_candidates <- function(A0) {
    # never-treated units: keep Inf (no candidates)
    if (!is.finite(A0)) return(Inf)

    cand <- seq(A0 - max_shift, A0)  # backward-only window
    if (clamp_years) cand <- cand[cand >= miny & cand <= maxy]

    if (snap_to_available) {
      cand <- vapply(cand, snap_year, numeric(1))
    }
    unique(cand)
  }


  # --- Block 5: weight builders
  # Sticky weights: Poisson pmf on backward distance d = A0 - cand (>= 0)
  sticky_weights <- function(A0, cand) {
    # Guard against Inf - Inf (never-treated shouldn't reach here, but be safe)
    if (!is.finite(A0) || any(!is.finite(cand))) {
      return(rep(1, length(cand)))
    }
    d <- A0 - cand
    w <- stats::dpois(d, lambda = shift)
    if (!is.finite(sum(w)) || sum(w) == 0) w <- rep(1, length(cand))
    w
  }

  # Hazard weights: h(i, cand)^hazard_power, with floor eps_haz
  hazard_weights <- function(uid, cand) {
    # Return uniform weights if no hazard column or non-finite candidates
    if (!has_haz || any(!is.finite(cand))) {
      return(rep(1, length(cand)))
    }
    cand <- as.integer(cand)
    # The hazard stored at year t is Pr(adopt at t + 1 | not yet adopted, X_t)
    # (prep_hazard_data), so adoption at candidate year c uses the row c - 1.
    # A candidate without a previous panel year gets the floor eps_haz.
    h <- haz_dt[J(uid, cand - 1L), h]
    h[is.na(h)] <- eps_haz
    w <- (pmax(h, eps_haz))^hazard_power
    if (!is.finite(sum(w)) || sum(w) == 0) w <- rep(1, length(cand))
    w
  }


  # --- Block 6: draw one simulated start year per unit
  draw_A_sim <- function(i, A0) {

    if (shift_type == "deterministic") {
      if (!is.finite(A0)) return(Inf)
      As <- A0 - shift
      As <- clamp_year(As)
      if (snap_to_available) As <- snap_year(As)
      return(As)
    }

    # Keep never-treated as never-treated
    if (!is.finite(A0)) return(Inf)

    # poisson (sticky) sampling over backward candidates
    cand <- make_candidates(A0)
    if (length(cand) == 1L) return(cand)

    w <- sticky_weights(A0, cand) * hazard_weights(i, cand)
    if (!is.finite(sum(w)) || sum(w) == 0) w <- rep(1, length(cand))

    sample(cand, size = 1L, prob = w)
  }

  #draw_A_sim(unit_dt$i[1], unit_dt$A_obs[1])

  # Apply per unit
  unit_dt[, A_sim := NA_real_]
  unit_dt[, A_sim := draw_A_sim(i, A_obs), by = i]

  # Return as named vector: A_sim[i]
  setNames(unit_dt$A_sim, unit_dt$i)
}


shuffle_a <- function(dt
  , id_col = "i"
  , year_col = "year"
  , hazard_col = "cpr"
  , start_year_col = "start_year"
) {
  # unit-level A
  unit_A <- unique(dt[, .(A = get(start_year_col)), by = id_col])

  treated_A <- unit_A[is.finite(A), A]

  n_by_year <- as.integer(table(treated_A))
  eligible_years <- as.integer(names(table(treated_A)))

  # weights table: one row per (unit, eligible year) with weight = hazard at that year
  w <- dt[get(year_col) %in% eligible_years,
          .(w = get(hazard_col)[1L]),
          by = c(id_col, year_col)]
  w[is.na(w), w := 0]  # keep it simple

  # initialize simulated assignment: Inf for everyone
  ids <- unit_A[[id_col]]
  A_sim <- setNames(rep(Inf, length(ids)), as.character(ids))

  # units still available to assign
  unassigned <- as.character(ids)

  for (k in seq_along(eligible_years)) {
    t <- eligible_years[k]
    need <- n_by_year[k]
    if (need == 0L) next

    cand <- w[get(year_col) == t & as.character(get(id_col)) %chin% unassigned,
              .(id = as.character(get(id_col)), w)]

    picked <- sample(cand$id, size = need, replace = FALSE, prob = cand$w)

    A_sim[picked] <- t
    unassigned <- setdiff(unassigned, picked)
  }

  #dt_out <- copy(dt)
  #dt_out[, start_year_sim := A_sim[as.character(get(id_col))]]
  #list(A_sim = A_sim, dt = dt_out)
  A_sim
}



# -------------------------------------------------------------------
# simulate_assignments()
#
# Purpose:
#   Simulate new treatment start times under staggered adoption,
#   preserving observed cohort sizes and assigning units using
#   hazard-weighted sampling.
#
# Data contract (dt must be long format with):
#   i          : unit identifier
#   year       : calendar year (integer, balanced panel)
#   start_year : observed treatment start year
#                - finite for treated units
#                - Inf for never-treated units
#   <haz_col>  : discrete-time hazard
#                Pr(A_i = year | A_i > year-1, X_it)
#
# Key assumptions:
#   - Treatment is absorbing: once treated, always treated.
#   - Hazard is evaluated at the target start year.
#   - Never-treated units (start_year == Inf) are right-censored.
#
# Options:
#   shift_type = "deterministic" | "replication" | "cohort"
#     deterministic : all cohorts shifted by a fixed amount
#     replication   : single random shift shared across cohorts
#     cohort        : cohort-specific random shifts
#
# Returns:
#   - A_sim : data.table with simulated start year per unit
#   - dt    : original dt with simulated treatment indicator W_sim
#
# Notes:
#   - The number of treated units is preserved by construction.
#   - The number of distinct cohort years is not preserved
#     when shift_type = "cohort".
# -------------------------------------------------------------------
#dt |> simulate_assignments(shift_type = "shift_only") |> names()
#dt |> simulate_assignments() |> names()

simulate_assignments <- function(
  dt,
  shift = 1,
  shift_type = c("replication", "cohort", "deterministic"),
  include_never_treated = TRUE,
  haz_col = "cpr",
  clamp_years = TRUE, # shift A is out-of-range
  snap_to_available = TRUE # shift A is in-range but missing
) {

  shift_type <- match.arg(shift_type)

  years_avail <- sort(unique(dt$year))
  miny <- min(years_avail)
  maxy <- max(years_avail)
  snap_year <- function(y) {
    if (!snap_to_available) return(y)
    years_avail[which.min(abs(years_avail - y))]
  }

  # One row per unit with observed start
  unit_dt <- unique(dt[, .(i, A_obs = start_year)])

  treated_units <- unit_dt[is.finite(A_obs)]
  never_units   <- unit_dt[!is.finite(A_obs)]

  # Observed cohort sizes among treated
  N_by <- treated_units[, .N, by = A_obs][order(A_obs)]
  cohorts <- N_by$A_obs

  # Build shift map for COHORT YEARS (finite only)
  if (shift_type == "deterministic") {
    shift_map <- data.table(A_obs = cohorts, A_shift = cohorts - shift)
  } else if (shift_type == "replication") {
    u <- sample(0:(2 * shift), 1)
    shift_map <- data.table(A_obs = cohorts, A_shift = cohorts - u)
  } else { # "cohort"
    u <- sample(0:(2 * shift), length(cohorts), replace = TRUE)
    shift_map <- data.table(A_obs = cohorts, A_shift = cohorts - u)
  }


  if (clamp_years) {
    shift_map[A_shift < miny, A_shift := miny]
    shift_map[A_shift > maxy, A_shift := maxy]
  }
  shift_map[, A_shift := vapply(A_shift, snap_year, numeric(1))]

  # Turn cohort shifts into counts needed per target year
  shift_map <- merge(shift_map, N_by, by = "A_obs", all.x = TRUE)
  target_counts <-
    shift_map[, .(N_target = sum(N)), by = A_shift][order(A_shift)]

  # Candidate pool
  pool_ids <- if (include_never_treated) unit_dt$i else treated_units$i

  # Simulated starts (default: Inf = never treated)
  A_sim <- data.table(i = pool_ids)
  A_sim[, A_sim := Inf]

  # Assign starts year-by-year (increasing target year)
  for (k in seq_len(nrow(target_counts))) {
    t0 <- target_counts$A_shift[k]
    n_need <- target_counts$N_target[k]
    if (n_need <= 0) next

    remaining_ids <- A_sim[is.infinite(A_sim), i]
    if (length(remaining_ids) == 0) break

    if (haz_col %in% colnames(dt)) {
      w_dt <- dt[year == t0 & i %in% remaining_ids
        , .(i, w = get(haz_col))
      ] |> unique()
    } else {
      w_dt <- data.table(i = remaining_ids, w = 1)
    }

    if (nrow(w_dt) == 0) next
    w_dt[is.na(w) | w < 0, w := 0]
    if (sum(w_dt$w) == 0) w_dt[, w := 1]

    take <- min(n_need, nrow(w_dt))
    chosen <- sample(w_dt$i, size = take, replace = FALSE, prob = w_dt$w)

    A_sim[i %in% chosen, A_sim := t0]
  }

  A_sim <- setNames(A_sim$A_sim, as.character(A_sim$i))
  A_sim
}

if (FALSE) {
  start_per_id <- rep(2006:2010, length.out = 10)
  start_per_id[start_per_id <= 2008] <- Inf
  dt <- data.table(
    i = rep(1:10, each = 10)
    , year = rep(2001:2010, times = 10)
    , start_year = rep(start_per_id, each = 10)
    #, cpr = rep(1, times = 100)
  )
  simulate_assignments(dt) |> print()
  simulate_assignments(dt, shift_type = "shift_only") |> print()

}




#' Identify pre-treatment rows under observed and (optionally) simulated adoption
#'
#' Returns a logical vector indicating which unit-time rows are pre-treatment
#' ("untreated") under the observed adoption time \code{start_year} and, if
#' provided, under a simulated adoption time \code{A_sim}.
#'
#' For each row, a row is considered pre-treatment if:
#' \enumerate{
#'   \item \code{year < start_year} (observed assignment), and
#'   \item if \code{A_sim} is provided, \code{year < A_sim[id]} (simulated assignment).
#' }
#'
#' @param dt A data.table containing at least \code{id_col}, \code{year_col},
#'   and \code{start_year_col}.
#' @param A_sim Optional simulated adoption time. Expected to be indexable by
#'   unit id coerced to character, i.e., accessed as \code{A_sim[as.character(id)]}.
#'   Typical use is a named vector with names equal to unit IDs.
#' @param id_col Name of the unit identifier column (default \code{"i"}).
#' @param year_col Name of the time column (default \code{"year"}).
#' @param start_year_col Name of the observed adoption-time column
#'   (default \code{"start_year"}). \code{Inf} values are allowed and imply
#'   "never treated"; comparison \code{year < Inf} evaluates to TRUE.
#'
#' @details
#' \code{A_sim} is treated as a unit-level vector of adoption years. It is
#' mapped to rows via \code{A_row <- unname(A_sim[as.character(id)])}.
#' If \code{A_sim} does not contain entries for some ids, the resulting
#' \code{A_row} will be \code{NA}, and the simulated pre-treatment condition
#' \code{year < A_row} will be \code{NA} for those rows, propagating to the
#' final mask.
#'
#' @return A logical vector of length \code{nrow(dt)}. If \code{A_sim} is provided
#'   but has missing ids, the result may contain \code{NA}.
#' @keywords internal
mask_pre <- function(dt
  , A_sim = NULL # vector: either length N indexable, or named by id
  , id_col = "i"
  , year_col = "year"
) {

  is_pre <- !is.na(dt$y0)          # Inf works naturally

  if (!is.null(A_sim)) {
    A_row <- unname(A_sim[as.character(dt$i)])
    is_pre_sim <- (dt[[year_col]] < A_row)
    is_pre <- is_pre & is_pre_sim
  }
  is_pre
}

#mask_pre(dt) |> table()
#mask_pre(dt, A_sim) |> table()
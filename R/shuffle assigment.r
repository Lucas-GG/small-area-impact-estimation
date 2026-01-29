library(data.table)
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

simulate_assignments <- function(
  dt,
  shift = 1,
  shift_type = c("replication", "cohort", "deterministic", "shift_only"),
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

  # SHIFT-ONLY MODE: in-time placebo (no reassignment)
  if (shift_type == "shift_only") {
    unit_dt[, A_sim := A_obs - shift]
    unit_dt[, A_obs := NULL]

    # Clamp to observed window (optional)
    if (clamp_years) {
      unit_dt[is.finite(A_sim) & A_sim < miny, A_sim := miny]
      unit_dt[is.finite(A_sim) & A_sim > maxy, A_sim := maxy]
    }

    # Snap to nearest observed year (optional; mostly useful if years have gaps)
    if (snap_to_available) {
      unit_dt[is.finite(A_sim),
        A_sim := vapply(A_sim, snap_year, numeric(1))
      ]
    }
    unit_dt <- setNames(unit_dt$A_sim, unit_dt$i)
    return(unit_dt)
  }


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
  , start_year_col = "start_year"
) {

  year <- dt[[year_col]]
  A_obs <- dt[[start_year_col]]
  is_pre_obs <- year < A_obs          # Inf works naturally

  if (is.null(A_sim)) {
    return(is_pre_obs)
  } else {
    idv <- dt[[id_col]]
    A_row <- unname(A_sim[as.character(idv)])
    is_pre_sim <- (year < A_row)
    is_pre_obs & is_pre_sim
  }
}

#mask_pre(dt) |> table()
#mask_pre(dt, A_sim) |> table()
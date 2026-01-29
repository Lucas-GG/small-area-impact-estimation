#' Choose history window lengths (Lset) based on available pre-treatment time span
#'
#' Chooses one or two window lengths used for within-unit history summaries
#' (e.g., \code{ybar_L}, \code{ydot_L}, \code{czero_L}), based on the amount of
#' pre-treatment history available before the earliest observed adoption time.
#'
#' @param dt A data.table containing at least \code{year_col} and \code{start_year_col}.
#' @param year_col Name of the time column (default \code{"year"}).
#' @param start_year_col Name of the observed adoption-time column
#'   (default \code{"start_year"}). May contain \code{Inf} for never-treated units.
#' @param min_L Minimum allowed history window length (default 2).
#' @param max_L Maximum allowed history window length (default 10).
#' @param shift_max Maximum amount (in years) by which adoption times may be
#'   shifted earlier in simulated assignments (default 2). Used to ensure
#'   the selected history window remains feasible under simulated shifts.
#'
#' @details
#' The function computes:
#' \itemize{
#'   \item \code{first_year}: earliest observed time in \code{dt}.
#'   \item \code{earliest_start}: earliest finite \code{start_year} minus \code{shift_max}.
#'   \item \code{pre_history_years = earliest_start - first_year}.
#' }
#' It then defines a "long" window:
#' \code{L_long = max(min_L, min(pre_history_years, max_L))}.
#'
#' If \code{L_long} is at least \code{min_L + 2}, a "short" window is also returned:
#' \code{L_short = max(min_L, floor(L_long / 2))}.
#' The output is the sorted unique set \code{c(L_short, L_long)}.
#'
#' If no finite \code{start_year} values exist (i.e., all units are never-treated),
#' the function returns \code{min_L}.
#'
#' @return An integer vector of window lengths. Either a scalar \code{min_L}
#'   or a length-2 vector \code{c(L_short, L_long)}.
#' @keywords internal
choose_lset <- function(dt
  , year_col = "year"
  , start_year_col = "start_year"
  , min_L = 2
  , max_L = 10
  , shift_max = 2
) {

  first_year <- min(dt[[year_col]], na.rm = TRUE)

  sy <- dt[[start_year_col]]
  if (!any(is.finite(sy))) return(min_L)

  earliest_start <- min(sy[is.finite(sy)], na.rm = TRUE)  - shift_max
  pre_history_years <- as.integer(earliest_start - first_year)

  L_long <- max(min_L, min(pre_history_years, max_L))

  if (L_long >= (min_L + 2L)) {
    L_short <- max(min_L, floor(L_long / 2))
    return(sort(unique(c(L_short, L_long))))
  }
  min_L
}


#' Drop data.table columns whose names match a regex.
#' Convenience for cleaning up previously-added posterior draw / summary columns in-place.
drop_by_pattern <- function(dt, pattern) {
  cols <- grep(pattern, names(dt), value = TRUE)
  if (length(cols)) dt[, (cols) := NULL]
  invisible(dt)
}

remove_names <- \(.mx) {
  dimnames(.mx) <- NULL
  .mx
}



rps <- \(y, mu, sigma, nsum = 100) {
  nsum <- sapply(seq_len(length(mu)), \(i) {
    sads:::qpoilog(1 - 1e-2, mu[i], sigma[i])
  })
  Fy <- \(qv) sapply(seq_len(length(qv)), \(q) sads:::ppoilog(q, mu, sigma))
  ysum <- 0:nsum
  indicator <- ifelse(ysum - y >= 0, 1, 0) # Heaviside
  score <- sum((indicator - Fy(ysum))^2)
  score
}


pi_poilog <- \(y, mu, sigma) {
  mapply(\(x, w, z) sads:::dpoilog(x, w, z)
    , y, mu, sigma
  )
}

ci_poilog <- \(mu, sigma, q = c(.025, .975)) {
  mapply(\(w, z) sads:::qpoilog(q, w, z)
    , mu, sigma
  ) %>% t
}

#minor modification to avoid an error
qfinder2 <- \(dist, want, coef, init = 0) {
  if (want < 0) return(0)
  if (want >= 1) return(Inf)
  q0 <- sum(do.call(dist, c(list(x = 0:init), coef)))

  if (is.nan(q0)) return(NaN)
  if (q0 >= want) return(init)
  step <- 1
  guess <- init + 2
  last <- init + 1
  cum <- q0
  repeat {
    my.q <- do.call(dist, c(list(x = last:guess), coef))
    my.sq <- sum(my.q)
    if (cum + my.sq > want) break
    if (my.sq < 10^(-10^3)) stop("quantile function did not converge!")
    last <- guess + 1
    step <- step * 2
    guess <- guess + step
    cum <- cum + my.sq
  }
  my.sq <- cum + cumsum(my.q)
  add <- which(my.sq < want - .Machine$double.eps)
  if (length(add)) last <- last + max(add)
  return(last)
}

environment(qfinder2) <- asNamespace("sads")
assignInNamespace("qfinder", qfinder2, ns = "sads")

#better neighborhood
W_est <- \(.data, W) {
  .data %>%
    filter(!is.na(y0)) %>%
    group_by(i) %>%
    reframe(n = sum(n), y = sum(.data$y0), r = y / n * 10^5) %>%
    select(r) %>%
    as.matrix %>%
    CARBayesST:::W.estimate(W, .)
}

## Function to get the structure matrix for the RW1
## Argument:
## n - length of the random walk
rw1 = function(n) {
    ## RW1
    L = matrix(0, n, n)
    for(i in 2:n) {
        L[c(i-1, i), i] = c(-1, 1)
    }
    Q = L %*% t(L)
    return (Q)
}

## ... and RW2
rw2 = function(n) {
    R = rw1(n)
    R = R[-c(1, n), ]
    Q = t(R) %*% R
    return (Q)
}

## get reference standard deviation derived by taking
## the geometric mean from the marginal standard deviations.
## Argument:
## Q - Structure matrix of the IGMRF
get.refSD <- function(Q){

    require(MASS)
    ## get the generalized inverse 
    gi <- ginv(Q)
    ## the marginal variances are on its diagonal
    mvar <- diag(gi)
    ## get the marginal standard deviations
    msd <- sqrt(mvar)
    ## compute the reference standard deviation as proposed
    ## by Sorbye and Rue 2013
    refSD <- exp(mean(log(msd)))

    return(refSD)
}

## correct aisolated areas
correct.island <- function(nb,poly){
  # 1) Identify islands
  no_neighs <- which(card(nb) == 0) # returns set cardinality.
  
  #2) For each island find nearest neighbour  
  k1nb <- knn2nb(knearneigh(coordinates(poly), k=1,longlat=!is.projected(poly)))
  
  #3) Assign a neighbour relationship between the island and its nearest neighbour
  
  is.symmetric.nb(nb, force=TRUE)
  nb1 <- nb
  nb1[no_neighs] <- k1nb[no_neighs]
  attr(nb1, "sym") <- is.symmetric.nb(nb1, force=TRUE) #! re-assign the now incorrect symmetry attribute
  nb2 <- make.sym.nb(nb1)
  return(nb2)
  }


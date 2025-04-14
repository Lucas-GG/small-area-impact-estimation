remove_names <- \(.mx) {
  dimnames(.mx) <- NULL
  .mx
}
#y <- 10 # actual value
#mu <- 4
#sigma <- 1




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
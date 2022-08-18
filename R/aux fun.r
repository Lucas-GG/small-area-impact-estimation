
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
#a matrix that contains per tree and node one subsampled observation
get_values_nodes <- \(.forest) {
  nodesX <- sapply(.forest, \(tree) tree$nodesX)
  rownames(nodesX) <- NULL

  #size of the largest tree
  nnodes <- max(nodesX)
  ntree <- ncol(nodesX)

  # matrix each row is a node, each column a tree
  valuesNodes  <- matrix(nrow = nnodes, ncol = ntree)

  # take outcome from the first tree (all the ame)
  y <- .forest[[1]]$y[, 2]
  n <- length(y)

  for (j in 1:ntree){
    #useNodes <- sort(unique(as.numeric(shuffledNodes)))
    useNodes <- sort(unique(as.numeric(nodesX[, j])))
    valuesNodes[useNodes, j] <-
      sapply(useNodes, \(k) {
        sample(y, 1, prob = as.numeric(nodesX[, j] == k) * .forest[[j]]$wts)
      })
    valuesNodes[useNodes, j]
    # match return position of (first) match
    #shuffledNodes <-  nodesX[rank(ind <- sample(1:n, n)), j]
    #valuesNodes[useNodes, j] <- y[ind[match(useNodes, shuffledNodes)]]
    #valuesNodes[useNodes, j]
  }
  return(valuesNodes)
}

# I think this is the BOP?
predict_interval_forest <- \(.forest, newdata) {
  predictNodes <- sapply(.forest, \(tree) prednode(tree, newdata))
  rownames(predictNodes) <- NULL
  valuesPredict <- 0 * predictNodes
  valuesNodes <- get_values_nodes(.forest)
  ntree <- ncol(valuesNodes)
  for (tree in 1:ntree){
    valuesPredict[, tree] <- valuesNodes[predictNodes[, tree], tree]
  }
  ci <- apply(valuesPredict, 1, quantile, c(.025, .975), na.rm = TRUE) %>% t
  colnames(ci) <- NULL
  sd <- apply(valuesPredict, 1, sd, na.rm = TRUE)
  cbind(sd, ci)
}

prednode <- \(tree, .data) {
  newdat <- model.frame(delete.response(terms(tree)), .data)
  newmat <- rpart:::rpart.matrix(newdat)
  rpart:::pred.rpart(tree, newmat)
}


#package randomForestSRC
#trtf
#quantregForest
#qrf quantile_forest
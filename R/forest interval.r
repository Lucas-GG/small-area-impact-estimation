# I think this is the BOP?
predict_interval_forest <- \(.forest, .dt) {
  if (!"wts" %in% names(.dt)) .dt[, wts := 1]
  predictNodes <- sapply(.forest, \(tree) prednode(tree, .dt))
  dim(predictNodes)
  rownames(predictNodes) <- NULL
  #dim(predictNodes)
  #summary(predictNodes)
  valuesPredict <- 0 * predictNodes
  #dim(valuesPredict) #ntrree by nrow(newdata)

  valuesNodes <- get_values_nodes(.forest)
  ntree <- length(valuesNodes)

  for (tree in 1:ntree){
    valuesPredict[, tree] <- valuesNodes[[tree]][paste0(predictNodes[, tree])]
  }

  #summary(valuesPredict)

  #ci <- apply(valuesPredict, 1, quantile, c(.025, .975), na.rm = TRUE) %>% t
  #colnames(ci) <- NULL
  #sd <- apply(valuesPredict, 1, sd, na.rm = TRUE)
  #cbind(sd, ci)
  valuesPredict
}

#This helper function determines which node a new observation falls into
#for a given tree
prednode <- \(tree, .dt) {
  model_frame <- model.frame(delete.response(terms(tree))
    , .dt, na.action = na.pass
  )
  # dim(model_frame)
  # Converts to the matrix format expected by `rpart
  model_mat <- rpart:::rpart.matrix(model_frame)
  # dim(model_mat)
  # Uses rpart's internal prediction function to find terminal nodes
  rpart:::pred.rpart(tree, model_mat)
}

# tree <- .forest[[1]]
# dt %>% add_prop %>% mutate(wts = 1) %>% summary
# prednode(.forest[[1]], dt %>% add_prop %>% mutate(wts = 1)) %>% length


#a sample of a single observation per node
get_values_nodes <- function(.forest) {
  nodesx <- lapply(.forest, function(tree) tree$where)
  ntree <- length(nodesx)

  # Create a list with ntree elements
  valuesNodes <- vector("list", ntree)

  for (j in 1:ntree) {
    # Take outcome from the j-th tree
    y_j <- .forest[[j]]$y[, 2]
    n_j <- length(y_j)

    # Create a random permutation of observations
    ind <- sample(1:n_j, n_j)

    # Shuffle node assignments according to permutation
    shuffledNodes <- nodesx[[j]][ind]

    # Get unique node IDs
    useNodes <- sort(unique(as.numeric(shuffledNodes)))

    # For each node, take the first observation in permutation 
    # that falls into it
    sampled_values <- y_j[ind[match(useNodes, shuffledNodes)]]

    # Create a named vector
    names(sampled_values) <- useNodes

    # Store in the list
    valuesNodes[[j]] <- sampled_values
  }

  valuesNodes
}


#newdata <- dt %>% add_prop %>% mutate(wts = 1)
#.forest <- forest(dt %>% add_prop)


#a matrix that contains per tree and node one subsampled observation
#get_values_nodes(.forest) %>% length
#prednode(.forest[[1]], dt %>% add_prop %>% mutate(wts = 1))

#predict_interval_forest(.forest, newdata) %>% head

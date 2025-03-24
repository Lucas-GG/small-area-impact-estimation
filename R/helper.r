# Helper function for building a single tree in parallel
build_single_tree <- function(i, params) {
  # Extract parameters
  data <- params$data
  k <- params$k
  min_bucket <- params$min_bucket
  max_depth <- params$max_depth
  mtry <- params$mtry
  exvar <- params$exvar
  wtype <- params$wtype
  shuffle_time <- params$shuffle_time
  
  # Process the tree
  data <- shuffle_start_cpp(data)
  data <- expandx_cpp(data)
  
  # Filter entries
  data <- subset(data, year < start_year)
  
  if (exvar) {
    data <- cluster_boot_cpp(data, "i", wtype, shuffle_time)
  }
  
  # Remove y and y25 if they exist
  data$y <- NULL
  data$y25 <- NULL
  
  data <- rselect_cpp(data, mtry)
  tree <- fit_single_tree_cpp(data, k, min_bucket, max_depth)
  
  return(tree)
}

#' Process a forest in parallel
#' @param data DataFrame to process
#' @param ntrees Number of trees to create
#' @param k Parameter for trees
#' @param min_bucket Minimum bucket size
#' @param max_depth Maximum tree depth
#' @param mtry Number of variables to try per split
#' @param ncores Number of cores to use
#' @param exvar Whether to use exvar
#' @param wtype Weight type
#' @param shuffle_time Whether to shuffle time
#' @return List of trees
process_forest_parallel <- function(data, ntrees, k, min_bucket, max_depth, mtry, ncores, 
                                   exvar, wtype, shuffle_time) {
  # Ensure we stay within system limits
  ncores <- min(ncores, parallel::detectCores())
  
  # Create a function to process a tree index
  process_tree <- function(i) {
    # Call the C++ function to process a single tree
    process_single_tree_cpp(data, k, min_bucket, max_depth, mtry, exvar, wtype, shuffle_time)
  }
  
  # Use mclapply to process in parallel
  parallel::mclapply(1:ntrees, process_tree, mc.cores = ncores)
}
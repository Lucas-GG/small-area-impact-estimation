// forest.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List forest_cpp(DataFrame data, 
               int ntrees = 20,
               double k = 10, 
               int min_bucket = 7, 
               int max_depth = 30,
               Nullable<int> mtry = R_NilValue,
               bool exvar = true,
               String wtype = "exp",
               bool shuffle_time = false) {
  
  // Create vector to store tree results
  List forest(ntrees);
  
  Environment globalEnv = Environment::global_env();
  Function shuffle_start = globalEnv["shuffle_start"];
  Function expandx = globalEnv["expandx"];
  Function cluster_boot = globalEnv["cluster_boot"];
  Function fit_single_tree = globalEnv["fit_single_tree"];
  Function rselect = globalEnv["rselect"];
  
  // Process trees sequentially (for stability)
  for (int i = 0; i < ntrees; i++) {
    if (i % 5 == 0) Rcpp::Rcout << "Processing tree " << (i+1) << " of " << ntrees << std::endl;
    
    // Make a copy of data for this tree
    DataFrame tree_data = clone(data);
    
    // Apply transformations
    tree_data = shuffle_start(tree_data);
    tree_data = expandx(tree_data);
    
    // Filter data
    NumericVector year = tree_data["year"];
    NumericVector start_year = tree_data["start_year"];
    LogicalVector keep(year.size());
    
    for (int j = 0; j < year.size(); j++) {
      keep[j] = year[j] < start_year[j];
    }
    
    Environment base = Environment::base_env();
    Function subset = base["subset"];
    tree_data = subset(tree_data, keep);
    
    // Apply bootstrap if needed
    if (exvar) {
      tree_data = cluster_boot(tree_data, "i", wtype, shuffle_time);
    }
    
    // Remove response variables 
    if (tree_data.containsElementNamed("y")) {
      // Use R's way of removing columns
      Function eval = base["eval"];
      Function parse = base["parse"];
      
      std::string removeCode = "function(df) { df$y <- NULL; return(df); }";
      Function remove_y = eval(parse(Named("text", removeCode)));
      tree_data = remove_y(tree_data);
    }
    
    if (tree_data.containsElementNamed("y25")) {
      // Use R's way of removing columns
      Function eval = base["eval"];
      Function parse = base["parse"];
      
      std::string removeCode = "function(df) { df$y25 <- NULL; return(df); }";
      Function remove_y25 = eval(parse(Named("text", removeCode)));
      tree_data = remove_y25(tree_data);
    }
    
    // Select variables for this tree
    tree_data = rselect(tree_data, mtry);
    
    // Fit tree
    forest[i] = fit_single_tree(tree_data, k, min_bucket, max_depth);
  }
  
  return forest;
}
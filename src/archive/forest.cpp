#include <Rcpp.h>
#include <omp.h>  // For OpenMP functions
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

// Forward declarations of functions from other files
// Add these declarations to make the compiler aware of these functions
DataFrame shuffle_start_cpp(DataFrame data, int a = 3);
DataFrame expandx_cpp(DataFrame data, int nlags = 0, bool y25 = true, bool recent = false);
DataFrame cluster_boot_cpp(DataFrame data, String cl = "i", String wtype = "mammen", bool shuffle_time = true);
DataFrame rselect_cpp(DataFrame data, Nullable<int> mtry = R_NilValue, Nullable<CharacterVector> yvar = R_NilValue);
List fit_single_tree_cpp(DataFrame data, double k = 10, int min_bucket = 7, int max_depth = 30);

// [[Rcpp::export]]
List forest_rcpp(DataFrame data, 
  int ntrees = 20,
  double k = 10, 
  int min_bucket = 7, 
  int max_depth = 30,
  Nullable<int> mtry = R_NilValue,
  int ncores = 1,
  bool exvar = true,
  String wtype = "exp",
  bool shuffle_time = false) {

// Create list to store trees
List trees(ntrees);

// Process one tree at a time (parallel processing disabled for debugging)
for (int i = 0; i < ntrees; i++) {
Rcpp::Rcout << "Processing tree " << (i+1) << " of " << ntrees << std::endl;

// Make a copy of data for this iteration
DataFrame thread_data = clone(data);

// Apply shuffle_start
thread_data = shuffle_start_cpp(thread_data);

// Apply expandx
thread_data = expandx_cpp(thread_data);

// Filter using R's subset function instead of dplyr::filter
Environment base = Environment::base_env();
Function subset = base["subset"];

// Extract columns with proper typing
NumericVector year = as<NumericVector>(thread_data["year"]);
NumericVector start_year = as<NumericVector>(thread_data["start_year"]);

// Create logical vector for filtering
LogicalVector keep_rows = year < start_year;

// Use the logical vector with subset
thread_data = subset(thread_data, keep_rows);

// Apply cluster_boot if exvar is true
if (exvar) {
thread_data = cluster_boot_cpp(thread_data, "i", wtype, shuffle_time);
}

// Remove y and y25 columns if they exist
if (thread_data.containsElementNamed("y")) {
thread_data["y"] = R_NilValue;
}
if (thread_data.containsElementNamed("y25")) {
thread_data["y25"] = R_NilValue;
}

// Apply feature selection
thread_data = rselect_cpp(thread_data, mtry);

// Fit a single tree
List tree = fit_single_tree_cpp(thread_data, k, min_bucket, max_depth);

// Store the tree in the list
trees[i] = tree;
}

return trees;
}

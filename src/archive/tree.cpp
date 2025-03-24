#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List fit_single_tree_cpp(DataFrame data, double k = 10, int min_bucket = 7, int max_depth = 30) {
  // Simply call the R version directly
  Function fit_tree = Environment::global_env()["fit_single_tree"];
  
  // Call the R function and return its result
  List result = fit_tree(data, k, min_bucket, max_depth);
  
  return result;
}
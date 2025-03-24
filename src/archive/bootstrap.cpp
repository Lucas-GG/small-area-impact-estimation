#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame shuffle_cl_cpp(DataFrame data, String cl = "i", String wtype = "mammen") {
  // Make a copy to avoid modifying the original
  DataFrame df = clone(data);
  
  // Add weights if they don't exist
  if (!df.containsElementNamed("wts")) {
    NumericVector ones(df.nrows(), 1.0);
    df["wts"] = ones;
  }
  
  // Get the cluster column
  if (!df.containsElementNamed(cl.get_cstring())) {
    stop("Clustering column '%s' not found", cl.get_cstring());
  }
  
  // Store cluster column in a temporary column named 'cl'
  df["cl"] = df[cl];
  
  // Call R functions directly for compatibility
  Environment base("package:base");
  Environment dplyr = Environment::namespace_env("dplyr");
  
  // Get unique clusters
  Function unique = base["unique"];
  SEXP unique_clusters = unique(df["cl"]);
  
  // Get count of distinct clusters
  Function n_distinct = dplyr["n_distinct"];
  int n_cl = as<int>(n_distinct(df["cl"]));
  
  // Get the weight generator from fwb
  Environment fwb = Environment::namespace_env("fwb");
  Function make_gen_weights = fwb["make_gen_weights"];
  Function gen_w = make_gen_weights(wtype);
  
  // Generate weights
  Function as_vector = base["as.vector"];
  NumericVector new_wts = as<NumericVector>(as_vector(gen_w(n_cl, 1)));
  
  // Create data frame with clusters and new weights
  DataFrame wdata = DataFrame::create(
    Named("cl") = unique_clusters,
    Named("new_wts") = new_wts
  );
  
  // Left join using R's merge function (more reliable than trying to do it in C++)
  Function merge = base["merge"];
  df = as<DataFrame>(merge(df, wdata, 
                          Named("by") = "cl", 
                          Named("all.x") = true,
                          Named("sort") = false));
  
  // Update weights: wts = wts * new_wts
  NumericVector weights = df["wts"];
  NumericVector new_weights = df["new_wts"];
  int n = weights.size();
  
  for (int i = 0; i < n; i++) {
    weights[i] = weights[i] * new_weights[i];
  }
  
  // Remove new_wts column
  CharacterVector names = df.names();
  CharacterVector new_names;
  std::vector<SEXP> new_cols;
  
  for (int i = 0; i < names.size(); i++) {
    if (std::string(names[i]) != "new_wts") {
      new_names.push_back(names[i]);
      new_cols.push_back(df[i]);
    }
  }
  
  // Reconstruct data frame without new_wts
  DataFrame result = DataFrame::create();
  for (int i = 0; i < new_names.size(); i++) {
    result.push_back(new_cols[i], std::string(new_names[i]));
  }
  
  // Normalize weights
  weights = result["wts"];
  double sum_wts = 0.0;
  for (int i = 0; i < n; i++) {
    sum_wts += weights[i];
  }
  
  for (int i = 0; i < n; i++) {
    weights[i] = weights[i] / sum_wts * n;
  }
  
  // Remove cl column
  names = result.names();
  new_names = CharacterVector();
  new_cols.clear();
  
  for (int i = 0; i < names.size(); i++) {
    if (std::string(names[i]) != "cl") {
      new_names.push_back(names[i]);
      new_cols.push_back(result[i]);
    }
  }
  
  // Final data frame
  DataFrame final_result = DataFrame::create();
  for (int i = 0; i < new_names.size(); i++) {
    final_result.push_back(new_cols[i], std::string(new_names[i]));
  }
  
  return final_result;
}

// [[Rcpp::export]]
DataFrame cluster_boot_cpp(DataFrame data, String cl = "i", 
                           String wtype = "mammen", bool shuffle_time = true) {
  DataFrame result = clone(data);
  
  if (shuffle_time) {
    result = shuffle_cl_cpp(result, "year", wtype);
  }
  
  result = shuffle_cl_cpp(result, cl, wtype);
  
  return result;
}
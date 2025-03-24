#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame rselect_cpp(DataFrame data, Nullable<int> mtry = R_NilValue, Nullable<CharacterVector> yvar = R_NilValue) {
  // If mtry is NULL, return the original data
  if (mtry.isNull()) {
    return data;
  }
  
  int mtry_value = as<int>(mtry);
  
  // If mtry is 0, set it to floor(sqrt(ncol(data)))
  if (mtry_value == 0) {
    mtry_value = floor(sqrt(data.ncol()));
  }
  
  // Set default yvar if NULL
  CharacterVector y_vars;
  if (yvar.isNull()) {
    y_vars = CharacterVector::create("n", "y0", "wts", "y0_hat");
  } else {
    y_vars = as<CharacterVector>(yvar);
  }
  
  // Get all column names
  CharacterVector all_names = data.names();
  
  // Identify which columns are y variables
  LogicalVector is_y_var(all_names.size(), false);
  for (int i = 0; i < all_names.size(); i++) {
    String col_name = all_names[i];
    for (int j = 0; j < y_vars.size(); j++) {
      if (col_name == y_vars[j]) {
        is_y_var[i] = true;
        break;
      }
    }
  }
  
  // Create vectors for y and x column names
  CharacterVector y_names;
  CharacterVector x_names;
  
  for (int i = 0; i < all_names.size(); i++) {
    if (is_y_var[i]) {
      y_names.push_back(all_names[i]);
    } else {
      x_names.push_back(all_names[i]);
    }
  }
  
  // Fix: Convert to same type before using std::min
  int x_names_size = static_cast<int>(x_names.size());
  int x_sample_size = std::min(mtry_value, x_names_size);
  
  // Randomly sample x_sample_size columns from x_names
  CharacterVector sampled_x_names;
  
  // Create a vector of indices and shuffle it
  IntegerVector indices = seq_len(x_names.size()) - 1; // 0-based
  indices = sample(indices, x_names.size(), false);
  
  // Take the first x_sample_size indices
  for (int i = 0; i < x_sample_size; i++) {
    sampled_x_names.push_back(x_names[indices[i]]);
  }
  
  // Combine y_names and sampled_x_names
  CharacterVector result_names;
  for (int i = 0; i < y_names.size(); i++) {
    result_names.push_back(y_names[i]);
  }
  for (int i = 0; i < sampled_x_names.size(); i++) {
    result_names.push_back(sampled_x_names[i]);
  }
  
  // Create a new DataFrame with only the selected columns
  DataFrame result = DataFrame::create();
  for (int i = 0; i < result_names.size(); i++) {
    String col_name = result_names[i];
    if (data.containsElementNamed(col_name.get_cstring())) {
      result.push_back(data[col_name], col_name);
    }
  }
  
  // Return the result
  return result;
}
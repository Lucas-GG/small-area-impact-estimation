#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
#include <numeric>
#include <algorithm>
#include <map>
#include <string>
#include <vector>


// Include all necessary Rcpp plugins
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//' @export
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

//' @export
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

//' @export
// [[Rcpp::export]]
DataFrame shuffle_start_cpp(DataFrame data, int a = 3) {
  // Ensure we have copies of the data to avoid modifying the original
  DataFrame df = clone(data);
  
  // Step 1: Modify start_year by subtracting a
  NumericVector start_year = df["start_year"];
  for (int i = 0; i < start_year.size(); i++) {
    if (!NumericVector::is_na(start_year[i]) && start_year[i] != R_PosInf) {
      start_year[i] -= a;
    }
  }
  
  // Get other columns we need
  NumericVector ids = df["i"];
  NumericVector years = df["year"];
  NumericVector st = df["st"];
  
  // Check if cpr column exists, if not add dummy values
  NumericVector cpr;
  if (df.containsElementNamed("cpr")) {
    cpr = df["cpr"];
  } else {
    cpr = NumericVector(ids.size(), 0.5);  // Default equal probability
  }
  
  // Create unique id-start_year-st combinations to mimic group_by & reframe
  std::map<std::pair<double, double>, double> uniqueStartYears;
  for (int i = 0; i < ids.size(); i++) {
    if (!NumericVector::is_na(start_year[i]) && start_year[i] != R_PosInf) {
      uniqueStartYears[std::make_pair(ids[i], st[i])] = start_year[i];
    }
  }
  
  // Convert unique combinations to a data frame
  NumericVector gdata_i;
  NumericVector gdata_start_year;
  
  for (auto const& entry : uniqueStartYears) {
    gdata_i.push_back(entry.first.first);
    gdata_start_year.push_back(entry.second);
  }
  
  // *** SAFETY CHECK ***
  // If we have too few unique IDs/start years, use a different approach
  if (gdata_i.size() < 5) {
    Rcpp::warning("Too few unique ID-start year combinations (%d). Using original start years.", gdata_i.size());
    return df; // Return original data with modified start years
  }
  
  // Count occurrences of each start_year
  std::map<double, int> startYearCounts;
  for (int i = 0; i < gdata_start_year.size(); i++) {
    startYearCounts[gdata_start_year[i]]++;
  }
  
  // Convert to R-style table and sort keys
  std::vector<double> start_year_values;
  std::vector<int> counts;
  
  for (auto const& entry : startYearCounts) {
    start_year_values.push_back(entry.first);
    counts.push_back(entry.second);
  }
  
  // Sort by start_year
  std::vector<size_t> indices(start_year_values.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&start_year_values](size_t a, size_t b) {
              return start_year_values[a] < start_year_values[b];
            });
  
  std::vector<double> sorted_start_years;
  std::vector<int> sorted_counts;
  
  for (size_t i = 0; i < indices.size(); i++) {
    sorted_start_years.push_back(start_year_values[indices[i]]);
    sorted_counts.push_back(counts[indices[i]]);
  }
  
  // Remove the last element as in R code: v0 <- v0[-length(v0)]
  if (!sorted_start_years.empty()) {
    sorted_start_years.pop_back();
    sorted_counts.pop_back();
  }
  
  // Years range 2006:2016
  std::vector<int> year_range;
  for (int yr = 2006; yr <= 2016; yr++) {
    year_range.push_back(yr);
  }
  
  // Equivalent to the lapply and accumulation of s
  std::set<double> selected_ids;
  std::vector<double> all_new_s;
  std::vector<double> all_new_start_years;
  
  for (size_t j = 0; j < sorted_counts.size(); j++) {
    int count_to_sample = sorted_counts[j];
    int year_j = year_range[j];
    
    // Filter data equivalent to: dt <- .data[(!(.data$i %in% s)) & .data$year == ((2006:2016))[j], ]
    std::vector<double> eligible_ids;
    std::vector<double> eligible_probs;
    
    for (int i = 0; i < ids.size(); i++) {
      if (selected_ids.find(ids[i]) == selected_ids.end() && years[i] == year_j) {
        eligible_ids.push_back(ids[i]);
        eligible_probs.push_back(cpr[i]);
      }
    }
    
    // Only proceed if we have eligible ids
    if (!eligible_ids.empty()) {
      // Sample from eligible ids with probability cpr
      NumericVector r_eligible_ids = wrap(eligible_ids);
      NumericVector r_eligible_probs = wrap(eligible_probs);
      
      // Handle case where we need more samples than available
      int n_samples = std::min(count_to_sample, (int)eligible_ids.size());
      
      if (n_samples > 0) {
        NumericVector sampled = sample(r_eligible_ids, n_samples, false, r_eligible_probs);
        
        // Add to our tracking variables
        for (int k = 0; k < sampled.size(); k++) {
          double id = sampled[k];
          selected_ids.insert(id);
          all_new_s.push_back(id);
          all_new_start_years.push_back(sorted_start_years[j]);
        }
      }
    }
  }
  
  // *** SAFETY CHECK ***
  // If we have too few sampled rows, use a different approach
  if (all_new_s.size() < 5) {
    Rcpp::warning("Too few samples generated (%d). Using original data with modified start years.", all_new_s.size());
    return df;
  }
  
  // Create new data frame with i and start_year columns
  DataFrame sl = DataFrame::create(
    Named("i") = wrap(all_new_s),
    Named("start_year") = wrap(all_new_start_years)
  );
  
  // Now perform the left join operation
  // First remove start_year from original data
  df["start_year"] = R_NilValue;
  
  // Create map from id to new start_year
  std::map<double, double> id_to_start_year;
  NumericVector sl_ids = sl["i"];
  NumericVector sl_start_years = sl["start_year"];
  
  for (int i = 0; i < sl_ids.size(); i++) {
    id_to_start_year[sl_ids[i]] = sl_start_years[i];
  }
  
  // Create new start_year vector for the result
  NumericVector new_start_years(ids.size());
  
  for (int i = 0; i < ids.size(); i++) {
    double id = ids[i];
    if (id_to_start_year.find(id) != id_to_start_year.end()) {
      new_start_years[i] = id_to_start_year[id];
    } else {
      new_start_years[i] = R_PosInf; // Inf in R
    }
  }
  
  // Make a clean copy to be safe
  DataFrame result = clone(df);
  
  // Add the new start_year to the result
  result["start_year"] = new_start_years;
  
  return result;
}

//' @export
// [[Rcpp::export]]
NumericVector gmean_cpp(DataFrame data, String x_name, String f_name, 
                       bool na_rm = true, bool fill = true) {
    NumericVector x = data[x_name];
    NumericVector f = data[f_name];
    
    // Create map of group means
    std::map<double, std::pair<double, int>> group_sums;
    
    // Calculate sum and count for each group
    for (int i = 0; i < x.size(); i++) {
        double f_val = f[i];
        double x_val = x[i];
        
        if (na_rm && NumericVector::is_na(x_val)) continue;
        
        if (group_sums.find(f_val) == group_sums.end()) {
            group_sums[f_val] = std::make_pair(0.0, 0);
        }
        
        group_sums[f_val].first += x_val;
        group_sums[f_val].second += 1;
    }
    
    // Calculate means
    std::map<double, double> group_means;
    for (auto& kv : group_sums) {
        if (kv.second.second > 0) {
            group_means[kv.first] = kv.second.first / kv.second.second;
        } else {
            group_means[kv.first] = NA_REAL;
        }
    }
    
    // Apply means to result vector
    NumericVector result(x.size());
    
    for (int i = 0; i < f.size(); i++) {
        double f_val = f[i];
        
        if (group_means.find(f_val) != group_means.end()) {
            result[i] = group_means[f_val];
        } else if (fill) {
            result[i] = NA_REAL;
        } else {
            result[i] = 0.0;
        }
    }
    
    return result;
}

//' @export
// [[Rcpp::export]]
DataFrame reset_y0_cpp(DataFrame data) {
    DataFrame result = clone(data);
    
    NumericVector year = result["year"];
    NumericVector start_year = result["start_year"];
    NumericVector y0 = result["y0"];
    
    // Calculate event_time
    NumericVector event_time(year.size());
    NumericVector y0_modified(y0.size());
    
    for (int i = 0; i < year.size(); i++) {
        event_time[i] = year[i] - start_year[i];
        
        // Set y0 to NA if year >= start_year
        if (year[i] < start_year[i]) {
            y0_modified[i] = y0[i];
        } else {
            y0_modified[i] = NA_REAL;
        }
    }
    
    // Update result dataframe
    result["event_time"] = event_time;
    result["y0"] = y0_modified;
    
    // Check if y0_25 exists
    if (result.containsElementNamed("y0_25")) {
        NumericVector y0_25 = result["y0_25"];
        NumericVector y0_25_modified(y0_25.size());
        
        for (int i = 0; i < y0_25.size(); i++) {
            if (year[i] < start_year[i]) {
                y0_25_modified[i] = y0_25[i];
            } else {
                y0_25_modified[i] = NA_REAL;
            }
        }
        
        result["y0_25"] = y0_25_modified;
    }
    
    return result;
}

//' @export
// [[Rcpp::export]]
DataFrame expandx_cpp(DataFrame data, 
                     int nlags = 0, 
                     bool y25 = true, 
                     bool recent = false) {
    // Reset y0 values
    DataFrame result = reset_y0_cpp(data);
    
    // Calculate various group means
    result["m_i"] = gmean_cpp(result, "y0", "i");
    result["m_st"] = gmean_cpp(result, "y0", "st");
    result["m_t"] = gmean_cpp(result, "y0", "year");
    result["m_u"] = gmean_cpp(result, "y0", "urb");
    result["m_c"] = gmean_cpp(result, "y0", "start_year");
    result["m_tr"] = gmean_cpp(result, "y0", "TribalStatus");
    
    // Handle lags if needed - calling back to R for add_hist
    if (nlags == R_PosInf) {
        // Fixed: Use R's length() function to get size
        Environment baseEnv = Environment::base_env();
        Function unique = baseEnv["unique"];
        Function length = baseEnv["length"];
        
        NumericVector years = result["year"];
        SEXP unique_years = unique(years);
        int n_years = as<int>(length(unique_years));
        nlags = n_years - 1;
        
        // Call R function add_hist
        Function add_hist = Environment::global_env()["add_hist"];
        result = add_hist(result, nlags);
    } else if (nlags > 0) {
        Function add_hist = Environment::global_env()["add_hist"];
        result = add_hist(result, nlags);
    }
    
    // Handle recent data if requested
    if (recent) {
        // Create y0_l10 column
        NumericVector year = result["year"];
        NumericVector start_year = result["start_year"];
        NumericVector y0 = result["y0"];
        NumericVector y0_l10(y0.size(), NA_REAL);
        
        for (int i = 0; i < y0.size(); i++) {
            if (year[i] > start_year[i] - 10) {
                y0_l10[i] = y0[i];
            }
        }
        
        result["y0_l10"] = y0_l10;
        
        // Calculate group means for recent data
        result["ml_i"] = gmean_cpp(result, "y0_l10", "i");
        result["ml_st"] = gmean_cpp(result, "y0_l10", "st");
        result["ml_u"] = gmean_cpp(result, "y0_l10", "urb");
        result["ml_c"] = gmean_cpp(result, "y0_l10", "start_year");
    }
    
    // Handle 25-year averages if requested
    if (y25 && result.containsElementNamed("y0_25")) {
        result["m25_i"] = gmean_cpp(result, "y0_25", "i");
        result["m25_st"] = gmean_cpp(result, "y0_25", "st");
        result["m25_t"] = gmean_cpp(result, "y0_25", "year");
        result["m25_u"] = gmean_cpp(result, "y0_25", "urb");
        result["m25_c"] = gmean_cpp(result, "y0_25", "start_year");
        result["m25_tr"] = gmean_cpp(result, "y0_25", "TribalStatus");
    }
    
    return as<DataFrame>(result);
}

//' @export
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


//' @export
// [[Rcpp::export]]
List fit_single_tree_cpp(DataFrame data, double k = 10, int min_bucket = 7, int max_depth = 30) {
  // Simply call the R version directly
  Function fit_tree = Environment::global_env()["fit_single_tree"];
  
  // Call the R function and return its result
  List result = fit_tree(data, k, min_bucket, max_depth);
  
  return result;
}

//' Process a single forest tree
//' @export
// [[Rcpp::export]]
List process_single_tree_cpp(DataFrame data, 
                           double k = 10, 
                           int min_bucket = 7, 
                           int max_depth = 30,
                           Nullable<int> mtry = R_NilValue,
                           bool exvar = true,
                           String wtype = "exp",
                           bool shuffle_time = false) {
  
  // Process a single tree
  DataFrame thread_data = clone(data);
  thread_data = shuffle_start_cpp(thread_data);
  thread_data = expandx_cpp(thread_data);
  
  // Filter entries
  NumericVector year = thread_data["year"];
  NumericVector start_year = thread_data["start_year"];
  LogicalVector keep_rows(year.size());
  
  for (int j = 0; j < year.size(); j++) {
    keep_rows[j] = year[j] < start_year[j];
  }
  
  // Use subset from base R
  Environment base = Environment::base_env();
  Function subset = base["subset"];
  thread_data = subset(thread_data, keep_rows);
  
  if (exvar) {
    thread_data = cluster_boot_cpp(thread_data, "i", wtype, shuffle_time);
  }
  
  // Remove y and y25 if they exist
  if (thread_data.containsElementNamed("y")) {
    thread_data["y"] = R_NilValue;
  }
  if (thread_data.containsElementNamed("y25")) {
    thread_data["y25"] = R_NilValue;
  }
  
  thread_data = rselect_cpp(thread_data, mtry);
  return fit_single_tree_cpp(thread_data, k, min_bucket, max_depth);
}

// TreeProcessor: Efficient tree processing in parallel
class TreeProcessor : public Worker {
private:
    const DataFrame& data;          // Reference to input data (no copying)
    double k;
    int min_bucket;
    int max_depth;
    Nullable<int> mtry;
    bool exvar;
    std::string wtype;
    bool shuffle_time;
    std::vector<List>& results;    // Reference to output results
    
public:
    // Constructor
    TreeProcessor(const DataFrame& data, 
                 double k, int min_bucket, int max_depth, 
                 Nullable<int> mtry, bool exvar, 
                 std::string wtype, bool shuffle_time,
                 std::vector<List>& results) 
        : data(data), k(k), min_bucket(min_bucket), max_depth(max_depth),
          mtry(mtry), exvar(exvar), wtype(wtype), shuffle_time(shuffle_time),
          results(results) {}
    
    // Process a range of indices
    void operator()(std::size_t begin, std::size_t end) {
        try {
            // Setup R function environment once per thread
            Environment base = Environment::base_env();
            Function subset = base["subset"];
            
            for (std::size_t i = begin; i < end; i++) {
                // Process tree i
                DataFrame thread_data = clone(data);
                
                // 1. Shuffle start years
                thread_data = shuffle_start_cpp(thread_data);
                
                // 2. Expand features
                thread_data = expandx_cpp(thread_data);
                
                // 3. Filter rows
                NumericVector year = thread_data["year"];
                NumericVector start_year = thread_data["start_year"];
                LogicalVector keep_rows(year.size(), false);
                
                for (int j = 0; j < year.size(); j++) {
                    if (year[j] < start_year[j]) {
                        keep_rows[j] = true;
                    }
                }
                
                // Use subset from R
                thread_data = subset(thread_data, keep_rows);
                
                // 4. Apply bootstrap if needed
                if (exvar) {
                    thread_data = cluster_boot_cpp(thread_data, "i", wtype, shuffle_time);
                }
                
                // 5. Remove target variables if present
                if (thread_data.containsElementNamed("y")) {
                    thread_data["y"] = R_NilValue;
                }
                if (thread_data.containsElementNamed("y25")) {
                    thread_data["y25"] = R_NilValue;
                }
                
                // 6. Select variables
                thread_data = rselect_cpp(thread_data, mtry);
                
                // 7. Fit tree
                results[i] = fit_single_tree_cpp(thread_data, k, min_bucket, max_depth);
            }
        } catch(std::exception& e) {
            Rcpp::stop(e.what());
        } catch(...) {
            Rcpp::stop("Unknown error occurred in parallel processing");
        }
    }
};

//' @export
// [[Rcpp::export]]
List forest_rcpp_parallel(DataFrame data, 
                         int ntrees = 20,
                         double k = 10, 
                         int min_bucket = 7, 
                         int max_depth = 30,
                         Nullable<int> mtry = R_NilValue,
                         int ncores = 1,
                         bool exvar = true,
                         String wtype = "exp",
                         bool shuffle_time = false) {
    
    // Prepare output vector
    std::vector<List> results(ntrees);
    
    // Prepare worker
    TreeProcessor processor(data, k, min_bucket, max_depth, 
                           mtry, exvar, wtype.get_cstring(), shuffle_time, 
                           results);
    
    // Execute - parallel or sequential based on ncores
    if (ncores <= 1) {
        processor(0, ntrees);  // Run sequentially
    } else {
        // For maximum safety, just use ncores as provided or query from R
        if (ncores > 1) {
            // Get available cores from R's parallel package
            Environment parallel = Environment::namespace_env("parallel");
            Function detectCores = parallel["detectCores"];
            int max_cores = as<int>(detectCores());
            
            if (ncores > max_cores) {
                ncores = max_cores;
            }
        }
        
        // Run in parallel with the determined number of cores
        parallelFor(0, ntrees, processor, ncores);
    }
    
    // Convert results to R list
    List r_results(ntrees);
    for (int i = 0; i < ntrees; i++) {
        r_results[i] = results[i];
    }
    
    return r_results;
}
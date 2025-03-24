#include <Rcpp.h>
#include <numeric>

using namespace Rcpp;

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
  
  // *** SAFETY CHECK ***
  // Make sure we have enough years for sampling
  if (sorted_start_years.empty() || sorted_counts.empty() || year_range.empty()) {
    Rcpp::warning("Insufficient data for sampling start years. Using original data with modified start years.");
    return df;
  }
  
  // Equivalent to the lapply and accumulation of s
  std::set<double> selected_ids;
  std::vector<double> all_new_s;
  std::vector<double> all_new_start_years;
  
  for (size_t j = 0; j < sorted_counts.size() && j < year_range.size(); j++) {
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
  // Create map from id to new start_year
  std::map<double, double> id_to_start_year;
  NumericVector sl_ids = sl["i"];
  NumericVector sl_start_years = sl["start_year"];
  
  // *** SAFETY CHECK ***
  // Ensure sl_ids and sl_start_years have the same length
  if (sl_ids.size() != sl_start_years.size()) {
    Rcpp::warning("ID and start_year vectors have different lengths. Using original data.");
    return df;
  }
  
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



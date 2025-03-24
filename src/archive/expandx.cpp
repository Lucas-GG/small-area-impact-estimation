#include <Rcpp.h>
#include <map>
#include <string>
#include <algorithm>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// Helper function for gmean - calculates group means
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

// Reset y0 values based on event time
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

// Main expandx function 
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
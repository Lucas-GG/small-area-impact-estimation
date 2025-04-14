#------------------------------------------------------------------------------
#' expand feautrues
#------------------------------------------------------------------------------
library(collapse)
library(ranger)

#-------------------------------------------------------------------
gmean <- \(.data, x, f) {
  x <- deparse(substitute(x))
  f <- deparse(substitute(f))
  fbetween(.data[[x]], .data[[f]], na.rm = TRUE, fill = TRUE)
}
#-------------------------------------------------------------------

#' this computes slope but it is too slow
s <- \(x) {
  x <- na.omit(x)
  if (length(x) < 2) return(NA)
  else coef(lm(x ~ I(-seq(x))))[2]
}

reset_y0 <- function(dt) {

  result <- copy(dt)
  # Perform all operations in-place with :=
  result[, event_time := year - start_year]
  result[year >= start_year, `:=`(y0 = NA_real_, y0_25 = NA_real_)]


  # Return the modified data.table
  result
}

add_hist <- \(.data, nlags = length(unique(.data$year)) - 1) {
  .data <- .data %>% reset_y0
  h_y0  <- .data %>% with(flag(y0, 1:nlags, i, year))

  for (i in 1:nlags) {
    set(.data, j = paste0("h_y0_l", i), value = h_y0[, i])
  }
  .data$h_y0 <- rowMeans(h_y0)
  .data$h_mi  <- .data %>% gmean(h_y0, i)
  .data$h_mst <- .data %>% gmean(h_y0, st)
  .data$h_mt  <- .data %>% gmean(h_y0, year)
  .data$h_mu  <- .data %>% gmean(h_y0, urb)
  .data$h_mc  <- .data %>% gmean(h_y0, start_year)
  .data$h_mtr <- .data %>% gmean(h_y0, TribalStatus)

  .data
}

expandx <- \(.data
  , nlags = 0
  , y25 = TRUE
#  , recent = 0
) {
  .data <- .data %>% reset_y0

  .data$m_i  <- .data %>% gmean(y0, i)
  .data$m_st <- .data %>% gmean(y0, st)
  .data$m_t  <- .data %>% gmean(y0, year)
  .data$m_u  <- .data %>% gmean(y0, urb)
  .data$m_c  <- .data %>% gmean(y0, start_year)
  .data$m_tr <- .data %>% gmean(y0, TribalStatus)

  if (nlags == Inf) {
    nlags <- length(unique(.data$year)) - 1
    .data <- .data %>% add_hist(nlags)
  } else if (nlags > 0) {
    .data <- .data %>% add_hist(nlags)
  }

  if (y25) {
    .data$m25_i <- .data %>% gmean(y0_25, i)
    .data$m25_st <- .data %>% gmean(y0_25, st)
    .data$m25_t <-  .data %>%  gmean(y0_25, year)
    .data$m25_u <-  .data %>% gmean(y0_25, urb)
    .data$m25_c <-  .data %>% gmean(y0_25, start_year)
    .data$m25_tr <-  .data %>% gmean(y0_25, TribalStatus)
  }

  .data
}


#dt %>% add_prop

library(xgboost)
#add propensity with optimized data.table operations
add_prop <- function(dt, nlags = 3) {

  # Create a shallow copy instead of deep copy (faster)
  result <- copy(dt) %>% expandx(nlags = nlags)

  sel_col <- grep(c("m_|h_"), names(result), value = TRUE)
  sel_col <- c(sel_col
    , "i", "year", "n", "n25", "uemp", "pov", "mhinc", "lat", "lon"
  )
  # Create the sdt subset using more efficient data.table syntax
  sdt <- result[
    year < start_year & year >= min(year) + nlags
    , c(.(w = as.numeric(year == start_year - 1))
        , .SD)
    , .SDcols = sel_col
  ]
  dim(result)
  dim(sdt)
  summary(sdt)
  with(sdt, tapply(is.na(h_y0_l1), year, sum))

  # Use data.table's efficient set operations for splitting
  set.seed(42)
  # Create watchlist
  watchlist <- get_watchlist(sdt, p = .5)
  # Train with optimized parameters for speed
  params <- list(
    objective = "binary:logistic"
    , eta = 0.1           # Faster learning rate
    , max_depth = 6
    , nthread = 1
   #, subsample = 0.8
   # , colsample_bytree = 0.8  # Use fewer features per tree (faster)
    , tree_method = "exact"#"hist"     # Faster histogram-based algorithm
  )

  mp_xgb <- xgb.train(
    params = params,
    data = watchlist$train,
    watchlist = watchlist,
    nrounds = 200,
    early_stopping_rounds = 5,
    verbose = 0
  )

  # Predict more efficiently - prepare predictors only once
  pred_cols <- names(sdt)[!names(sdt) %in% c("w", "i", "is_train")]

  # Create matrix directly without using model.matrix (more efficient)
  # But we need to ensure the column order matches the training data
  pred_matrix <- model.matrix(~ . - 1
    , data = model.frame(~ ., result[, ..pred_cols], na.action = na.pass)
    , na.action = "na.pass"
  )

  # Use best iteration from early stopping
  best_iter <- mp_xgb$best_iteration

  # Add predictions directly with :=
  result[, cpr := predict(mp_xgb, pred_matrix,
           iteration_range = c(0, best_iter)
         )]

}

# Get watchlist
get_watchlist <- \(sdt, p = .5) {
  # Get unique IDs using data.table's uniqueN (much faster)
  id_dt <- unique(sdt[, .(i)])
  id_dt[, rn := .I]  # Add row numbers

  # Sample row numbers instead of values (more efficient)
  n_train <- as.integer(p * nrow(id_dt))
  train_rows <- sample.int(nrow(id_dt), n_train)

  # Mark training set
  id_dt[, is_train := rn %in% train_rows]

  # Join back efficiently with data.table
  sdt[id_dt, is_train := i.is_train, on = "i"]

  # Split data without creating copies
  # Just use logical indexing when needed rather than creating separate objects

  # XGBoost operation (not data.table specific, but optimized)
  # Convert to matrix format - can't avoid this step with XGBoost
  train_matrix <- model.matrix(w ~ . - 1 - i - is_train
    , data = sdt[is_train == TRUE]
  )
  train_label <- sdt[is_train == TRUE, w]

  valid_matrix <- model.matrix(w ~ . - 1 - i - is_train
    , data = sdt[is_train == FALSE]
  )
  valid_label <- sdt[is_train == FALSE, w]

  # Clean up large objects early
  rm(id_dt)
  gc()  # Force garbage collection to free memory

  dim(train_matrix)
  dim(valid_matrix)
  length(train_label)
  length(valid_label)
  # Create DMatrix objects
  dtrain <- xgb.DMatrix(train_matrix, label = train_label)
  dvalid <- xgb.DMatrix(valid_matrix, label = valid_label)

  # Create watchlist
  watchlist <- list(train = dtrain, valid = dvalid)

  # Clean up more objects
  rm(train_matrix, valid_matrix)
  gc()

  watchlist
}




#add propensity
add_prop_forest <- function(dt) {
  # Create a copy to avoid modifying the original unintentionally
  result <- copy(dt) %>% expandx()
  # Create the sdt subset using data.table operations
  sdt <- result[year < start_year, .(
    w = as.numeric(year == start_year - 1),
    n, n25,
    m_i, m25_i,
    m_st, m25_st,
    m_t, m25_t,
    m_u, m25_u,
    m_tr, m25_tr,
    uemp, pov, mhinc,
    lat, lon
  )]
  summary(sdt)

  # Fit random forest model
  mp <- ranger(w ~ ., data = sdt, probability = TRUE, num.threads = 1)

  # Predict probabilities and add to result
  pred <- predict(mp, result, num.threads = 1)$predictions[, 2]
  result[, cpr := pred]

  # Print summary to make output visible
  cat("add_prop_dt: Added propensity scores using RANGER\n")
  # Return the modified data.table
  result
}
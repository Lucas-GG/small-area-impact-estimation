
#' Shuffle start year
shuffle_start <- function(.dt, lag = 2, rnd = TRUE) {
  if (rnd) a <- sample(0:lag, 1)

  .dt <- copy(.dt)


  # Adjust start_year (in-place modification)
  .dt[, start_year := start_year - a]

  # Get count of unique i values per start_year directly
  v0 <- .dt[, .(count = uniqueN(i)), by = start_year][order(start_year)]

  # If you need it as a named vector like table() would produce
  v0 <- setNames(v0$count, v0$start_year)

  # get rid of last cohort
  v0 <- v0[-length(v0)]

  # Initialize empty result and sampling vector
  s <- NULL
  sl_list <- vector("list", length(v0))

  # Process each year
  years <- as.numeric(names(v0))
  year_range <- 2006:2016

  for (j in seq_len(length(v0))) {
    # Filter data for the current year, excluding already sampled IDs
    current_dt <- .dt[!(i %in% s) & year == year_range[j]]

    # Sample IDs with probability proportional to cpr
    if (nrow(current_dt) > 0) {
      # Get unique i values and their corresponding cpr values
      sample_dt <- unique(current_dt[, .(i, cpr)])

      # Sample IDs
      new_s <- sample(sample_dt$i, min(v0[j], nrow(sample_dt))
                      , prob = sample_dt$cpr, replace = FALSE)

      # Update the overall sample list
      s <- c(s, new_s)

      # Create result for this iteration
      sl_list[[j]] <- data.table(i = new_s, start_year_new = years[j])
    }
  }

  # Combine the results
  sl <- rbindlist(sl_list)

#  # Remove old start_year column and join with the new values
#  .dt[, start_year := NULL]

  # Join using data.table's efficient merge
  .dt <- merge(.dt, sl, by = "i", all.x = TRUE)
  .dt[is.na(start_year_new), start_year_new := Inf]
  
#  # Set missing start_year values to Inf
#  .dt[is.na(start_year), start_year := Inf]

  .dt[, start_year := do.call(pmin, .SD)
    , .SD = c("start_year", "start_year_new")
  ]
  .dt[, start_year_new := NULL]

  return(.dt)
}
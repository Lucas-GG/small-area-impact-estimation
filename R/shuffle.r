
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

#' Defaults to MCAR (completly at random) but MAR can be selected instead
shuffle_start_old <- \(.data, type = "mar", a = 3) {
  .gdata <- .data %>%
    # first things first, recursively change the start year
    mutate(start_year = start_year - a) %>%
    # a historical for before the earliest intervention
    mutate(
      tots = ifelse(year < min(start_year), y + y25, NA) #
      , totn = ifelse(year < min(start_year), n + n25, NA) #
    ) %>%
    group_by(i, start_year, st) %>%
    reframe(
      my = sum(tots, na.rm = TRUE) / sum(totn, na.rm = TRUE) * 10^5
      , start_year = first(start_year)
      , n = mean(totn, na.rm = TRUE)
    )
  #------------------#
  # shuffle start year
  if (type == "mar") {
    .gdata <-  .gdata %>% shuffle_start_mar
  } else if (type == "mcar") {
    .gdata <- .gdata %>% mutate(start_year = sample(start_year))
  } else if (type == "mnar") {
    .gdata <- .gdata %>%
      group_by(st) %>%
      mutate(start_year = sample(start_year)) %>%
      ungroup
  }
  .data <- left_join(
    select(.data, - start_year)
    , select(.gdata, c(i, start_year))
    , by = "i"
  ) %>%
    arrange(i, year) %>%
    mutate(
      y0 = ifelse(year < start_year, y0, NA)
      , y25 = ifelse(year < start_year, y25, NA)
    )
  return(.data)
}

# MAR means missing depends on observables
#' Most dengerous would be dependence on past outcome
#' More relsitic is dependence on size (i.e.,
#' the larger the county, the larger the chances to get at least one training)
#' #' to create a mar selection
#' the procedure computes hstorical suicide rate
#' from 1999 to min(start_year)-1 , i.e., <= 2006-3-1=2002
#' this variable is not saved, so it cannot introduce noise anywhere
#' further, forest define subsamples recursively laggin the start
#' but does not use mar (it uses mcar)
shuffle_start_mar <- \(.gdata) {
  v0  <- table(.gdata$start_year)
  v0 <- v0[-length(v0)]
  s <- NULL
  .gdata$size <- .gdata$n
  s <<- NULL
  sl <- lapply(seq_len(length(v0)), \(j) {
    a <- with(.gdata[!.gdata$i %in% s, ], sample(i, v0[j], prob = size))
    s <<- c(s, a)
    #print(paste(j, length(s)))
    data.frame(i = a, start_year = as.numeric(names(v0))[j])
  }) %>% bind_rows
  .gdata$start_year <- NULL
  .gdata <- left_join(.gdata, sl, by = "i") %>%
    mutate(start_year = ifelse(is.na(start_year), Inf, start_year))
  return(.gdata)
}




#dt0 %>% shuffle_start("mar") %>% with(., cor(start_year, n, method = "spearman"))

#shuffle_start_mar <- \(.gdata) {
#  p0 <- cumsum(prop.table(table(.gdata$start_year)))
#  p0 <- p0[- length(p0)]
#  z0 <- qlogis(p0)
#  q <- length(z0)
#  n <- nrow(.gdata)
#  mz0 <- matrix(z0, n, q, byrow = TRUE)
#  eta <- drop(scale(ntile(.gdata$my, 5)))
#  cumpr <- plogis(mz0 + eta)
#  p <- t(apply(cumpr, 1, \(x) diff(c(0, x, 1))))
#  start_year <- apply(p, 1, \(pi) rmultinom(1, 1, pi)) %>%
#    apply(., 2, \(v) {
#      levels(ordered(.gdata$start_year))[which(v == 1)] %>%
#        as.numeric
#    })
#  .gdata$start_year <- start_year
#  return(.gdata)
#}

shuffle_start_old2 <- \(.data, a = 3) {
  # first things first, recursively change the start year
  .data <- .data %>% mutate(start_year = start_year - a)

  .gdata <- .data %>%
    group_by(i, start_year, st) %>%
    reframe(start_year = first(start_year))

  v0  <- table(.gdata$start_year)
  cumsum(v0)
  v0 <- v0[-length(v0)]
  s <- NULL
  s <<- NULL
#  j =10
  sl <- lapply(seq_len(length(v0)), \(j) {
    dt <- .data[(!(.data$i %in% s)) & .data$year == ((2006:2016))[j], ]
    nrow(dt)
    new_s <- with(dt, sample(i, v0[j], prob = cpr))
    s <<- c(s, new_s)
    sdt <- data.frame(i = new_s, start_year = as.numeric(names(v0))[j])
    sdt
  }) %>% bind_rows
  .data$start_year <- NULL
  .data <- left_join(.data, sl, by = "i") %>%
    mutate(start_year = ifelse(is.na(start_year), Inf, start_year))
  #  .data$cpr <- NULL
  return(.data)
}

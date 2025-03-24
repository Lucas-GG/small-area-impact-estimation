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

reset_y0 <- \(.data) {
  .data <- .data %>%
    mutate(event_time = year - start_year
      , y0 = ifelse(year < start_year, y0, NA)
      , y0_25 = ifelse(year < start_year, y0_25, NA)
    )
}

add_hist <- \(.data, nlags = length(unique(.data$year)) - 1) {
  .data <- .data %>% reset_y0
  .data$h_y0 <- .data %>% with(flag(y0, 1:nlags, i, year))
  .data
}

expandx <- \(.data
  , nlags = 0
  , y25 = TRUE
  , recent = FALSE
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

  if (recent) {
    .data <- .data %>%
      mutate(
        y0_l10 = ifelse(year > start_year - 10, y0, NA) # no so far ago
      )
    .data$ml_i <- .data %>% gmean(y0_l10, i)
    .data$ml_st <- .data %>% gmean(y0_l10, st)
    .data$ml_u <- .data %>% gmean(y0_l10, urb)
    .data$ml_c <- .data %>% gmean(y0_l10, start_year)
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
#.data %>%  exp




#dt %>% expandx %>% select(m_i:m_c) %>% summary
#dt %>% add_avg %>% select(mi_y0:mc_y0) %>% summary


colnames(.data)

add_prop <- \(.data) {
  sdt <- .data %>%
    mutate(w = as.numeric(year == start_year - 1)) %>%
    filter(year < start_year) %>%
    select(w
      , n, n25
      , m_i, m25_i
      , m_st, m25_st
      , m_t, m25_t
      , m_u, m25_u
      , m_tr, m25_tr
      , uemp, pov, mhinc
      , lat, lon
    )

  mp <- ranger(w ~ ., data = sdt, probability = TRUE, num.threads = 1)
  .data$cpr <- predict(mp, .data, num.threads = 1)$predictions[, 2]
  .data
}
#.data %>% expandx %>% with(tapply(cpr, year, \(x) mean(x>0)))

#system.time(filter(placebo_dt, replication == 2) %>% expandx)







add_avg <- \(.data) {
  .data <- .data %>% reset_y0

  .data <- .data %>%
    #    group_by(i) %>%
    #    mutate(
    #      m_slp = ifelse(sum(!is.na(y0)) >1, coef(lm(y0 ~ year))[2], 0)
    #    ) %>%
    #    ungroup %>%
    mutate(
      mi_y0 = ave(y0, i, FUN = \(x) mean(x, na.rm = TRUE))
      , mt_y0 = ave(y0, year, FUN = \(x) mean(x, na.rm = TRUE))
      , mst_y0 = ave(y0, st, FUN = \(x) mean(x, na.rm = TRUE))
      , mc_y0 = ave(y0, start_year, FUN = \(x) mean(x, na.rm = TRUE)) #cohort
      #, mi_y25 = ave(y25, i, FUN = \(x) mean(x, na.rm = TRUE))
      #, mt_y25 = ave(y25, year, FUN = \(x) mean(x, na.rm = TRUE))
      #, mst_y25 = ave(y25, st, FUN = \(x) mean(x, na.rm = TRUE))
      #, mc_y25 = ave(y25, start_year, FUN = \(x) mean(x, na.rm = TRUE)) #cohort
      #, mtr_y0 = ave(y0, paste0(year, region), FUN = \(x) mean(x, na.rm = TRUE))
    )
  return(.data)
}
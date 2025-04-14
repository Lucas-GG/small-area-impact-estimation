library(lme4)
fm_me <- formula(y0 ~ offset(log(x)) - 1 + (1 | i))

impute_me <- \(.dt, fm) {
  if (!"wts" %in% names(.dt)) .dt[, wts := 1]
  glmer(fm, family = "poisson", .dt, weights = wts) %>%
    predict(.dt, type = "response")
}
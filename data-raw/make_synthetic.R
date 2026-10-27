# Synthetic county-year panel (data/synthetic_dt.rds) with the structure of the private panel, so
# that the code can be run without the restricted data. It does not reproduce the paper's results.
#
# Real values, from public sources: the 516 counties' adjacency (data/W, row = county `i`), years
# 1999-2016, NCHS urbanicity class `urb`, youth (10-24) population `n`, and the county covariates
# `uemp`, `pov`, `mhinc` (unemployment, poverty, median household income) and `AIAN` (American
# Indian/Alaska Native share), from Census Bureau and BLS releases. County identifiers are dropped.
#
# Generated:
#   start_year  adoption year, 2006-2016 or Inf (never adopts); about 70% of counties adopt, more
#               likely when larger; not the GLS adoption years.
#   y           youth suicide deaths, Poisson with rate n * exp(eta):
#                 eta = log(9e-5) + urbanicity effect + 0.015 (year - 2007) + u_i + s_i (year - 2007) / 10
#                       + log(0.85) after adoption,
#               u_i a spatially smoothed county effect (sd about 0.3), s_i a county trend (sd 0.1).
#               The known program effect is a rate ratio of 0.85 in every adopting county.
#   long        share of long trainings in a post-adoption year: missing before adoption and in about
#               half of the years after the adoption year, zero in about 60% of adopting counties, otherwise Beta(2, 5).
#
# Usage (with the private panel, to copy the public columns): Rscript data-raw/make_synthetic.R
library(data.table)
set.seed(20261027)
real <- as.data.table(readRDS("data/notpublic/compact_dt.rds"))
dt <- real[, .(i, year, urb, n, uemp, pov, mhinc, AIAN)]
setorder(dt, i, year)
load("data/W")
counties <- dt[, .(n = mean(n), urb = urb[1]), by = i]
stopifnot(nrow(counties) == nrow(W), all(counties$i == seq_len(nrow(W))))

# Adoption: 70% of counties on average, more likely when larger; year uniform on 2006-2016.
p <- plogis(qlogis(0.7) + 0.4 * as.numeric(scale(log(counties$n))))
counties[, start_year := fifelse(runif(.N) < p, as.numeric(sample(2006:2016, .N, replace = TRUE)), Inf)]

# County effects: iid noise smoothed over neighbours, rescaled; county trends.
z <- rnorm(nrow(W))
deg <- pmax(rowSums(W), 1)
u <- as.numeric(z + (W %*% z) / deg)
counties[, `:=`(u = 0.3 * u / sd(u), slope = rnorm(.N, 0, 0.1))]
urb_effect <- c("1.LCMetro" = -0.3, "2.LFMetro" = -0.2, "3.MMetro" = -0.1, "4.SMetro" = 0, "5.Micro" = 0.15, "6.NonCore" = 0.3)
counties[, b_urb := urb_effect[as.character(urb)]]

dt <- merge(dt, counties[, .(i, start_year, u, slope, b_urb)], by = "i")
dt[, post := year >= start_year]
dt[, eta := log(9e-5) + b_urb + 0.015 * (year - 2007) + u + slope * (year - 2007) / 10 + fifelse(post, log(0.85), 0)]
dt[, y := rpois(.N, n * exp(eta))]

# Long-training share, post-adoption years only.
share <- counties[, .(i, s = fifelse(runif(.N) < 0.6, 0, rbeta(.N, 2, 5)))]
dt <- merge(dt, share, by = "i")
# Recorded in the adoption year (adoption is the first year with a recorded event) and in about half the later years.
dt[, long := fifelse(post & (year == start_year | runif(.N) > 0.5), pmin(pmax(s + fifelse(s > 0, rnorm(.N, 0, 0.05), 0), 0), 1), NA_real_)]

out <- dt[, .(i, year, y, n, start_year, urb, uemp, pov, mhinc, AIAN, long)]
setorder(out, i, year)
saveRDS(out, "data/synthetic_dt.rds")
print(out[, .(counties = uniqueN(i), rows = .N, deaths_per_county_year = mean(y), share_zero = mean(y == 0),
  rate_per_100k = 1e5 * sum(y) / sum(n), adopting = uniqueN(i[is.finite(start_year)]))])

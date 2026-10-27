# Supplement S8 and Section 2.4: equivalence of TWFE and the cohort x group ETWFE, without covariates,
# on the case-study data (actual adoption) and one placebo dataset (scenario 1, replication 1):
#   (a) year effects gamma_t (relative to the first year: only differences are identified);
#   (b) predictions summed over each cohort x group cell, all counties, all fitted (untreated) years;
#   (c) noncore predictions over the evaluated county-years (placebo: recoverable event years 0-4;
#       case study: post-adoption event years 0-4), by county-year and summed, over the rows both
#       models predict and over all evaluated rows with TWFE's missing predictions as 0;
#   (d) predictions for the imputed years (case study: post-adoption; placebo: masked years), summed
#       over whole cohort x group cells (the out-of-sample counterpart of (b));
# each with and without the population offset log(n). Then, over all 600 placebo datasets, the
# relative difference in the noncore summed prediction (c) of the cohort x group ETWFE and of the
# pooled-cohort ETWFE against TWFE (point fits, offset, no covariates). Counties with no untreated
# deaths have no TWFE fit (their county effect diverges to -Inf); their TWFE predictions count as 0.
# Writes results/numbers/etwfe_twfe_equivalence{,_placebo}.csv (numbers quoted in Section 2.4) and the
# per-dataset sums to data/derived/pooled_cohorts/. Usage: Rscript analysis/S8_pooled_cohorts/s8_03_equivalence.R
library(data.table)
invisible(lapply(list.files("R", "[.]r$", full.names = TRUE), source))
panel <- as.data.table(readRDS(panel_file()))[, .(i, year, start_year, urb, y, n)]
assignments <- readRDS("data/simulation/simulated_assignments_with_meta.rds")$assignments

fit <- function(dt, fe, offset) {
  fml <- stats::as.formula(paste("y0 ~ 1 |", fe))
  if (offset) fixest::fepois(fml, data = dt, offset = ~log(n), notes = FALSE, warn = FALSE)
  else fixest::fepois(fml, data = dt, notes = FALSE, warn = FALSE)
}
diffs <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  c(max_abs = max(abs(x[ok] - y[ok])), max_rel = max(abs(x[ok] / y[ok] - 1)[y[ok] != 0]), n = sum(ok))
}

compare <- function(dt, label, eval_rows) {
  dt[, cell := as.vector(cohort_cell(cohort, ncf, y0))]
  train <- !is.na(dt$y0)
  out <- list()
  for (offset in c(TRUE, FALSE)) {
    tw <- fit(dt, "i + year", offset); et <- fit(dt, "cell + year", offset)
    # (a) year effects, relative to the first year
    ga <- fixest::fixef(tw)$year; gb <- fixest::fixef(et)$year
    yrs <- intersect(names(ga), names(gb))
    a <- diffs(ga[yrs] - ga[yrs][1], gb[yrs] - gb[yrs][1])
    # (b) predictions summed over each cell, fitted (untreated) rows, TWFE NA (no fit) as 0
    p_tw <- predict(tw, newdata = dt, type = "response"); p_et <- predict(et, newdata = dt, type = "response")
    p_tw0 <- fifelse(is.na(p_tw), 0, p_tw)
    cs <- data.table(cell = dt$cell[train], tw = p_tw0[train], et = p_et[train], obs = dt$y0[train])[,
      .(tw = sum(tw), et = sum(et), obs = sum(obs)), by = cell]
    b <- diffs(cs$tw, cs$et)
    out[[length(out) + 1]] <- data.table(data = label, offset, comparison = c("(a) year effects", "(b) cell sums, fitted rows"),
      max_abs = c(a["max_abs"], b["max_abs"]), max_rel = c(a["max_rel"], b["max_rel"]), n = c(a["n"], b["n"]))
    # (c) noncore evaluated county-years
    ev_all <- eval_rows(dt)
    ev <- ev_all & is.finite(p_tw) & is.finite(p_et)
    cr <- diffs(p_tw[ev], p_et[ev])
    tot <- c(sum(p_tw[ev]), sum(p_et[ev]))
    tot0 <- c(sum(p_tw0[ev_all]), sum(p_et[ev_all]))
    # (d) imputed years (no untreated outcome), summed over whole cells, TWFE NA as 0
    imp <- is.na(dt$y0)
    dc <- data.table(cell = dt$cell[imp], tw = p_tw0[imp], et = p_et[imp])[, .(tw = sum(tw), et = sum(et)), by = cell]
    dd <- diffs(dc$tw, dc$et)
    out[[length(out) + 1]] <- data.table(data = label, offset,
      comparison = c("(c) evaluated, by county-year", "(c) evaluated, summed (both predict)",
        "(c) evaluated, summed (TWFE NA as 0)", "(d) imputed years, cell sums"),
      max_abs = c(cr["max_abs"], abs(diff(tot)), abs(diff(tot0)), dd["max_abs"]),
      max_rel = c(cr["max_rel"], abs(tot[2] / tot[1] - 1), abs(tot0[2] / tot0[1] - 1), dd["max_rel"]),
      n = c(cr["n"], cr["n"], sum(ev_all), dd["n"]))
  }
  out[[length(out) + 1]] <- data.table(data = label, offset = NA, comparison = "cells merged (zero deaths)",
    max_abs = NA, max_rel = NA, n = length(attr(cohort_cell(dt$cohort, dt$ncf, dt$y0), "merged")))
  rbindlist(out)
}

# Case study: actual adoption; evaluated rows are noncore adopters' event years 0-4.
cs_dt <- copy(panel); set_y0(cs_dt)
res_cs <- compare(cs_dt, "case study", function(d) d$ncf == 1 & is.finite(d$cohort) & d$year >= d$cohort & d$year <= d$cohort + 4)
# Placebo: scenario 1, replication 1; evaluated rows are recoverable noncore event years 0-4.
pl_dt <- copy(panel); set_y0(pl_dt, A_sim = assignments[, 1, 1])
res_pl <- compare(pl_dt, "placebo (MAR, r = 1)", function(d) d$ncf == 1 & is.finite(d$cohort) & d$year >= d$cohort &
  d$year <= d$cohort + 4 & d$year < d$start_year)
res <- rbind(res_cs, res_pl)
setorder(res, data, comparison, -offset)
fwrite(res, "results/numbers/etwfe_twfe_equivalence.csv")
print(res, digits = 3)

# All placebo datasets: noncore summed prediction over the evaluated rows every model predicts.
one <- function(k, r) {
  dt <- copy(panel); set_y0(dt, A_sim = assignments[, r, k])
  dt[, cell := as.vector(cohort_cell(cohort, ncf, y0))]
  p <- lapply(c(twfe = "i + year", cell = "cell + year", pooled = "cohort + year"),
    function(fe) predict(fit(dt, fe, TRUE), newdata = dt, type = "response"))
  ev <- dt$ncf == 1 & is.finite(dt$cohort) & dt$year >= dt$cohort & dt$year <= dt$cohort + 4 & dt$year < dt$start_year &
    Reduce(`&`, lapply(p, is.finite))
  data.table(scenario = k, r = r, rows = sum(ev), twfe = sum(p$twfe[ev]), cell = sum(p$cell[ev]), pooled = sum(p$pooled[ev]))
}
all <- rbindlist(parallel::mclapply(seq_len(600), function(j) one((j - 1) %/% 200 + 1, (j - 1) %% 200 + 1), mc.cores = 18))
all[, `:=`(rel_cell = cell / twfe - 1, rel_pooled = pooled / twfe - 1)]
dir.create("data/derived/pooled_cohorts", recursive = TRUE, showWarnings = FALSE)
fwrite(all, "data/derived/pooled_cohorts/equivalence_placebo_all.csv")
summ <- melt(all, id.vars = "scenario", measure.vars = c("rel_cell", "rel_pooled"), variable.name = "etwfe")[,
  .(mean = mean(value), p10 = quantile(value, 0.1), p90 = quantile(value, 0.9), mean_abs = mean(abs(value))), by = .(scenario, etwfe)]
summ[, scenario := c("MAR", "Latent level", "Latent trend")[scenario]]
fwrite(summ, "results/numbers/etwfe_twfe_equivalence_placebo.csv")
print(summ, digits = 3)

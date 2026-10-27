# Supplement S8: placebo scores and case-study effects of the models with adoption-cohort effects
# pooled across noncore and other counties, and of the intermediate RW (cohort x group intercepts,
# one deviation walk per cohort). Fits: 02_fit_placebo.qmd (`*_pooled`, `rw_cohort_walks`).
#   Placebo: year-0-4 ledger on the evaluation county-years of Table 2 (03a samples), and the paired
#     change in year-4 noncore RMSE from pooled cohorts to cohort x group cells for every model.
#   Case study (with covariates, the 82 counties of Table 3): ETWFE and RW with pooled cohorts and
#     the intermediate RW, with the panel, draws and seeds of 05_case_study.qmd.
# Writes data/derived/pooled_cohorts/{ledger.rds,case.csv} and results/numbers/s8_pooled_vs_cells.csv.
# Usage: Rscript analysis/S8_pooled_cohorts/s8_01_score.R
suppressMessages({library(data.table); invisible(lapply(list.files("R", "[.]r$", full.names = TRUE), source))})
load("data/W"); rownames(W) <- NULL
options(xpn.conf = list(Lset = c(2, 4), graph = Matrix::Matrix(W)))
Sys.setenv(OMP_NUM_THREADS = 1); INLA::inla.setOption(num.threads = 1)
out <- "data/derived/pooled_cohorts"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dgps <- c("MAR", "Latent level", "Latent trend")
specs <- c(wocov = "Without covariates", wcov = "With covariates")
labels <- c(rw_cohort_walks = "RW, cohort walks", etwfe_pooled = "ETWFE, pooled", rw_pooled = "RW, pooled",
  rwi_pooled = "RW-IID, pooled", rwst_pooled = "RW-ST, pooled", xgb_pooled = "GB-Ens, pooled")

# Placebo: each model on the county-years of the Table 2 sample.
ledger <- rbindlist(lapply(names(specs), function(sp) {
  keys <- unique(readRDS(sprintf("data/simulation/placebo_sample_%s.rds", sp))[, .(dgp = as.character(dgp), r, i, s)])
  rbindlist(lapply(names(labels), function(code) {
    fits <- rbindlist(lapply(1:3, function(k) {
      x <- readRDS(sprintf("data/simulation/placebo/%s/%s_%s.rds", sp, code, k))
      x[, `:=`(dgp = dgps[k], s = year - start_year_sim, nc = ordered(urb != "6.NonCore", labels = c("NonCore", "Core")), model = labels[[code]])]
      x[s >= 0 & s <= 4][keys[dgp == dgps[k]], on = .(dgp, r, i, s), nomatch = 0L]
    }))
    stopifnot(all(is.finite(fits$y0_hat)), nrow(fits) == nrow(keys))
    fits[, dgp := ordered(dgp, levels = dgps)]
    l <- aggregate_placebo(fits, seed = 0203)$aggregate[, specification := specs[[sp]]]
    rm(fits); gc()
    l
  }))
}))
saveRDS(ledger, file.path(out, "ledger.rds"))

# Paired change in year-4 noncore RMSE, pooled cohorts -> cohort x group cells (S8 text).
main <- as.data.table(readRDS("_output/analysis/output/03b/placebo_scores_by_replication.rds"))[nc == "NonCore" & s == 4]
pooled <- ledger[nc == "NonCore" & s == 4 & grepl(", pooled$", model)][, model := sub(", pooled$", "", model)]
chg <- rbindlist(lapply(c("ETWFE", "RW", "RW-IID", "RW-ST", "GB-Ens"), function(m) {
  z <- merge(main[model == m, .(specification, dgp = as.character(dgp), r, cell = sq_diff)],
    pooled[model == m, .(specification, dgp = as.character(dgp), r, pool = sq_diff)], by = c("specification", "dgp", "r"))
  z[, c(list(model = m), paired_change(cell, pool, rmse = TRUE)), by = .(specification, dgp)]
}))
fwrite(chg, "results/numbers/s8_pooled_vs_cells.csv")

# Case study: the 05 panel and sample; seeds of 05 for ETWFE and RW, 20261103 for the intermediate.
covariates <- c("uemp", "pov", "mhinc", "aian")
panel <- backdate_covariate_panel(as.data.table(readRDS(panel_file())))[, c("i", "year", "start_year", "urb", "y", "n", covariates), with = FALSE]
seeds <- fread("data/derived/case_study/case_models_and_seeds.csv")[spec == "wcov"]
case_fit <- function(label, seed, pooled, impute, ...) {
  options(etwfe.group = !pooled)
  dt <- copy(panel); set_y0(dt); set.seed(seed)
  fit_blind(dt, impute, nboot = 200L, ...)
  options(etwfe.group = TRUE)
  dt[, c("i", "year", "start_year", "cohort", "urb", "y", "n", "y0_hat", paste0("y0_post_", 1:200)), with = FALSE][,
    `:=`(spec = "wcov", model = label, s = year - start_year, nc = fifelse(urb == "6.NonCore", "NonCore", "Core"))]
}
f <- rbind(
  case_fit("ETWFE, pooled", seeds[model == "ETWFE", seed], TRUE, impute_fe, fm = y0 ~ uemp + pov + mhinc + aian | cell + year),
  case_fit("RW, pooled", seeds[model == "RW", seed], TRUE, impute_inla,
    fit_fun = function(x) fit_nested(x, structural_only = TRUE, iid = FALSE, covariates = covariates)),
  case_fit("RW, cohort walks", 20261103L, FALSE, impute_inla,
    fit_fun = function(x) fit_nested(x, structural_only = TRUE, iid = FALSE, covariates = covariates, walks = "cohort")))
cur <- rbindlist(lapply(list.files("data/derived/case_study", "^fit_w(o)?cov_.*[.]rds$", full.names = TRUE),
  function(p) readRDS(p)$fit[, .(i, year, start_year, urb, y0_hat)]))
cur <- cur[urb == "6.NonCore" & is.finite(start_year) & year - start_year >= 0 & year - start_year <= 4]
common <- cur[, .(ok = all(is.finite(y0_hat))), by = .(i, s = year - start_year)][ok == TRUE, .(i, s)]
e <- case_effects(case_cumulative(f[nc == "NonCore"][common, on = .(i, s), nomatch = 0L]), balanced = TRUE)[s == 4]
stopifnot(all(e$counties == 82))
fwrite(e, file.path(out, "case.csv"))

y4 <- performance_summary(ledger)[nc == "NonCore" & s == 4 & metric %in% c("RMSE", "coverage")]
print(dcast(y4, specification + model ~ dgp + metric, value.var = "estimate"), digits = 3)
print(chg[, .(specification, dgp, model, change = round(percent_change), se = round(mc_se))])
print(e[, .(model, rr = round(rr, 2), lo = round(rr_lo, 2), hi = round(rr_hi, 2))])

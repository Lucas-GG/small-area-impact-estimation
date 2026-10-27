# Supplement S8, Table S10: ETWFE and RW with adoption-cohort effects pooled across noncore and
# other counties ("Cohort"), against the main cohort x noncore specification, and for RW the
# intermediate with cohort x noncore intercepts and one walk per cohort.
# Placebo, noncore, event year 4, both specifications; case-study year-4 rate ratios.
#   Cohort x noncore and TWFE rows: the 03b ledger behind Table 2 (checked against Table 2).
#   Cohort and intermediate rows: s8_01_score.R (same evaluation county-years).
# RMSE comparisons are paired by replication on the same squared errors (no Poisson noise).
# Usage: Rscript analysis/S8_pooled_cohorts/s8_02_table.R
library(data.table)
library(gt)
invisible(lapply(list.files("R", "[.]r$", full.names = TRUE), source))
stem <- "results/tables/tableS10_cohort_noncore"
specs <- c("Without covariates", "With covariates")
dgps <- c("MAR", "Latent level", "Latent trend")
cnc <- c(ETWFE = "ETWFE, cohort x noncore", RW = "RW, cohort x noncore")
mid <- "RW, cohort walks"
y4 <- function(x) x[nc == "NonCore" & s == 4]

main <- y4(as.data.table(readRDS("_output/analysis/output/03b/placebo_scores_by_replication.rds")))[
  specification %in% specs & model %in% c("TWFE", "ETWFE", "RW")]
main[, model := fcase(model == "ETWFE", cnc[["ETWFE"]], model == "RW", cnc[["RW"]], default = "TWFE")]
s8 <- y4(readRDS("data/derived/pooled_cohorts/ledger.rds"))[model %in% c("ETWFE, pooled", "RW, pooled", mid)]
s8[, model := sub(", pooled$", "", model)]
ledger <- rbind(main, s8, fill = TRUE)
ledger[, dgp := as.character(dgp)]

perf <- performance_summary(ledger)
vs_twfe <- relative_to_model(ledger, "TWFE")[metric == "RMSE"]
chg <- rbind(relative_to_model(ledger[model %in% c("ETWFE", cnc[["ETWFE"]])], "ETWFE")[model == cnc[["ETWFE"]]],
  relative_to_model(ledger[model %in% c("RW", cnc[["RW"]], mid)], "RW")[model %in% c(cnc[["RW"]], mid)])[metric == "RMSE"]
num <- function(x, d) mt_minus(sprintf(paste0("%.", d, "f"), x))
long <- merge(dcast(perf[metric %in% c("RMSE", "coverage")], specification + dgp + model ~ metric, value.var = c("estimate", "mc_se")),
  vs_twfe[, .(specification, dgp, model, pct = percent_change, pct_se = mc_se)], by = c("specification", "dgp", "model"), all.x = TRUE)
long <- merge(long, chg[, .(specification, dgp, model, chg = percent_change, chg_se = mc_se)], by = c("specification", "dgp", "model"), all.x = TRUE)
def <- function(m) fcase(m == "TWFE", "", m %in% cnc, "Cohort × noncore", m == mid, "Cohort × noncore intercepts, cohort walks", default = "Cohort")
long[, `:=`(
  rmse = sprintf("%s (%s)", num(estimate_RMSE, 2), num(mc_se_RMSE, 2)),
  vs = fifelse(model == "TWFE", "—", sprintf("%s (%s)", num(pct, 0), num(pct_se, 0))),
  cov = num(estimate_coverage, 2),
  change = fifelse(is.na(chg), "", sprintf("%s (%s)", mt_minus(sprintf("%+.0f", chg)), num(chg_se, 0))),
  Model = sub(",.*", "", model), Definition = def(model))]

# Check: the cohort x noncore and TWFE rows reproduce Table 2.
t2 <- fread("results/tables/table2_model_comparison.csv")
for (g in dgps) {
  a <- long[model %in% c("TWFE", cnc) & dgp == g, .(specification, model = sub(",.*", "", model), rmse2 = num(estimate_RMSE, 2),
    cov2 = num(estimate_coverage, 2), pct1 = round(pct, 1))]
  b <- t2[model %in% c("TWFE", "ETWFE", "RW"), .(specification, model, rmse2 = sub(" .*", "", get(paste0("RMSE_", g))),
    cov2 = sub(" .*", "", get(paste0("Coverage_", g))), pct1 = suppressWarnings(as.numeric(sub(" .*", "", get(paste0("vs TWFE (%)_", g))))))]
  m <- merge(a, b, by = c("specification", "model"))
  stopifnot(nrow(m) == 6, all(m$rmse2.x == m$rmse2.y), all(m$cov2.x == m$cov2.y), all(is.na(m$pct1.y) | abs(m$pct1.x - m$pct1.y) < 1e-9))
}
cat("Cohort x noncore and TWFE rows reproduce Table 2 (RMSE, coverage, % vs TWFE).\n")

w <- dcast(long, specification + Model + Definition ~ dgp, value.var = c("rmse", "vs", "cov", "change"), sep = "__")
w[, `:=`(specification = factor(specification, specs), Model = factor(Model, c("TWFE", "ETWFE", "RW")),
  Definition = factor(Definition, c("", "Cohort", "Cohort × noncore intercepts, cohort walks", "Cohort × noncore")))]
setorder(w, specification, Model, Definition)
w[, Model := fifelse(rowid(specification, Model) == 1, as.character(Model), "")]
w[, `:=`(group = as.character(specification), Definition = as.character(Definition))]
cols <- unlist(lapply(dgps, function(g) paste0(c("rmse", "vs", "cov", "change"), "__", g)))
w <- w[, c("group", "Model", "Definition", cols), with = FALSE]

# Case study: year-4 rate ratio (95% interval), balanced sample, with covariates.
ce <- function(path, models) fread(path)[s == 4 & spec == "wcov" & grepl("Balanced", sample) & model %in% models]
cs <- rbind(ce("data/derived/case_study/case_effects.csv", c("TWFE", "ETWFE", "RW"))[, model := fifelse(model == "TWFE", "TWFE", cnc[model])],
  fread("data/derived/pooled_cohorts/case.csv")[, model := sub(", pooled$", "", model)], fill = TRUE)
stopifnot(nrow(cs) == 6, all(cs$counties == 82))
cs[, `:=`(group = "Case study (with covariates)", Model = sub(",.*", "", model), Definition = def(model),
  rr = mt_minus(sprintf("%.2f (%.2f, %.2f)", rr, rr_lo, rr_hi)))]
cs[, ord := match(model, c("TWFE", "ETWFE", cnc[["ETWFE"]], "RW", mid, cnc[["RW"]]))]
setorder(cs, ord)
cs[, Model := fifelse(rowid(Model) == 1, Model, "")]
case_rows <- cs[, c(list(group = group, Model = Model, Definition = Definition), setNames(rep(list(""), length(cols)), cols))]
case_rows[[cols[1]]] <- cs$rr
tab <- rbind(w, case_rows)
fwrite(tab, paste0(stem, ".csv"))

build <- function(fmt) {
  g <- gt(as.data.frame(tab), groupname_col = "group", rowname_col = "Model")
  for (d in dgps) g <- tab_spanner(g, d, columns = paste0(c("rmse", "vs", "cov", "change"), "__", d), id = d)
  lab <- function(d) setNames(list(mt_label("RMSE\n(SE)", fmt), mt_label("% vs TWFE\n(SE)", fmt), "Coverage", mt_label("Change,\n% (SE)", fmt)),
    paste0(c("rmse", "vs", "cov", "change"), "__", d))
  g <- cols_label(g, Definition = "", .list = do.call(c, lapply(dgps, lab)))
  mt_style(g, stub = c("Model", "Definition"))
}
# 12 scenario columns do not fit 6.5 in at 9 pt: set for a landscape page (9 in text width).
width <- save_manuscript_table(build, stem, textwidth_in = 9, tabcolsep = "2pt")
cat(sprintf("LaTeX width at 9 pt: %.2f in; overfull boxes: %d\n", width$width_in, width$overfull))
print(tab[, c("group", "Model", "Definition", grep("__MAR$", names(tab), value = TRUE)), with = FALSE])

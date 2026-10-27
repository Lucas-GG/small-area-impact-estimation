# Lead section of the half-A report: the fixed rule and its headline results.
# Reads saved aggregate summaries only; never recalculates screening statistics.

load_decision <- function(root) {
  # The decision file is fixed; read it as is (it has no trailing newline).
  path <- file.path(root, 'analysis/S6_screening/s3_decision.yml')
  d <- yaml::yaml.load(paste(readLines(path, warn = FALSE), collapse = '\n'))
  stopifnot(d$diagnostic %in% c('pre_rmse', 'pre_rms_pearson', 'ho_raw', 'ho_standardized', 'ho_uncertainty', 'p_two'),
    is.numeric(d$cutoff), d$missing_diagnostic_handling %in% c('keep', 'exclude'))
  d
}

# Rows of an aggregate summary for the fixed rule and the unscreened baseline.
rule_rows <- function(s, decision, subgroup_name = 'noncore') {
  policy <- decision$missing_diagnostic_handling
  s[subgroup == subgroup_name & (diagnostic == 'unscreened' |
    (diagnostic == decision$diagnostic & cutoff == decision$cutoff & missing_handling == policy))]
}

# One row per scenario x model: unscreened versus the fixed rule.
fixed_rule_table <- function(s, decision, subgroup_name = 'noncore') {
  v <- rule_rows(s, decision, subgroup_name)
  base <- v[diagnostic == 'unscreened', .(scenario, model, rmse_base = rmse, relative_base = relative_rmse,
    coverage_base = coverage, n_base = n_retained)]
  rule <- v[diagnostic == decision$diagnostic, .(scenario, model, rmse_ratio, rmse_ratio_mc_se,
    relative_rule = relative_rmse, coverage_rule = coverage, width_rule = width, excluded_share, n_rule = n_retained)]
  out <- merge(base, rule, by = c('scenario', 'model'))
  stopifnot(nrow(out) == nrow(base), nrow(out) > 0)
  out[order(match(scenario, scenario_order), match(model, model_order))]
}

fixed_rule_gt <- function(tab) {
  v <- tab[, .(Scenario = scenario, Model = model,
    `Counties (unscreened)` = round(n_base, 1), `Excluded` = excluded_share,
    `RMSE ratio` = sprintf('%.2f (%.2f)', rmse_ratio, rmse_ratio_mc_se),
    `Relative RMSE: unscreened` = relative_base, `Relative RMSE: screened` = relative_rule,
    `Coverage: unscreened` = coverage_base, `Coverage: screened` = coverage_rule)]
  gt::gt(as.data.frame(v), groupname_col = 'Scenario') |>
    gt::fmt_percent(columns = c(`Excluded`, `Relative RMSE: unscreened`, `Relative RMSE: screened`,
      `Coverage: unscreened`, `Coverage: screened`), decimals = 1) |>
    gt::tab_spanner('Relative aggregate RMSE', columns = gt::starts_with('Relative')) |>
    gt::tab_spanner('95% aggregate coverage', columns = gt::starts_with('Coverage')) |>
    gt::cols_label(`Relative RMSE: unscreened` = 'Unscreened', `Relative RMSE: screened` = 'Screened',
      `Coverage: unscreened` = 'Unscreened', `Coverage: screened` = 'Screened') |>
    gt::tab_options(table.font.size = 12, data_row.padding = 3)
}

# Values the decision file recorded from half A next to what the saved results give now.
decision_reference_table <- function(s, s_common, decision) {
  ref <- decision$half_a_reference
  tab <- fixed_rule_table(s, decision)
  non_twfe <- tab[model != 'TWFE']
  e <- fixed_rule_table(s_common, decision)[model != 'TWFE']
  top10 <- s[subgroup == 'noncore' & diagnostic == 'pre_rmse' & cutoff == 10 &
    missing_handling == decision$missing_diagnostic_handling, mean(excluded_share)]
  data.table::data.table(
    Quantity = c('Noncore coverage, all models: unscreened', 'Noncore coverage, all models: screened',
      'Noncore coverage, non-TWFE on E: unscreened', 'Noncore coverage, non-TWFE on E: screened',
      'Noncore RMSE ratio, non-TWFE: lowest', 'Noncore RMSE ratio, non-TWFE: highest',
      'Pre-RMSE top 10%: noncore share excluded'),
    Recorded = c(ref$noncore_coverage_mean_all_models$unscreened, ref$noncore_coverage_mean_all_models$screened,
      ref$noncore_coverage_mean_non_twfe_on_E$unscreened, ref$noncore_coverage_mean_non_twfe_on_E$screened,
      ref$noncore_abs_rmse_ratio_non_twfe_range[[1]], ref$noncore_abs_rmse_ratio_non_twfe_range[[2]],
      ref$pre_rmse_top10_noncore_excluded_share),
    Current = c(mean(tab$coverage_base), mean(tab$coverage_rule), mean(e$coverage_base), mean(e$coverage_rule),
      min(non_twfe$rmse_ratio), max(non_twfe$rmse_ratio), top10))
}

pct1 <- function(x) sprintf('%.1f%%', 100 * x)

# Plain-language statement of the fixed rule.
rule_label <- function(decision) {
  what <- switch(decision$diagnostic,
    ho_standardized = '|held-out standardized error|', ho_uncertainty = '|held-out uncertainty-scaled error|',
    ho_raw = '|held-out raw error|', pre_rmse = 'pre-period RMSE', pre_rms_pearson = 'pre-period Pearson RMS',
    p_two = 'two-sided predictive mid-p')
  op <- if (decision$diagnostic == 'p_two') '<' else '>'
  sprintf('%s %s %s', what, op, format(decision$cutoff))
}

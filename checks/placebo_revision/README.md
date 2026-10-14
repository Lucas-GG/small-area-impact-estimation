# Placebo revision: Phase 0 only

New parameterized report: `analysis/3.2_revision_placebo_summary.qmd`.
`specification` accepts `no_covariates`, `with_covariates`, or `both`.
`recompute: true` reruns the computations; the default reads the new audit CSVs.
The report and CSVs are new; legacy scripts and results are preserved.

The numerical reference is the whole-sample historical algorithm in git
`18b3ae2`, paths `analyses/3.2 summary placebo.rmd` and
`analyses/with covariates/3.2.1 summary placebo.rmd`. Source snapshots are in
`_output/placebo_revision/phase0/source_reference/`.
The current covariate summary is not the historical reference: it evaluates B
and uses a different restriction. All four saved historical HTML tables were
additionally checked against replayed numbers within their three-decimal rounding.

Original seed 0203, draw-column order, row order, Poisson noise, cumulative
sums, and shared scoring functions are preserved. The historical screen takes
the median pre-RMSE by scenario/replication/county across model rows, then uses
the strict scenario/subgroup 90th-percentile threshold. No new percentile
definition, complete-window filter, target, or rule is substituted in the replay.
The new aggregation is compared with a literal reference expression using
the same predictive draws. The reproduction CSV contains each metric/model/
scenario/subgroup/horizon cell, including the historical standard deviations.

The diagnostic audit loads both physical halves of each specification, checks
their assignment hashes and county IDs, and joins their own diagnostics onto
the whole main-simulation sample. Keys and label translations are exported.
Missing shares are computed over distinct county/model/scenario/simulation keys;
the same diagnostic applies at every event-year row for that county key.
Unmatched keys and unexplained unavailable diagnostics stop the audit.

Outputs: `_output/placebo_revision/phase0/`, including the output index and input
path logs. HTML: `_output/analysis/3.2_revision_placebo_summary.html`.
Phase 1 and Phase 2 have not been run.

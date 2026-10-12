# Half-A available-window revision — 2026-10-11

No models were fit. The split seed (2032026), assignment, original main/backdated checkpoints, legacy A input and all B inputs are unchanged. Exploration loads only the physically A-only `half_A_available.rds` and asserts `half == "A"`; preparation deserializes the unsplit saved RDS objects and immediately filters A before calculations. RDS cannot selectively deserialize county rows.

## Changed source files

- `R/heldout_exploration.R`: new available-window summaries, predictive cumulative intervals, A-only preparation/cache, diagnostics, county averages, k-stratified/pooled associations, Monte Carlo summaries, retention/random screens, size-stratum flag rates, common-sample sensitivity and ranks.
- `R/heldout_prepare.R`: route no-covariate preparation to the revised A-only workflow; preserve the covariate legacy route.
- `R/heldout_analysis.R`: validate six decision fields and reject the blank template before confirmation.
- `analysis/3.1a.1_heldout_compute.qmd`: revised self-contained exploration and exports, calibration, cutoff curves, sensitivity targets, coverage, ranks and scenario robustness.
- `analysis/3.1a.4_heldout_decision.yml`: blank six-field template; no choice recorded.
- `analysis/3.1b heldout confirmation.qmd` and `analysis/with_covariates/3c.1b heldout confirmation.qmd`: intersection-sample placeholders only; neither rendered. Confirmation needs the new target/k/missing-policy implementation before use.
- `checks/tests/check_available_windows.R`: hand-built target, availability, predictive-draw, missing-policy and association acceptance checks.
- `analysis/heldout README.md`: updated entry points and output locations.
- This revision note.

## Saved schema and assumptions

Main checkpoint columns are `r` (simulation), `i` (county), `year`, `start_year_sim`, `y`, `y0_hat`, and `y0_post_1` through `y0_post_200`. Event time is `year - start_year_sim`. The draw columns contain latent expected counts, not predictive observed counts. There are no saved yearly predictive quantiles or post-period variances. We generate seeded Poisson counts conditional on those saved latent draws, preserving draw indices across years, then sum matching predictive draws over the selected window. This assumes independent Poisson observation noise conditional on the model's latent draw. It does not refit any model. Cumulative mean uncertainty is the variance of the corresponding latent draw sums. Invalid draws make coverage unavailable; coverage denominators and missing shares are exported separately.

Backdated summaries use `sim`, `model`, `id`, `x`, `scenario`, `y_ho`, `mu_ho_mean`, `mu_ho_var`, `ppp_ho_upper`, `ppp_ho_lower`, `y_post`, `mu_post_mean`, `pre_rmse`, `pre_rms_pearson`. Paired legacy inclusive tails recover mid-p exactly; an unpaired tail stays unavailable.

New primary targets sum the finite observable years at event times 0–4 with positive expected counts. Complete-five-year sensitivities preserve the original nonnegative-count eligibility. County averages use k >= 2 by default, include the current simulation and are explanatory. Top-decile boundary ties count as positives. Pre-RMSE top-percent screens use the nearest requested count and county ID to resolve ties. Associations exclude strata with fewer than five paired rows or a single outcome class. Finite statistics are averaged equally over simulations; MC SE uses their between-simulation SD / sqrt(number finite). Random-screen bands average the within-simulation 2.5/97.5 percentiles; they are reference bands, not confidence intervals for the mean. The operational screen stays per model; the intersection is a sensitivity with target availability paired across all six models.

## Generated outputs and checks

New data: `data/derived/heldout/half_A_available.rds`, 18 A-only summary caches under `post_A/`, and `available_window_counts.csv`. Revised CSVs/plots are under `_output/heldout/exploration/available_window/`; rendered HTML is under `_output/analysis/`. Existing outputs outside the explicitly replaced exploration are preserved.

The available raw target increases from 187,602 complete-window rows to 367,885 available-window rows across 18 model/scenario cells (639,210 saved A rows total). All complete-window point targets reproduce the saved old targets within 1e-8; the document also checks simulation-level eligibility and baseline errors to three decimals. The cancellation and single-year hand tests pass. Cached preparation completed in 7.36 seconds without errors; the uncached remaining six files took 153.317 seconds after resuming the initial cache pass. Rendering runtime is exported in `runtime.csv` when rendering finishes.

Rendering completed: all computation chunks took 494.024 seconds (about 8.2 minutes); the complete HTML was generated roughly nine minutes after launch. The 12.4 MB self-contained HTML was verified to contain the final sections and closing HTML tag, then moved to `_output/analysis/3.1a.1_heldout_compute.html`. Because the exploration is excluded from the project-wide render list, an explicit render writes beside its source by default; the later presentation revision adds two included project documents that render centrally. The excluded legacy computation document renders beside its source; `--output-dir` is only accepted for project documents. There are 22 CSV exports and four PNGs. All 43,200 complete-window eligibility/baseline checks passed, as did identical retained-count checks for all 2,520 common-sample rule/target/scenario groups. Size-stratum exports include all strata 1–6. No render error was observed. The acceptance tests also verify common-sample pairing and finite-simulation MC SE.

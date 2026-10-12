# Analysis workflow

These notebooks compare six imputation models for county youth suicide counts, then apply them to the observed program-adoption case. Each R Markdown file now starts with its purpose; each R chunk begins with one or two explanatory comment lines. The documentation describes the code as it currently stands, including exploratory chunks, rather than changing its computations.

## Where to start

For the main placebo workflow, read **1 → 2 → 3.0 → 3.1 → 3.2**. Step 1 generates assignments; step 2 fits models; step 3.0 prepares saved predictions; step 3.1 examines pre-fit reliability; step 3.2 produces cumulative performance summaries. Step 3.3 is an exploratory branch, not a prerequisite. Step 4 is the observed-adoption case study, independent of the placebo assignment branch. The `with_covariates/` files repeat fitting, preparation and evaluation with additional predictors; `archive/` holds earlier experiments.

## Backdated held-out diagnostic

[2.1.refit placebo for heldout error.qmd](2.1.refit%20placebo%20for%20heldout%20error.qmd) is an optional branch after the main fits. It backdates finite placebo starts by three years (configurable), reuses the six model specifications, rebuilds GB-Ens features under the new mask, and joins county-level held-out summaries to the unchanged saved main-fit predictions. It verifies assignment dates and observed input values before fitting each model/scenario. It does not test associations between held-out and post-period errors yet.

Rendering defaults to a preflight and writes `data/fits/backdate/eligibility.csv`. Set `run_fits: true` to run the backdated fits and save `heldout.rds`, `heldout.csv`, summary checkpoints and failure counts in `data/fits/backdate/`. `simulations`, `scenarios` and `models` can restrict a first run. The helper file is `../R/backdate_helpers.R`; synthetic checks are in `checks/tests/check_helpers.R`. The shared fixed-effects bootstrap now consistently fits its `.dt` argument, without changing any saved main outputs.

## Data and quantities

- `i` identifies counties, `r` replications, `s` years since observed or simulated adoption, `dgp` the assignment scenario, and `nc` the county group. Scenarios 1–3 are MAR (history-based hazard), latent level and latent trend.
- Model labels are TWFE, ETWFE, RW, RW-IID, RW-ST and GB-Ens. Adjusted models add standardized `uemp`, `pov`, `mhinc` and `aian` predictors.
- `y0_hat` is the fitted untreated mean; `y0_post_*` columns initially hold latent mean draws. Performance and case-study scripts add Poisson noise for predictive count draws. Cumulative sums use existing row order.
- Main performance tables first average cumulative outcomes and draws across counties within each replication, then summarize errors across replications. Pre-fit plots instead compute county-error RMSE within percentile bins and average those RMSEs across replications; these are different summaries.
- Reliability selection pools median pre-treatment RMSE across model rows within county/replication, then applies a scenario/county-group 90th-percentile threshold. It is not a separate top-decile exclusion within each model and replication.
- `sq_diff` is MSE after `summ_score()`; **`se_diff` is the correct RMSE**. Use `se_diff` directly in publication tables. Taking the square root of a displayed three-decimal MSE recreates the Excel rounding bug. `intv` is interval score, not interval width.
- In the case study, the later transformation replaces untreated draws with cumulative effect draws `y - y0`; negative effects correspond to fewer observed events than predicted untreated events.

## Current execution details

The notebooks use repository-root-relative paths. Quarto executes from the project root. Saved assignments/main fits are in `data/simulation/`, backdated fits in `data/fits/`, and prepared A/B inputs in `data/derived/heldout/`. Paper-ready exports go to `results/figures/` and `results/tables/`; reports and diagnostic exports go to `_output/`. Shared support is directly in `R/`. The private input `data/notpublic/compact_dt.rds` is required by fitting scripts and is not part of the public illustrative data.

These are working analysis notebooks, not all standalone renderable pipelines. In particular, no-covariate step 3.0 does not save `sim_fits.rds`; adjusted TWFE fitting explicitly runs only worlds 2–3; step 3.3 and archived notebooks use legacy objects/helpers; the case-study table chunk references undefined `atb`, and its training chunks require a `long` column in the loaded results. Some older notes in step 1 describe a 20-replication pilot or an older filename, but the active code uses 200 replications and saves `simulated_assignments_with_meta.rds`.

Rendering fitting or assignment notebooks executes expensive work and writes saved results. The summary notebooks also generate random predictive draws and export files. During directory cleanup, 2.1 and 2c.1 were rendered with fitting disabled; no models or half-B confirmation analyses were run.

## File and chunk guide

The numbered items below follow each file's chunk order; they summarize the comments added inside the chunks.

### [1. generate placebo data.rmd](1.%20generate%20placebo%20data.rmd)

Generate 200 placebo assignments for each of three assignment mechanisms while keeping observed outcomes fixed. Save assignments and simulation metadata for the fitting scripts.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load the county-year panel, reconstruct untreated observations and choose history lags. Store the lag configuration in options(xpn.conf) for downstream helper functions.
3. Estimate observed-history hazards and latent county levels/trends, then sample 200 assignments per world. Combine MAR, latent-level and latent-trend assignments into a county x replication x world array.
4. Package assignment arrays with seeds, model specifications and data summaries. Save data/simulation/simulated_assignments_with_meta.rds for downstream model fitting.

### [2. fit models placebo.rmd](2.%20fit%20models%20placebo.rmd)

Fit TWFE, ETWFE, RW, RW-IID, RW-ST and GB-Ens under the saved placebo assignments without covariate adjustment. Save county/year predictions and draws for later evaluation.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load the original panel, rebuild untreated labels, configure history/spatial inputs and load placebo assignments.
3. Define fit_sim(): apply each saved assignment, fit replications in parallel, and save model/scenario RDS files. Retain only years before observed adoption, so withheld placebo outcomes remain verifiable.
4. Fit and save GB-Ens predictions for all three assignment scenarios using 200 replications and 200 draws.
5. Limit nested numerical threads, then fit and save RW-ST predictions for all three scenarios.
6. Fit and save the RW model for all three scenarios.
7. Fit and save RW-IID predictions for all three scenarios.
8. Define TWFE/ETWFE formulas, then fit and save TWFE predictions for all three scenarios.
9. Fit and save ETWFE predictions using the adoption-year and calendar-year fixed-effects formula.

### [3.0 summary placebo.rmd](3.0%20summary%20placebo.rmd)

Assemble the no-covariate replication results, attach pre-fit diagnostics and retain a common post-placebo sample for all six models. This is preparation, not the final performance-table calculation.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Combine saved model/scenario results and label models, scenarios, county groups and event time.
3. Compute county/model pre-treatment RMSE and join it back to the prediction rows. These unrounded diagnostics support the later reliability restriction.
4. Keep event years 0-4 with nonmissing predictions and six model rows per county/year/replication. This prepares the common comparison sample; this file has no final saveRDS() call.

### [3.1 prefit.rmd](3.1%20prefit.rmd)

Examine how cumulative prediction error varies with pre-treatment fit quality, especially at year 4 among noncore counties. Use this diagnostic to assess the motivation for the reliability restriction.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load the combined saved predictions and pre-treatment RMSE diagnostics.
3. Accumulate outcomes, fitted means and posterior draws within each county/model/replication. Use the existing row order; this chunk does not sort event times.
4. Bin counties by pre-fit RMSE percentile and summarize post-fit RMSE across replications. Also explore a GAM-based knee for TWFE at year 4 in the noncore latent-level scenario.
5. Plot noncore year-4 post-fit RMSE against pre-fit percentiles, with replication quantile bands. Export the diagnostic figure to results/figures/pre_fit_plot.png.

### [3.2 summary placebo.rmd](3.2%20summary%20placebo.rmd)

Evaluate cumulative untreated outcomes in the no-covariate placebo study for the reliability-restricted and full samples. Produce replication-aggregated performance tables and diagnostic figures for core and noncore counties.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load combined saved predictions and inspect a few county/replication rows.
3. Convert latent mean draws into predictive count draws by adding Poisson observation noise. This changes posterior columns in memory; it is needed for predictive coverage and scoring.
4. Accumulate outcomes, fitted means and posterior draws within each county/model/replication. Use the existing row order; this chunk does not sort event times.
5. Retain a commented-out option for conversion to rates per 100,000 population. No rate conversion is executed here; subsequent results use counts.
6. Build a shared county reliability flag using median pre-fit RMSE across model rows. Keep rows below the scenario/county-group 90th percentile of that pooled diagnostic.
7. Average cumulative outcomes and draws across retained counties within each replication, then score and export. The correct RMSE is se_diff; formatting to three decimals happens after summ_score().
8. Repeat county aggregation and scoring for the full sample, without the reliability restriction. Export the full-sample performance table; se_diff is already the correctly computed RMSE.
9. Define model colors, reshape replication scores, and build performance boxplots by scenario/horizon. The main plot helper selects absolute error, CRPS or interval score, rather than table-level RMSE.
10. Compare full and restricted samples graphically and write the selected performance figures. These calls reuse scored replication means rather than fitting any new model.
11. Define coverage summaries and plots using the binomial Monte Carlo SE sqrt(p*(1-p)/R). Display nominal 95% coverage and clipped normal-approximation Monte Carlo uncertainty bands.
12. Draw and export coverage figures for county groups and sample variants using the helpers above.

### [3.3 other analysis.rmd](3.3%20other%20analysis.rmd)

Explore placebo performance by pre-fit difficulty strata and alternative summaries. Later chunks retain older plotting/table experiments and depend on objects not created by the earlier chunks.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load saved predictions that already include pre-fit difficulty strata and inspect posterior columns.
3. Convert latent mean draws into predictive count draws by adding Poisson observation noise. This changes posterior columns in memory; it is needed for predictive coverage and scoring.
4. Accumulate outcomes, fitted means and posterior draws within each county/model/replication. Use the existing row order; this chunk does not sort event times.
5. Retain a commented-out option for conversion to rates per 100,000 population. No rate conversion is executed here; subsequent results use counts.
6. Inspect the distribution of the existing rmse_q difficulty labels before stratified aggregation.
7. Flag extreme pre-fit errors and score county means separately within each difficulty stratum. The proposed flag-based exclusion is commented out; the active path keeps all rows.
8. Plot replication-level absolute errors (sqrt of one squared discrepancy), CRPS and interval scores by difficulty. Although labelled RMSE in these plots, a single replication value is an absolute error.
9. Summarize coverage and Monte Carlo uncertainty by scenario, horizon, model and county group. This older helper pools difficulty-stratum rows, so its n is not necessarily the replication count.
10. Define an older four-metric plot helper with TWFE reference bands and a nominal-coverage line. Its expected se column and lowercase model names differ from the current scoring/label conventions.
11. Inspect and plot older cumulative results at years 0, 2 and 4. Requires sdt_cum and legacy column names from an earlier session; they are not created above.
12. Format older replication-level performance medians and IQRs into a summary table. Requires sdt_cum with legacy grouping columns; this is not the main appendix-table route.

### [4. the case.rmd](4.%20the%20case.rmd)

Apply the six imputation models to observed program adoption and summarize cumulative effects, county heterogeneity and associations with long-training intensity. Compare full and reliability-restricted samples.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load the case-study panel, rebuild untreated labels, and configure history lags and the county adjacency graph.
3. Define a wrapper that fits one case-study model with 400 draws and saves its prediction rows.
4. Fit and save all six case-study models using observed adoption rather than placebo assignments.
5. Load saved case-study predictions, combine models and inspect missing posterior draws.
6. Compute pre-fit RMSE, add Poisson predictive noise and cumulate post-adoption outcomes/draws for years 0-4. Apply the pooled reliability rule, label models and save data/all_case_study.rds.
7. Reload cumulative case-study predictions and replace untreated predictions/draws with effects y - y0. After this transformation, posterior columns represent cumulative effects rather than untreated outcomes.
8. Summarize draw-wise county-mean cumulative effects and 95% quantiles for full/restricted samples.
9. Plot and export noncore mean cumulative effects across horizons, models and sample variants.
10. Export an older case-study table from atb. The atb object and event_time column are not created in this file; this chunk needs prior session state.
11. Compute draw-wise interquartile ranges of effects across counties as a robust heterogeneity measure. Summarize their posterior means and 95% quantiles for full/restricted samples.
12. Plot and export noncore county-effect heterogeneity using the IQR summaries.
13. Compute county mean long-training share and posterior Spearman correlations with cumulative effects. Summarize these associations for full/restricted samples and define a >0.75 training-share group.
14. Plot noncore training-effect correlations and export the correlation summary table.
15. Compare posterior mean cumulative effects for counties above/below the 0.75 long-training-share threshold. Plot the noncore comparison separately for full and reliability-restricted samples.

### [with_covariates/2.1 fit models placebo with covariates.rmd](with_covariates/2.1%20fit%20models%20placebo%20with%20covariates.rmd)

Repeat placebo model fitting with standardized socioeconomic/demographic covariates. Save covariate-adjusted prediction rows for comparison with the outcome-history-only models.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load the panel, spatial graph and placebo assignments; standardize uemp, pov, mhinc and aian. Configure the same history lags used by the no-covariate branch.
3. Define fit_sim(): apply each saved assignment, fit replications in parallel, and save model/scenario wcov_imputation RDS files. Retain only years before observed adoption, so withheld placebo outcomes remain verifiable.
4. Fit and save GB-Ens for all three scenarios, adding the four standardized covariates as predictors.
5. Limit nested numerical threads, then fit and save the covariate-adjusted RW-ST model.
6. Fit and save the covariate-adjusted RW model for all three scenarios.
7. Fit and save covariate-adjusted RW-IID predictions for all three scenarios.
8. Define covariate-adjusted TWFE/ETWFE formulas and fit TWFE for scenarios 2 and 3 only. The active worlds=2:3 call does not generate MAR TWFE results in this chunk.
9. Fit and save covariate-adjusted ETWFE predictions for all three assignment scenarios.

### [with_covariates/3.0.1 summary placebo wcov.rmd](with_covariates/3.0.1%20summary%20placebo%20wcov.rmd)

Combine covariate-adjusted replication predictions, attach pre-fit diagnostics and build the common six-model comparison sample. Save the prepared results for adjusted diagnostics and performance summaries.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Combine saved model/scenario results and label models, scenarios, county groups and event time.
3. Compute county/model pre-treatment RMSE and join it back to the prediction rows. These unrounded diagnostics support the later reliability restriction.
4. Keep event years 0-4 with nonmissing predictions and six model rows per county/year/replication.
5. Save the prepared covariate-adjusted predictions and pre-fit diagnostics to data/simulation/sim_fits_wcov.rds.

### [with_covariates/3.1.1 prefit.rmd](with_covariates/3.1.1%20prefit.rmd)

Repeat the pre-fit versus post-fit error diagnostic for covariate-adjusted models. Compare year-4 noncore performance across pre-fit percentile bins and assignment scenarios.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load the combined covariate-adjusted saved predictions and pre-treatment RMSE diagnostics.
3. Accumulate outcomes, fitted means and posterior draws within each county/model/replication. Use the existing row order; this chunk does not sort event times.
4. Bin counties by pre-fit RMSE percentile and summarize post-fit RMSE across replications. Also explore a GAM-based knee for TWFE at year 4 in the noncore latent-level scenario.
5. Plot noncore year-4 post-fit RMSE against pre-fit percentiles, with replication quantile bands. Export the diagnostic figure to results/figures/pre_fit_plot.png.

### [with_covariates/3.2.1 summary placebo.rmd](with_covariates/3.2.1%20summary%20placebo.rmd)

Evaluate cumulative untreated outcomes for covariate-adjusted placebo models. Produce full/restricted sample performance tables and figures analogous to the no-covariate analysis.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load combined covariate-adjusted predictions and inspect a few county/replication rows.
3. Convert latent mean draws into predictive count draws by adding Poisson observation noise. This changes posterior columns in memory; it is needed for predictive coverage and scoring.
4. Accumulate outcomes, fitted means and posterior draws within each county/model/replication. Use the existing row order; this chunk does not sort event times.
5. Retain a commented-out option for conversion to rates per 100,000 population. No rate conversion is executed here; subsequent results use counts.
6. Build a shared county reliability flag using median pre-fit RMSE across model rows. Keep rows below the scenario/county-group 90th percentile of that pooled diagnostic.
7. Average cumulative outcomes and draws across retained counties within each replication, then score and export. The correct RMSE is se_diff; formatting to three decimals happens after summ_score().
8. Repeat county aggregation and scoring for the full sample, without the reliability restriction. Export the full-sample performance table; se_diff is already the correctly computed RMSE.
9. Define the covariate-analysis colors, score reshaping and performance boxplot helper. Use these helpers to compare cumulative prediction errors across sample variants.
10. Compare full and restricted samples graphically and write the selected performance figures. These calls reuse scored replication means rather than fitting any new model.
11. Define coverage summaries and plots using the binomial Monte Carlo SE sqrt(p*(1-p)/R). Display nominal 95% coverage and clipped normal-approximation Monte Carlo uncertainty bands.
12. Draw and export coverage figures for county groups and sample variants using the helpers above.

### [archive/0.1 test.rmd](archive/0.1%20test.rmd)

Prototype assignment simulation, history expansion and hazard estimation on a small synthetic county panel. This archived scratchpad mixes toy and private-data experiments and is not the production simulation pipeline.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Construct a ten-county, ten-year synthetic panel with covariates, ring adjacency and Poisson outcomes. Configure history/spatial inputs and generate untreated-history features for testing.
3. Simulate and visualize placebo assignments on the toy panel, excluding never-treated units from assignment.
4. Try private-panel loading, history expansion and alternative hazard learners. This scratch chunk mixes toy/private inputs and requires id_col, start_year_col and Lset from session state.

### [archive/0.2. testing .rmd](archive/0.2.%20testing%20.rmd)

Explore covariate specifications and alternative imputation/inference methods before the main workflow. This archived notebook mixes single-assignment diagnostics with expensive legacy replication fits.

1. Load packages and shared functions from R/, and set the random seed and display options.
2. Load the private panel, define a held-out untreated test set and apply the first saved placebo assignment.
3. Standardize socioeconomic/demographic covariates and label core versus noncore counties.
4. Compare TWFE and ETWFE on withheld outcomes using the shared scoring functions. Round only the resulting displayed summary, after RMSE has been computed.
5. Compare TWFE specifications with covariates and additional state/county-type time interactions.
6. Compare ETWFE covariate specifications and summarize withheld-outcome performance.
7. Explore history-based XGBoost predictions with/without extra covariates and sticky assignment generation. The final comparison references xgb2, which is not defined in this notebook.
8. Compare RW, RW-IID and RW-ST, each with and without covariates, on the held-out test set.
9. Define two alternative spatial-temporal INLA specifications for experimentation. The second fit_spt definition replaces the first; graph/time inputs rely on session configuration.
10. Run and time 200 legacy forest/imputation replications, then save the fitted results. Requires placebo_dt and legacy forest helpers not established in the preceding chunks.
11. Run and save a legacy TWFE inference experiment over 200 replications. Expects a replication-indexed fit object and the bbx_inference_nox helper.
12. Load forest results and run/save the analogous ETWFE inference experiment.
13. Run and save a legacy fit_ime INLA experiment over 200 replications.
14. Run and save a legacy autoregressive INLA experiment over 200 replications.
15. Estimate replication-specific adjacency inputs, run the legacy spatial INLA model and save results. Requires W_est(), the legacy fit object and graph configuration from prior session state.


## Held-out exploration and confirmation

`3.0 prepare placebo results.rmd` now has a separate preparation chunk combining backdated summaries without changing its original prediction object. It saves the fixed county split and physically separate A/B inputs. The same preparation can run alone via `Rscript R/heldout_prepare.R`.

`3.1a.1_heldout_compute.qmd` compares all diagnostics on A, computes within-simulation correlations/AUC for realized and county-averaged targets, plots the primary targets with county-bootstrap bands, and exports cutoff tables. It never loads B diagnostic rows or writes a decision.

`3.1b heldout confirmation.qmd` validates the hand-written decision before reading B, then evaluates that rule alone. Do not render it until the choice is recorded. Details and decision syntax: [held-out README](../heldout/README.md).

## Covariate counterpart of the held-out workflow

[with_covariates/2c.1. refit placebo for held out error.qmd](with_covariates/2c.1.%20refit%20placebo%20for%20held%20out%20error.qmd) reuses 2.1's shared implementation, changing only covariate preprocessing/specifications and input/output branches. `3c.0` prepares separate covariate inputs with the original county split. `3c.2` applies the hand-written no-covariate rule to each covariate model separately, on B by default, without retuning. `3c.1` can be skipped. See the [covariate workflow README](with_covariates/README.md) for commands and dependencies.


`with_covariates/3c.1b heldout confirmation.qmd` evaluates the unchanged no-covariate decision on the original half-B counties using covariate results. It reuses the original assignment, loads only covariate half B, and writes separate confirmation tables. No covariate half-A exploration is run. Render only after the covariate prepared inputs and the hand-written decision file exist.

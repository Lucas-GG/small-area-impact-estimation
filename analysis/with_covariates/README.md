# Covariate branch with the fixed no-covariate reliability rule

The main model specifications and historical fit outputs stay unchanged. The new backdated notebook reuses shared helpers; its differences from 2.1 are covariate preprocessing, model specifications, main-fit paths and output paths. Covariates `uemp`, `pov`, `mhinc` and `aian` are scaled over the original panel exactly as in 2c. GB-Ens rebuilds masked-outcome features and uses these four predictors in fitting, propensity estimation and the leakage audit.

1. `2c. fit models placebo with covariates.rmd` is the existing main fitting notebook. **Do not rerun it for this addition.** All 18 saved main covariate model/scenario files currently exist in `data/simulation/wcov/`, including TWFE MAR (although the current main notebook's active TWFE call lists only scenarios 2 and 3).
2. `2c.1. refit placebo for held out error.qmd` runs only second, backdated fits, using the same saved assignment array, default three-year window, seeds and draw settings as 2.1. Its default `run_fits: false` runs a preflight. Enable fitting explicitly with:

   ```sh
   quarto render "analysis/with_covariates/2c.1. refit placebo for held out error.qmd" -P run_fits:true
   ```

   Summaries and both mid-p tails go to `data/fits/backdate_wcov/`. No full draws are persisted. Original main fits and no-covariate backdated results are never overwritten.
3. `3c.0 prepare placebo results.rmd` reads the saved main covariate files from `data/simulation/wcov/` (historical `data/` fallback). Newly prepared prediction rows go to `data/derived/heldout/with_covariates/sim_fits_wcov.rds`, preserving the historical `data/simulation/sim_fits_wcov.rds`. Its added final chunk combines backdated covariate checkpoints and reuses the existing `data/derived/heldout/county_assignment.rds`. It never creates a second county split. If checkpoints are not available yet, it reports the dependency.

   You can prepare only the small held-out tables, using the already saved main combined prediction file, without rerunning the large 3c.0 assembly:

   ```sh
   Rscript -e 'source("R/backdate_helpers.R"); source("R/heldout_prepare.R"); prepare_heldout(find_backdate_root(), with_covariates=TRUE)'
   ```

4. Skip `3c.1 prefit.rmd` when using the transferred rule. Keep it only as historical exploration; do not use it to retune the rule.
5. `3c.1b heldout confirmation.qmd` retains its existing filename and is now the **covariate half-B estimate** report. It reuses `R/heldout_estimate.R` with `with_covariates=TRUE`, the fixed no-covariate decision and the original county assignment. It evaluates the recorded strict |z| > 2 rule and only the pre-specified mid-p < .10 and pre-RMSE top-10% sensitivities. Available event years 0–4 with k ≥ 1, k-stratified county associations, missing diagnostics kept, and Poisson-noised predictive draws use the same shared formulas as 3.1b. Tables cover noncore own populations, all-half populations, TWFE-eligible E, common screened samples, ranks, unique flag contributions and model-averaged gain estimates/MC SEs. Covariate A is the matching descriptive reference; original no-covariate A is reported separately. No rule is retuned and no model is fitted.

   Outputs are isolated under `_output/heldout/with_covariates/confirmation/`. Only the requested half's snapshot is prepared if needed; the other half's snapshot is not created. Saved mixed checkpoints/main fits are deserialized and immediately restricted to the requested assigned IDs before computation. Every input path and metadata access is logged. The A dry run uses an isolated render project:

   ```sh
   quarto render _output/heldout/with_covariates/confirmation/render_A/covariate_half_A_dry_run.qmd -P half:A
   ```

   The report defaults to B, but B execution is gated separately from the completed no-covariate B run: the user must explicitly authorize covariate B, and its branch-specific authorization marker must then exist. Before B is opened, the committed decision bytes and assignment hash must match the covariate A dry run. The evaluated decision, its MD5/commit provenance, and input logs are saved beside the half-specific estimates.

6. `3c.2 summary placebo.rmd` requires the hand-written **no-covariate** `analysis/3.1a.4_heldout_decision.yml` and all 18 prepared covariate checkpoints. It applies that diagnostic and **the same numeric cutoff** to each covariate model's own diagnostics, separately by county/simulation/model. It does not recompute a 90th-percentile cutoff. Missing diagnostics cannot qualify for the restricted sample and are counted separately.

The default evaluation half in 3c.2 is **B**, reusing the original fixed split. This makes the covariate analysis an evaluation of the transferred rule, without new exploration on B. Full-sample tables mean all available counties in B; restricted tables contain each model's passing counties. Their county composition can differ across models, so differences reflect both predictions and retained samples. This is the explicitly selected model-specific application, rather than a common intersection across models. `fixed_rule_counts.csv` and `applied_decision.yml` record retention and the decision file's hash. For p-values, exclusion is `p_two < cutoff`; other diagnostics use `abs(diagnostic) > cutoff`, with equality retained.

Validation: 2c.1 was rendered with fitting disabled. No new model fits or simulations were run. The configuration/scaling and fixed-rule tests passed. 3c.0, 3c.1b and 3c.2 were syntax-checked, not executed; they need the new fits and your decision file. The preflight emitted the existing INLA/macOS `/bin/kstat` notices, but rendered successfully.

Rendering uses Quarto HTML settings and restores the source directory before the final renderer reopens relative source paths. To render without fitting, pass `-P run_fits:false`: completed saved summaries are displayed if available; otherwise it runs only preflight. This prevents a failed HTML render from requiring another fit run.

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
5. `3c.1b heldout confirmation.qmd` is the covariate counterpart of 3.1b. It reads only `data/derived/heldout/with_covariates/half_B.rds`, checks the original saved assignment/hash, and uses the same no-covariate `analysis/3.1a.4_heldout_decision.yml`. It evaluates only the recorded diagnostic/cutoff against covariate realized and county-averaged targets, with within-simulation Spearman correlations, AUC, replication counts and exclusion-performance tables. Outputs go separately to `_output/heldout/with_covariates/confirmation/`. It does not load half A, rerun exploration, choose a new rule or fit any model. Missing decisions or inputs stop execution; partial checkpoints are explicitly provisional.

   Once the covariate inputs and the fixed decision exist, render it with:

   ```sh
   quarto render "analysis/with_covariates/3c.1b heldout confirmation.qmd"
   ```

6. `3c.2 summary placebo.rmd` requires the hand-written **no-covariate** `analysis/3.1a.4_heldout_decision.yml` and all 18 prepared covariate checkpoints. It applies that diagnostic and **the same numeric cutoff** to each covariate model's own diagnostics, separately by county/simulation/model. It does not recompute a 90th-percentile cutoff. Missing diagnostics cannot qualify for the restricted sample and are counted separately.

The default evaluation half in 3c.2 is **B**, reusing the original fixed split. This makes the covariate analysis an evaluation of the transferred rule, without new exploration on B. Full-sample tables mean all available counties in B; restricted tables contain each model's passing counties. Their county composition can differ across models, so differences reflect both predictions and retained samples. This is the explicitly selected model-specific application, rather than a common intersection across models. `fixed_rule_counts.csv` and `applied_decision.yml` record retention and the decision file's hash. For p-values, exclusion is `p_two < cutoff`; other diagnostics use `abs(diagnostic) > cutoff`, with equality retained.

Validation: 2c.1 was rendered with fitting disabled. No new model fits or simulations were run. The configuration/scaling and fixed-rule tests passed. 3c.0, 3c.1b and 3c.2 were syntax-checked, not executed; they need the new fits and your decision file. The preflight emitted the existing INLA/macOS `/bin/kstat` notices, but rendered successfully.

Rendering uses Quarto HTML settings and restores the source directory before the final renderer reopens relative source paths. To render without fitting, pass `-P run_fits:false`: completed saved summaries are displayed if available; otherwise it runs only preflight. This prevents a failed HTML render from requiring another fit run.

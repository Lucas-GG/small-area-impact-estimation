# Held-out evaluation

1. Run `Rscript R/heldout_prepare.R` from the repository root after backdated checkpoints are saved. This never fits a model. For the revised no-covariate exploration it requires the existing assignment and all 18 checkpoints, immediately selects A, converts paired legacy inclusive tails exactly to mid-p, and writes `half_A_available.rds`. It preserves the legacy A/B inputs and caches A-only post summaries in `data/derived/heldout/post_A/`. The covariate preparation retains its existing separate workflow. A legacy upper tail alone cannot be converted and is marked unavailable. New fits save both mid-p tails with a convention marker.
2. The original preparation created `county_assignment.rds` and `.csv` with seed 2032026. Stratification uses quintiles of saved MAIN TWFE MAR mean pre-period expected counts (equal weight to each fit's pre-period mean). County IDs break size ties; a seeded shuffle balances A/B within each stratum. Never delete/recreate this assignment after exploration. Refreshing inputs reuses it.
3. Render `analysis/3.1a.1_heldout_compute.qmd` only when intentionally refreshing half-A computations. It now includes county-level statistics, aggregate screening (10), the TWFE-eligible sensitivity (10.1), and common screened sample (10.2). Then render `analysis/3.1a.2_heldout_exploration.qmd` and `analysis/3.1a.3_heldout_details.qmd` for the summary and details. HTML lives in `_output/analysis/`; CSV downloads stay under `_output/analysis/output/3.1a/`.
4. After choosing a rule using A, write `analysis/3.1a.4_heldout_decision.yml` **by hand**. Example syntax (this is only an example, not a chosen rule):

```yaml
diagnostic: ho_uncertainty
cutoff: 2
target_definition: "Cumulative available-window standardized error; k >= 1"
missing_diagnostic_handling: keep
date: "2026-10-09"
rationale: "Write the rationale for the choice made using half A here."
```

Allowed diagnostics: `pre_rmse`, `pre_rms_pearson`, `ho_raw`, `ho_standardized`, `ho_uncertainty`, `p_two`. The exclusion direction is absolute diagnostic > cutoff, except `p_two < cutoff`.
5. Before confirmation, align its targets, k-stratified evaluation and missing policy with the revised exploration. Only after that implementation and the choice are fixed, render `analysis/3.1b heldout confirmation.qmd`. It checks the decision before reading B and evaluates only the chosen rule. Missing/invalid decisions cause an error. This document was deliberately not rendered during implementation.

CSV outputs include within-simulation Spearman correlations and AUC for realized and county-averaged absolute targets, their summaries, per-county average-target replication counts, cutoff performance and availability. Top-decile boundaries use all finite target rows within the relevant county half; boundary ties are retained. Pairwise diagnostic availability affects the evaluated rows, not the target-class definition. Comparisons with differing denominators require care. County-averaged targets include the current simulation and are explanatory, not independent predictive validation.

When checkpoint files are produced by an already-running R process, editing the helper on disk does not change functions loaded in that process. Its paired inclusive tails are still recoverable during preparation. A checkpoint containing only one inclusive tail needs the missing predictive information (or a new backdated fit) to supply the mid-p diagnostic; the analysis never guesses it and never refits automatically.

The 2026-10-11 revision is documented in [the revision notes](checks/heldout_revision/README.md). It keeps complete-window and county-average complete-window targets as sensitivities.

## Readable presentation pages

`3.1a.2_heldout_exploration.qmd` is the short visual half-A report; `3.1a.3_heldout_details.qmd` contains complete interactive report tables and all CSV downloads. They read existing saved statistical outputs; the computation step remains the computation/refresh document. The only numerical correction in the presentation workflow is the random band (quantiles after simulation averaging), implemented separately in `R/heldout_random_band.R`. Main targets, ranks, screens and other statistics are unchanged. Both new pages write/use byte-identical CSV copies under `_output/analysis/output/3.1a/`. `DT` is installed in the project-local `_output/R-library` and loaded by the presentation helper. Render both pages explicitly with `quarto render analysis/3.1a.2_heldout_exploration.qmd` and `quarto render analysis/3.1a.3_heldout_details.qmd`; the project directs them to `_output/analysis/`. Keep their HTML together and their relative `output/3.1a/` CSV links alongside them. HTML and downloadable CSVs are under `_output/analysis/`. No decision is selected.


## File order and generated artifacts

- `3.1a.1_heldout_compute.qmd`: all half-A computations, including sections 10, 10.1 and 10.2.
- `3.1a.2_heldout_exploration.qmd`: render the visual summary.
- `3.1a.3_heldout_details.qmd`: render full tables and downloads.
- `3.1a.4_heldout_decision.yml`: record the choice manually before 3.1b confirmation.

Quarto HTML and document support files go to `_output/`; shared libraries go to `_output/libs/` and caches go to `.quarto/`. Quarto manages temporary render files and removes them from the source tree after successful rendering. Verification source files are in `checks/`, with generated audit artifacts in `_output/checks/`. Previous source-adjacent render artifacts are preserved in `_output/source_artifacts/`.


The additive aggregate analysis is in `R/heldout_aggregate.R`; acceptance checks are `Rscript checks/tests/check_aggregate.R`. It preserves the existing CSVs and decision file, writes only new `aggregate_*.csv` exports, reuses A-only seeded predictive draws, and stops if noncore metadata or A-only/hash checks fail. Figures and scorecards are section 10 of the main report. The full replication ledger, paired MC SEs, random-exclusion bands, oracle availability, population ranges and matched-window manuscript comparisons are in the details report and downloadable CSVs. The oracle uses the other simulations only and is unavailable when a county lacks at least two such observations. Pointwise exploratory minimum flags are not selection-adjusted.

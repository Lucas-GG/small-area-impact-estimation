# Supplement S3: held-out screening of poorly predicted counties

Asks whether pre-treatment diagnostics can flag counties whose counterfactuals are predicted poorly, and whether excluding them improves the noncore aggregate estimate. A screening rule was chosen on a random half of the counties (A) and its gain estimated on the other half (B). Helpers are in `R/screening/`. Run from the repository root, after the main placebo fits (`02_fit_placebo.qmd`, both specifications).

| Notebook | Does | Writes |
|---|---|---|
| `s3_01_backdate_refit.qmd` | Moves every placebo start three years earlier, refits the four S3 models (TWFE, ETWFE, GB-Ens, nested RW-IID; RW and RW-ST dropped 2026-10-21) with those years masked, and summarizes the held-out predictions per county. Then rebuilds the half-A/B inputs. `-P covariates:true` for the covariate models. Fits only with `-P run_fits:true` (~1 day for both). | `data/fits/backdate[_wcov]/`, `data/derived/heldout/` |
| `s3_02_half_A.qmd` | Half A: the fixed rule, its results and sensitivities, and the exploration that chose it. `-P refresh:true` recomputes the saved half-A results first. | `_output/heldout/`, `_output/analysis/output/3.1a/` |
| `s3_decision.yml` | The rule recorded by hand from half A before half B was opened (|held-out standardized error| > 2, missing diagnostics kept). Never edited. | — |
| `s3_03_half_B.qmd` | Half B: the gain from the fixed rule without and with covariates, then the county-influence decomposition. Rendered explicitly only. | `_output/heldout/confirmation/`, `_output/heldout/with_covariates/confirmation/`, `_output/heldout/influence/` |

```sh
quarto render analysis/S6_screening/s3_01_backdate_refit.qmd -P run_fits:true
quarto render analysis/S6_screening/s3_01_backdate_refit.qmd -P run_fits:true -P covariates:true -o s3_01_backdate_refit_covariates.html
quarto render analysis/S6_screening/s3_02_half_A.qmd -P refresh:true
quarto render analysis/S6_screening/s3_03_half_B.qmd
```

## Design notes

- **Split.** `data/derived/heldout/county_assignment.rds` (seed 2032026) stratifies counties by quintiles of mean pre-period expected count from the main TWFE MAR fits, with a separate missing-size stratum. Models are fit on all counties; the split governs only evaluation. Never recreate it.
- **Halves.** Half-A and half-B inputs are physically separate files; the half-A report asserts it reads A only. Half B could not be read until the decision file was committed (checked by MD5 and git).
- **Correction after the decision.** The placebo fits were corrected (models now see the simulated adoption year as their cohort) after the rule was recorded. Both halves are evaluated on the corrected fits; the rule was not revised.
- **Targets.** County level: absolute error of the cumulative count over available event years 0–4, divided by the square root of the predicted total. Aggregate: signed mean cumulative error over retained noncore counties, with RMSE across replications and coverage of the 95% aggregate predictive interval.
- **Mid-p tails.** `p_upper = P(Y > y) + 0.5 P(Y = y)`, `p_lower = P(Y < y) + 0.5 P(Y = y)`, `p_two = min(1, 2 min(p_upper, p_lower))`; paired legacy inclusive tails are converted exactly, a lone tail is unavailable.

Tests for this branch are in `checks/tests/` (`check_aggregate.R`, `check_twfe_common.R`, `check_common_screen.R`, `check_half_estimate.R`, `check_heldout_influence.R`, …).

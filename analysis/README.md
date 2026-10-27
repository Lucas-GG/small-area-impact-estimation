# Replication guide

Run the notebooks from the repository root, in this order. Each one names its inputs and outputs at the top. All notebooks load the helpers in `R/` (top level only; supplement S6 also loads `R/screening/`, S7 `R/priors/`).

| Step | Notebook | Does | Reads | Writes | Runtime (18 cores) |
|---|---|---|---|---|---|
| 1 | `01_generate_assignments.qmd` | 200 placebo adoption schedules for each of three assignment mechanisms (MAR, latent level, latent trend) | panel | `data/simulation/simulated_assignments_with_meta.rds` | minutes |
| 2 | `02_fit_placebo.qmd` | Placebo fits of every reported model: the six main models, the S5 covariate variants and the S8 pooled-cohort and intermediate models (`-P covariates:true` for the covariate specification) | panel, assignments, `data/W` | `data/simulation/placebo/{wocov,wcov}/` | ~15 h per specification (RW-ST dominates; the S8 models double it) |
| 3a | `03a_prepare_placebo.qmd` | Common six-model sample, event years 0–4, both specifications; the S5 samples | step 2 | `data/simulation/placebo_sample_*.rds` | ~15 min |
| 3b | `03b_placebo_results.qmd` | Placebo performance: Table 2, Figure 1, Tables S1–S3, Figures S1–S5, S8 | step 3a | `results/`, `_output/analysis/output/03b/` | minutes |
| 4 | `04_injected_effects.qmd` | Injected effects and the heterogeneity test (no refits): Figure D1, Table S4, Figure S6 | step 3a | `results/`, `_output/analysis/output/04/` | ~1 h |
| 5 | `05_case_study.qmd` | Case study: Table 3, Figure 2, Tables S5–S7, Figure S7 | panel, `04` null distributions | `results/`, `data/derived/case_study/` | ~1 h |
| S6 | `S6_screening/` | Screening counties by pre-adoption fit (Table S8); see [S6_screening/README.md](S6_screening/README.md) | step 2 fits | `data/fits/`, `data/derived/heldout/` | ~1.5 days |
| S7 | `S7_priors/` | Sensitivity to the variance-component priors (Table S9): `s4_01_fit_candidates.qmd` fits the ten prior settings (`-P run_fits:true`), `s4_05_prior_sensitivity.qmd` tabulates them | assignments | `data/derived/prior_selection/` | ~13 h per specification |
| S8 | `S8_pooled_cohorts/` | Pooled vs cohort × group cells (Table S10), ETWFE–TWFE equivalence; see [S8_pooled_cohorts/README.md](S8_pooled_cohorts/README.md) | steps 2–5 | `results/`, `data/derived/pooled_cohorts/` | ~30 min |

File prefixes in `S6_screening/` (`s3_`) and `S7_priors/` (`s4_`) keep the supplement numbers of an earlier draft.

## Conventions

- `i` is the county, `r` the placebo replication, `year` the calendar year, `s` years since (simulated) adoption, `dgp` the assignment scenario, `nc` the county group (noncore or core).
- `start_year` is the real adoption year and `start_year_sim` the placebo one. Models never see either: `set_y0()` builds the masked outcome `y0` and the analyst-visible `cohort` (the simulated year in the placebo, the real year in the case study), and `fit_blind()` fits each model on a copy without the real years or raw outcomes.
- `y0_hat` is the fitted untreated mean and `y0_post_1`…`y0_post_200` are latent mean draws; summaries add Poisson noise to obtain predictive counts.
- Saved fits contain only years before the real adoption, so every placebo prediction can be checked against an untreated outcome.

## Data

The notebooks read the county-year panel from `panel_file()` (`R/aux.r`): the restricted panel `data/notpublic/compact_dt.rds` by default, or the synthetic panel with `SAI_PANEL=data/synthetic_dt.rds` (see the top-level README). `data/W` is the county adjacency matrix. `01` writes the placebo schedules to `data/simulation/simulated_assignments_with_meta.rds`; every later step reads them, so the placebo fits reproduce only from the same schedules.

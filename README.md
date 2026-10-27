# small-area-impact-estimation

Code for "Estimating intervention effects in small areas: comparing statistical approaches in rural youth suicide prevention". The paper compares six models for predicting untreated county outcomes from sparse counts (TWFE, ETWFE, RW, RW-IID, RW-ST and GB-Ens) in a placebo study with reassigned adoption times, injects known effects to assess power and a test of effect heterogeneity, and applies the models to the adoption of the Garrett Lee Smith (GLS) youth suicide prevention program in noncore counties.

## Data

The analysis panel holds annual youth (10–24) suicide deaths, youth population, GLS adoption year, program activity and four county covariates for 516 counties, 1999–2016. **It is not distributed**: county death counts below 10 are restricted by NCHS small-cell rules, and the GLS program data may not be released. [TODO: how to request access.]

`data/synthetic_dt.rds` lets the code run without the restricted data. It has the same structure; population, urbanicity and the covariates are the public Census values; adoption years, program activity and deaths are generated, with a known program effect (rate ratio 0.85 after adoption). It does not reproduce the paper's results. `data-raw/make_synthetic.R` documents exactly what is real and what is generated. `data/W` is the county adjacency matrix.

## Layout

- `analysis/`: the pipeline, in run order (see [analysis/README.md](analysis/README.md)); supplements in `analysis/S6_screening/`, `analysis/S7_priors/`, `analysis/S8_pooled_cohorts/`.
- `R/`: model, feature, scoring and table helpers (`R/screening/`, `R/priors/`: helpers of supplements S6 and S7).
- `results/`: the exhibits of the paper and supplement, named by their number, and the numbers quoted in the text (see [results/README.md](results/README.md)).
- `data-raw/`: the script that builds the synthetic panel.

Intermediate files go to `data/simulation/`, `data/derived/` and `_output/` (not tracked).

## Running

Requirements: R 4.5 with data.table, fixest, INLA (from the R-INLA repository, <https://www.r-inla.org>), xgboost, lme4, Matrix, matrixStats, MASS, ggplot2, gt, tidyverse, knitr, rmarkdown, yaml, xml2 and zip; Quarto 1.9. The results were produced with R 4.5.3, INLA 25.10.19, fixest 0.14.2, xgboost 3.2.1.1, data.table 1.18.4, lme4 2.0.1 and Quarto 1.9.36 on macOS (20 cores, 128 GB).

Run every notebook from the repository root, in the order of [analysis/README.md](analysis/README.md); for example `quarto render analysis/02_fit_placebo.qmd -P covariates:true`. The full placebo study takes about two days on 18 cores (the RW-ST fits dominate).

With the restricted panel at `data/notpublic/compact_dt.rds`, the pipeline reproduces the reported numbers. Every fit is seeded per replication, with one exception: the saved fits of TWFE, pooled-cohort ETWFE and the two covariate variants of Supplement S5 come from an earlier, unseeded version of `02_fit_placebo.qmd`, so a refit reproduces their point predictions (and RMSE) exactly but their coverage, interval score and CRPS only within Monte Carlo error.

On the synthetic panel, set `SAI_PANEL` and, for a quick check, reduce the replications and draws (about an hour on 18 cores):

```
export SAI_PANEL=data/synthetic_dt.rds
quarto render analysis/01_generate_assignments.qmd -P nrep:10
quarto render analysis/02_fit_placebo.qmd -P nrep:10 -P nboot:20 -P 'models:["twfe","etwfe","rw","rwi","rwst","xgb"]'
quarto render analysis/02_fit_placebo.qmd -P covariates:true -P nrep:10 -P nboot:20 -P 'models:["twfe","etwfe","rw","rwi","rwst","xgb","twfe_v","etwfe_v","etwfe_pooled"]'
quarto render analysis/03a_prepare_placebo.qmd
quarto render analysis/03b_placebo_results.qmd -P refresh:true
quarto render analysis/04_injected_effects.qmd -P refresh:true
quarto render analysis/05_case_study.qmd -P refit:true -P nboot:20
```

These runs overwrite `results/`; `git checkout results` restores the paper's exhibits.

## Licence

[TODO: licence. Proposed: MIT for the code, CC BY 4.0 for `results/` and the synthetic data.]

# Supplement S8: cohort fixed effects pooled across county groups

The main analyses define the adoption-cohort effects of ETWFE, RW, RW-IID, RW-ST and the GB-Ens baseline separately for noncore and other counties (cohort × group cells). This supplement compares them with cohort effects pooled across the two groups, and checks how far ETWFE agrees with TWFE. Run from the repository root, after the main pipeline (`02`, `03a`, `03b`, `05`).

| Script | Does | Writes |
|---|---|---|
| `s8_01_score.R` | Scores the pooled-cohort models (`*_pooled` in `02_fit_placebo.qmd`) and the intermediate RW (cohort × group intercepts, one deviation walk per cohort; `rw_cohort_walks`) on the Table 2 evaluation county-years; paired change in RMSE from pooled cohorts to cells for every model; fits pooled ETWFE, pooled RW and the intermediate to the case study | `data/derived/pooled_cohorts/{ledger.rds,case.csv}`, `results/numbers/s8_pooled_vs_cells.csv` |
| `s8_02_table.R` | Table S10: pooled cohorts vs cohort × noncore cells for ETWFE and RW (and the intermediate), placebo and case study; checks its cohort × noncore rows against Table 2 | `results/tables/tableS10_cohort_noncore.*` |
| `s8_03_equivalence.R` | ETWFE vs TWFE without covariates: year effects, cell sums over fitted and imputed years, noncore evaluated predictions, with and without the population offset; then all 600 placebo datasets (numbers quoted in Section 2.4) | `results/numbers/etwfe_twfe_equivalence{,_placebo}.csv` |

Run time: `s8_01` about 30 minutes; `s8_02` and `s8_03` a few minutes.

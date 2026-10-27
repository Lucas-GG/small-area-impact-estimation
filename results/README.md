# Manuscript and supplement exhibits

`figures/` and `tables/` hold the exhibits of the paper and its supplement, named after their number in the manuscript. `numbers/` holds the CSVs behind numbers quoted in the text that have no exhibit of their own. Each file is written by one notebook or script in `analysis/`. Tables are written as CSV (the numbers) and, where built with `R/tables.r`, as HTML, LaTeX and Word; the Word layout in the manuscript was finished by hand, so only the numbers are meant to reproduce exactly.

| Exhibit | File | Written by |
|---|---|---|
| Table 2 | `tables/table2_model_comparison.*` | `analysis/03b_placebo_results.qmd` |
| Figure 1 | `figures/fig1_placebo_horizon_noncore.png` | `03b` |
| Table 3 | `tables/table3_case_study.*` | `analysis/05_case_study.qmd` |
| Figure 2 | `figures/fig2_case_effects.png` | `05` |
| Figure D1 | `figures/figD1_heterogeneity_power.png` | `analysis/04_injected_effects.qmd` |
| Table S1 | `tables/tableS1_performance_noncore.*` | `03b` |
| Table S2 | `tables/tableS2_performance_core.*` | `03b` |
| Table S3 | `tables/tableS3_county_effects.*` | `03b` |
| Figure S1 | `figures/figS1_placebo_errors_noncore.png` | `03b` |
| Figure S2 | `figures/figS2_placebo_errors_core.png` | `03b` |
| Figure S3 | `figures/figS3_placebo_horizon_balanced_noncore.png` | `03b` |
| Figure S4 | `figures/figS4_placebo_horizon_core.png` | `03b` |
| Figure S5 | `figures/figS5_county_effects_noncore.png` | `03b` |
| Table S4 | `tables/tableS4_injection_mde.csv` | `04` |
| Figure S6 | `figures/figS6_injection_power.png` | `04` |
| Figure S7 | `figures/figS7_case_effects_sensitivity.png` | `05` |
| Table S5 | `tables/tableS5_case_effects.*` | `05` |
| Table S6 | `tables/tableS6_case_heterogeneity.*` | `05` |
| Table S7 | `tables/tableS7_case_effect_modification.*` | `05` |
| Figure S8 | `figures/figS8_covariate_variants.png` | `03b` |
| Table S8 | `tables/tableS8_screening.*` | `analysis/S6_screening/s3_04_exhibits.qmd` |
| Table S9 | `tables/tableS9_prior_sensitivity.*` | `analysis/S7_priors/s4_05_prior_sensitivity.qmd` |
| Table S10 | `tables/tableS10_cohort_noncore.*` | `analysis/S8_pooled_cohorts/s8_02_table.R` |

Numbers quoted in the text:

| File | Used in | Written by |
|---|---|---|
| `numbers/covariate_change.csv` | Section 3.3.1: change in RMSE and interval score from adding covariates | `03b` |
| `numbers/covariate_variants.csv` | Supplement S5: TWFE and ETWFE with common vs varying covariate slopes (data of Figure S8) | `03b` |
| `numbers/appendixD_heterogeneity_size_power.csv`, `numbers/appendixD_heterogeneity_null_diagnostics.csv` | Appendix D.3: size and power of the heterogeneity test; placebo null diagnostics | `04` |
| `numbers/etwfe_twfe_equivalence.csv`, `numbers/etwfe_twfe_equivalence_placebo.csv` | Section 2.4: agreement of ETWFE and TWFE without covariates | `analysis/S8_pooled_cohorts/s8_03_equivalence.R` |

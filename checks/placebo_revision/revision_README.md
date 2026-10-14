# Whole-sample placebo revision

Existing scripts and their results are preserved. New outputs are under
`_output/placebo_revision/revision/`; caches are under `.quarto/`.

Run from the project root:

```r
source('R/placebo_revision_gap_audit.R')
gate <- run_revision_gap_audit()
stopifnot(gate$proceed)
source('R/placebo_revision.R')
run_placebo_revision(specification='no_covariates')
run_placebo_revision(specification='with_covariates')
source('R/placebo_revision_verify.R')
verify_placebo_revision(specification='no_covariates')
verify_placebo_revision(specification='with_covariates')
source('R/placebo_revision_presentation.R')
revision_export_paper_outputs()
```

Then render `analysis/3.2_revision_placebo_performance.qmd`. Its specification
parameter accepts `both`, `no_covariates`, and `with_covariates`.

The exact cumulative-year audit found zero gaps at every observed horizon in
both inputs. Every missing five-year window is trailing truncation: a county
observed at s has all years 0 through s. The legacy horizon-specific primary
samples therefore remain unchanged. The gate conservatively stops if any cell
has more than 0.1% shortened cumulative totals or any ordering/duplicate issue.
The synthetic audit check distinguishes trailing truncation, missing starts,
and missing middle years.

S1 is full. S2 is omitted because all prepared county-years already have six
finite model predictions, including TWFE. S3 preserves the original pooled
median pre-RMSE and strict 90th-percentile threshold. S4 uses the fixed own-model
strict |z| > 2 screen, retaining missing diagnostics. S5 intersects all six
S4 retained sets within S1. Each specification uses its own diagnostics.

B1/B3/B4/B5 are separately labelled balanced counterparts: county/simulations
with all years 0–4 are used at every horizon. S3/S4/S5 flags and percentile bins
are frozen from the full legacy reference, rather than recalibrated. At year 4
all primary/balanced counterpart results must coincide; counts in balanced
arms must stay constant across horizons. These sensitivities do not replace S1.

Counts and excluded shares are reported by horizon. Missing diagnostic shares
are shown separately for every model, including TWFE. Paired gain MC SEs use
simulation-level ratio influence functions. Upper-tail RMSE gives equal weight
to each simulation with a nonempty bin; its exported counts show missing bins.
Top 5%/1% bins overlap the highest decile.

Pre-period Poisson summaries use exact saved main-fit county-year values; no
refits occur. S1 uses the full prepared county population; B1 restricts those
totals to balanced counties. Dispersion and zero shares report their valid
denominators. The raw Pearson sums and counts reproduce the physically separated
held-out input summaries within floating-point precision.

The null IQR executes the existing `iqr_all` expression from
`analysis/4. the case.rmd`: posterior **mean** of draw-wise county IQRs. Null
aggregate effects use observed minus Poisson-noised cumulative predictions;
posterior summaries are formed within each simulation, then their quantiles
across simulations are exported. Predictive draws reuse the seed/column/row
order of the legacy placebo summary.

Verification includes Phase 0 S1/S3 reproduction, A/B assignment and physical
half guards, shared-screen agreement on overlapping counties, S5 subset and
identical observed aggregates, balanced constant counts and year-4 equality.
Every substantive input path is logged. The output index includes intended
paper placements and descriptive filenames for all exported tables/figures.

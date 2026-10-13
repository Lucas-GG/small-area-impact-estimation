# Half-A aggregate-estimand screening

Run `Rscript R/heldout_aggregate.R` after the existing county-level compute step, then render the two current presentation documents. The optional Quarto computation wrapper is `analysis/3.1a.1_heldout_compute.qmd`. New exports use only `aggregate_*.csv` names in `_output/analysis/output/3.1a/`; predictive caches are A-only and live in `.quarto/heldout-aggregate/`.

The analysis preserves original CSVs and the decision file byte-for-byte. It checks A county IDs against the fixed assignment, rejects any non-A row, verifies the assignment hash, and immediately restricts deserialized main fits to A before any calculation. Half-B evaluation inputs are never opened. Noncore comes from `urb == "6.NonCore"` in the original county panel. Screens are computed in the original full-A pool and then restricted to each subgroup, preserving top-percent denominators.

Acceptance checks:

- `Rscript checks/tests/check_aggregate.R`: signed-error cancellation, proportional paired RMSE-ratio uncertainty, missing oracle pairs, exact matched-window manuscript reference, B-row/hash rejection, and seeded predictive reconstruction.
- All 18 saved model/scenario draw groups reproduce their existing valid county predictive intervals exactly.
- All retained counts match the existing replication ledger and retention summary.
- Random-band quantiles use simulation-averaged squared errors and then RMSE ratios, with 100 seeded draw indices; unscreened random ratios equal 1.
- Existing CSV and decision checksums are verified before/after computation and rendering.

The available-window manuscript-formula check agrees within 5.11e-15 per replication. **The complete-five-year manuscript reference fails the paired within-MC-error equality check in all 36 subgroup/model/scenario cells.** The new available-window estimand includes partial histories and different eligible county denominators, so it cannot numerically reproduce the complete-horizon RMSE. Both comparisons and their paired MC SEs are reported; no target window was changed to force agreement.

The leave-one-simulation-out variance oracle has 168–195 valid paired replications out of 200, depending on subgroup/scenario. When a current eligible county has fewer than two other usable simulations, its variance and the oracle replication are unavailable. The independent, unbiased variance criterion does not guarantee the lowest actual RMSE in the presence of bias/dependence, and some oracle ratios exceed 1.

At the explored minima with missing diagnostics kept, all 18 model/scenario comparisons in each subgroup have a pointwise 1.96-MC-SE upper bound below 1. These are exploratory, without cutoff-selection or multiple-comparison adjustment; no decision is recorded. Both subgroups, both missing policies, signed and standardized metrics, exact minima, prediction-interval availability, size-population ranges, and resolution checks are in the current main/details pages and CSVs. The three mid-p grid cutoffs 0.001, 0.0025, 0.005 have identical retained sets/errors. Population ranges describe mean pre-adoption outcome-age population within strata defined by expected counts; stratum 5 is largest by expected count.

The former standalone wrapper was consolidated into the single compute document. Its source is preserved in `.quarto/source-history/3.1a-compute/`; the R computation module and outputs are unchanged.

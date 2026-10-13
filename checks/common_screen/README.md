# Common-screen aggregate sensitivity

Run `Rscript R/heldout_common_screen.R` after the original aggregate analysis and TWFE-eligible common-population sensitivity (10.1), then render the exploration and details reports. The optional compute wrapper is `analysis/3.1a.1_heldout_compute.qmd`.

New CSVs live under `_output/analysis/output/3.1a/common_screen/` and use new `common_screen_` filenames. No previous CSV, computation, figure, table or decision-template value is modified. All intermediate/render artifacts follow the project's `.quarto`/`_output` conventions.

For each simulation/scenario, noncore E county IDs are aligned across the six saved A-only model caches. The existing `keep_matrix_A` evaluates fixed diagnostic thresholds for z=2 and predictive mid-p=.05/.10, for both missing policies. Retained masks use the same model-aligned `Reduce(&, screens)` intersection as section 8. No target-specific extra exclusion is needed because section 10.1 already requires all six targets to be finite, on the same available window. Aggregate formulas and matching predictive draw indices reuse the existing aggregate helpers.

The comparison table preserves unscreened and own-screen E results from 10.1 and adds common-screen results. RMSE ranks are ascending and coverage ranks descending, with shared ranks for ties. Relative RMSE uses mean retained observed aggregate and a covariance-aware MC SE; RMSE ratios compare to each model's unscreened E baseline with paired influence-function MC SEs. Sole-flag contributions count counties flagged by exactly one model; overlaps are counted once in the shared contribution. Pooled shares divide total flag counts by total common exclusions. Mean per-simulation shares and MC SEs are also exported.

`Rscript checks/tests/check_common_screen.R` independently checks exported membership IDs against the per-model subset ledger, identical counts and observed aggregates, own-screen reproduction against 10.1, exclusive/shared accounting, every unscreened numerical metric/MC SE, A-only/hash guards and prior-file checksums. Also rerun `checks/tests/check_aggregate.R` and `checks/tests/check_twfe_common.R`. Original output hashes are snapshotted in `.quarto/common-screen/original-md5.rds`; render/source snapshots are in that same folder.

The latent-level reversal persists on the identical common retained sample: TWFE goes from first unscreened to sixth for all three rules. Rank changes versus own-screen depend on scenario/rule and are exported, rather than treated as a general dominance claim. No screening rule is selected.

The former standalone wrapper was consolidated into the single compute document. Its source is preserved in `.quarto/source-history/3.1a-compute/`; the R computation module and outputs are unchanged.

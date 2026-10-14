# Post-hoc held-out county influence

This is a descriptive appendix. It does not select a rule, refit any model, reconstruct post-period targets, change the decision or overwrite any prior evaluation output. It reads the physically separate A/B snapshots and validated saved cumulative prediction caches for the two branches. The cached `y_check - mu_check` is the existing available-window signed county error; raw main-fit contents are prohibited. Cache provenance uses assignment MD5 and saved fit/input metadata. The shared `derive_heldout` gained `diagnostics_only=TRUE`, which returns the same original diagnostic scores before deriving any old targets/averages; its default behavior is unchanged. Original masks reuse `keep_matrix_A`, and all original half/hash guards and `aggregate_summary` formulas remain in use.

Run from the repository root:

```
Rscript -e 'source("R/heldout_influence.R"); run_heldout_influence(normalizePath("."))'
quarto render analysis/3.1d_heldout_influence_appendix.qmd --output heldout_screening_influence_appendix.html
mv _output/analysis/heldout_screening_influence_appendix.html _output/heldout/influence/heldout_screening_influence_appendix.html
Rscript checks/tests/check_heldout_influence.R
```

New `influence_*.csv` files and the self-contained short report live in `_output/heldout/influence/`. Initial temporary render sources are under `.quarto/heldout-influence-render/`. All prior statistical exports, A reports, A/B estimate outputs, the decision and assignment are checksum-protected. The input path log distinguishes cached reads from fit metadata checks.

Each rule's county contribution is the mean of E × e_j/n on the same finite aggregate simulations, with zero for unretained/absent counties. Contributions are signed; negative values describe cancellation, so top positive shares can exceed 100%. Concentration ranks contributions, while top-five profiles and top-k deletion rank the decrease in unscreened RMSE when deleting a county (county ID breaks ties). Ratio influence ranks absolute changes, with signed changes also exported. Every original screen, including the full-half percentile reference, is frozen during deletion. A county is removed from both arms where it appears, without replacing another flagged county or recomputing a percentile threshold. Empty arms are unavailable and finite/paired counts are exported.

Profiles retain the original size stratum and its main-TWFE-MAR pre-period expected-count reference. Population is saved outcome-age exposure n averaged across all finite panel years. Pseudo-treated frequency uses all pseudo-treated appearances; errors and conditional flag rates use eligible target appearances. All-simulation flag rates are exported too. Recurrence is within the same branch/half/scenario across models; A and B have disjoint counties.

Flag capture weights top-five county/simulation appearances by the positive part of their unscreened contribution. It describes coincidence of flags and influential realizations, not a causal allocation of the gain. Summary gap statistics average the five non-TWFE models across scenarios; full tables include all six. Both signed and absolute A/B gaps and the original full-sample ratios are exported. The existing estimand changes with deletion, including its denominator, so influence differs from contribution.

Checks independently re-sum all 288 MSE decompositions, check 432 reproduction comparisons, verify A/B IDs and absence of raw-fit reads, test signed cancellation and absent-county deletion, preserve all original checksums, and verify k=0 deletion reproduces the full estimate. The report ends with seven descriptive lines and proposes no rule.

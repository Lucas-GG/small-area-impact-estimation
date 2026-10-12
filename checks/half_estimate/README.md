# Fixed-rule half-B estimate: A-only dry run

The maintained report is `analysis/3.1b heldout confirmation.qmd` (the existing filename is retained; its title and framing are now half-B estimate). `params.half` defaults to B. `R/heldout_estimate.R` evaluates only the recorded strict |z| > 2 rule and the two fixed sensitivities: mid-p < .10 and pre-RMSE top 10%. The YAML requires `estimation_plan` with four nonempty entries; no mechanical decision table is evaluated.

The Phase-1 render uses an isolated project in `_output/heldout/confirmation/render_A/`, copying the maintained report to `half_A_dry_run.qmd`, with output-dir `..`. This avoids changing the main project configuration or rendering A reports. Reproduce with:

```
quarto render _output/heldout/confirmation/render_A/half_A_dry_run.qmd -P half:A
Rscript checks/tests/check_half_estimate.R
```

Refresh the copied qmd from the maintained source whenever it changes. New outputs have `half_A_` prefixes inside `_output/heldout/confirmation/`. They include the rendered dry run, decision copy/MD5, commit and assignment provenance, every analysis input path and access type, per-simulation signed-error ledger, aggregate summaries, county-level associations/retention, common membership, flag accounting, ranks, numerical reproduction and model-averaged estimates. Original A CSV/report hashes are checked within the engine. File metadata checks of complete-panel fit paths verify cached provenance; their contents are never opened during the A dry run.

Shared functions gained optional `expected_half` arguments, defaulting to A: aggregate guard and predictive reconstruction, available-window derivation, screen masks/retention, and aggregate predictive-interval grouping. The estimator calls these same functions for either half. No separate implementation of aggregate RMSE/relative RMSE/paired-ratio formulas, percentile screening or predictive intervals is introduced. MC SEs for the model-averaged summaries average influence functions within simulation indices, preserving model/scenario dependence. The summary relative ratio includes the covariance of both screened and unscreened observed means. A full-precision reference is computed from the original A per-simulation exports; the fixed YAML's rounded A reference values are displayed separately.

No half-B file has been opened. The default B path stops before loading B unless the explicit user instruction `run B` has been received and the session's authorization marker `_output/heldout/confirmation/run_B_authorized` exists. Do not create that marker before the user instruction. Before B is loaded, the decision must be committed with identical bytes, and its MD5 and the assignment hash must match the completed A dry run. B then loads its physically separate half-only snapshot and reconstructs the available-window targets from saved draws restricted to B, using the shared A formulas and seeds. An incomplete or incompatible saved input stops rather than inventing a target or diagnostic. B outputs use new `half_B_` filenames; the final descriptive summary is ten lines.

The historical `checks/tests/check_common_screen.R` preservation snapshot includes the old blank decision MD5. Its only mismatch after the user replaced the decision is that intentional YAML change; its original CSVs remain unchanged. The updated Phase-1 engine separately verifies all original A CSVs and reports, with the revised decision MD5 pinned independently.

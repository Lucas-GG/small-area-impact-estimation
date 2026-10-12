# Half-A report presentation revision

The visual summary is `analysis/3.1a.2_heldout_exploration.qmd`; full report tables and methods are in `analysis/3.1a.3_heldout_details.qmd`. Both render into `_output/analysis/` and link to each other. All 22 original CSVs are copied byte-for-byte (hard links where possible) to `_output/analysis/output/3.1a/`, with an equal-weight scenario candidate scorecard and three separately labelled corrected-band CSVs. CSV mirrors beside the central HTML preserve local download links. The very large replication-level ledgers are downloadable in full; report-level tables are displayed in full using DT, 25 rows per page, column filters and filtered CSV buttons. No full table is printed statically on the details page except the one-row band-width example.

## Computation scope

All actual screening statistics, targets, diagnostic associations, coverage, ranks, MC SEs, county assignment, seeds and fitted models are unchanged. The pages read the existing saved summaries, assert physical half-A-only input and verify its assignment hash. Half B is never opened and no rule is selected. The six-candidate object is defined once in the shared presentation helper; defaults are err_z, missing keep, k-stratified associations. The original spaced-name exploration remains the computation/refresh document. Presentation fails if its saved count summaries predate the A input.

The one numerical correction is the random band: regenerate the same 100 seeded permutations per simulation, verify every per-simulation mean against its saved value, average each random draw index across simulations, then take quantiles across those averages. All 604,800 regenerated per-rule/simulation means (raw and standardized cumulative targets) match to 1e-10 with matching finite status. The corrected per-rule bands use the actual exclusions per simulation; the single grey plotting band uses an independent excluded-fraction grid rounded to county counts. Old random-bound columns remain unchanged as an audit trail. For RW / MAR / err_z / mid-p < 0.05 / missing keep, width is 0.1155267 before and 0.008332636 after (11.55 and 0.83 percentage points). These are random-reference bounds, not a confidence interval for the observed benefit.

## Changed source files

- `analysis/3.1a.2_heldout_exploration.qmd`: new concise main report, inline computed takeaways, gt scorecards, requested visual sections and links.
- `analysis/3.1a.3_heldout_details.qmd`: new paginated full tables, sensitivities, methods, correction audit and CSV catalogue.
- `R/heldout_presentation.R`: saved-result loading, shared shortlist/order/colours, formatting, presentation plots, DT tables and CSV mirrors.
- `R/heldout_random_band.R`: isolated random-band correction and cache; original random means are asserted equal.
- `checks/tests/check_presentation.R`: averaging-before-quantiles counterexample and mean-preservation checks.
- `_quarto.yml`: narrow the legacy exploration exclusion so the new pages belong to the project and render centrally.
- `.gitignore`: ignore generated `analysis/output/` CSVs and cache.
- `analysis/heldout README.md`: entry points and presentation workflow.
- `checks/heldout_revision/README.md`: correct the older explicit-render advice.
- This revision note and generated `acceptance.json`.

`DT` was installed into the project-local `_output/R-library`, without changing the shared R library. Initial rendering attempts exposed a project-exclusion flag issue, chunk-option syntax, and a ribbon mapping inherited from the main curve layer; all were corrected. Final successful render logs are `_output/heldout/presentation-main-render.log` and `presentation-details-render.log`. The main page is below 2 MB and has at most 18 rows per displayed table. Exact size, link and inline-value checks are recorded after final validation.

## Final acceptance results

Both final sources rendered successfully (exit code 0). Main HTML: 1,978,111 bytes; largest table: 18 rows and five columns. Details: 20 full report-level DT tables, all with 25-row pages, filters and CSV buttons. All 76 cross-page/anchor/CSV links resolve. Both copies of every candidate scorecard cell match the saved scorecard CSV to two decimals; summary old/new row counts, added rows and increase match the count CSV. The 22 original CSV exports remain byte-identical. There are 26 downloadable CSVs. Figure inspection confirms consistent diagnostic colours, fixed model order, labelled candidate points and one grey random band. Focused averaging-order and mean-preservation tests pass. Structured results are in `acceptance.json`.

The legacy available-window revision note was also corrected to explain that Quarto's output-directory flag applies to project documents, rather than excluded standalone documents. No original simulation/heldout inputs, decisions, or statistical helpers were changed.

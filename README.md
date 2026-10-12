# small_area_impact

Code for estimating small-area impacts and comparing placebo prediction methods for county youth suicide counts. Private outcome data remain under `data/notpublic/` and are not distributed.

## Project layout

- `analysis/`: numbered notebooks; `with_covariates/` is the covariate branch, `archive/` contains historical analyses.
- `R/`: all shared model, feature and analysis helpers, directly in this folder.
- `data/simulation/`: saved assignments and original main simulation fits.
- `data/fits/`: separate backdated fits and replication/county summaries.
- `data/derived/heldout/`: fixed county assignment and prepared A/B inputs, with a separate covariate branch.
- `results/figures/` and `results/tables/`: paper-ready figures and tables, including existing saved results.
- `docs/`: manuscripts and supporting documents.
- `checks/`: verification scripts and audit report sources; generated audit files go to `_output/checks/`.
- `_output/`: rendered reports, diagnostic exports and preserved legacy generated artifacts.
- `.quarto/`: Quarto project state/cache; the former nested cache is preserved inside it as `legacy_analysis_cache/`.

The manually chosen held-out rule belongs in `analysis/3.1a.4_heldout_decision.yml`, alphabetically after exploration and before confirmation. No decision file is created automatically. Both confirmation branches and 3c.2 use this same file.

## Running notebooks

Run commands from the project root. The root `_quarto.yml` sets project-root execution and sends rendered HTML to `_output/`. The backdated notebooks default to `run_fits:false`, displaying existing saved summaries or running preflight. To start new backdated fits, explicitly pass `-P run_fits:true`.

```sh
quarto render "analysis/2.1.refit placebo for heldout error.qmd" -P run_fits:false
quarto render "analysis/with_covariates/2c.1. refit placebo for held out error.qmd" -P run_fits:false
Rscript R/heldout_prepare.R
```

Preparation reuses the existing split. Covariate preparation is available through `prepare_heldout(root, with_covariates=TRUE)` in `R/heldout_prepare.R`. This reads saved summaries and never fits a model.

Exploration and confirmation documents must be rendered explicitly; they are excluded from a whole-project render. Render confirmation only after the decision is fixed. Audit documents are also explicit opt-ins. See `analysis/README.md` and `analysis/with_covariates/README.md` for details. Do not rerun the original main fitting notebooks just to produce reports.

Saved historical RDS files may contain provenance paths using the former folder layout/name. Their contents were preserved during cleanup; those provenance strings are not used to locate current inputs. The directory migration log is `_output/checks/directory_cleanup/moves.json`.

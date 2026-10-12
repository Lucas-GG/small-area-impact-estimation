# Directory cleanup

The project is organized with numbered notebooks in `analysis/`, shared code directly in `R/`, saved main simulations in `data/simulation/`, backdated summaries in `data/fits/`, prepared evaluation inputs in `data/derived/`, and paper figures/tables in `results/`. The fixed rule is read from `analysis/3.1a.4_heldout_decision.yml`; the `.1` suffix places it after exploration and before confirmation. It is not created automatically.

`moves.json` records moves performed before the final project-root rename. Legacy generated files were preserved in `_output/legacy/`, the nested cache was preserved under root `.quarto/legacy_analysis_cache/`, and old root archive sources were preserved under `analysis/archive/legacy_root/`.

Validation: all shared R code and notebook chunks parsed; synthetic helper, covariate and diagnostic checks passed; the moved assignment's hash still matched the saved half-A snapshot. Both backdated notebooks rendered successfully with `run_fits:false`, displaying saved results. No model, simulation, exploration or half-B confirmation analysis was rerun. INLA emitted its existing macOS `/bin/kstat` notices; no render errors occurred.

Original main fits, backdated fit summaries and the fixed assignment were moved without editing their serialized contents. Old absolute provenance strings within RDS files remain historical. No GitHub repository rename or Git commit was performed.

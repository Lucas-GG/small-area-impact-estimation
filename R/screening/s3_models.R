# Screening supplement (S6) models and their main placebo fits (02_fit_placebo.qmd).
s3_models <- c("TWFE", "ETWFE", "GB-Ens", "RW-IID")
s3_codes <- c(TWFE = "twfe", ETWFE = "etwfe", `GB-Ens` = "xgb", `RW-IID` = "rwi")

# Saved main placebo fit of one model and scenario (1-3).
s3_main_path <- function(root, model, scenario, with_covariates = FALSE) {
  code <- s3_codes[[model]]
  file.path(root, "data/simulation/placebo", if (with_covariates) "wcov" else "wocov", sprintf("%s_%s.rds", code, scenario))
}

# Fixed per-model seed index from the original six-model order (TWFE, ETWFE, RW, RW-IID,
# RW-ST, GB-Ens), so predictive-draw seeds do not depend on which models S3 includes.
s3_seed_index <- c(TWFE = 1L, ETWFE = 2L, RW = 3L, `RW-IID` = 4L, `RW-ST` = 5L, `GB-Ens` = 6L)

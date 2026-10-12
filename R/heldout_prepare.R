# Preparation only: split before analysis; never estimate a model.
prepare_heldout_legacy <- function(root, split_seed = 2032026L, with_covariates = FALSE) {
  library(data.table)
  shared <- file.path(root, 'data/derived/heldout')
  out <- if(with_covariates) file.path(shared,'with_covariates') else shared
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  assignment_file <- file.path(shared,'county_assignment.rds')
  if(with_covariates && !file.exists(assignment_file))
    stop('Create the no-covariate county assignment first; the covariate branch reuses it.')
  if (!file.exists(assignment_file)) {
    # A single reference model/scenario gives one county-size measure shared by all analyses.
    # Pool county-specific pre-period expected counts across saved MAIN TWFE MAR fits.
    p <- file.path(root,'data/simulation/wocov/new_imputation_twfe_1.rds')
    raw <- readRDS(p)
    z <- as.data.table(raw)[,.(i,r,year,start_year_sim,y0_hat)]
    rm(raw); gc()
    per_fit <- z[year < start_year_sim & is.finite(y0_hat) & y0_hat >= 0,
      .(pre_mean=mean(y0_hat)),by=.(id=i,sim=r)]
    sizes <- per_fit[,.(mean_pre_expected_count=mean(pre_mean),n_reference_fits=.N),by=id]
    ids <- as.integer(dimnames(readRDS(file.path(root,'data/simulation/simulated_assignments_with_meta.rds'))$assignments)[[1]])
    sizes <- merge(data.table(id=ids),sizes,by='id',all.x=TRUE)
    setorder(sizes,mean_pre_expected_count,id,na.last=TRUE)
    sizes[, stratum := 6L]
    sizes[is.finite(mean_pre_expected_count), stratum := pmin(5L,1L+floor((seq_len(.N)-1L)*5L/.N))]
    # Deterministic ID tie-breaking for quintiles; odd strata differ by at most one county.
    set.seed(split_seed)
    sizes[, half := sample(rep(c('A','B'),length.out=.N)),by=stratum]
    sizes[, `:=`(split_seed=split_seed,reference='Main TWFE MAR: mean of per-fit pre-period expected counts')]
    setorder(sizes,id); saveRDS(sizes,assignment_file); fwrite(sizes,file.path(shared,'county_assignment.csv'))
    rm(z,per_fit); gc()
  }
  assignment <- readRDS(assignment_file)
  stopifnot(all(assignment$split_seed==split_seed),!anyDuplicated(assignment$id))
  fit_dir <- if(with_covariates) 'data/fits/backdate_wcov' else 'data/fits/backdate'
  paths <- list.files(file.path(root,fit_dir),'^summary_.*\\.rds$',full.names=TRUE)
  if (!length(paths)) stop('No backdated checkpoints available')
  # Prefer group checkpoints, including a partial run; never alter the original summaries.
  x <- rbindlist(lapply(paths,readRDS),fill=TRUE)
  stopifnot(!anyDuplicated(x[,.(sim,model,scenario,id)]))
  if (!'tail_convention' %in% names(x)) x[,tail_convention:=NA_character_]
  if (!'ppp_ho_upper' %in% names(x)) x[,ppp_ho_upper:=NA_real_]
  if (!'ppp_ho_lower' %in% names(x)) x[,ppp_ho_lower:=NA_real_]
  legacy <- is.na(x$tail_convention) | !grepl('^mid-p', x$tail_convention)
  both <- is.finite(x$ppp_ho_upper) & is.finite(x$ppp_ho_lower)
  recover <- legacy & both
  equality <- x$ppp_ho_upper[recover]+x$ppp_ho_lower[recover]-1
  stopifnot(all(equality >= -1e-12 & equality <= 1+1e-12))
  x[recover, `:=`(ppp_ho_upper=ppp_ho_upper-.5*pmax(0,equality),
                  ppp_ho_lower=ppp_ho_lower-.5*pmax(0,equality),tail_convention='mid-p recovered from inclusive pair')]
  x[legacy & !both, `:=`(ppp_ho_upper=NA_real_,ppp_ho_lower=NA_real_,tail_convention='unavailable: legacy pair absent')]
  # Remove stale inclusive aliases/two-sided values; analysis derives two-sided mid-p afresh.
  x[, intersect(c('ppp_ho','ppp_ho_two_sided'),names(x)) := NULL]
  x <- merge(x,assignment,by='id',all.x=TRUE)
  stopifnot(!anyNA(x$half))
  manifest <- data.table(file=basename(paths),bytes=file.info(paths)$size,modified=as.character(file.info(paths)$mtime))
  fwrite(manifest,file.path(out,'input_manifest.csv'))
  for (h in c('A','B')) {
    rows <- x[half==h]
    saveRDS(list(rows=rows,manifest=manifest,prepared=Sys.time(),expected_groups=18L,
      assignment_md5=unname(tools::md5sum(assignment_file))),file.path(out,paste0('half_',h,'.rds')))
  }
  invisible(manifest)
}
# The no-covariate default updates only the revised A input; existing B files are untouched.
prepare_heldout <- function(root, split_seed=2032026L, with_covariates=FALSE) {
  if (with_covariates) return(prepare_heldout_legacy(root, split_seed, with_covariates=TRUE))
  source(file.path(root,'R/heldout_exploration.R'))
  prepare_exploration_A(root,split_seed)
}
if (sys.nframe()==0L) {
  source('R/backdate_helpers.R')
  started<-proc.time()[['elapsed']]
  prepare_heldout(find_backdate_root())
  cat('Preparation elapsed seconds:',proc.time()[['elapsed']]-started,'\n')
}

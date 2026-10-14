# Phase 0 only. Existing summaries and outputs are never written by this script.
library(data.table)
run_placebo_revision_phase0 <- function(root=getwd(), specification=c('no_covariates','with_covariates')) {
  specification <- match.arg(specification)
  out <- file.path(root,'_output/placebo_revision/phase0',specification)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  log <- data.table(path=character(),read_at=character(),purpose=character())
  record <- function(path,purpose) {
    log <<- rbind(log,data.table(path=normalizePath(path),read_at=as.character(Sys.time()),purpose=purpose))
    fwrite(log,file.path(out,'revision_input_paths.csv'))
  }
  read <- function(path,purpose) {record(path,purpose); readRDS(path)}
  source_path <- file.path(root,'R/heldout_analysis.R')
  record(source_path,'shared diagnostic derivation'); source(source_path)
  cov <- specification=='with_covariates'
  main <- read(file.path(root,'data/simulation',if(cov)'sim_fits_wcov.rds' else 'sim_fits.rds'),'existing main-simulation input, whole sample')
  main <- as.data.table(main)
  # Release predictive draw columns: the join audit uses original county keys only.
  keys <- unique(main[,.(model=as.character(model),scenario=as.character(dgp),sim=r,id=i,subgroup=as.character(nc))])
  rm(main); gc()
  model_map <- c(twfe='TWFE',etwfe='ETWFE',rw='RW',rwi='RW-IID',rwst='RW-ST',xgb='GB-Ens')
  scenario_map <- c('1'='MAR','2'='Latent level','3'='Latent trend')
  translate <- function(x,map) {hit<-x %in% names(map); x[hit]<-unname(map[x[hit]]); x}
  keys[,model:=translate(model,model_map)]
  keys[,scenario:=translate(scenario,scenario_map)]
  stopifnot(all(keys$model %in% unname(model_map)),all(keys$scenario %in% unname(scenario_map)))
  fwrite(rbindlist(list(data.table(key_type='model',input=names(model_map),output=unname(model_map)),
    data.table(key_type='scenario',input=names(scenario_map),output=unname(scenario_map)),
    data.table(key_type='column',input=c('r','i','dgp'),output=c('sim','id','scenario')))),file.path(out,'revision_key_translation.csv'))
  assignment_path <- file.path(root,'data/derived/heldout/county_assignment.rds')
  assignment <- read(assignment_path,'fixed assignment and hash validation')
  input_dir <- file.path(root,'data/derived/heldout',if(cov)'with_covariates' else '')
  diagnostics <- rbindlist(lapply(c('A','B'),function(h) {
    obj <- read(file.path(input_dir,paste0('half_',h,'.rds')),paste('own-specification diagnostic half',h))
    stopifnot(identical(obj$assignment_md5,unname(tools::md5sum(assignment_path))),
      all(obj$rows$half==h),all(obj$rows$id %in% assignment[half==h,id]))
    derive_heldout(obj$rows,diagnostics_only=TRUE)
  }),fill=TRUE)
  stopifnot(!anyDuplicated(diagnostics[,.(model,scenario,sim,id)]))
  d <- diagnostics[,.(model,scenario,sim,id,ho_standardized,ppp_ho_upper,ppp_ho_lower,p_two,
    n_pre,n_ho,heldout_status,fit_error,diagnostic_present=TRUE)]
  joined <- merge(keys,d,by=c('model','scenario','sim','id'),all.x=TRUE)
  joined[,missing_reason:=fcase(is.na(diagnostic_present),'unmatched county key',
    is.finite(ho_standardized)&is.finite(p_two),'available',
    !is.na(n_pre)&n_pre<2,'fewer than two backdated pre-period years',
    !is.na(heldout_status)&heldout_status!='ok',paste0('saved heldout_status: ',heldout_status),
    default='unexplained unavailable diagnostic')]
  joined[,specification:=specification]
  fwrite(joined,file.path(out,'revision_diagnostic_join.csv'))
  audit <- joined[,.(counties=.N,missing_z=sum(!is.finite(ho_standardized)),
    missing_mid_p=sum(!is.finite(p_two)),missing_z_share=mean(!is.finite(ho_standardized)),
    missing_mid_p_share=mean(!is.finite(p_two))),by=.(specification,model,scenario,subgroup)]
  fwrite(audit,file.path(out,'revision_missing_diagnostics.csv'))
  fwrite(joined[,.(counties=.N),by=.(specification,model,scenario,subgroup,missing_reason)],file.path(out,'revision_missing_reasons.csv'))
  # No rounded HTML comparison can establish numerical-precision reproduction.
  writeLines(c('STOP: Phase 0 reproduction is not established.',
    'Saved HTML tables round to three decimal places; an unrounded reference has not been identified.',
    'Current 3c.2 uses half B and fixed-rule screening, not whole-sample historical top-10% screening.',
    'Identify the published covariate historical source/reference before declaring reproduction.'),file.path(out,'revision_phase0_status.txt'))
  if(any(joined$missing_reason %in% c('unmatched county key','unexplained unavailable diagnostic')))
    stop('Unexplained diagnostic missingness: inspect revision_missing_reasons.csv. Phase 0 stopped.')
  invisible(audit)
}

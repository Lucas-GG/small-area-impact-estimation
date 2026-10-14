library(data.table)
verify_placebo_revision <- function(root=getwd(),specification=c('no_covariates','with_covariates')) {
  specification<-match.arg(specification);source(file.path(root,'R/placebo_revision.R'));ctx<-revision_context(root)
  read<-function(name)ctx$csv(file.path(ctx$out,paste0('revision_',specification,'_',name,'.csv')),'revision verification')
  flags<-read('county_flags');ledger<-read('scores_by_simulation');retention<-read('retention_by_simulation')
  stopifnot(all(!flags$S5_common_z2|flags$S4_own_z2),all(!flags$B5_balanced_common_z2|flags$B4_balanced_own_z2))
  balanced_counts<-retention[grepl('^B',sample),.(horizons=uniqueN(s),base_counts=uniqueN(n_base),retained_counts=uniqueN(n_retained)),by=.(sample,dgp,nc,model,r)]
  stopifnot(all(balanced_counts$horizons==5L),all(balanced_counts$base_counts==1L),all(balanced_counts$retained_counts==1L))
  common<-ledger[sample%in%c('S5_common_z2','B5_balanced_common_z2'),.(n_models=uniqueN(model),observed_range=max(y)-min(y)),by=.(sample,nc,dgp,s,r)]
  stopifnot(all(common$n_models==6L),all(common$observed_range<1e-12))
  fwrite(common,file.path(ctx$out,paste0('revision_',specification,'_S5_observed_aggregate_check.csv')))
  map<-data.table(sample=c('S1_full','S3_historical','S4_own_z2','S5_common_z2'),balanced=c('B1_balanced_full','B3_balanced_historical','B4_balanced_own_z2','B5_balanced_common_z2'))
  a<-merge(ledger[s==4&sample%in%map$sample],map,by='sample');b<-copy(ledger[s==4&sample%in%map$balanced]);b[,balanced:=sample]
  metrics<-c('y','y0_hat','diff','sq_diff','crps','cvg','intv');keys<-c('balanced','nc','dgp','model','s','r')
  a[,(metrics):=lapply(.SD,as.numeric),.SDcols=metrics];b[,(metrics):=lapply(.SD,as.numeric),.SDcols=metrics]
  check<-merge(melt(a,id.vars=keys,measure.vars=metrics,value.name='primary'),melt(b,id.vars=keys,measure.vars=metrics,value.name='balanced_value'),by=c(keys,'variable'))
  check[,absolute_difference:=abs(primary-balanced_value)]
  stopifnot(all(check$absolute_difference<1e-12|is.na(check$absolute_difference)))
  fwrite(check,file.path(ctx$out,paste0('revision_',specification,'_balanced_year4_check.csv')))
  source_path<-file.path(root,'R/heldout_exploration.R');ctx$record(source_path,'verify shared screen on held-out overlap');source(source_path)
  source_path<-file.path(root,'R/heldout_analysis.R');ctx$record(source_path,'verify original diagnostics and assignment hash');source(source_path)
  assignment_path<-file.path(root,'data/derived/heldout/county_assignment.rds');assignment<-ctx$read(assignment_path,'assignment hash guard')
  physical<-rbindlist(lapply(c('A','B'),function(half_label) {
    path<-file.path(root,'data/derived/heldout',if(specification=='with_covariates')'with_covariates' else '',paste0('half_',half_label,'.rds'))
    obj<-ctx$read(path,'own-specification diagnostic original')
    stopifnot(identical(obj$assignment_md5,unname(tools::md5sum(assignment_path))),all(obj$rows$half==half_label),
      all(obj$rows$id%in%assignment[half==half_label,id]))
    derive_heldout(obj$rows,diagnostics_only=TRUE)[,.(dgp=scenario,model,r=sim,i=id,half,z_original=ho_standardized)]
  }))
  overlap<-merge(flags,physical,by=c('dgp','model','r','i'))
  stopifnot(nrow(overlap)==nrow(flags))
  rule<-data.table(diagnostic='ho_standardized',cutoff=2,rule_type='absolute',rule_id=1L)
  checked<-overlap[,{
    expected<-keep_matrix_A(data.table(half=half[1],id=i,ho_standardized=z_original),rule,expected_half=half[1])$keep[,1]
    list(overlap_counties=.N,revision_retained=sum(S4_own_z2),heldout_function_retained=sum(expected),flags_identical=all(expected==S4_own_z2))
  },by=.(dgp,nc,model,r,half)]
  stopifnot(all(checked$flags_identical))
  fwrite(checked,file.path(ctx$out,paste0('revision_',specification,'_S4_heldout_overlap_check.csv')))
  phase0<-read('phase0_reproduction');stopifnot(all(phase0$absolute_difference<1e-10|is.na(phase0$absolute_difference)))
  balanced_coverage<-retention[grepl('^B',sample),.(minimum_horizons=min(s),maximum_horizons=max(s),all_counts_constant=uniqueN(n_base)==1L&&uniqueN(n_retained)==1L),by=.(sample,dgp,nc,model,r)]
  fwrite(balanced_coverage,file.path(ctx$out,paste0('revision_',specification,'_balanced_counts_check.csv')))
  message(specification,': overlap, subset, common observed aggregates, balanced counts, year-4 identity, and Phase 0 checks passed')
  invisible(TRUE)
}

# Replay the published whole-sample algorithms; never render legacy scripts.
placebo_revision_reproduce <- function(root=getwd(),specification=c('no_covariates','with_covariates')) {
  specification <- match.arg(specification)
  out <- file.path(root,'_output/placebo_revision/phase0',specification)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  record <- function(p,purpose) {
    lp<-file.path(out,'revision_input_paths.csv')
    log<-if(file.exists(lp))fread(lp) else data.table(path=character(),read_at=character(),purpose=character())
    fwrite(rbind(log,data.table(path=normalizePath(p),read_at=as.character(Sys.time()),purpose=purpose)),lp)
  }
  score_path<-file.path(root,'R/scoring.r'); record(score_path,'shared original scoring functions'); source(score_path)
  p<-file.path(root,'data/simulation',if(specification=='with_covariates')'sim_fits_wcov.rds' else 'sim_fits.rds')
  record(p,'whole-sample original performance replay'); fits<-as.data.table(readRDS(p))
  reference_path<-file.path(root,'_output/placebo_revision/phase0/source_reference',
    if(specification=='with_covariates')'published_covariates_summary.rmd' else 'published_no_covariates_summary.rmd')
  record(reference_path,'historical source replay reference, git 18b3ae2')
  historical_source<-readLines(reference_path,warn=FALSE)
  stopifnot(any(grepl('median(rmse_pre)',historical_source,fixed=TRUE)),
    any(grepl('quantile(rmse_pre_agg, 0.9',historical_source,fixed=TRUE)))
  window<-fits[,.(n_years=uniqueN(s),complete_five_year=identical(sort(unique(as.integer(s))),0:4)),
    by=.(dgp,nc,model,r,i)]
  fwrite(window[,.(county_simulations=.N,incomplete=sum(!complete_five_year)),by=.(dgp,nc,model)],
    file.path(out,'revision_window_audit.csv'))
  if(any(!window$complete_five_year))stop('Existing input has incomplete five-year windows; Phase 0 stopped.')
  # Both published source versions use seed 0203, column order, and original row order.
  set.seed(0203)
  draw_cols <- grep('^y0_post_',names(fits),value=TRUE)
  fits[, (draw_cols) := lapply(.SD,function(lam) rpois(.N,lam)),.SDcols=draw_cols]
  fits[, c('n','y','y0_hat',draw_cols) := lapply(.SD,cumsum),
       by=.(dgp,nc,model,r,i),.SDcols=c('n','y','y0_hat',draw_cols)]
  fits[, rmse_pre_agg := median(rmse_pre),by=.(dgp,r,i)]
  fits[, keep := rmse_pre_agg < quantile(rmse_pre_agg,.9,na.rm=TRUE),by=.(dgp,nc)]
  # Literal reference aggregation from the historical source at 18b3ae2.
  reference <- function(x) {
    x[,c(lapply(.SD,mean),list(pcor=cor(y,y0_hat,use='complete.obs'),
      scor=cor(y,y0_hat,use='complete.obs',method='spearman'))),
      by=.(nc,dgp,model,s,r),.SDcols=c('y','y0_hat',draw_cols)]
  }
  # New parameterized implementation, with independently expressed arithmetic means.
  revised <- function(x) {
    x[,c(lapply(.SD,function(v)sum(v)/length(v)),list(pcor=cor(y,y0_hat,use='complete.obs'),
      scor=cor(y,y0_hat,use='complete.obs',method='spearman'))),
      by=.(nc,dgp,model,s,r),.SDcols=c('y','y0_hat',draw_cols)]
  }
  summarize <- function(x,historical) {
    score_avg(x); setkey(x,nc,dgp,s,model)
    x[,summ_score(.SD,sd=historical),by=.(nc,dgp,model,s)]
  }
  checks<-list()
  for(sample in c('S1_full','S3_historical')) {
    historical<-sample=='S3_historical'
    x<-if(historical)fits[keep==TRUE] else fits
    ref<-summarize(reference(x),historical)
    rev<-summarize(revised(x),historical)
    ids<-c('nc','dgp','model','s')
    metrics<-setdiff(names(ref),ids)
    a<-melt(ref,id.vars=ids,measure.vars=metrics,variable.name='metric',value.name='reference')
    b<-melt(rev,id.vars=ids,measure.vars=metrics,variable.name='metric',value.name='revised')
    c<-merge(a,b,by=c(ids,'metric'))
    c[,max_absolute_difference:=abs(revised-reference)]
    c[,reproduced:=fifelse(is.finite(reference)&is.finite(revised),
      max_absolute_difference<=1e-10*pmax(1,abs(reference)),
      (is.na(reference)&is.na(revised))|(reference==revised))]
    c[,`:=`(specification=specification,sample=sample)]
    checks[[sample]]<-c
    fwrite(rev,file.path(out,paste0('revision_',sample,'_performance.csv')))
  }
  result<-rbindlist(checks)
  setnames(result,c('nc','dgp','s'),c('subgroup','scenario','horizon'))
  fwrite(result,file.path(out,'revision_reproduction_checks.csv'))
  saved_checks<-rbindlist(lapply(c('S1_full','S3_historical'),function(sample_label) {
    legacy<-file.path(root,'results/tables',paste0('performance_table_avg',
      if(sample_label=='S1_full')'_full' else '',
      if(specification=='with_covariates')'_wcov' else '', '.html'))
    record(legacy,'rounded saved-table cross-check')
    a<-as.data.table(rvest::html_table(rvest::read_html(legacy))[[1]])
    ids<-c('nc','dgp','model','s'); metrics<-setdiff(names(a),ids)
    a[,(metrics):=lapply(.SD,function(v)as.numeric(gsub(',','',gsub('−','-',as.character(v))))),.SDcols=metrics]
    ref<-result[sample==sample_label,.(nc=subgroup,dgp=scenario,model,s=horizon,metric,reference)]
    z<-merge(melt(a,id.vars=ids,measure.vars=metrics,variable.name='metric',value.name='saved'),ref,by=c(ids,'metric'))
    z[,`:=`(sample=sample_label,specification=specification)]
    z[,within_rounding:=abs(saved-reference)<=.00050000001 | (is.na(saved)&is.na(reference))]
    z
  }))
  fwrite(saved_checks,file.path(out,'revision_saved_table_crosscheck.csv'))
  if(any(!saved_checks$within_rounding,na.rm=TRUE))stop('Saved-table rounding mismatch; Phase 0 stopped.')
  writeLines(c('Reference: published source at git commit 18b3ae2.',
    'Whole-sample input; original column/row order; seed 0203; pooled pre-RMSE strict percentile.',
    'Comparison is against numerical source replay, not rounded saved HTML.',
    'STOP after Phase 0. No Phase 1 or Phase 2 results generated.'),file.path(out,'revision_phase0_status.txt'))
  if(any(!result$reproduced,na.rm=TRUE))stop('Numerical reproduction mismatch; inspect revision_reproduction_checks.csv')
  invisible(result)
}

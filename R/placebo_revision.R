# Whole-sample placebo revision. Legacy scripts, rule, inputs, and outputs are read-only.
library(data.table)
revision_context <- function(root=getwd()) {
  out<-file.path(root,'_output/placebo_revision/revision'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  record<-function(p,purpose) {
    p<-normalizePath(p,mustWork=TRUE); lp<-file.path(out,'revision_input_paths.csv')
    z<-if(file.exists(lp))fread(lp) else data.table(path=character(),read_at=character(),purpose=character())
    fwrite(rbind(z,data.table(path=p,read_at=as.character(Sys.time()),purpose=purpose)),lp)
  }
  read<-function(p,purpose){record(p,purpose);readRDS(p)}
  csv<-function(p,purpose){record(p,purpose);fread(p)}
  list(root=root,out=out,record=record,read=read,csv=csv)
}
revision_case_iqr <- function(context) {
  p<-file.path(context$root,'analysis/4. the case.rmd'); context$record(p,'reuse Figure 4 IQR expression')
  lines<-readLines(p,warn=FALSE); a<-grep('^iqr_all <- fits',lines)[1]; b<-grep('^iqr_all\\[, sample',lines)[1]-1L
  expr<-parse(text=paste(lines[a:b],collapse='\n'))[[1]][[3]]
  blocks<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name('{')),as.list(expr))
  stopifnot(length(blocks)==1L)
  block<-blocks[[1]]
  function(draws)eval(block,envir=list(.SD=as.data.table(draws)),enclos=environment())
}
revision_mc_mean <- function(x) {
  x<-x[is.finite(x)]; c(estimate=mean(x),mc_se=if(length(x)>1)sd(x)/sqrt(length(x)) else NA_real_,n=length(x))
}
revision_rmse <- function(x) {
  x<-x[is.finite(x)]; r<-sqrt(mean(x^2))
  c(estimate=r,mc_se=if(r>0&&length(x)>1)sd(x^2)/(2*r*sqrt(length(x))) else 0,n=length(x))
}
revision_protect <- function(root) {
  p<-c(list.files(file.path(root,'results'),recursive=TRUE,full.names=TRUE),
    file.path(root,c('analysis/3.2 summary placebo.qmd','analysis/with_covariates/3c.2 summary placebo.rmd','analysis/3.1a.4_heldout_decision.yml')))
  p<-p[file.exists(p)&!dir.exists(p)]; tools::md5sum(p)
}
run_placebo_revision <- function(root=getwd(),specification=c('no_covariates','with_covariates'),window_policy='legacy_horizons') {
  specification<-match.arg(specification); ctx<-revision_context(root); before<-revision_protect(root)
  for(p in c('R/scoring.r','R/heldout_analysis.R')){ctx$record(file.path(root,p),'shared scoring / fixed decision validation');source(file.path(root,p))}
  decision_path<-file.path(root,'analysis/3.1a.4_heldout_decision.yml');ctx$record(decision_path,'fixed rule'); decision<-read_decision(root)
  stopifnot(decision$diagnostic=='ho_standardized',decision$cutoff==2,decision$missing_diagnostic_handling=='keep')
  p<-file.path(root,'data/simulation',if(specification=='with_covariates')'sim_fits_wcov.rds' else 'sim_fits.rds')
  fits<-as.data.table(ctx$read(p,'whole-sample main predictions'))
  gate<-ctx$csv(file.path(ctx$out,'revision_gap_gate.csv'),'required cumulative-year gap gate')
  stopifnot(all(grepl('^PROCEED',gate$gap_gate)))
  draws<-grep('^y0_post_',names(fits),value=TRUE)
  output<-function(z,name)fwrite(z,file.path(ctx$out,paste0('revision_',specification,'_',name,'.csv')))
  window<-fits[,.(n_years=uniqueN(s),complete=setequal(s,0:4)&&.N==5L),by=.(dgp,nc,model,r,i)]
  output(window[,.(county_simulations=.N,incomplete=sum(!complete)),by=.(dgp,nc,model)],'window_check')
  if(any(!window$complete)&&is.null(window_policy))stop('Incomplete five-year windows: user clarification required before changing S1.')
  if(!is.null(window_policy))stopifnot(window_policy %in% c('legacy_horizons','complete_five_year'))
  eligible<-fits[,.(n_models=uniqueN(model),rows=.N,twfe_finite=any(as.character(model)=='TWFE'&is.finite(y0_hat))),by=.(dgp,nc,r,i,s)]
  identical_population<-all(eligible$n_models==6L & eligible$rows==6L & eligible$twfe_finite)
  output(eligible[,.(county_years=.N,not_six_models=sum(n_models!=6L|rows!=6L),missing_twfe=sum(!twfe_finite),S2_identical_to_S1=identical_population),by=.(dgp,nc)],'population_check')
  if(!identical_population)stop('Prepared input does not establish S2 identical to S1; review population_check.')
  if(!identical(window_policy,'legacy_horizons'))stop('Primary results must preserve legacy horizons; use the separately labelled balanced sensitivity.')
  # Read specification-specific audit joins; retain original key translations explicitly.
  diagnostics<-ctx$csv(file.path(root,'_output/placebo_revision/phase0',specification,'revision_diagnostic_join.csv'),'own-specification full-sample diagnostics')
  assignment_path<-file.path(root,'data/derived/heldout/county_assignment.rds')
  assignment<-ctx$read(assignment_path,'assignment and physical half guards')
  physical<-rbindlist(lapply(c('A','B'),function(half_label) {
    half_path<-file.path(root,'data/derived/heldout',if(specification=='with_covariates')'with_covariates' else '',paste0('half_',half_label,'.rds'))
    obj<-ctx$read(half_path,'verify own-specification diagnostic provenance')
    stopifnot(identical(obj$assignment_md5,unname(tools::md5sum(assignment_path))),
      all(obj$rows$half==half_label),all(obj$rows$id %in% assignment[half==half_label,id]))
    derive_heldout(obj$rows,diagnostics_only=TRUE)[,.(model,scenario,sim,id,half,z_original=ho_standardized,p_original=p_two)]
  }))
  provenance<-merge(diagnostics,physical,by=c('model','scenario','sim','id'),all.x=TRUE)
  stopifnot(all(abs(provenance$ho_standardized-provenance$z_original)<1e-12|is.na(provenance$ho_standardized)&is.na(provenance$z_original)),
    all(abs(provenance$p_two-provenance$p_original)<1e-12|is.na(provenance$p_two)&is.na(provenance$p_original)))
  lookup<-diagnostics[,.(dgp=scenario,model,r=sim,i=id,ho_standardized,p_two,ppp_ho_upper,ppp_ho_lower,missing_reason)]
  stopifnot(!anyDuplicated(lookup[,.(dgp,model,r,i)]))
  fits[,`:=`(dgp=as.character(dgp),model=as.character(model))]
  fits[lookup,`:=`(ho_standardized=i.ho_standardized,p_two=i.p_two,missing_reason=i.missing_reason),on=.(dgp,model,r,i)]
  stopifnot(!anyNA(fits$missing_reason))
  # Historical pooling, threshold, and strict inequality exactly match Phase 0.
  fits[,pooled_pre_rmse:=median(rmse_pre),by=.(dgp,r,i)]
  fits[,S3_historical:=pooled_pre_rmse<quantile(pooled_pre_rmse,.9,na.rm=TRUE),by=.(dgp,nc)]
  fits[,S4_own_z2:=!is.finite(ho_standardized)|abs(ho_standardized)<=2]
  fits[,S5_common_z2:=all(S4_own_z2),by=.(dgp,r,i,s)]
  stopifnot(all(!fits$S5_common_z2|fits$S4_own_z2))
  window[,`:=`(dgp=as.character(dgp),model=as.character(model))]
  fits[window,balanced_0_4:=i.complete,on=.(dgp,nc,model,r,i)]
  fits[,`:=`(S1_full=TRUE,B1_balanced_full=balanced_0_4,
    B3_balanced_historical=balanced_0_4 & S3_historical,
    B4_balanced_own_z2=balanced_0_4 & S4_own_z2,B5_balanced_common_z2=balanced_0_4 & S5_common_z2)]
  samples<-c('S1_full','S3_historical','S4_own_z2','S5_common_z2','B1_balanced_full','B3_balanced_historical','B4_balanced_own_z2','B5_balanced_common_z2')
  county<-unique(fits[,c('dgp','nc','model','r','i','rmse_pre','ho_standardized','p_two','balanced_0_4',samples),with=FALSE])
  horizons<-unique(fits[,c('dgp','nc','model','r','i','s','ho_standardized','balanced_0_4',samples),with=FALSE])
  retention<-rbindlist(lapply(samples,function(sample_label){
    base<-if(grepl('^B',sample_label))horizons[balanced_0_4==TRUE] else horizons
    base[,.(sample=sample_label,n_base=.N,
    n_retained=sum(get(sample_label)%in%TRUE),excluded_share=1-mean(get(sample_label)%in%TRUE),
    missing_diagnostic_share=mean(!is.finite(ho_standardized))),by=.(dgp,nc,model,s,r)]}))
  output(retention,'retention_by_simulation')
  output(retention[,.(n_base=mean(n_base),n_retained=mean(n_retained),excluded_share=mean(excluded_share),
    excluded_share_mc_se=sd(excluded_share)/sqrt(.N),missing_diagnostic_share=mean(missing_diagnostic_share),
    missing_diagnostic_share_mc_se=sd(missing_diagnostic_share)/sqrt(.N)),by=.(sample,dgp,nc,model,s)],'retention')
  output(county,'county_flags')
  # Reapply the shared held-out screen on the overlapping county keys.
  exploration_path<-file.path(root,'R/heldout_exploration.R');ctx$record(exploration_path,'shared held-out screening implementation');source(exploration_path)
  overlap<-merge(county,physical,by.x=c('dgp','model','r','i'),by.y=c('scenario','model','sim','id'))
  rule<-data.table(diagnostic='ho_standardized',cutoff=2,rule_type='absolute',rule_id=1L)
  count_check<-overlap[,{
    screen<-keep_matrix_A(data.table(half=half[1L],id=i,ho_standardized=z_original),rule,expected_half=half[1L])$keep[,1L]
    list(overlap_counties=.N,revision_retained=sum(S4_own_z2),heldout_function_retained=sum(screen),flags_identical=all(screen==S4_own_z2))
  },by=.(dgp,nc,model,r,half)]
  output(count_check,'S4_heldout_overlap_check');stopifnot(all(count_check$flags_identical))
  # Same Poisson seed, column ordering and cumulative row order as the historical summary.
  set.seed(0203)
  fits[,(draws):=lapply(.SD,function(lam)rpois(.N,lam)),.SDcols=draws]
  fits[,c('n','y','y0_hat',draws):=lapply(.SD,cumsum),by=.(dgp,nc,model,r,i),.SDcols=c('n','y','y0_hat',draws)]
  ledger<-list();checks<-list()
  for(sample_label in samples) {
    x<-if(sample_label=='S1_full')fits else fits[get(sample_label)==TRUE]
    agg<-x[,c(lapply(.SD,mean),list(pcor=cor(y,y0_hat,use='complete.obs'),scor=cor(y,y0_hat,use='complete.obs',method='spearman'))),
      by=.(nc,dgp,model,s,r),.SDcols=c('y','y0_hat',draws)]
    score_avg(agg);agg[,sample:=sample_label]
    ledger[[sample_label]]<-agg[,setdiff(names(agg),draws),with=FALSE]
    if(sample_label %in% c('S1_full','S3_historical')) {
      ref<-ctx$csv(file.path(root,'_output/placebo_revision/phase0',specification,paste0('revision_',sample_label,'_performance.csv')),'Phase 0 numerical reference')
      setkey(agg,nc,dgp,s,model)
      current<-agg[,summ_score(.SD,sd=sample_label=='S3_historical'),by=.(nc,dgp,model,s)]
      ids<-c('nc','dgp','model','s');metrics<-setdiff(names(ref),ids)
      z<-merge(melt(ref,id.vars=ids,measure.vars=metrics,value.name='reference'),melt(current,id.vars=ids,measure.vars=metrics,value.name='current'),by=c(ids,'variable'))
      z[,`:=`(sample=sample_label,absolute_difference=abs(current-reference))]
      checks[[sample_label]]<-z
      stopifnot(all(z$absolute_difference<1e-10|is.na(z$absolute_difference)))
    }
    message(specification,' ',sample_label,' aggregation complete')
  }
  score_ledger<-rbindlist(ledger);output(score_ledger,'scores_by_simulation');output(rbindlist(checks),'phase0_reproduction')
  common_check<-score_ledger[sample%in%c('S5_common_z2','B5_balanced_common_z2'),.(n_models=uniqueN(model),
    observed_range=max(y)-min(y)),by=.(sample,nc,dgp,s,r)]
  output(common_check,'S5_observed_aggregate_check')
  stopifnot(all(common_check$n_models==6L),all(common_check$observed_range<1e-12))
  group<-c('sample','nc','dgp','model','s')
  performance<-score_ledger[,{
    metrics<-intersect(c('diff','abs_diff','sq_diff','rel_diff','abs_rel_diff','sq_rel_diff','crps','deviance','cvg','intv','pcor','scor'),names(.SD))
    a<-rbindlist(lapply(metrics,function(m){v<-revision_mc_mean(get(m));data.table(metric=m,estimate=v[1],mc_se=v[2],n_simulations=v[3])}))
    v<-revision_rmse(diff);vr<-revision_rmse(rel_diff)
    rbind(a,data.table(metric='RMSE',estimate=v[1],mc_se=v[2],n_simulations=v[3]),
      data.table(metric='relative_RMSE',estimate=vr[1],mc_se=vr[2],n_simulations=vr[3]))
  },by=group]
  output(performance,'performance')
  # Paired gain vs TWFE on the identical S1 base; S2 omitted because it is identical.
  metric_map<-c(RMSE='sq_diff',relative_RMSE='sq_rel_diff',crps='crps',MAE='abs_diff',relative_MAE='abs_rel_diff',
    interval_score='intv',deviance='deviance',coverage='cvg',bias='diff',relative_bias='rel_diff',Pearson='pcor',Spearman='scor')
  baseline<-copy(score_ledger[sample%in%c('S1_full','B1_balanced_full')&model=='TWFE',c('sample','nc','dgp','s','r',unname(metric_map)),with=FALSE])
  setnames(baseline,unname(metric_map),paste0('b_',unname(metric_map)))
  paired<-merge(score_ledger[sample%in%c('S1_full','B1_balanced_full')],baseline,by=c('sample','nc','dgp','s','r'))
  gains<-paired[,rbindlist(lapply(names(metric_map),function(m) {
    x<-get(metric_map[[m]]);b<-get(paste0('b_',metric_map[[m]]))
    good<-is.finite(x)&is.finite(b);x<-x[good];b<-b[good];power<-if(m%in%c('RMSE','relative_RMSE')).5 else 1
    ratio<-(mean(x)/mean(b))^power
    influence<-if(power==1)(x-mean(x))/mean(b)-mean(x)*(b-mean(b))/mean(b)^2 else .5*ratio*(x/mean(x)-b/mean(b))
    data.table(metric=m,percent_difference=100*(ratio-1),mc_se=100*sd(influence)/sqrt(length(x)),n_simulations=length(x))
  })),by=.(sample,nc,dgp,model,s)]
  output(gains,'gain_vs_TWFE')
  # County upper tails: bins formed separately within each replication/model/scenario/subgroup.
  for(diag_label in c('pre_rmse','heldout_abs_z')) {
    value<-if(diag_label=='pre_rmse')fits$rmse_pre else abs(fits$ho_standardized)
    fits[,bin_value:=value]
    ranks<-unique(fits[,.(dgp,nc,model,r,i,bin_value)])
    ranks[,percentile:=fifelse(is.finite(bin_value),(frank(bin_value,ties.method='average',na.last='keep')-.5)/sum(is.finite(bin_value)),NA_real_),by=.(dgp,nc,model,r)]
    fits[ranks,percentile:=i.percentile,on=.(dgp,nc,model,r,i)]
    upper<-rbindlist(lapply(c(paste0('decile_',1:10),'top_5_percent','top_1_percent'),function(bin_label) {
      take<-if(bin_label=='top_5_percent')fits$percentile>=.95 else if(bin_label=='top_1_percent')fits$percentile>=.99 else pmin(10L,floor(fits$percentile*10)+1L)==as.integer(sub('decile_','',bin_label))
      rbindlist(lapply(c('S1_full','B1_balanced_full'),function(sample_label) {
        fits[s %in% c(0,4)&take==TRUE&get(sample_label)==TRUE,.(mse=mean((y-y0_hat)^2),counties=.N),by=.(dgp,nc,model,s,r)][,`:=`(bin=bin_label,sample=sample_label)]
      }))
    }))
    output(upper,paste0('upper_tail_',diag_label,'_by_simulation'))
    summary<-upper[,{
      root_mse<-sqrt(mean(mse));list(RMSE=root_mse,mc_se=sd(mse)/(2*root_mse*sqrt(.N)),n_simulations=.N,mean_counties=mean(counties))
    },by=.(sample,dgp,nc,model,s,bin)]
    output(summary,paste0('upper_tail_',diag_label))
  }
  # Case-study effect IQR: execute its existing draw-wise IQR body without rewriting it.
  iqr_function<-revision_case_iqr(ctx)
  effects<-rbindlist(lapply(c('S1_full','B1_balanced_full'),function(sample_label)fits[get(sample_label)==TRUE,{
    mat<-as.matrix(.SD);effect_draws<-y-mat
    case_iqr<-iqr_function(effect_draws)
    list(aggregate_effect=mean(vapply(as.data.table(effect_draws),mean,numeric(1),na.rm=TRUE),na.rm=TRUE),
      point_prediction_effect=mean(y-y0_hat),county_effect_IQR=case_iqr$diff_y0_hat,
      IQR_posterior_q025=case_iqr$q025,IQR_posterior_q975=case_iqr$q975)
  },by=.(dgp,nc,model,s,r),.SDcols=draws][,sample:=sample_label]))
  output(effects,'null_benchmarks_by_simulation')
  null<-effects[,rbindlist(lapply(c('aggregate_effect','county_effect_IQR'),function(m) {
    qs<-quantile(get(m),c(.025,.5,.975),na.rm=TRUE,names=FALSE)
    data.table(metric=m,q025=qs[1],median=qs[2],q975=qs[3],n_simulations=.N)
  })),by=.(sample,dgp,nc,model,s)]
  output(null,'null_benchmarks')
  calibration<-rbindlist(lapply(c('S1_full','B1_balanced_full'),function(sample_label)county[get(sample_label)==TRUE,.(
    heldout_mid_p_001=mean(p_two<.01,na.rm=TRUE),heldout_mid_p_005=mean(p_two<.05,na.rm=TRUE),
    heldout_mid_p_010=mean(p_two<.10,na.rm=TRUE),available=sum(is.finite(p_two)),missing_share=mean(!is.finite(p_two))),by=.(dgp,nc,model)][,sample:=sample_label]))
  output(calibration,'heldout_calibration')
  pre<-ctx$csv(file.path(ctx$out,paste0('revision_',specification,'_preperiod_county_totals.csv')),'cached exact pre-period totals for both sample bases')
  pre<-merge(pre,county[,.(dgp,nc,model,r,i,balanced_0_4)],by=c('dgp','nc','model','r','i'))
  poisson<-rbindlist(lapply(c('S1_full','B1_balanced_full'),function(sample_label) {
    z<-if(sample_label=='B1_balanced_full')pre[balanced_0_4==TRUE] else pre
    z[,.(dispersion_ratio=sum(pearson_sum)/sum(pearson_n),pearson_county_years=sum(pearson_n),
      observed_zero_share=sum(observed_zeros)/sum(zero_n),expected_zero_share=sum(expected_zeros)/sum(zero_n),
      zero_county_years=sum(zero_n),nonpositive_mean=sum(nonpositive_mean),positive_y_zero_mean=sum(positive_y_zero_mean)),by=.(nc,dgp,model)][,sample:=sample_label]
  }))
  output(poisson,'poisson_diagnostics')
  stopifnot(identical(before,revision_protect(root)))
  message(specification,' revision computation complete; protected files unchanged')
  invisible(performance)
}

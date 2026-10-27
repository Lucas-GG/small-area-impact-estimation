# Additive, A-only common-population sensitivity. Existing section-10 exports are read-only.
twfe_common_prediction_cache <- function(root, obj, model, scenario) {
  source(file.path(root,'R/screening/s3_models.R'),local=TRUE);codes<-s3_codes
  sc<-match(scenario,c('MAR','Latent level','Latent trend'))
  main<-s3_main_path(root,model,sc)
  path<-file.path(root,'.quarto/heldout-aggregate',sprintf('%s_%s.rds',codes[[model]],sc))
  if(!file.exists(path)) stop('Missing A-only saved-fit predictive cache: ',path,'. No complete-panel or half-B file is opened by this sensitivity.')
  cached<-readRDS(path)
  expected<-list(assignment=obj$assignment_md5,predictive_version=1L,main_path=main,
    bytes=file.info(main)$size,mtime=as.numeric(file.info(main)$mtime),
    input=as.numeric(file.info(file.path(root,'data/derived/heldout/half_A_available.rds'))$mtime))
  if(!identical(cached$signature,expected)) stop('Saved-fit provenance is stale or ambiguous: ',path)
  p<-as.data.table(cached$rows)
  stopifnot(all(p$half=='A'),all(p$id %in% unique(obj$rows$id)),!anyDuplicated(p[,.(sim,id)]))
  p
}

twfe_common_membership <- function(z, caches) {
  out<-list()
  for(sc in names(caches)) {
    p<-caches[[sc]]
    e<-p[k_check>=1 & is.finite(y_check) & is.finite(mu_check) & mu_check>0,
      .(scenario=sc,sim,id,k_twfe=k_check,sum_y_twfe=y_check,sum_mu_twfe=mu_check)]
    prepared<-z[model=='TWFE' & scenario==sc & k>=1 & is.finite(sum_y-sum_mu),.(scenario,sim,id)]
    if(!fsetequal(e[,.(scenario,sim,id)],prepared)) stop('TWFE eligibility cannot be identified consistently from saved fit cache and prepared targets: ',sc)
    stopifnot(!anyDuplicated(e[,.(scenario,sim,id)]))
    out[[sc]]<-e
  }
  e<-rbindlist(out)
  groups<-z[model=='TWFE',unique(.SD),.SDcols=c('scenario','sim')]
  stopifnot(nrow(unique(e[,.(scenario,sim)]))==nrow(groups))
  e
}

twfe_common_group <- function(u, keys, keep, expected_half='A') {
  stopifnot(expected_half %in% c('A','B'),all(u$half==expected_half),nrow(keep)==nrow(u),ncol(keep)==nrow(keys))
  nr<-colSums(keep)
  y<-as.numeric(crossprod(keep,u$sum_y))/nr
  mu<-as.numeric(crossprod(keep,u$sum_mu))/nr
  error<-y-mu;error[nr==0]<-NA_real_
  z<-error/sqrt(as.numeric(crossprod(keep,u$sum_mu))/nr^2)
  bounds<-matrix(NA_real_,ncol(keep),2)
  invalid<-rep(0,ncol(keep))
  if(nrow(u)) {
    mat<-do.call(rbind,u$draws);good<-apply(mat,1,function(v)all(is.finite(v)))
    mat[!is.finite(mat)]<-0
    invalid<-as.numeric(crossprod(keep,!good))
    valid<-invalid==0 & nr>0
    totals<-crossprod(keep,mat)/nr
    if(any(valid)) bounds[valid,]<-t(apply(totals[valid,,drop=FALSE],1,quantile,probs=c(.025,.975),names=FALSE,type=7))
  }
  a<-copy(keys)
  a[,`:=`(n_eligible=nrow(u),n_retained=nr,excluded_share=if(nrow(u))1-nr/nrow(u) else NA_real_,
    error=error,error_z=z,observed_aggregate=y,predicted_aggregate=mu,
    lower=bounds[,1],upper=bounds[,2],coverage=as.numeric(y>=bounds[,1]&y<=bounds[,2]),width=bounds[,2]-bounds[,1],
    n_invalid_predictive_counties=invalid)]
  a
}

twfe_common_random_accumulate <- function(acc, u, keys, nr, base_error, seed) {
  random_error<-with_backdate_seed(seed,replicate(100,{idx<-sample.int(nrow(u));cumsum(u$signed_error[idx])/seq_len(nrow(u))}))
  squared<-matrix(NA_real_,nrow(keys),100);ok<-nr>0
  squared[ok,]<-random_error[nr[ok],,drop=FALSE]^2
  if(is.null(acc)) acc<-list(sq=matrix(0,nrow(keys),100),base=matrix(0,nrow(keys),100),n=matrix(0,nrow(keys),100),keys=copy(keys))
  valid<-is.finite(squared);squared[!valid]<-0
  acc$sq<-acc$sq+squared;acc$base<-acc$base+valid*base_error^2;acc$n<-acc$n+valid
  acc
}

run_twfe_common_population <- function(root) {
  library(data.table)
  source(file.path(root,'R/screening/heldout_aggregate.R'))
  source(file.path(root,'R/screening/heldout_analysis.R'));source(file.path(root,'R/backdate.r'));source(file.path(root,'R/screening/heldout_exploration.R'))
  obj<-readRDS(file.path(root,'data/derived/heldout/half_A_available.rds'));aggregate_guard(obj,root)
  base_out<-file.path(root,'_output/analysis/output/3.1a')
  original_paths<-c(list.files(base_out,'[.]csv$',full.names=TRUE),file.path(root,'analysis/S6_screening/s3_decision.yml'))
  before<-tools::md5sum(original_paths)
  out<-file.path(base_out,'twfe_common');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  z<-derive_exploration_A(obj$rows)$rows
  meta<-aggregate_county_metadata(root,unique(z$id));z<-merge(z,meta,by='id',all.x=TRUE)
  stopifnot(all(z$half=='A'),!anyNA(z$noncore))
  scenarios<-c('MAR','Latent level','Latent trend')
  twfe_caches<-setNames(lapply(scenarios,function(sc)twfe_common_prediction_cache(root,obj,'TWFE',sc)),scenarios)
  e<-twfe_common_membership(z,twfe_caches)
  e<-merge(e,meta,by='id',all.x=TRUE)
  source(file.path(root,'R/screening/s3_models.R'),local=TRUE);models<-s3_models
  # Every model must evaluate exactly E, with matching county windows and observations.
  aligned<-merge(e,z,by=c('scenario','sim','id'),all.x=TRUE,allow.cartesian=TRUE,suffixes=c('_E',''))
  stopifnot(nrow(aligned)==nrow(e)*length(models),all(aligned$k>=1),all(is.finite(aligned$sum_y-aligned$sum_mu)),
    all(aligned$k==aligned$k_twfe),all(abs(aligned$sum_y-aligned$sum_y_twfe)<1e-8))
  rules<-unique(fread(file.path(base_out,'rule_grid.csv'))[,.(diagnostic,cutoff,rule_type)])
  stopifnot(fsetequal(rules,exploration_rules()[,.(diagnostic,cutoff,rule_type)]));rules[,rule_id:=seq_len(.N)]
  baseline_key<-data.table(diagnostic='unscreened',cutoff=NA_real_,rule_type='baseline',missing_handling='none')
  ledgers<-list();randoms<-list();decompositions<-list();memberships<-list();checks<-list()
  group_index<-0L
  for(sc in scenarios) for(m in models) {
    group_index<-group_index+1L
    pred<-if(m=='TWFE') copy(twfe_caches[[sc]]) else twfe_common_prediction_cache(root,obj,m,sc)
    pred[,half:=NULL]
    v<-merge(z[scenario==sc & model==m],pred,by=c('sim','id'),all.x=TRUE)
    stopifnot(all(v$k[is.na(v$k_check)]==0));v[is.na(k_check),k_check:=0L]
    stopifnot(all(v$half=='A'),all(v$k==v$k_check),all(abs(v$sum_y-v$y_check)<1e-8,na.rm=TRUE),all(abs(v$sum_mu-v$mu_check)<1e-8,na.rm=TRUE))
    valid<-v$k>=1 & vapply(v$draws,function(a)all(is.finite(a)),logical(1))
    pi<-t(vapply(v$draws[valid],function(a)quantile(a,c(.025,.975),names=FALSE,type=7),numeric(2)))
    stopifnot(all(abs(pi[,1]-v$pi_lower[valid])<1e-8),all(abs(pi[,2]-v$pi_upper[valid])<1e-8))
    v[,signed_error:=sum_y-sum_mu]
    acc<-list()
    for(simulation in sort(unique(v$sim))) {
      reference<-v[sim==simulation]
      screen<-keep_matrix_A(reference,rules)
      full<-reference[k>=1 & is.finite(signed_error)]
      ids<-e[scenario==sc & sim==simulation,id]
      full[,in_E:=id %in% ids]
      stopifnot(sum(full$in_E)==length(ids))
      for(sg in c('noncore','all_A')) {
        original<-full[if(sg=='noncore')noncore else rep(TRUE,.N)]
        u<-original[in_E==TRUE];stopifnot(nrow(u)>0)
        # Preserve the original model's calibrated screen; evaluate only its E rows.
        take<-match(u$id,reference$id)
        keep<-cbind(unscreened=TRUE,screen$keep[take,,drop=FALSE])
        keys<-rbindlist(list(baseline_key,screen$keys[,.(diagnostic,cutoff,rule_type,missing_handling)]))
        a<-twfe_common_group(u,keys,keep)
        a[,`:=`(subgroup=sg,scenario=sc,model=m,sim=simulation)]
        ledgers[[length(ledgers)+1L]]<-a
        seed<-83000000L+group_index*10000L+simulation*10L+match(sg,c('noncore','all_A'))
        acc[[sg]]<-twfe_common_random_accumulate(acc[[sg]],u,keys,a$n_retained,a$error[1],seed)
        checks[[length(checks)+1L]]<-data.table(subgroup=sg,scenario=sc,model=m,sim=simulation,n_E=nrow(u),
          observed_E=a$observed_aggregate[1],count_all_screens_valid=all(a$n_retained==colSums(keep)))
        if(sg=='noncore' && m!='TWFE') {
          c_rows<-original[in_E==FALSE]
          if(nrow(c_rows)) memberships[[length(memberships)+1L]]<-c_rows[,.(scenario,model,sim,id,half,positive_pre,y_ho,k)]
          pop<-list(full=original,E=u,C=c_rows)
          shares<-nrow(c_rows)/nrow(original)
          for(nm in names(pop)) {
            q<-twfe_common_group(pop[[nm]],baseline_key,matrix(TRUE,nrow(pop[[nm]]),1))
            q[,`:=`(subgroup=paste0('noncore/',nm),scenario=sc,model=m,sim=simulation,population=nm,
              C_share=shares,n_full=nrow(original),n_E=nrow(u),n_C=nrow(c_rows))]
            decompositions[[length(decompositions)+1L]]<-q
          }
        }
      }
    }
    for(sg in names(acc)) {
      r<-acc[[sg]];q<-copy(r$keys);q[,`:=`(subgroup=sg,scenario=sc,model=m)]
      q<-q[rep(seq_len(.N),each=100)]
      q[,`:=`(draw=rep(1:100,times=nrow(r$keys)),rmse_ratio=as.vector(t(sqrt(r$sq/r$base))),n_sims=as.vector(t(r$n)))]
      randoms[[length(randoms)+1L]]<-q
    }
    message('TWFE common population complete ',group_index,'/18: ',m,' / ',sc)
  }
  d<-rbindlist(ledgers);summary<-aggregate_summary(d)
  draws<-rbindlist(randoms)
  # E is TWFE's original population: preserve its seeded random reference exactly too.
  original_twfe_draws<-fread(file.path(base_out,'aggregate_random_by_draw.csv'))[model=='TWFE' & diagnostic!='oracle']
  draws<-rbindlist(list(draws[model!='TWFE'],original_twfe_draws),use.names=TRUE,fill=TRUE)
  group_keys<-c('subgroup','scenario','model','diagnostic','cutoff','rule_type','missing_handling')
  bands<-draws[,.(random_lo=quantile(rmse_ratio,.025,na.rm=TRUE),random_hi=quantile(rmse_ratio,.975,na.rm=TRUE),random_mean=mean(rmse_ratio,na.rm=TRUE)),by=group_keys]
  summary<-merge(summary,bands,by=group_keys,all.x=TRUE)
  summary<-summary[order(subgroup!='noncore',scenario,model,diagnostic,cutoff,missing_handling)]
  shortlist<-summary[diagnostic=='unscreened' | (diagnostic=='p_two' & cutoff %in% c(.05,.1)) | (diagnostic=='ho_standardized' & cutoff %in% c(2,2.5)) | (diagnostic=='pre_rmse' & cutoff==10)]
  minima<-summary[diagnostic!='unscreened' & is.finite(rmse_ratio)][order(rmse_ratio,cutoff),.SD[1],by=.(subgroup,scenario,model,diagnostic,missing_handling)]
  minima[,pointwise_significantly_below_one:=rmse_ratio+1.96*rmse_ratio_mc_se<1]
  decomposition<-rbindlist(decompositions);decomp_summary<-aggregate_summary(decomposition)
  decomp_summary[,`:=`(population=sub('noncore/','',subgroup),subgroup='noncore')]
  shares<-decomposition[,.(C_share=mean(C_share),C_share_mc_se=sd(C_share)/sqrt(.N),mean_n_full=mean(n_full),mean_n_E=mean(n_E),mean_n_C=mean(n_C)),by=.(scenario,model,population)]
  decomp_summary<-merge(decomp_summary,shares,by=c('scenario','model','population'))
  base_checks<-rbindlist(checks)
  equal_population<-base_checks[,.(n_models=uniqueN(model),same_n=uniqueN(n_E)==1,same_observed=max(observed_E)-min(observed_E)<1e-8,valid_counts=all(count_all_screens_valid)),by=.(subgroup,scenario,sim)]
  stopifnot(all(equal_population$n_models==length(models)),all(equal_population$same_n),all(equal_population$same_observed),all(equal_population$valid_counts))
  original<-fread(file.path(base_out,'aggregate_summary.csv'))
  join_keys<-c('subgroup','scenario','model','diagnostic','cutoff','rule_type','missing_handling')
  metrics<-c('rmse','bias','mae','relative_rmse','coverage','width','n_retained','mean_observed_aggregate')
  check_twfe<-merge(original[model=='TWFE' & diagnostic!='oracle'],summary[model=='TWFE'],by=join_keys,suffixes=c('_original','_common'))
  stopifnot(nrow(check_twfe)==nrow(summary[model=='TWFE']))
  for(v in metrics) stopifnot(all(abs(check_twfe[[paste0(v,'_original')]]-check_twfe[[paste0(v,'_common')]])<1e-10,na.rm=TRUE))
  old_d<-fread(file.path(base_out,'aggregate_per_simulation.csv'))
  count_checks<-merge(d,old_d,by=c(join_keys,'sim'),suffixes=c('_E','_own'))
  stopifnot(nrow(count_checks)==nrow(d),all(count_checks$n_retained_E<=count_checks$n_retained_own))
  count_checks[,counts_within_original:=n_retained_E<=n_retained_own]
  rep_metrics<-c('n_eligible','n_retained','excluded_share','error','error_z','observed_aggregate','predicted_aggregate','lower','upper','coverage','width')
  twfe_rep<-count_checks[model=='TWFE']
  for(v in rep_metrics) {
    x<-twfe_rep[[paste0(v,'_E')]];y<-twfe_rep[[paste0(v,'_own')]]
    stopifnot(identical(is.na(x),is.na(y)),all(abs(x-y)<1e-10,na.rm=TRUE))
  }
  count_checks<-count_checks[,c(join_keys,'sim','n_retained_E','n_retained_own','counts_within_original'),with=FALSE]
  gain<-d[,.(baseline_coverage=coverage[diagnostic=='unscreened'][1]),by=.(subgroup,scenario,model,sim)]
  gd<-merge(d,gain,by=c('subgroup','scenario','model','sim'))
  gains<-gd[,.(baseline_coverage=mean(baseline_coverage,na.rm=TRUE),coverage_gain=mean(coverage-baseline_coverage,na.rm=TRUE),coverage_gain_mc_se=sd(coverage-baseline_coverage,na.rm=TRUE)/sqrt(sum(is.finite(coverage-baseline_coverage)))),by=join_keys]
  summary<-merge(summary,gains,by=join_keys)
  shortlist<-merge(shortlist,gains,by=join_keys)
  side_by_side<-merge(summary[(diagnostic=='ho_standardized' & cutoff==2)|(diagnostic=='p_two' & cutoff==.05)],original,by=join_keys,suffixes=c('_E','_own'))
  twfe_numeric<-setdiff(names(original)[vapply(original,is.numeric,logical(1))],c('cutoff','n_oracle_unavailable'))
  for(v in twfe_numeric) {
    x<-check_twfe[[paste0(v,'_original')]];y<-check_twfe[[paste0(v,'_common')]]
    stopifnot(identical(is.na(x),is.na(y)),all(abs(x-y)<1e-10,na.rm=TRUE))
  }
  count_check_summary<-count_checks[,.(n_simulations=.N,max_retained_excess=max(n_retained_E-n_retained_own),all_counts_within_original=all(counts_within_original)),by=join_keys]
  full_check<-merge(original[subgroup=='noncore' & model!='TWFE' & diagnostic=='unscreened'],decomp_summary[population=='full'],by=c('subgroup','scenario','model'),suffixes=c('_original','_decomposition'))
  for(v in metrics) stopifnot(all(abs(full_check[[paste0(v,'_original')]]-full_check[[paste0(v,'_decomposition')]])<1e-8))
  e_check<-merge(summary[subgroup=='noncore' & model!='TWFE' & diagnostic=='unscreened'],decomp_summary[population=='E'],by=c('subgroup','scenario','model'),suffixes=c('_common','_decomposition'))
  for(v in metrics) stopifnot(all(abs(e_check[[paste0(v,'_common')]]-e_check[[paste0(v,'_decomposition')]])<1e-8))
  identity<-dcast(decomposition,scenario+model+sim~population,value.var=c('error','observed_aggregate','n_retained'))
  reconstruct<-function(nm) (identity[[paste0('n_retained_E')]]*identity[[paste0(nm,'_E')]]+
    fifelse(identity$n_retained_C>0,identity$n_retained_C*identity[[paste0(nm,'_C')]],0))/identity$n_retained_full
  stopifnot(all(abs(reconstruct('error')-identity$error_full)<1e-8),all(abs(reconstruct('observed_aggregate')-identity$observed_aggregate_full)<1e-8))
  complement<-rbindlist(memberships)
  history<-complement[,.(rows=.N,positive_pre_rows=sum(positive_pre,na.rm=TRUE),positive_heldout_rows=sum(y_ho>0,na.rm=TRUE),missing_history_rows=sum(is.na(positive_pre)|!is.finite(y_ho))),by=.(scenario,model)]
  rank_keys<-c('subgroup','scenario','diagnostic','cutoff','missing_handling')
  rankings<-copy(summary[diagnostic=='unscreened' | (diagnostic=='ho_standardized' & cutoff==2)])
  rankings[,`:=`(absolute_rank=frank(rmse,ties.method='min'),relative_rank=frank(relative_rmse,ties.method='min'),coverage_rank=frank(-coverage,ties.method='min')),by=rank_keys]
  old_ranks<-copy(original[diagnostic=='unscreened' | (diagnostic=='ho_standardized' & cutoff==2)])
  old_ranks[,`:=`(original_absolute_rank=frank(rmse,ties.method='min'),original_relative_rank=frank(relative_rmse,ties.method='min'),original_coverage_rank=frank(-coverage,ties.method='min')),by=rank_keys]
  rankings<-merge(rankings,old_ranks[,c(rank_keys,'model','original_absolute_rank','original_relative_rank','original_coverage_rank'),with=FALSE],by=c(rank_keys,'model'))
  old_scorecard<-fread(file.path(base_out,'aggregate_scorecard.csv'))
  old_scorecard_with_baseline<-rbindlist(list(original[diagnostic=='unscreened'],old_scorecard),use.names=TRUE,fill=TRUE)
  exports<-list(twfe_common_membership=e,twfe_common_per_simulation=d,twfe_common_summary=summary,twfe_common_scorecard=shortlist,
    twfe_common_curve_minima=minima,twfe_common_random_by_draw=draws,twfe_common_random_bands=bands,
    twfe_common_decomposition_per_simulation=decomposition,twfe_common_decomposition=decomp_summary,
    twfe_common_complement_membership=complement,twfe_common_complement_history=history,
    twfe_common_population_checks=base_checks,twfe_common_equal_population_checks=equal_population,
    twfe_common_twfe_baseline_check=check_twfe,twfe_common_rankings=rankings,
    twfe_common_side_by_side=side_by_side,twfe_common_retained_count_checks=count_check_summary,
    twfe_common_retained_count_checks_per_simulation=count_checks,
    twfe_common_original_scorecard_with_unscreened=old_scorecard_with_baseline)
  for(nm in names(exports)) fwrite(exports[[nm]],file.path(out,paste0(nm,'.csv')))
  stopifnot(identical(before,tools::md5sum(names(before))))
  aggregate_guard(obj,root)
  invisible(exports)
}

if(sys.nframe()==0L) run_twfe_common_population(normalizePath(getwd()))

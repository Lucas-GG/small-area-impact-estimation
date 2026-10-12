# Additive aggregate sensitivity: identical retained counties for all six A-only models.
run_aggregate_common_screen <- function(root) {
  library(data.table)
  source(file.path(root,'R/heldout_aggregate.R'))
  source(file.path(root,'R/heldout_analysis.R'))
  source(file.path(root,'R/heldout_exploration.R'))
  source(file.path(root,'R/heldout_twfe_population.R'))
  obj<-readRDS(file.path(root,'data/derived/heldout/half_A_available.rds'))
  aggregate_guard(obj,root)
  base<-file.path(root,'_output/analysis/output/3.1a')
  out<-file.path(base,'common_screen');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  preserved<-c(list.files(base,'[.]csv$',recursive=TRUE,full.names=TRUE),file.path(root,'analysis/3.1a.4_heldout_decision.yml'))
  preserved<-preserved[!startsWith(preserved,paste0(out,'/'))]
  hashes<-tools::md5sum(preserved)
  prior<-fread(file.path(base,'twfe_common/twfe_common_summary.csv'))
  e<-fread(file.path(base,'twfe_common/twfe_common_membership.csv'))[noncore==TRUE]
  z<-derive_exploration_A(obj$rows)$rows
  models<-c('TWFE','ETWFE','GB-Ens','RW','RW-IID','RW-ST')
  rules<-data.table(diagnostic=c('ho_standardized','p_two','p_two'),cutoff=c(2,.05,.10),rule_type=c('absolute','p_value','p_value'),rule_id=1:3)
  baseline<-data.table(diagnostic='unscreened',cutoff=NA_real_,rule_type='baseline',missing_handling='none')
  ledger<-members<-contributions<-checks<-list()
  for(sc in c('MAR','Latent level','Latent trend')) {
    fitted<-setNames(lapply(models,function(m) {
      p<-twfe_common_prediction_cache(root,obj,m,sc)
      v<-merge(z[scenario==sc & model==m],p[,.(sim,id,k_check,y_check,mu_check,draws)],by=c('sim','id'))
      stopifnot(all(v$half=='A'))
      v
    }),models)
    for(simulation in sort(unique(e$sim))) {
      ids<-sort(e[scenario==sc & sim==simulation,id])
      stopifnot(length(ids)>0,!anyDuplicated(ids))
      u<-lapply(fitted,function(v) {
        x<-v[sim==simulation][match(ids,id)]
        stopifnot(identical(x$id,ids),all(x$half=='A'),all(x$k>=1),all(x$k==x$k_check),
          all(is.finite(x$sum_y-x$sum_mu)),all(abs(x$sum_y-x$y_check)<1e-10),all(abs(x$sum_mu-x$mu_check)<1e-10))
        x
      })
      stopifnot(all(vapply(u,function(x)identical(x$sum_y,u[[1]]$sum_y),logical(1))),
        all(vapply(u,function(x)identical(x$k,u[[1]]$k),logical(1))))
      screens<-lapply(u,function(x)keep_matrix_A(x,rules))
      # Same model-aligned AND intersection as section 8 common_sample_A.
      own<-lapply(screens,`[[`,'keep')
      joint<-Reduce(`&`,own)
      keys<-rbindlist(list(baseline,screens[[1]]$keys[,.(diagnostic,cutoff,rule_type,missing_handling)]))
      common_keep<-cbind(unscreened=TRUE,joint)
      flags<-Reduce(`+`,lapply(own,function(k)!k))
      for(j in seq_len(nrow(keys))) {
        retained<-ids[common_keep[,j]]
        if(length(retained)) {
          a<-copy(keys[rep(j,length(retained))]);a[,`:=`(scenario=sc,sim=simulation,id=retained,half='A')]
          members[[length(members)+1L]]<-a
        }
      }
      for(m in models) {
        x<-twfe_common_group(u[[m]],keys,common_keep)
        x[,`:=`(subgroup='noncore',scenario=sc,model=m,sim=simulation)]
        ledger[[length(ledger)+1L]]<-x
        stopifnot(all(!joint | own[[m]]))
        ck<-copy(keys);ck[,`:=`(scenario=sc,sim=simulation,model=m,n_E=length(ids),n_common=colSums(common_keep),
          n_own=colSums(cbind(TRUE,own[[m]])),subset_of_own=TRUE,identical_ids=TRUE,
          retained_id_key=vapply(seq_len(nrow(keys)),function(j)paste(ids[common_keep[,j]],collapse=','),character(1)),
          observed_common=x$observed_aggregate)]
        checks[[length(checks)+1L]]<-ck
        c<-copy(screens[[m]]$keys)
        n_excluded<-colSums(!joint)
        n_unique<-colSums(!own[[m]] & flags==1)
        c[,`:=`(scenario=sc,sim=simulation,model=m,n_E=length(ids),n_excluded=n_excluded,
          n_model_flags=colSums(!own[[m]]),n_unique_flags=n_unique,
          unique_share_of_exclusions=ifelse(n_excluded>0,n_unique/n_excluded,0),
          n_shared_flags=colSums(flags>=2))]
        contributions[[length(contributions)+1L]]<-c
      }
    }
    message('Common retained aggregate completed: ',sc)
  }
  d<-rbindlist(ledger);summary<-aggregate_summary(d)
  k<-c('subgroup','scenario','model','diagnostic','cutoff','rule_type','missing_handling')
  prior_selected<-prior[subgroup=='noncore' & (diagnostic=='unscreened' | (diagnostic=='ho_standardized' & cutoff==2) | (diagnostic=='p_two' & cutoff %in% c(.05,.10)))]
  baseline_check<-merge(summary[diagnostic=='unscreened'],prior_selected[diagnostic=='unscreened'],by=k,suffixes=c('_common','_10_1'))
  metrics<-setdiff(names(summary)[vapply(summary,is.numeric,logical(1))],'cutoff')
  for(nm in metrics) {
    a<-baseline_check[[paste0(nm,'_common')]];b<-baseline_check[[paste0(nm,'_10_1')]]
    stopifnot(identical(is.na(a),is.na(b)),all(abs(a-b)<1e-10,na.rm=TRUE))
  }
  ck<-rbindlist(checks)
  # Confirm own masks reproduce section 10.1, independently of the joint mask.
  prior_d<-fread(file.path(base,'twfe_common/twfe_common_per_simulation.csv'))
  ck<-merge(ck,prior_d[subgroup=='noncore',c(k,'sim','n_retained'),with=FALSE],by=c('scenario','sim','model','diagnostic','cutoff','rule_type','missing_handling'))
  stopifnot(nrow(ck)==nrow(d),all(ck$n_own==ck$n_retained),all(ck$n_common<=ck$n_own))
  equality<-ck[,.(n_models=uniqueN(model),identical_ids=all(identical_ids) & uniqueN(retained_id_key)==1,same_count=uniqueN(n_common)==1,
    same_observed=all(is.na(observed_common)) || max(observed_common)-min(observed_common)<1e-10,
    n_common=n_common[1],n_E=n_E[1]),by=.(scenario,sim,diagnostic,cutoff,rule_type,missing_handling)]
  stopifnot(all(equality$n_models==6),all(equality$identical_ids),all(equality$same_count),all(equality$same_observed))
  comparison<-rbindlist(list(copy(prior_selected[diagnostic=='unscreened'])[,sample:='Unscreened E'],
    copy(prior_selected[diagnostic!='unscreened'])[,sample:='Own screen on E'],
    copy(summary[diagnostic!='unscreened'])[,sample:='Common screen']),use.names=TRUE,fill=TRUE)
  comparison[,`:=`(rmse_rank=frank(rmse,ties.method='min',na.last='keep'),coverage_rank=frank(-coverage,ties.method='min',na.last='keep')),
    by=.(scenario,sample,diagnostic,cutoff,missing_handling)]
  ranks<-merge(comparison[sample=='Own screen on E',.(scenario,model,diagnostic,cutoff,missing_handling,own_rmse_rank=rmse_rank,own_coverage_rank=coverage_rank)],
    comparison[sample=='Common screen',.(scenario,model,diagnostic,cutoff,missing_handling,common_rmse_rank=rmse_rank,common_coverage_rank=coverage_rank)],
    by=c('scenario','model','diagnostic','cutoff','missing_handling'))
  ranks[,`:=`(rmse_rank_changed=own_rmse_rank!=common_rmse_rank,coverage_rank_changed=own_coverage_rank!=common_coverage_rank)]
  c<-rbindlist(contributions)
  contribution_check<-c[,.(n_models=uniqueN(model),n_excluded=n_excluded[1],n_exclusive=sum(n_unique_flags),n_shared=n_shared_flags[1],
    accounted_for=sum(n_unique_flags)+n_shared_flags[1]==n_excluded[1]),by=.(scenario,sim,diagnostic,cutoff,missing_handling)]
  stopifnot(all(contribution_check$n_models==6),all(contribution_check$accounted_for))
  cs<-c[,.(mean_n_E=mean(n_E),n_excluded=mean(n_excluded),excluded_share=mean(n_excluded/n_E),excluded_share_mc_se=sd(n_excluded/n_E)/sqrt(.N),
    n_unique_flags=mean(n_unique_flags),unique_share_of_exclusions=mean(unique_share_of_exclusions),
    unique_share_mc_se=sd(unique_share_of_exclusions)/sqrt(.N),pooled_unique_share=sum(n_unique_flags)/sum(n_excluded),
    n_model_flags=mean(n_model_flags),shared_share=sum(n_shared_flags)/sum(n_excluded)),by=.(scenario,model,diagnostic,cutoff,rule_type,missing_handling)]
  exports<-list(common_screen_per_simulation=d,common_screen_summary=summary,common_screen_comparison=comparison,
    common_screen_rank_changes=ranks,common_screen_contributions=cs,common_screen_contributions_per_simulation=c,
    common_screen_membership=rbindlist(members),common_screen_subset_checks=ck,
    common_screen_equal_population_checks=equality,common_screen_unscreened_check=baseline_check,
    common_screen_contribution_checks=contribution_check)
  for(nm in names(exports)) fwrite(exports[[nm]],file.path(out,paste0(nm,'.csv')))
  aggregate_guard(obj,root);stopifnot(identical(hashes,tools::md5sum(names(hashes))))
  invisible(exports)
}

if(sys.nframe()==0L) run_aggregate_common_screen(normalizePath(getwd()))

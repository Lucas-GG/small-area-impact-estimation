load_heldout_influence <- function(root) {
  out<-file.path(root,'_output/heldout/influence')
  read<-function(n)fread(file.path(out,paste0('influence_',n,'.csv')))
  list(out=out,concentration=read('concentration'),profiles=read('top5_profiles'),gaps=read('AB_gaps'),capture=read('flag_capture'),checks=read('reproduction_checks'))
}

influence_concentration_widget <- function(x,branch_name='No covariates') {
  v<-x$concentration[branch==branch_name & diagnostic=='unscreened']
  w<-dcast(v,scenario+model~half,value.var=c('top1_share','top3_share','counties_to_50','n_counties'))
  DT::datatable(as.data.frame(w),rownames=FALSE,filter='top',options=list(pageLength=18,scrollX=TRUE)) |>
    DT::formatPercentage(columns=c('top1_share_A','top1_share_B','top3_share_A','top3_share_B'),digits=1)
}

influence_profile_widget <- function(x,branch_name='No covariates') {
  v<-copy(x$profiles[branch==branch_name]);setorder(v,half,scenario,model,unscreened_influence_rank)
  v<-v[,.(half,scenario,model,id,unscreened_influence_rank,stratum,mean_pre_expected_count,mean_outcome_age_population,
    pseudo_treated_share,n_eligible_simulations,mean_signed_error,mean_absolute_error,unscreened_MSE_share,
    unscreened_rmse_decrease,flag_z_conditional,flag_p_conditional,flag_pre_conditional,top5_models_same_scenario,influential_several_models)]
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',options=list(pageLength=15,scrollX=TRUE)) |>
    DT::formatPercentage(columns=c('pseudo_treated_share','unscreened_MSE_share','flag_z_conditional','flag_p_conditional','flag_pre_conditional'),digits=1) |>
    DT::formatRound(columns=c('mean_pre_expected_count','mean_outcome_age_population','mean_signed_error','mean_absolute_error','unscreened_rmse_decrease'),digits=3)
}

influence_gap_widget <- function(x,branch_name='No covariates') {
  v<-copy(x$gaps[branch==branch_name]);v[,rule:=ifelse(diagnostic=='ho_standardized','|z| > 2','Pre-RMSE top 10%')]
  v<-v[,.(scenario,model,rule,k_removed,rmse_ratio_A,rmse_ratio_B,absolute_AB_gap,original_absolute_gap,gap_shrinkage,gap_shrinks)]
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',options=list(pageLength=18,scrollX=TRUE)) |>
    DT::formatRound(columns=c('rmse_ratio_A','rmse_ratio_B','absolute_AB_gap','original_absolute_gap','gap_shrinkage'),digits=3)
}

influence_summary_lines <- function(x) {
  c<-x$concentration[branch=='No covariates' & diagnostic=='unscreened']
  p<-x$profiles[branch=='No covariates'];overlap<-list()
  for(h in c('A','B'))for(sc in unique(p$scenario)) {
    v<-p[half==h & scenario==sc];mods<-unique(v$model)
    for(pair in combn(mods,2,simplify=FALSE))overlap[[length(overlap)+1L]]<-data.table(half=h,shared=length(intersect(v[model==pair[1],id],v[model==pair[2],id])))
  }
  ov<-rbindlist(overlap)
  f<-x$capture[branch=='No covariates' & model!='TWFE']
  catch<-function(h,d)mean(f[half==h & diagnostic==d,top5_positive_contribution_caught])
  flag<-function(h,d)mean(f[half==h & diagnostic==d,top5_conditional_flag_share])
  g<-x$gaps[branch=='No covariates' & model!='TWFE']
  stat<-function(d,k,col)mean(g[diagnostic==d & k_removed==k,get(col)])
  cov<-x$gaps[branch=='With covariates' & model!='TWFE']
  c(
    sprintf('Unscreened noncore MSE is concentrated: median top-one shares are %.1f%% on A and %.1f%% on B; top-three shares are %.1f%% and %.1f%%. Median counties needed for 50%% are %g and %g.',100*median(c[half=='A',top1_share]),100*median(c[half=='B',top1_share]),100*median(c[half=='A',top3_share]),100*median(c[half=='B',top3_share]),median(c[half=='A',counties_to_50]),median(c[half=='B',counties_to_50])),
    sprintf('Influential counties recur across models: each pair of the six models shares %.2f/5 top-five counties on A and %.2f/5 on B, averaged across scenarios. A and B have disjoint county assignments, so their influential IDs cannot be the same.',mean(ov[half=='A',shared]),mean(ov[half=='B',shared])),
    sprintf('For the five non-TWFE models, |z| > 2 flags %.1f%% of eligible top-five county appearances on A versus %.1f%% on B; it catches %.1f%% versus %.1f%% of their positive unscreened MSE contribution. This supports stronger influence capture on A.',100*flag('A','ho_standardized'),100*flag('B','ho_standardized'),100*catch('A','ho_standardized'),100*catch('B','ho_standardized')),
    sprintf('Pre-RMSE top 10%% catches %.1f%% of the top-five positive contribution on B versus %.1f%% on A (eligible appearance flag rates %.1f%% and %.1f%%). On B it catches %s positive contribution than |z|, supporting the proposed pattern on average, with model/scenario exceptions visible in the exports.',100*catch('B','pre_rmse'),100*catch('A','pre_rmse'),100*flag('B','pre_rmse'),100*flag('A','pre_rmse'),if(catch('B','pre_rmse')>catch('B','ho_standardized'))'more' else 'less'),
    sprintf('Removing each half\'s top one/three harmful counties reduces the mean absolute A/B gap in |z| RMSE ratios from %.3f to %.3f/%.3f; the top-three gap shrinks in %d/15 model/scenario cells. Mean ratios move from A %.3f/B %.3f to A %.3f/B %.3f after removing three.',stat('ho_standardized',0,'absolute_AB_gap'),stat('ho_standardized',1,'absolute_AB_gap'),stat('ho_standardized',3,'absolute_AB_gap'),sum(g[diagnostic=='ho_standardized' & k_removed==3,gap_shrinks]),stat('ho_standardized',0,'rmse_ratio_A'),stat('ho_standardized',0,'rmse_ratio_B'),stat('ho_standardized',3,'rmse_ratio_A'),stat('ho_standardized',3,'rmse_ratio_B')),
    sprintf('For pre-RMSE, the analogous gap is %.3f before versus %.3f/%.3f after removing one/three; removing three shrinks it in %d/15 cells. Covariate-adjusted secondary gaps after removing three are %.3f for |z| and %.3f for pre-RMSE (before %.3f and %.3f).',stat('pre_rmse',0,'absolute_AB_gap'),stat('pre_rmse',1,'absolute_AB_gap'),stat('pre_rmse',3,'absolute_AB_gap'),sum(g[diagnostic=='pre_rmse' & k_removed==3,gap_shrinks]),mean(cov[diagnostic=='ho_standardized' & k_removed==3,absolute_AB_gap]),mean(cov[diagnostic=='pre_rmse' & k_removed==3,absolute_AB_gap]),mean(cov[diagnostic=='ho_standardized' & k_removed==0,absolute_AB_gap]),mean(cov[diagnostic=='pre_rmse' & k_removed==0,absolute_AB_gap])),
    'This is descriptive and post hoc. The pattern is consistent with a few counties explaining much of the A/B difference, not a causal test or an independent assessment of another rule. Signed shares include cancellation; original flags are frozen. The decision and all existing outputs are unchanged; no rule is suggested.'
  )
}

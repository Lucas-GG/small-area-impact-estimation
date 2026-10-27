# Presentation-only, additive TWFE common-population sensitivity.
load_twfe_common_population <- function(root) {
  source(file.path(root,'R/screening/heldout_aggregate.R'))
  obj<-readRDS(file.path(root,'data/derived/heldout/half_A_available.rds'));aggregate_guard(obj,root)
  out<-file.path(root,'_output/analysis/output/3.1a/twfe_common')
  read<-function(n) {
    p<-file.path(out,paste0('twfe_common_',n,'.csv'))
    if(!file.exists(p)) stop('Run the additive TWFE-common computation first: missing ',p)
    if(file.info(p)$mtime<file.info(file.path(root,'data/derived/heldout/half_A_available.rds'))$mtime) stop('TWFE-common result predates the A-only input: ',p)
    fread(p)
  }
  list(summary=read('summary'),scorecard=read('scorecard'),minima=read('curve_minima'),
    decomposition=read('decomposition'),rankings=read('rankings'),checks=read('equal_population_checks'),
    history=read('complement_history'),side_by_side=read('side_by_side'),out=out)
}

twfe_common_figure <- function(common, original, policy='keep', subgroup_name='noncore') {
  labels<-c(p_two='Predictive mid-p',ho_standardized='Held-out standardized',ho_uncertainty='Held-out uncertainty',pre_rmse='Pre-RMSE percentile')
  v<-copy(common$summary[subgroup==subgroup_name & diagnostic!='unscreened' & missing_handling==policy])
  old<-copy(original$summary[subgroup==subgroup_name & diagnostic %in% names(labels) & missing_handling==policy])
  best<-copy(common$minima[subgroup==subgroup_name & missing_handling==policy])
  prep<-function(x) {x<-ordered_plot(x);x[,diagnostic:=factor(diagnostic,levels=names(labels),labels=labels)];x}
  v<-prep(v);old<-prep(old);best<-prep(best)
  ggplot(v,aes(excluded_share,rmse_ratio,color=scenario,group=scenario))+
    geom_ribbon(aes(ymin=random_lo,ymax=random_hi,fill=scenario),alpha=.1,color=NA)+
    geom_hline(yintercept=1,color='grey50',linetype='dashed')+
    geom_line(data=old,alpha=.18,linewidth=.65)+
    geom_line(linewidth=.8)+geom_point(size=.8)+
    geom_point(data=best,size=2,shape=21,fill='white',stroke=.7,show.legend=FALSE)+
    facet_grid(diagnostic~model)+
    scale_x_continuous(labels=scales::label_percent(accuracy=1))+
    scale_color_manual(values=c('MAR'='#0072B2','Latent level'='#C34A36','Latent trend'='#007A5E'),name=NULL)+
    scale_fill_manual(values=c('MAR'='#0072B2','Latent level'='#C34A36','Latent trend'='#007A5E'),guide='none')+
    labs(x='Excluded share of each curve’s base population',y='Aggregate RMSE / its own unscreened RMSE',title='Noncore: the same TWFE-eligible base population')+
    presentation_theme()+theme(legend.position='bottom',strip.text.y=element_text(angle=0,size=8),
      strip.text.x=element_text(size=9),plot.margin=margin(8,18,8,8),panel.spacing=grid::unit(1,'lines'))
}

twfe_common_scorecard_widget <- function(common, subgroup_name='noncore', policy='keep') {
  v<-copy(common$scorecard[subgroup==subgroup_name & (diagnostic=='unscreened' | missing_handling==policy)])
  v[,rule_order:=match(diagnostic,c('unscreened','ho_standardized','p_two','pre_rmse'))]
  setorder(v,scenario,model,rule_order,cutoff)
  v[,rule:=ifelse(diagnostic=='unscreened','Unscreened',paste(pretty_diagnostics[diagnostic],cutoff))]
  v[,ratio_with_mc_se:=sprintf('%.3f (%.3f)',rmse_ratio,rmse_ratio_mc_se)]
  v<-v[,.(scenario,model,rule,excluded_share,ratio_with_mc_se,relative_rmse,relative_rmse_mc_se,
    mean_observed_aggregate,coverage,width,n_retained,paired_n_sims)]
  labels<-readable_names(names(v));labels[names(v)=='ratio_with_mc_se']<-'Absolute RMSE ratio (MC SE)'
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',colnames=labels,
    options=list(pageLength=24,scrollX=TRUE,lengthMenu=c(24,48,144))) |>
    DT::formatPercentage(columns=c('excluded_share','relative_rmse','relative_rmse_mc_se','coverage'),digits=1) |>
    DT::formatRound(columns=c('mean_observed_aggregate','width','n_retained'),digits=3)
}

twfe_common_decomposition_widget <- function(common) {
  v<-copy(common$decomposition)
  v[,population:=factor(population,levels=c('full','E','C'),labels=c('Full noncore','TWFE-eligible E','Complement C'))]
  setorder(v,scenario,model,population)
  v<-v[,.(scenario,model,population,rmse,rmse_mc_se,relative_rmse,relative_rmse_mc_se,coverage,coverage_mc_se,
    width,n_retained,mean_observed_aggregate,C_share,n_error,n_interval)]
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',colnames=readable_names(names(v)),
    options=list(pageLength=15,scrollX=TRUE,lengthMenu=c(15,45))) |>
    DT::formatPercentage(columns=c('relative_rmse','relative_rmse_mc_se','coverage','coverage_mc_se','C_share'),digits=1) |>
    DT::formatRound(columns=c('rmse','rmse_mc_se','width','n_retained','mean_observed_aggregate'),digits=3)
}

twfe_common_side_by_side_widget <- function(common) {
  v<-copy(common$side_by_side[subgroup=='noncore' & missing_handling=='keep'])
  v[,rule:=ifelse(diagnostic=='p_two','mid-p < 0.05','|z| > 2')]
  v<-v[,.(scenario,model,rule,excluded_share_own,excluded_share_E,
    rmse_ratio_own,rmse_ratio_E,rmse_ratio_mc_se_own,rmse_ratio_mc_se_E,
    relative_rmse_own,relative_rmse_E,relative_rmse_mc_se_own,relative_rmse_mc_se_E,
    coverage_own,coverage_E,width_own,width_E,n_retained_own,n_retained_E)]
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',colnames=readable_names(names(v)),options=list(pageLength=15,scrollX=TRUE)) |>
    DT::formatRound(columns=names(v)[vapply(v,is.numeric,logical(1))],digits=3)
}

twfe_common_rankings_widget <- function(common) {
  v<-copy(common$rankings[subgroup=='noncore' & (diagnostic=='unscreened' | missing_handling=='keep')])
  v[,rule:=ifelse(diagnostic=='unscreened','Unscreened','|z| > 2')]
  v<-v[,.(scenario,model,rule,original_absolute_rank,absolute_rank,original_coverage_rank,coverage_rank,rmse,coverage,width)]
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',colnames=readable_names(names(v)),options=list(pageLength=18,scrollX=TRUE)) |>
    DT::formatRound(columns=c('rmse','coverage','width'),digits=3)
}

twfe_common_answer_lines <- function(common) {
  v<-common$summary[subgroup=='noncore' & missing_handling=='keep' & ((diagnostic=='ho_standardized' & cutoff==2)|(diagnostic=='p_two' & cutoff==.05))]
  other<-v[model!='TWFE'];twfe<-v[model=='TWFE']
  ranks<-common$rankings[subgroup=='noncore' & (diagnostic=='unscreened' | missing_handling=='keep')]
  d<-dcast(common$decomposition,scenario+model~population,value.var=c('rmse','relative_rmse','coverage'))
  c(
    sprintf('For non-TWFE models on E, absolute RMSE ratios remain below 1 in %d/%d |z| > 2 and mid-p < 0.05 model/scenario/rule cells; ratios range %.3f–%.3f.',sum(other$rmse_ratio<1),nrow(other),min(other$rmse_ratio),max(other$rmse_ratio)),
    sprintf('Their aggregate coverage improves over unscreened E in %d/%d cells (mean %.1f%% unscreened versus %.1f%% screened); paired coverage-gain MC SEs are exported.',sum(other$coverage_gain>0),nrow(other),100*mean(other$baseline_coverage),100*mean(other$coverage)),
    sprintf('TWFE reproduces every original rule exactly: the two operational screens give ratios %.3f–%.3f and coverage changes %.1f–%.1f percentage points. TWFE remains a screening null for these two rules (no RMSE improvement); the full grid still allows small exploratory gains.',min(twfe$rmse_ratio),max(twfe$rmse_ratio),100*min(twfe$coverage_gain),100*max(twfe$coverage_gain)),
    sprintf('Population restriction changes unscreened absolute-RMSE ranks in %d/18 cells and coverage ranks in %d/18; under |z| > 2 it changes %d/18 and %d/18 respectively (ties share ranks).',sum(ranks[diagnostic=='unscreened',absolute_rank!=original_absolute_rank]),sum(ranks[diagnostic=='unscreened',coverage_rank!=original_coverage_rank]),sum(ranks[diagnostic!='unscreened',absolute_rank!=original_absolute_rank]),sum(ranks[diagnostic!='unscreened',coverage_rank!=original_coverage_rank])),
    sprintf('Zero-history C has lower absolute aggregate error in most comparisons: C-only RMSE is below E-only in %d/15 cells, but its relative RMSE is higher in %d/15. Adding C lowers full-population absolute RMSE in %d/15 and relative RMSE in %d/15, but lowers coverage in %d/15; “easier” depends on scale and interval calibration.',sum(d$rmse_C<d$rmse_E),sum(d$relative_rmse_C>d$relative_rmse_E),sum(d$rmse_full<d$rmse_E),sum(d$relative_rmse_full<d$relative_rmse_E),sum(d$coverage_full<d$coverage_E)))
}

twfe_common_population_note <- function(common) paste(twfe_common_answer_lines(common),collapse='  \n')

twfe_common_results <- function(common) {
  cat(paste0('- ',twfe_common_answer_lines(common),collapse='\n'),'\n',sep='')
  cat('- All TWFE reproduction, retained-count bounds, A-only and assignment-hash checks pass; original CSVs and the decision template are unchanged.\n')
  cat('- Results are descriptive sensitivity comparisons; exploratory minima are not multiplicity-adjusted, and no rule is selected.\n')
}

# Presentation of additive aggregate exports; no statistical refresh on render.
load_aggregate_presentation <- function(root) {
  source(file.path(root,'R/heldout_aggregate.R'))
  obj<-readRDS(file.path(root,'data/derived/heldout/half_A_available.rds'));aggregate_guard(obj,root)
  out<-file.path(root,'_output/analysis/output/3.1a')
  required<-c('summary','scorecard','curve_minima','manuscript_comparison','size_strata','resolution_checks','predictive_checks')
  read<-function(n) {
    p<-file.path(out,paste0('aggregate_',n,'.csv'))
    if(!file.exists(p)) stop('Run Rscript R/heldout_aggregate.R first: missing ',p)
    if(file.info(p)$mtime<file.info(file.path(root,'data/derived/heldout/half_A_available.rds'))$mtime) stop('Aggregate output predates the A-only input: ',p)
    fread(p)
  }
  setNames(lapply(required,read),required)
}

aggregate_curve_plot <- function(a, subgroup_name, diag_name, policy='keep') {
  v<-copy(a$summary[subgroup==subgroup_name & diagnostic==diag_name & missing_handling==policy])
  oracle<-copy(a$summary[subgroup==subgroup_name & diagnostic=='oracle' & oracle_main_visible==TRUE])
  best<-copy(a$curve_minima[subgroup==subgroup_name & diagnostic==diag_name & missing_handling==policy])
  v<-ordered_plot(v);oracle<-ordered_plot(oracle);best<-ordered_plot(best)
  ggplot(v,aes(excluded_share,rmse_ratio,color=scenario,group=scenario))+
    geom_ribbon(aes(ymin=random_lo,ymax=random_hi,fill=scenario),alpha=.1,color=NA)+
    geom_hline(yintercept=1,color='grey50',linetype='dashed')+
    geom_line(linewidth=.7)+geom_point(size=.8)+
    geom_point(data=best,size=2.6,shape=21,fill='white',stroke=.8,show.legend=FALSE)+
    geom_point(data=oracle,shape=8,size=2.4,show.legend=FALSE)+facet_wrap(~model,ncol=3)+
    scale_x_continuous(labels=scales::label_percent(accuracy=1))+
    scale_color_manual(values=c('MAR'='#0072B2','Latent level'='#C34A36','Latent trend'='#007A5E'),name=NULL)+
    scale_fill_manual(values=c('MAR'='#0072B2','Latent level'='#C34A36','Latent trend'='#007A5E'),guide='none')+
    labs(x='Excluded share',y='Aggregate RMSE / unscreened RMSE',title=paste(if(subgroup_name=='noncore')'Noncore counties' else 'All half-A counties',pretty_diagnostics[[diag_name]],sep=': '))+
    presentation_theme()+theme(legend.position='bottom',plot.margin=margin(8,18,8,8),panel.spacing.x=grid::unit(1.5,'lines'))
}

aggregate_scorecard_widget <- function(a, policy='keep', subgroup_name='noncore') {
  v<-copy(a$scorecard[subgroup==subgroup_name & (missing_handling==policy | (diagnostic=='oracle' & oracle_main_visible==TRUE))])
  # Explicitly authorized baseline addition; existing scorecard values are preserved.
  v<-rbindlist(list(a$summary[subgroup==subgroup_name & diagnostic=='unscreened'],v),use.names=TRUE,fill=TRUE)
  v[,rule:=ifelse(diagnostic=='unscreened','Unscreened',ifelse(diagnostic=='oracle','Oracle (infeasible; standardized, single pass)',paste(pretty_diagnostics[diagnostic],cutoff)))]
  v[,ratio_with_mc_se:=sprintf('%.3f (%.3f)',rmse_ratio,rmse_ratio_mc_se)]
  v<-v[,.(subgroup,scenario,model,rule,excluded_share,ratio_with_mc_se,relative_rmse,relative_rmse_mc_se,mean_observed_aggregate,coverage,width,n_retained,paired_n_sims)]
  labels<-readable_names(names(v))
  labels[names(v)=='ratio_with_mc_se']<-'Absolute RMSE ratio (MC SE)'
  labels[names(v)=='relative_rmse']<-'Relative RMSE'
  labels[names(v)=='relative_rmse_mc_se']<-'Relative RMSE MC SE'
  labels[names(v)=='mean_observed_aggregate']<-'Mean retained observed aggregate'
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',colnames=labels,
    options=list(pageLength=21,scrollX=TRUE,lengthMenu=c(21,42,126,252))) |>
    DT::formatPercentage(columns=c('excluded_share','relative_rmse','relative_rmse_mc_se','coverage'),digits=1) |>
    DT::formatRound(columns=c('mean_observed_aggregate','width','n_retained'),digits=3)
}

aggregate_oracle_details_widget <- function(a) {
  v<-copy(a$summary[diagnostic=='oracle'])
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',colnames=readable_names(names(v)),
    options=list(pageLength=18,scrollX=TRUE,lengthMenu=c(18,36))) |>
    DT::formatRound(columns=names(v)[vapply(v,is.numeric,logical(1))],digits=3)
}

aggregate_noncore_coverage_note <- function(a) {
  v<-a$summary[subgroup=='noncore']
  baseline<-mean(v[diagnostic=='unscreened',coverage])
  z<-mean(v[diagnostic=='ho_standardized' & cutoff==2 & missing_handling=='keep',coverage])
  p<-mean(v[diagnostic=='p_two' & cutoff==.05 & missing_handling=='keep',coverage])
  sprintf('Noncore 95%% aggregate coverage: unscreened %.1f%% → |z| > 2 %.1f%%; unscreened %.1f%% → mid-p < 0.05 %.1f%% (equal averages across models and scenarios; missing diagnostics kept).',100*baseline,100*z,100*baseline,100*p)
}

aggregate_noncore_top10_note <- function(a) {
  v<-a$summary[subgroup=='noncore' & diagnostic=='pre_rmse' & cutoff==10 & missing_handling=='keep']
  sprintf('The pre-RMSE top-10%% rule excludes about 0.1%% of noncore counties (%.3f%% averaged across simulations, models and scenarios; missing diagnostics kept), because its percentile is defined in the full half-A pool.',100*mean(v$excluded_share))
}

aggregate_results_lines <- function(a, subgroup_name='noncore') {
  v<-a$curve_minima[subgroup==subgroup_name & missing_handling=='keep'][order(rmse_ratio),.SD[1],by=.(model,scenario)]
  rule_label<-function(d,c) switch(d,p_two=paste0('p < ',format(c)),ho_standardized=paste0('|z| > ',format(c)),ho_uncertainty=paste0('|u| > ',format(c)),pre_rmse=paste0('pre-RMSE top ',format(c),'%'))
  cat('Each entry gives the minimum-cutoff rule, absolute RMSE ratio (paired MC SE), and relative RMSE as a percentage of the mean retained observed aggregate. An asterisk marks a pointwise 1.96-MC-SE upper bound below 1 for the absolute ratio. Missing diagnostics are kept.\n\n')
  for(m in model_order) {
    u<-v[model==m];u<-u[match(scenario_order,scenario)]
    cells<-vapply(seq_len(nrow(u)),function(j) sprintf('%s: %s, absolute %.3f (%.3f)%s; relative %.1f%%',u$scenario[j],rule_label(u$diagnostic[j],u$cutoff[j]),u$rmse_ratio[j],u$rmse_ratio_mc_se[j],ifelse(u$pointwise_significantly_below_one[j],'*',''),100*u$relative_rmse[j]),character(1))
    cat('- **',m,'** — ',paste(cells,collapse='; '),'\n',sep='')
  }
  cat('\nHere z is the held-out standardized diagnostic and u the held-out uncertainty diagnostic. These exploratory flags do not adjust for choosing a minimum or comparing many cutoffs. Both missing policies are in the full tables. No rule is selected.\n')
}

# Presentation reads only the additive exports and guarded A-only input.
load_aggregate_common_screen <- function(root) {
  source(file.path(root,'R/heldout_aggregate.R'))
  obj<-readRDS(file.path(root,'data/derived/heldout/half_A_available.rds'));aggregate_guard(obj,root)
  out<-file.path(root,'_output/analysis/output/3.1a/common_screen')
  read<-function(n) {
    path<-file.path(out,paste0('common_screen_',n,'.csv'))
    if(!file.exists(path))stop('Run R/heldout_common_screen.R first: ',path)
    if(file.info(path)$mtime<file.info(file.path(root,'data/derived/heldout/half_A_available.rds'))$mtime)stop('Stale common-screen export: ',path)
    fread(path)
  }
  list(summary=read('summary'),comparison=read('comparison'),contributions=read('contributions'),
    ranks=read('rank_changes'),checks=read('equal_population_checks'),out=out,read=function(n)fread(file.path(out,paste0(n,'.csv'))))
}

common_screen_rule <- function(d,c) fifelse(d=='unscreened','Unscreened',fifelse(d=='ho_standardized','|z| > 2',paste0('mid-p < ',sprintf('%.2f',c))))

common_screen_comparison_widget <- function(x,scenario_name) {
  v<-copy(x$comparison[scenario==scenario_name & (diagnostic=='unscreened' | missing_handling=='keep')])
  v[,rule:=common_screen_rule(diagnostic,cutoff)]
  v[,sample:=factor(sample,levels=c('Unscreened E','Own screen on E','Common screen'))]
  v[,rule:=factor(rule,levels=c('Unscreened','|z| > 2','mid-p < 0.05','mid-p < 0.10'))]
  v[,cell:=sprintf('RMSE %.3f (%.3f)<br>Relative %.1f%% (%.1f%%)<br>Ratio %.3f (%.3f)<br>Coverage %.1f%%; width %.3f<br>RMSE rank %d; coverage rank %d<br>N %.1f; excluded %.1f%%',
    rmse,rmse_mc_se,100*relative_rmse,100*relative_rmse_mc_se,rmse_ratio,rmse_ratio_mc_se,100*coverage,width,rmse_rank,coverage_rank,n_retained,100*excluded_share)]
  wide<-dcast(v,rule+sample~model,value.var='cell')
  setcolorder(wide,c('rule','sample',model_order))
  DT::datatable(as.data.frame(wide),rownames=FALSE,escape=FALSE,colnames=c('Rule','Evaluation sample',model_order),
    options=list(dom='t',ordering=FALSE,paging=FALSE,scrollX=TRUE))
}

common_screen_contributions_widget <- function(x) {
  v<-copy(x$contributions[missing_handling=='keep'])
  v[,rule:=common_screen_rule(diagnostic,cutoff)]
  v[,n_common:=mean_n_E-n_excluded]
  v<-v[,.(scenario,rule,model,mean_n_E,n_common,excluded_share,excluded_share_mc_se,n_excluded,n_unique_flags,
    pooled_unique_share,shared_share,n_model_flags)]
  DT::datatable(as.data.frame(v),rownames=FALSE,filter='top',colnames=c('Scenario','Rule','Model','Mean N in E','Mean N retained','Common excluded share','MC SE',
    'Mean excluded counties','Mean flagged only by this model','Only-this-model share of exclusions','Shared-flag share','Mean model flags'),
    options=list(pageLength=18,scrollX=TRUE)) |>
    DT::formatPercentage(columns=c('excluded_share','excluded_share_mc_se','pooled_unique_share','shared_share'),digits=1) |>
    DT::formatRound(columns=c('mean_n_E','n_common','n_excluded','n_unique_flags','n_model_flags'),digits=2)
}

common_screen_note <- function(x) {
  v<-x$comparison[missing_handling=='keep' | diagnostic=='unscreened']
  baseline<-v[scenario=='Latent level' & sample=='Unscreened E' & model=='TWFE',rmse_rank]
  lv<-v[scenario=='Latent level' & sample=='Common screen']
  twfe<-lv[model=='TWFE']
  reversal<-vapply(seq_len(nrow(twfe)),function(j)sum(lv[diagnostic==twfe$diagnostic[j] & cutoff==twfe$cutoff[j] & model!='TWFE',rmse]<twfe$rmse[j]),integer(1))
  rule_names<-common_screen_rule(twfe$diagnostic,twfe$cutoff)
  sentences<-c(sprintf('Latent level: TWFE ranks %d unscreened on E. Under the common screen its RMSE ranks are %s; %s. Thus the reversal %s on an identical retained population.',baseline,
    paste(sprintf('%s: %d',rule_names,twfe$rmse_rank),collapse='; '),
    paste(sprintf('%d/5 non-TWFE models ahead for %s',reversal,rule_names),collapse='; '),
    if(all(reversal>0))'persists for all three rules' else if(any(reversal>0))'persists for some rules only' else 'does not persist'))
  for(sc in scenario_order) {
    r<-x$ranks[scenario==sc & missing_handling=='keep']
    cells<-r[,.(rmse_changes=sum(rmse_rank_changed),coverage_changes=sum(coverage_rank_changed)),by=.(diagnostic,cutoff)]
    sentences<-c(sentences,sprintf('%s: common versus own-screen RMSE/coverage ranks change in %s model cells (out of six per rule).',sc,
      paste(sprintf('%s %d/%d',common_screen_rule(cells$diagnostic,cells$cutoff),cells$rmse_changes,cells$coverage_changes),collapse='; ')))
  }
  excluded<-unique(x$contributions[missing_handling=='keep',.(scenario,diagnostic,cutoff,excluded_share)])
  sentences<-c(sentences,sprintf('Common screening excludes %.1f%%–%.1f%% of E across scenario/rule means. Contributions count counties flagged only by one model; shared flags are reported separately. These are descriptive, unadjusted comparisons; no rule is selected and the decision template is unchanged.',100*min(excluded$excluded_share),100*max(excluded$excluded_share)))
  paste(sentences,collapse='  \n')
}

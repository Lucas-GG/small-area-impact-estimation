# Presentation correction only: quantiles of simulation-averaged random reductions.
# Existing point estimates, target definitions and screening decisions are untouched.
simulation_averaged_band <- function(draws_by_sim) {
  averages <- colMeans(draws_by_sim, na.rm=TRUE)
  averages <- averages[is.finite(averages)]
  if(!length(averages)) return(c(mean=NA_real_,lo=NA_real_,hi=NA_real_))
  c(mean=mean(averages),lo=quantile(averages,.025,names=FALSE),hi=quantile(averages,.975,names=FALSE))
}
random_draw_reductions <- function(ref, excluded) {
  if(is.null(ref)||excluded>=ref$n) return(rep(NA_real_,100L))
  removed<-if(excluded) ref$prefix[excluded,] else rep(0,ncol(ref$prefix))
  1-(ref$sum-removed)/(ref$n-excluded)/ref$baseline
}
refresh_presentation_bands <- function(root, obj, destination) {
  source(file.path(root,'R/backdate.r'))
  source(file.path(root,'R/screening/heldout_exploration.R'))
  stopifnot(all(obj$rows$half=='A'))
  saved<-file.path(root,'_output/heldout/exploration/available_window')
  input<-file.path(root,'data/derived/heldout/half_A_available.rds')
  signature<-list(version=1L,a_md5=unname(tools::md5sum(input)),
    source_info=file.info(file.path(saved,'retention_per_simulation.csv'))[,c('size','mtime')])
  cache<-file.path(destination,'random_band_cache.rds')
  if(file.exists(cache) && identical(readRDS(cache)$signature,signature)) return(invisible(NULL))
  trace<-fread(file.path(saved,'retention_per_simulation.csv'),select=c('scenario','model','sim','diagnostic','cutoff',
    'missing_handling','target','n_eligible','n_retained','random_mean'))[target %in% c('err_raw','err_z')]
  rules<-exploration_rules()
  keys<-rbindlist(list(copy(rules)[,missing_handling:='keep'],copy(rules)[,missing_handling:='exclude']))
  token<-function(d,c,p) paste(d,c,p,sep='|')
  key_tokens<-token(keys$diagnostic,keys$cutoff,keys$missing_handling)
  all_rules<-list();all_grid<-list();grid<-seq(0,.95,length.out=96L)
  source(file.path(root,'R/screening/s3_models.R'),local=TRUE);models<-s3_models;scenarios<-c('MAR','Latent level','Latent trend')
  for(sc in scenarios) for(m in models) for(target_name in c('err_raw','err_z')) {
    message('Correcting random reference: ',m,' / ',sc,' / ',target_name)
    # Match the legacy merge's within-simulation county order exactly.
    rows<-obj$rows[scenario==sc & model==m][order(id,sim)]
    counts<-trace[scenario==sc & model==m & target==target_name]
    rule_sum<-rule_n<-matrix(0,nrow(keys),100L)
    grid_sum<-grid_n<-matrix(0,length(grid),100L)
    for(s in sort(unique(rows$sim))) {
      y<-rows[sim==s,get(target_name)];y<-y[is.finite(y)]
      ref<-random_reference(y,2032032L+match(sc,scenarios)*100000L+s*100L+match(target_name,exploration_targets))
      v<-counts[sim==s];v<-v[match(key_tokens,token(diagnostic,cutoff,missing_handling))]
      stopifnot(nrow(v)==nrow(keys),!anyNA(v$sim),all(v$n_eligible==length(y)))
      actual<-t(vapply(v$n_eligible-v$n_retained,function(e) random_draw_reductions(ref,e),numeric(100)))
      # Regenerated per-simulation means must reproduce the existing exported results.
      check<-rowMeans(actual,na.rm=TRUE)
      stopifnot(all(abs(check-v$random_mean)<1e-10,na.rm=TRUE),identical(is.finite(check),is.finite(v$random_mean)))
      good<-is.finite(actual);rule_sum<-rule_sum+replace(actual,!good,0);rule_n<-rule_n+good
      g<-t(vapply(floor(length(y)*grid+.5),function(e) random_draw_reductions(ref,e),numeric(100)))
      good<-is.finite(g);grid_sum<-grid_sum+replace(g,!good,0);grid_n<-grid_n+good
    }
    summarize<-function(sums,ns) {
      avg<-sums/ns
      data.table(random_mean_corrected=rowMeans(avg,na.rm=TRUE),
        random_lo_corrected=apply(avg,1,quantile,.025,na.rm=TRUE,names=FALSE),
        random_hi_corrected=apply(avg,1,quantile,.975,na.rm=TRUE,names=FALSE),
        n_sims=apply(ns,1,min),n_random_draws=100L)
    }
    a<-cbind(copy(keys),summarize(rule_sum,rule_n));a[,`:=`(scenario=sc,model=m,target=target_name)]
    b<-cbind(data.table(excluded_share=grid),summarize(grid_sum,grid_n));b[,`:=`(scenario=sc,model=m,target=target_name)]
    all_rules[[length(all_rules)+1L]]<-a;all_grid[[length(all_grid)+1L]]<-b
  }
  a<-rbindlist(all_rules);b<-rbindlist(all_grid)
  fwrite(a,file.path(destination,'random_bands_by_rule.csv'))
  fwrite(b,file.path(destination,'random_band_plot_grid.csv'))
  old<-fread(file.path(saved,'retention_summary.csv'))
  comparison<-merge(a,old[,.(scenario,model,target,diagnostic,cutoff,missing_handling,random_lo,random_hi)],
    by=c('scenario','model','target','diagnostic','cutoff','missing_handling'))
  comparison[,`:=`(old_width=random_hi-random_lo,new_width=random_hi_corrected-random_lo_corrected)]
  fwrite(comparison,file.path(destination,'random_band_width_comparison.csv'))
  saveRDS(list(signature=signature,prepared=Sys.time()),cache)
}

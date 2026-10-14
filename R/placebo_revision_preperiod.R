library(data.table)
run_placebo_revision_preperiod <- function(root=getwd(),specification=c('no_covariates','with_covariates')) {
  specification<-match.arg(specification);source(file.path(root,'R/placebo_revision.R'));ctx<-revision_context(root)
  lookup<-ctx$csv(file.path(root,'_output/placebo_revision/phase0',specification,'revision_diagnostic_join.csv'),'whole-sample pre-period county keys')
  mapping<-c(TWFE='twfe',ETWFE='etwfe','GB-Ens'='xgb',RW='rw','RW-IID'='rwi','RW-ST'='rwst')
  scenarios<-c('MAR','Latent level','Latent trend');rows<-list()
  cache_dir<-file.path(root,'.quarto/placebo-revision-preperiod');dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  for(model_label in names(mapping))for(scenario_index in 1:3) {
    scenario_label<-scenarios[scenario_index]
    path<-file.path(root,'data/simulation',if(specification=='with_covariates')'wcov' else 'wocov',
      paste0(if(specification=='with_covariates')'wcov_imputation_' else 'new_imputation_',mapping[[model_label]],'_',scenario_index,'.rds'))
    cache<-file.path(cache_dir,paste0(specification,'_',mapping[[model_label]],'_',scenario_index,'.rds'))
    signature<-list(path=normalizePath(path),bytes=file.info(path)$size,mtime=as.character(file.info(path)$mtime))
    if(file.exists(cache)) {
      saved<-ctx$read(cache,'cached exact pre-period totals')
      if(identical(saved$signature,signature)){rows[[length(rows)+1]]<-saved$rows;next}
    }
    message('Pre-period read: ',specification,' ',model_label,' ',scenario_label)
    raw<-as.data.table(ctx$read(path,'saved main pre-period predictions; no refit'))
    pre<-raw[year<start_year_sim & is.finite(y)&is.finite(y0_hat),.(r,i,y,mu=y0_hat)]
    rm(raw);gc()
    key<-lookup[model==model_label & scenario==scenario_label,.(r=sim,i=id,nc=subgroup)]
    pre<-merge(pre,key,by=c('r','i'))
    result<-pre[,.(pearson_sum=sum((y[mu>0]-mu[mu>0])^2/mu[mu>0]),pearson_n=sum(mu>0),
      observed_zeros=sum(y[mu>=0]==0),expected_zeros=sum(exp(-mu[mu>=0])),zero_n=sum(mu>=0),
      nonpositive_mean=sum(mu<=0),positive_y_zero_mean=sum(mu==0&y>0)),by=.(nc,r,i)]
    result[,`:=`(model=model_label,dgp=scenario_label)]
    saveRDS(list(signature=signature,rows=result),cache);rows[[length(rows)+1]]<-result
    rm(pre);gc();message('Pre-period totals saved: ',specification,' ',model_label,' ',scenario_label)
  }
  all<-rbindlist(rows)
  fwrite(all,file.path(ctx$out,paste0('revision_',specification,'_preperiod_county_totals.csv')))
  summary<-all[,.(dispersion_ratio=sum(pearson_sum)/sum(pearson_n),pearson_county_years=sum(pearson_n),
    observed_zero_share=sum(observed_zeros)/sum(zero_n),expected_zero_share=sum(expected_zeros)/sum(zero_n),
    zero_county_years=sum(zero_n),nonpositive_mean=sum(nonpositive_mean),positive_y_zero_mean=sum(positive_y_zero_mean)),by=.(nc,dgp,model)]
  fwrite(summary,file.path(ctx$out,paste0('revision_',specification,'_poisson_diagnostics.csv')))
  invisible(summary)
}

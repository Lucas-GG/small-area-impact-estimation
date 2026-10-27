# Fixed-decision evaluation; B execution additionally needs explicit session authorization.
estimate_gains <- function(d) {
  own<-d[subgroup=='noncore' & population=='own' & (diagnostic=='unscreened' | missing_handling=='keep')]
  b<-own[diagnostic=='unscreened',.(scenario,model,sim,eb=error,yb=observed_aggregate,cb=coverage)]
  x<-merge(own[diagnostic=='ho_standardized' & cutoff==2],b,by=c('scenario','model','sim'))
  influence<-x[, {
    ok<-is.finite(error)&is.finite(eb)&is.finite(observed_aggregate)&is.finite(yb)&is.finite(coverage)&is.finite(cb)
    s<-error[ok];b<-eb[ok];ys<-observed_aggregate[ok];yb<-yb[ok];r<-sqrt(mean(s^2)/mean(b^2));rr<-r*mean(yb)/mean(ys)
    data.table(sim=sim[ok],absolute_ratio=r,relative_ratio=rr,coverage_change=mean(coverage[ok]-cb[ok]),
      absolute_if=(s^2-r^2*b^2)/(2*r*mean(b^2)),
      relative_if=rr*((s^2-mean(s^2))/(2*mean(s^2))-(b^2-mean(b^2))/(2*mean(b^2))+(yb-mean(yb))/mean(yb)-(ys-mean(ys))/mean(ys)),
      coverage_if=coverage[ok]-cb[ok]-mean(coverage[ok]-cb[ok]))
  },by=.(scenario,model)]
  influence[,group:=ifelse(model=='TWFE','TWFE','Non-TWFE models')]
  summarize<-function(v,overall=FALSE) {
    if(overall)v[,scenario:='Overall']
    cells<-unique(v[,.(scenario,group,model,absolute_ratio,relative_ratio,coverage_change)])
    if(overall)cells<-unique(v[,.(source_scenario,scenario,group,model,absolute_ratio,relative_ratio,coverage_change)])
    point<-cells[,.(absolute_ratio=mean(absolute_ratio),relative_ratio=mean(relative_ratio),coverage_change=mean(coverage_change)),by=.(scenario,group)]
    clustered<-v[,.(absolute_if=mean(absolute_if),relative_if=mean(relative_if),coverage_if=mean(coverage_if)),by=.(scenario,group,sim)]
    se<-clustered[,.(absolute_ratio_mc_se=sd(absolute_if)/sqrt(.N),relative_ratio_mc_se=sd(relative_if)/sqrt(.N),coverage_change_mc_se=sd(coverage_if)/sqrt(.N),n_simulations=.N),by=.(scenario,group)]
    merge(point,se,by=c('scenario','group'))
  }
  influence[,source_scenario:=scenario]
  rbindlist(list(summarize(copy(influence)),summarize(copy(influence),TRUE)))
}

run_half_estimate <- function(root, half='B', report_path=file.path(root,'analysis/S6_screening/s3_03_half_B.qmd'), with_covariates=FALSE, authorize_B=FALSE) {
  library(data.table)
  stopifnot(half %in% c('A','B'),is.logical(authorize_B),length(authorize_B)==1L,!is.na(authorize_B))
  out<-if(with_covariates)file.path(root,'_output/heldout/with_covariates/confirmation') else file.path(root,'_output/heldout/confirmation');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  prefix<-paste0('half_',half,'_')
  log<-data.table(path=character(),access=character(),timestamp_utc=character())
  record<-function(path,access='read') {
    log<<-rbind(log,data.table(path=normalizePath(path,mustWork=FALSE),access=access,timestamp_utc=format(Sys.time(),tz='UTC',usetz=TRUE)))
    fwrite(log,file.path(out,paste0(prefix,'input_paths.csv')))
  }
  read_rds<-function(file,...) {
    if(half=='A' && (grepl('half_B|/post_B/|/estimate_B/',file)||(!with_covariates && grepl('/data/simulation/',file))))stop('A dry run prohibits B or complete-panel fit access: ',file)
    record(file);base::readRDS(file,...)
  }
  read_csv<-function(file){record(file);fread(file)}
  record(file.path(root,'R/screening/heldout_estimate.R'),'source')
  record(report_path,'report source')
  if(file.exists(file.path(dirname(report_path),'_quarto.yml')))record(file.path(dirname(report_path),'_quarto.yml'),'isolated render configuration')
  record(file.path(root,'_quarto.yml'),'render configuration')
  env<-new.env(parent=globalenv());env$readRDS<-read_rds
  env$file.info<-function(...){paths<-unlist(list(...));for(p in paths)record(p,'metadata only');base::file.info(...)}
  for(f in c('heldout_analysis.R','../backdate.r','heldout_exploration.R','heldout_aggregate.R','heldout_twfe_population.R')) {
    path<-file.path(root,'R','screening',f);record(path,'source');sys.source(path,envir=env)
  }
  decision_path<-file.path(root,'analysis/S6_screening/s3_decision.yml');record(decision_path,'yaml')
  decision<-env$read_decision(root)
  record(file.path(root,'.git'),'git metadata')
  tracked<-length(system2('git',c('-C',shQuote(root),'ls-files','--',shQuote('analysis/S6_screening/s3_decision.yml')),stdout=TRUE))>0
  commit<-NA_character_
  if(tracked) {
    candidates<-system2('git',c('-C',shQuote(root),'log','-1','--format=%H','--','analysis/S6_screening/s3_decision.yml'),stdout=TRUE)
    if(length(candidates)) {
      committed_path<-tempfile('decision-provenance-',tmpdir=out)
      status<-system2('git',c('-C',shQuote(root),'show',paste0(candidates[1],':analysis/S6_screening/s3_decision.yml')),stdout=committed_path)
      if(status==0 && identical(unname(tools::md5sum(committed_path)),decision$md5))commit<-candidates[1]
    }
  }
  assignment_path<-file.path(root,'data/derived/heldout/county_assignment.rds');assignment<-read_rds(assignment_path)
  assignment_md5<-unname(tools::md5sum(assignment_path))
  if(half=='B') {
    if(!authorize_B && !file.exists(file.path(out,'run_B_authorized')))stop("B remains unopened: explicitly render the B estimate or set authorize_B=TRUE; use half='A' for a dry run.")
    if(!tracked || !length(commit) || is.na(commit))stop('The evaluated decision must be committed before B is opened.')
    blob_path<-tempfile('decision-committed-',tmpdir=out)
    status<-system2('git',c('-C',shQuote(root),'show',paste0(commit,':analysis/S6_screening/s3_decision.yml')),stdout=blob_path)
    if(status!=0||!identical(unname(tools::md5sum(blob_path)),decision$md5))stop('Decision differs from its committed version; B remains unopened.')
    phase1<-read_rds(file.path(out,'half_A_provenance.rds'))
    stopifnot(identical(phase1$assignment_md5,assignment_md5),identical(phase1$decision_md5,decision$md5),half=='B')
  }
  preserved<-c(list.files(file.path(root,'_output/analysis/output/3.1a'),'[.]csv$',recursive=TRUE,full.names=TRUE),list.files(file.path(root,'_output/analysis'),'^3[.]1a.*[.]html$',full.names=TRUE))
  before<-tools::md5sum(preserved)
  input<-if(with_covariates)file.path(root,'data/derived/heldout/with_covariates',paste0('half_',half,'.rds')) else file.path(root,'data/derived/heldout',if(half=='A')'half_A_available.rds' else 'half_B.rds')
  if(with_covariates && !file.exists(input)) {
    path<-file.path(root,'R/screening/heldout_covariate_input.R');record(path,'source');sys.source(path,envir=env)
    env$prepare_covariate_selected_half(root,half,input,read_input=read_rds)
  }
  obj<-read_rds(input);env$aggregate_guard(obj,root,half)
  source(file.path(root,'R/screening/s3_models.R'),local=TRUE);models<-s3_models;scenarios<-c('MAR','Latent level','Latent trend')
  rules<-data.table(diagnostic=c(decision$diagnostic,vapply(decision$sensitivity_rules,`[[`,character(1),'diagnostic')),
    cutoff=c(decision$cutoff,vapply(decision$sensitivity_rules,`[[`,numeric(1),'cutoff')),rule_type=c('absolute','p_value','top_percent'),rule_id=1:3)
  rules[diagnostic=='pre_rmse',rule_type:='top_percent'];rules[diagnostic=='p_two',rule_type:='p_value']
  baseline<-data.table(diagnostic='unscreened',cutoff=NA_real_,rule_type='baseline',missing_handling='none')
  keys<-rbindlist(list(baseline,rules[,.(diagnostic,cutoff,rule_type)][,missing_handling:='keep']))
  fits<-list();rows<-list()
  cache_dir<-file.path(out,paste0('estimate_',half));if(half=='B' || with_covariates)dir.create(cache_dir,showWarnings=FALSE)
  for(sc in scenarios)for(m in models) {
    p<-if(half=='A' && !with_covariates)env$twfe_common_prediction_cache(root,obj,m,sc) else env$aggregate_predictive_A(root,m,sc,obj,cache_dir,expected_half=half,input_path=input,main_dir=file.path(root,'data/simulation',if(with_covariates)'wcov' else 'wocov'),file_prefix=if(with_covariates)'wcov_imputation' else 'new_imputation')
    v<-merge(obj$rows[scenario==sc & model==m],p[,.(sim,id,k_check,y_check,mu_check,draws,legacy_y,legacy_mu)],by=c('sim','id'),all.x=TRUE)
    v[is.na(k_check),k_check:=0L]
    if(half=='B' || with_covariates) {
      v[,`:=`(k=k_check,sum_y=ifelse(k_check>=1,y_check,NA_real_),sum_mu=ifelse(k_check>=1,mu_check,NA_real_),
        err_raw=ifelse(k_check>=1,abs(y_check-mu_check),NA_real_),err_z=ifelse(k_check>=1,abs(y_check-mu_check)/sqrt(mu_check),NA_real_),
        k5_raw=abs(legacy_y-legacy_mu),k5_z=abs(legacy_y-legacy_mu)/sqrt(legacy_mu))]
      v[,c('pi_lower','pi_upper','coverage'):=list(NA_real_,NA_real_,NA_real_)]
      for(j in which(v$k>=1)) {
        a<-v$draws[[j]];if(all(is.finite(a))) {
          b<-quantile(a,c(.025,.975),names=FALSE,type=7)
          set(v,j,'pi_lower',b[1]);set(v,j,'pi_upper',b[2]);set(v,j,'coverage',as.numeric(v$sum_y[j]>=b[1]&&v$sum_y[j]<=b[2]))
        }
      }
    }
    stopifnot(all(v$k==v$k_check),all(abs(v$sum_y-v$y_check)<1e-10,na.rm=TRUE),all(abs(v$sum_mu-v$mu_check)<1e-10,na.rm=TRUE))
    fits[[paste(sc,m)]]<-v;rows[[length(rows)+1L]]<-v
  }
  z<-env$derive_exploration_A(rbindlist(rows,fill=TRUE),expected_half=half)$rows
  if(half=='A' && !with_covariates) {
    meta_e<-read_csv(file.path(root,'_output/analysis/output/3.1a/twfe_common/twfe_common_membership.csv'))
    meta_c<-read_csv(file.path(root,'_output/analysis/output/3.1a/twfe_common/twfe_common_complement_membership.csv'))
    noncore_ids<-unique(c(meta_e[noncore==TRUE,id],meta_c$id))
    z[,noncore:=id %in% noncore_ids]
  } else {
    meta<-env$aggregate_county_metadata(root,unique(z$id));z<-merge(z,meta,by='id',all.x=TRUE)
  }
  ledgers<-retention<-assoc<-strata<-contrib<-membership<-checks<-list()
  for(sc in scenarios)for(simulation in sort(unique(z$sim))) {
    full<-setNames(lapply(models,function(m)z[scenario==sc & sim==simulation & model==m][order(id)]),models)
    screens<-lapply(full,function(v)env$keep_matrix_A(v,rules,expected_half=half)$keep[,1:3,drop=FALSE])
    idsE<-full[['TWFE']][k>=1 & is.finite(sum_y-sum_mu),id]
    for(m in models) {
      v<-full[[m]]
      a<-env$retention_group_A(v,rules,targets=c('err_z','err_raw'),random_baseline=FALSE,expected_half=half)[missing_handling=='keep']
      a[,`:=`(scenario=sc,model=m,sim=simulation)];retention[[length(retention)+1L]]<-a
    }
    for(sg in c('noncore','all_half')) {
      eids<-if(sg=='noncore')full[['TWFE']][id %in% idsE & noncore==TRUE,id] else idsE
      e<-lapply(full,function(v)v[match(eids,id)])
      stopifnot(all(vapply(e,function(v)!anyNA(v$id)&&all(v$k>=1)&&all(is.finite(v$sum_y-v$sum_mu)),logical(1))))
      stopifnot(all(vapply(e,function(v)identical(v$sum_y,e[[1]]$sum_y)&&identical(v$k,e[[1]]$k),logical(1))))
      ownE<-lapply(models,function(m)screens[[m]][match(eids,full[[m]]$id),,drop=FALSE]);names(ownE)<-models
      joint<-Reduce(`&`,ownE);nflags<-Reduce(`+`,lapply(ownE,function(a)!a))
      for(m in models) {
        v<-full[[m]];eligible<-v$k>=1 & is.finite(v$sum_y-v$sum_mu) & if(sg=='noncore')v$noncore else TRUE
        u<-v[eligible];mask<-screens[[m]][eligible,,drop=FALSE]
        for(pop in c('own','E','common')) {
          counties<-if(pop=='own')u else e[[m]]
          keep<-cbind(TRUE,switch(pop,own=mask,E=ownE[[m]],common=joint))
          a<-env$twfe_common_group(counties,keys,keep,expected_half=half)
          a[,`:=`(subgroup=paste(sg,pop,sep='/'),scenario=sc,model=m,sim=simulation)]
          ledgers[[length(ledgers)+1L]]<-a
        }
        c<-copy(keys[-1]);nex<-colSums(!joint)
        c[,`:=`(subgroup=sg,scenario=sc,sim=simulation,model=m,n_E=length(eids),n_excluded=nex,n_unique=colSums(!ownE[[m]] & nflags==1),n_shared=colSums(nflags>=2))]
        contrib[[length(contrib)+1L]]<-c
        stopifnot(all(!joint|ownE[[m]]))
      }
      checks[[length(checks)+1L]]<-data.table(subgroup=sg,scenario=sc,sim=simulation,identical_ids=all(vapply(e,function(v)identical(v$id,eids),logical(1))),identical_observations=TRUE,subset_of_every_model=TRUE)
      for(j in 1:3) {
        take<-eids[joint[,j]];if(length(take))membership[[length(membership)+1L]]<-data.table(subgroup=sg,scenario=sc,sim=simulation,diagnostic=rules$diagnostic[j],cutoff=rules$cutoff[j],id=take,half=half)
      }
    }
  }
  for(t in c('err_z','err_raw')) {
    v<-z[k>=1,.(scenario,model,sim,k,score=abs(ho_standardized),error=get(t))]
    a<-env$association_stats(v,c('scenario','model','sim','k'));a[,target:=t];strata[[length(strata)+1L]]<-a
    b<-a[,.(n=sum(n[usable]),spearman=if(any(is.finite(spearman)))weighted.mean(spearman,n,na.rm=TRUE) else NA_real_,auc=if(any(is.finite(auc)))weighted.mean(auc,n,na.rm=TRUE) else NA_real_),by=.(scenario,model,sim)]
    b[,target:=t];assoc[[length(assoc)+1L]]<-b
  }
  d<-rbindlist(ledgers);summary<-env$aggregate_summary(d)
  split_population<-function(v){v[,population:=sub('.*/','',subgroup)];v[,subgroup:=sub('/.*','',subgroup)];v}
  d<-split_population(d);summary<-split_population(summary)
  summary[,rule_label:=ifelse(diagnostic=='unscreened','Unscreened',ifelse(diagnostic==decision$diagnostic,'Recorded rule','Sensitivity'))]
  summary[,`:=`(rmse_rank=frank(rmse,ties.method='min'),coverage_rank=frank(-coverage,ties.method='min')),by=.(subgroup,population,scenario,diagnostic,cutoff)]
  r<-rbindlist(retention);a<-rbindlist(assoc);c<-rbindlist(contrib)
  gain<-estimate_gains(d)
  original_dir<-file.path(root,'_output/analysis/output/3.1a')
  old_own<-read_csv(file.path(original_dir,'aggregate_per_simulation.csv'));old_own[,`:=`(population='own',subgroup=ifelse(subgroup=='all_A','all_half',subgroup))]
  old_E<-read_csv(file.path(original_dir,'twfe_common/twfe_common_per_simulation.csv'));old_E[,`:=`(population='E',subgroup=ifelse(subgroup=='all_A','all_half',subgroup))]
  old_common<-read_csv(file.path(original_dir,'common_screen/common_screen_per_simulation.csv'));old_common[,population:='common']
  old_d<-rbindlist(list(old_own,old_E,old_common),fill=TRUE)
  gain_no_cov_A<-estimate_gains(old_d)
  gain_A<-if(with_covariates) {
    if(half=='A')estimate_gains(d) else estimate_gains(read_csv(file.path(out,'half_A_aggregate_per_simulation.csv')))
  } else gain_no_cov_A
  gain<-merge(gain,gain_A,by=c('scenario','group'),suffixes=c('_estimate','_A_reference'))
  for(v in c('absolute_ratio','relative_ratio','coverage_change'))gain[,paste0(v,'_difference'):=get(paste0(v,'_estimate'))-get(paste0(v,'_A_reference'))]
  if(with_covariates) {
    ref<-copy(gain_no_cov_A);setnames(ref,setdiff(names(ref),c('scenario','group')),paste0(setdiff(names(ref),c('scenario','group')),'_no_cov_A'))
    gain<-merge(gain,ref,by=c('scenario','group'))
  }
  reproduce<-data.table()
  if(half=='A' && !with_covariates) {
    join<-c('subgroup','population','scenario','model','sim','diagnostic','cutoff','rule_type','missing_handling')
    selected<-d[diagnostic=='unscreened' | diagnostic=='ho_standardized']
    match<-merge(selected,old_d,by=join,suffixes=c('_new','_A'))
    # The existing 10.2 exports noncore only; compare every available reference row.
    expected<-selected[!(population=='common' & subgroup=='all_half')]
    stopifnot(nrow(match)==nrow(expected))
    quantities<-c('error','error_z','observed_aggregate','predicted_aggregate','coverage','width','lower','upper','n_retained','excluded_share')
    reproduce<-rbindlist(lapply(quantities,function(v)match[,.(quantity=v,max_abs_difference=max(abs(get(paste0(v,'_new'))-get(paste0(v,'_A'))),na.rm=TRUE),missing_agrees=identical(is.na(get(paste0(v,'_new'))),is.na(get(paste0(v,'_A'))))),by=.(population)]))
    refs<-list(own=read_csv(file.path(original_dir,'aggregate_summary.csv')),E=read_csv(file.path(original_dir,'twfe_common/twfe_common_summary.csv')),common=read_csv(file.path(original_dir,'common_screen/common_screen_summary.csv')))
    for(pop in names(refs)) {
      old<-refs[[pop]];old[,subgroup:=ifelse(subgroup=='all_A','all_half',subgroup)]
      j<-merge(summary[population==pop & diagnostic %in% c('unscreened','ho_standardized')],old,by=setdiff(join,c('population','sim')),suffixes=c('_new','_A'))
      for(v in c('rmse','rmse_mc_se','rmse_ratio','rmse_ratio_mc_se','relative_rmse','relative_rmse_mc_se','coverage','coverage_mc_se','width','width_mc_se','n_retained','excluded_share','rmse_rank')) {
        if(!paste0(v,'_A') %in% names(j))next
        reproduce<-rbind(reproduce,data.table(population=pop,quantity=paste0('summary_',v),max_abs_difference=max(abs(j[[paste0(v,'_new')]]-j[[paste0(v,'_A')]]),na.rm=TRUE),missing_agrees=identical(is.na(j[[paste0(v,'_new')]]),is.na(j[[paste0(v,'_A')]]))))
      }
    }
    for(v in c('absolute_ratio','relative_ratio','coverage_change','absolute_ratio_mc_se','relative_ratio_mc_se','coverage_change_mc_se')) {
      reproduce<-rbind(reproduce,data.table(population='estimation_plan',quantity=v,max_abs_difference=max(abs(gain[[paste0(v,'_estimate')]]-gain[[paste0(v,'_A_reference')]])),missing_agrees=TRUE))
    }
    rank_ref<-read_csv(file.path(original_dir,'common_screen/common_screen_comparison.csv'))[subgroup=='noncore' & (diagnostic=='unscreened' | diagnostic=='ho_standardized' & cutoff==2 & missing_handling=='keep')]
    rank_ref[,population:=ifelse(sample=='Common screen','common','E')]
    ranked<-merge(summary[subgroup=='noncore' & population %in% c('E','common') & (diagnostic=='unscreened' | diagnostic=='ho_standardized')],rank_ref,by=c('population','scenario','model','diagnostic','cutoff'),suffixes=c('_new','_A'))
    reproduce<-rbind(reproduce,data.table(population='rankings',quantity='RMSE ranks',max_abs_difference=max(abs(ranked$rmse_rank_new-ranked$rmse_rank_A)),missing_agrees=TRUE))
    reproduce[,within_numerical_precision:=missing_agrees & max_abs_difference<1e-10]
    if(!all(reproduce$within_numerical_precision))stop('A numerical reproduction differs: ',paste(reproduce[within_numerical_precision==FALSE,quantity],collapse=', '))
  }
  twfe<-merge(summary[model=='TWFE' & population=='own'],summary[model=='TWFE' & population=='E'],by=c('subgroup','scenario','model','diagnostic','cutoff','missing_handling'),suffixes=c('_own','_E'))
  for(v in c('rmse','rmse_mc_se','relative_rmse','relative_rmse_mc_se','rmse_ratio','rmse_ratio_mc_se','coverage','coverage_mc_se','width','width_mc_se','n_retained','excluded_share'))stopifnot(all(abs(twfe[[paste0(v,'_own')]]-twfe[[paste0(v,'_E')]])<1e-10,na.rm=TRUE))
  csummary<-c[,.(n_E=mean(n_E),n_excluded=mean(n_excluded),excluded_share=mean(n_excluded/n_E),n_unique=mean(n_unique),unique_share=sum(n_unique)/sum(n_excluded),shared_share=sum(n_shared)/sum(n_excluded)),by=.(subgroup,scenario,model,diagnostic,cutoff)]
  account<-c[,.(accounted_for=sum(n_unique)+n_shared[1]==n_excluded[1]),by=.(subgroup,scenario,sim,diagnostic,cutoff)]
  stopifnot(all(account$accounted_for))
  availability<-z[,.(rows=.N,counties=uniqueN(id),simulations=uniqueN(sim),available_target_rows=sum(k>=1),missing_recorded_diagnostic=sum(!is.finite(ho_standardized)),minimum_k=min(k),maximum_k=max(k)),by=.(scenario,model)]
  exports<-list(availability=availability,replication_counts=z[,.(counties=.N,available_counties=sum(k>=1)),by=.(scenario,model,sim)],aggregate_per_simulation=d,aggregate_summary=summary,
    associations_by_k=rbindlist(strata),associations_per_simulation=a,associations_summary=env$mc_summary(a,c('spearman','auc'),c('scenario','model','target')),
    retention_per_simulation=r,retention_summary=env$mc_summary(r,c('excluded_share','reduction','retained_coverage','n_retained'),c('scenario','model','diagnostic','cutoff','target')),
    estimates=gain,unique_flags=csummary,common_membership=rbindlist(membership),population_checks=rbindlist(checks),flag_accounting=account,twfe_own_E_checks=twfe,
    rankings=summary[subgroup=='noncore' & population %in% c('E','common') & scenario=='Latent level',.(population,model,diagnostic,cutoff,rule_label,rmse,relative_rmse,coverage,width,rmse_rank,coverage_rank)],reproduction_check=reproduce,input_manifest=as.data.table(obj$manifest))
  for(nm in names(exports))fwrite(exports[[nm]],file.path(out,paste0(prefix,nm,'.csv')))
  provenance<-list(half=half,with_covariates=with_covariates,A_reference=if(with_covariates)'covariate A dry run (original no-covariate A additionally reported)' else 'original no-covariate A',decision_md5=decision$md5,decision_commit=commit,repository_head=system2('git',c('-C',shQuote(root),'rev-parse','HEAD'),stdout=TRUE),assignment_md5=assignment_md5,
    framing=if(half=='A')'dry run on A — not a B estimate' else 'half-B estimate',half_a_reference=decision$half_a_reference,estimation_plan=decision$estimation_plan)
  saveRDS(provenance,file.path(out,paste0(prefix,'provenance.rds')));yaml::write_yaml(provenance,file.path(out,paste0(prefix,'provenance.yml')))
  file.copy(decision_path,file.path(out,paste0(prefix,'evaluated_decision.yml')),overwrite=TRUE)
  writeLines(decision$md5,file.path(out,paste0(prefix,'evaluated_decision.md5')))
  env$aggregate_guard(obj,root,half);stopifnot(identical(before,tools::md5sum(names(before))),identical(decision$md5,unname(tools::md5sum(decision_path))))
  list(half=half,decision=decision,provenance=provenance,tables=exports,out=out,prefix=prefix)
}

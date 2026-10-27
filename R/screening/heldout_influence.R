# Descriptive/post-hoc influence of existing cached county errors; flags are frozen.
influence_aggregate_summary <- function(E,N,Y,M,baseN,keys,sim,group='full') {
  k<-ncol(E);ns<-length(sim)
  d<-copy(keys[rep(seq_len(k),each=ns)])
  d[,`:=`(subgroup=group,scenario='current',model='current',sim=rep(sim,k),
    error=as.vector(E),error_z=as.vector(E/sqrt(M/N^2)),observed_aggregate=as.vector(Y/N),
    coverage=NA_real_,width=NA_real_,n_retained=as.vector(N),excluded_share=as.vector(1-N/baseN))]
  d[which(as.vector(N)==0),c('error','error_z','observed_aggregate'):=NA_real_]
  aggregate_summary(d)
}

run_heldout_influence <- function(root) {
  library(data.table)
  out<-file.path(root,'_output/heldout/influence');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  log<-data.table(path=character(),access=character(),timestamp_utc=character())
  record<-function(path,access='read') {
    log<<-rbind(log,data.table(path=normalizePath(path,mustWork=FALSE),access=access,timestamp_utc=format(Sys.time(),tz='UTC',usetz=TRUE)))
    fwrite(log,file.path(out,'influence_input_paths.csv'))
  }
  read_rds<-function(path) {if(grepl('/data/simulation/',path))stop('Influence cannot open raw fits or rederive targets: ',path);record(path);base::readRDS(path)}
  readRDS<-read_rds
  read_csv<-function(path){record(path);fread(path)}
  for(f in c('heldout_analysis.R','heldout_aggregate.R','heldout_exploration.R')) {
    p<-file.path(root,'R','screening',f);record(p,'source');source(p,local=TRUE)
  }
  # Bind the unchanged shared aggregate helper for the pure summary adapter.
  environment(influence_aggregate_summary)<-environment()
  original<-c(list.files(file.path(root,'_output/analysis/output/3.1a'),'[.]csv$',recursive=TRUE,full.names=TRUE),
    list.files(file.path(root,'_output/heldout/confirmation'),'[.](csv|html|yml)$',recursive=TRUE,full.names=TRUE),
    list.files(file.path(root,'_output/heldout/with_covariates/confirmation'),'[.](csv|html|yml)$',recursive=TRUE,full.names=TRUE),
    list.files(file.path(root,'_output/analysis'),'^3[.]1a.*[.]html$',full.names=TRUE),
    file.path(root,'analysis/S6_screening/s3_decision.yml'),file.path(root,'data/derived/heldout/county_assignment.rds'))
  before<-tools::md5sum(original);saveRDS(before,file.path(out,'influence_original_hashes.rds'))
  assignment_path<-file.path(root,'data/derived/heldout/county_assignment.rds');assignment<-as.data.table(read_rds(assignment_path))
  panel<-as.data.table(read_rds(file.path(root,Sys.getenv('SAI_PANEL','data/notpublic/compact_dt.rds'))))
  meta<-panel[,.(noncore=all(as.character(urb)=='6.NonCore'),mean_outcome_age_population=mean(n[is.finite(n)])),by=.(id=i)]
  stopifnot(!anyDuplicated(meta$id))
  source(file.path(root,'R/screening/s3_models.R'),local=TRUE);models<-s3_models;scenarios<-c('MAR','Latent level','Latent trend')
  codes<-s3_codes
  rules<-data.table(diagnostic=c('ho_standardized','p_two','pre_rmse'),cutoff=c(2,.10,10),rule_type=c('absolute','p_value','top_percent'),rule_id=1:3)
  keys<-rbindlist(list(data.table(diagnostic='unscreened',cutoff=NA_real_,rule_type='baseline',missing_handling='none'),rules[,.(diagnostic,cutoff,rule_type)][,missing_handling:='keep']))
  all_counties<-concentrations<-county_contributions<-loos<-profiles<-top_deletions<-reproductions<-captures<-baseline_summaries<-list()
  group_index<-0L
  for(branch in c('No covariates','With covariates'))for(half in c('A','B')) {
    cov<-branch=='With covariates'
    inp<-if(cov)file.path(root,'data/derived/heldout/with_covariates',paste0('half_',half,'.rds')) else file.path(root,'data/derived/heldout',if(half=='A')'half_A_available.rds' else 'half_B.rds')
    obj<-read_rds(inp);aggregate_guard(obj,root,half)
    # Derive ONLY the same diagnostic scores; no post-period target is derived here.
    z<-derive_heldout(obj$rows,diagnostics_only=TRUE)
    source_out<-if(cov)file.path(root,'_output/heldout/with_covariates/confirmation') else file.path(root,'_output/heldout/confirmation')
    reference_d<-read_csv(file.path(source_out,paste0('half_',half,'_aggregate_per_simulation.csv')))[subgroup=='noncore' & population=='own']
    reference_s<-read_csv(file.path(source_out,paste0('half_',half,'_aggregate_summary.csv')))[subgroup=='noncore' & population=='own']
    for(sc in scenarios)for(m in models) {
      group_index<-group_index+1L
      input_v<-z[scenario==sc & model==m][order(sim,id)]
      cache<-if(!cov && half=='A')file.path(root,'.quarto/heldout-aggregate',sprintf('%s_%s.rds',codes[[m]],match(sc,scenarios))) else file.path(source_out,paste0('estimate_',half),sprintf('%s_%s.rds',codes[[m]],match(sc,scenarios)))
      saved<-read_rds(cache)
      main<-s3_main_path(root,m,match(sc,scenarios),cov)
      record(main,'metadata only');record(inp,'metadata only')
      expected<-list(assignment=obj$assignment_md5,predictive_version=1L,main_path=main,bytes=file.info(main)$size,mtime=as.numeric(file.info(main)$mtime),input=as.numeric(file.info(inp)$mtime))
      if(!identical(saved$signature,expected))stop('Stale saved county targets; no regeneration permitted: ',cache)
      cached<-as.data.table(saved$rows)[,.(sim,id,k=k_check,sum_y=y_check,sum_mu=mu_check)]
      stopifnot(all(saved$rows$half==half),all(saved$rows$id %in% assignment$id[assignment$half==half]),!anyDuplicated(cached[,.(sim,id)]))
      rm(saved)
      v<-merge(input_v,cached,by=c('sim','id'),all.x=TRUE,suffixes=c('_input',''))
      if('k_input' %in% names(v))stopifnot(all(v$k_input==v$k,na.rm=TRUE))
      v<-merge(v,meta,by='id',all.x=TRUE);stopifnot(!anyNA(v$noncore),all(v$half==half))
      # Evaluate the original full-half percentile reference, including k=0 rows.
      v[,`:=`(keep_z=FALSE,keep_p=FALSE,keep_pre=FALSE)]
      for(r in unique(v$sim)) {
        at<-which(v$sim==r);screen<-keep_matrix_A(v[at],rules,expected_half=half)$keep[,1:3,drop=FALSE]
        for(j in 1:3)set(v,at,c('keep_z','keep_p','keep_pre')[j],screen[,j])
      }
      sims<-sort(unique(v$sim));ns<-length(sims)
      pseudo<-v[noncore==TRUE,.(n_pseudo_treated=uniqueN(sim[is.finite(x)])),by=id]
      v<-v[noncore==TRUE & k>=1 & is.finite(sum_y-sum_mu)]
      v[,e:=sum_y-sum_mu];stopifnot(!anyDuplicated(v[,.(sim,id)]))
      ids<-sort(unique(v$id));nc<-length(ids)
      e<-y<-mu<-matrix(0,ns,nc);present<-matrix(FALSE,ns,nc)
      idx<-cbind(match(v$sim,sims),match(v$id,ids));e[idx]<-v$e;y[idx]<-v$sum_y;mu[idx]<-v$sum_mu;present[idx]<-TRUE
      masks<-list(present,present,present,present)
      for(j in 1:3) {masks[[j+1L]][]<-FALSE;masks[[j+1L]][idx]<-v[[c('keep_z','keep_p','keep_pre')[j]]]}
      N<-vapply(masks,rowSums,numeric(ns));S<-vapply(masks,function(a)rowSums(e*a),numeric(ns))
      Y<-vapply(masks,function(a)rowSums(y*a),numeric(ns));M<-vapply(masks,function(a)rowSums(mu*a),numeric(ns))
      E<-S/N;baseN<-rowSums(present)
      full<-influence_aggregate_summary(E,N,Y,M,baseN,keys,sims)
      full[,`:=`(branch=branch,half=half,scenario=sc,model=m)]
      baseline_summaries[[length(baseline_summaries)+1L]]<-full
      rd<-reference_d[scenario==sc & model==m]
      rs<-reference_s[scenario==sc & model==m]
      point<-copy(keys[rep(1:4,each=ns)]);point[,`:=`(sim=rep(sims,4),error=as.vector(E),n_retained=as.vector(N))]
      join<-c('sim','diagnostic','cutoff','rule_type','missing_handling')
      check<-merge(point,rd,by=join,suffixes=c('_new','_saved'));stopifnot(nrow(check)==nrow(point),all(check$n_retained_new==check$n_retained_saved))
      repcheck<-data.table(quantity=c('per_simulation_error','retained_count'),max_abs_difference=c(max(abs(check$error_new-check$error_saved),na.rm=TRUE),max(abs(check$n_retained_new-check$n_retained_saved))))
      j<-merge(full,rs,by=c('diagnostic','cutoff','rule_type','missing_handling'),suffixes=c('_new','_saved'))
      stopifnot(nrow(j)==4)
      for(q in c('rmse','rmse_mc_se','rmse_ratio','rmse_ratio_mc_se'))repcheck<-rbind(repcheck,data.table(quantity=q,max_abs_difference=max(abs(j[[paste0(q,'_new')]]-j[[paste0(q,'_saved')]]),na.rm=TRUE)))
      repcheck[,`:=`(branch=branch,half=half,scenario=sc,model=m,within_precision=max_abs_difference<1e-10)];stopifnot(all(repcheck$within_precision))
      reproductions[[length(reproductions)+1L]]<-repcheck
      group_contrib<-list()
      for(rule in 1:4) {
        valid<-is.finite(E[,rule]);weight<-matrix(0,ns,nc)
        weight[valid,]<-(e*masks[[rule]])[valid,,drop=FALSE]*(E[valid,rule]/N[valid,rule])
        C<-colSums(weight)/sum(valid);mse<-mean(E[valid,rule]^2)
        stopifnot(abs(sum(C)-mse)<1e-10)
        cc<-data.table(id=ids,C=C,MSE=mse,share=C/mse,retained_simulations=colSums(masks[[rule]]),n_simulations=ns,n_usable_simulations=sum(valid))
        setorder(cc,-C,id);cc[,`:=`(contribution_rank=seq_len(.N),cumulative_share=cumsum(share))]
        cc[,`:=`(branch=branch,half=half,scenario=sc,model=m,diagnostic=keys$diagnostic[rule],cutoff=keys$cutoff[rule])]
        county_contributions[[length(county_contributions)+1L]]<-cc;group_contrib[[rule]]<-cc
        concentrations[[length(concentrations)+1L]]<-data.table(branch=branch,half=half,scenario=sc,model=m,diagnostic=keys$diagnostic[rule],cutoff=keys$cutoff[rule],MSE=mse,RMSE=sqrt(mse),
          top1_share=sum(head(cc$share,1)),top3_share=sum(head(cc$share,3)),top5_share=sum(head(cc$share,5)),top10_share=sum(head(cc$share,10)),
          counties_to_50=which(cc$cumulative_share>=.5)[1],n_counties=nc,n_negative_contributions=sum(C<0),negative_share=sum(C[C<0])/mse,sum_contributions=sum(C),decomposition_difference=sum(C)-mse,n_usable_simulations=sum(valid))
      }
      loo_group<-list()
      for(jj in seq_along(ids)) {
        Nj<-N;Sj<-S;Yj<-Y;Mj<-M
        for(rule in 1:4) {take<-masks[[rule]][,jj];Nj[,rule]<-N[,rule]-take;Sj[,rule]<-S[,rule]-e[,jj]*take;Yj[,rule]<-Y[,rule]-y[,jj]*take;Mj[,rule]<-M[,rule]-mu[,jj]*take}
        a<-influence_aggregate_summary(Sj/Nj,Nj,Yj,Mj,baseN-present[,jj],keys,sims)
        a<-merge(a,full[,.(diagnostic,cutoff,rule_type,missing_handling,full_rmse=rmse,full_ratio=rmse_ratio)],by=c('diagnostic','cutoff','rule_type','missing_handling'))
        a[,`:=`(id=ids[jj],branch=branch,half=half,scenario=sc,model=m,rmse_change=rmse-full_rmse,ratio_change=rmse_ratio-full_ratio)]
        loo_group[[jj]]<-a
      }
      lo<-rbindlist(loo_group)
      influence<-lo[diagnostic=='unscreened',.(id,unscreened_rmse=full_rmse,loo_unscreened_rmse=rmse,unscreened_rmse_decrease=full_rmse-rmse)]
      setorder(influence,-unscreened_rmse_decrease,id);influence[,unscreened_influence_rank:=seq_len(.N)]
      lo<-merge(lo,influence,by='id');lo[,ratio_influence_rank:=frank(-abs(ratio_change),ties.method='min'),by=.(diagnostic,cutoff)]
      loos[[length(loos)+1L]]<-lo
      topids<-head(influence$id,5)
      profile<-v[id %in% topids,.(n_eligible_simulations=.N,eligible_simulation_share=.N/ns,mean_signed_error=mean(e),mean_absolute_error=mean(abs(e)),
        flag_z_conditional=mean(!keep_z),flag_p_conditional=mean(!keep_p),flag_pre_conditional=mean(!keep_pre),
        flag_z_all_simulations=sum(!keep_z)/ns,flag_p_all_simulations=sum(!keep_p)/ns,flag_pre_all_simulations=sum(!keep_pre)/ns),by=id]
      profile<-merge(profile,pseudo,by='id');profile[,pseudo_treated_share:=n_pseudo_treated/ns]
      profile<-merge(profile,influence,by='id');profile<-merge(profile,assignment[,.(id,stratum,mean_pre_expected_count)],by='id');profile<-merge(profile,meta[,.(id,mean_outcome_age_population)],by='id')
      profile<-merge(profile,group_contrib[[1]][,.(id,unscreened_MSE_share=share,contribution_rank)],by='id')
      profile[,`:=`(branch=branch,half=half,scenario=sc,model=m)];profiles[[length(profiles)+1L]]<-profile
      base_weight<-(e*present)*(E[,1]/N[,1]);positive<-pmax(base_weight,0)
      for(rule in 2:4) {
        flagged<-present & !masks[[rule]];top<-ids %in% topids
        captures[[length(captures)+1L]]<-data.table(branch=branch,half=half,scenario=sc,model=m,diagnostic=keys$diagnostic[rule],cutoff=keys$cutoff[rule],
          top5_conditional_flag_share=sum(flagged[,top,drop=FALSE])/sum(present[,top,drop=FALSE]),
          top5_positive_contribution_caught=sum(positive[,top,drop=FALSE]*flagged[,top,drop=FALSE])/sum(positive[,top,drop=FALSE]),
          all_positive_contribution_caught=sum(positive*flagged)/sum(positive))
      }
      for(k in c(0L,1L,3L)) {
        remove<-if(k)head(influence$id,k) else integer()
        Nk<-N;Sk<-S;Yk<-Y;Mk<-M
        for(rule in 1:4)if(k) {take<-masks[[rule]][,ids %in% remove,drop=FALSE];Nk[,rule]<-N[,rule]-rowSums(take);Sk[,rule]<-S[,rule]-rowSums(e[,ids %in% remove,drop=FALSE]*take);Yk[,rule]<-Y[,rule]-rowSums(y[,ids %in% remove,drop=FALSE]*take);Mk[,rule]<-M[,rule]-rowSums(mu[,ids %in% remove,drop=FALSE]*take)}
        a<-influence_aggregate_summary(Sk/Nk,Nk,Yk,Mk,baseN-if(k)rowSums(present[,ids %in% remove,drop=FALSE]) else 0,keys,sims)
        a[,`:=`(branch=branch,half=half,scenario=sc,model=m,k_removed=k,removed_ids=paste(remove,collapse=','))]
        top_deletions[[length(top_deletions)+1L]]<-a
      }
      v[,`:=`(branch=branch,half=half,scenario=sc,model=m)]
      all_counties[[length(all_counties)+1L]]<-v[,.(branch,half,scenario,model,sim,id,x,k,sum_y,sum_mu,e,keep_z,keep_p,keep_pre)]
      message('Influence complete ',group_index,'/72: ',branch,' / ',half,' / ',sc,' / ',m)
    }
    aggregate_guard(obj,root,half)
  }
  profile<-rbindlist(profiles)
  mult<-profile[,.(top5_models_same_scenario=uniqueN(model),models=paste(sort(unique(model)),collapse=', ')),by=.(branch,half,scenario,id)]
  profile<-merge(profile,mult,by=c('branch','half','scenario','id'));profile[,influential_several_models:=top5_models_same_scenario>=2]
  top<-rbindlist(top_deletions)
  gaps<-dcast(top[diagnostic %in% c('ho_standardized','pre_rmse')],branch+scenario+model+diagnostic+cutoff+k_removed~half,value.var=c('rmse_ratio','rmse_ratio_mc_se'))
  gaps[,`:=`(B_minus_A=rmse_ratio_B-rmse_ratio_A,absolute_AB_gap=abs(rmse_ratio_B-rmse_ratio_A))]
  ref<-gaps[k_removed==0,.(branch,scenario,model,diagnostic,cutoff,original_absolute_gap=absolute_AB_gap)]
  gaps<-merge(gaps,ref,by=c('branch','scenario','model','diagnostic','cutoff'));gaps[,`:=`(gap_shrinkage=original_absolute_gap-absolute_AB_gap,gap_shrinks=absolute_AB_gap<original_absolute_gap-1e-12)]
  exports<-list(influence_county_errors_flags=rbindlist(all_counties),influence_county_MSE_contributions=rbindlist(county_contributions),influence_concentration=rbindlist(concentrations),
    influence_leave_one_out=rbindlist(loos),influence_top5_profiles=profile,influence_topk_out=top,influence_AB_gaps=gaps,influence_flag_capture=rbindlist(captures),
    influence_full_sample=rbindlist(baseline_summaries),influence_reproduction_checks=rbindlist(reproductions))
  for(name in names(exports))fwrite(exports[[name]],file.path(out,paste0(name,'.csv')))
  stopifnot(identical(before,tools::md5sum(names(before))))
  writeLines('Descriptive and post hoc; original flags and percentile references frozen. No refit, target reconstruction, decision change or rule suggestion.',file.path(out,'influence_scope.txt'))
  invisible(exports)
}

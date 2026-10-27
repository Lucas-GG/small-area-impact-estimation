# Additive aggregate-estimand analysis. Existing checkpoints/exports are read-only.
aggregate_guard <- function(obj, root, expected_half='A') {
  stopifnot(expected_half %in% c('A','B'))
  assignment<-readRDS(file.path(root,'data/derived/heldout/county_assignment.rds'))
  stopifnot(all(obj$rows$id %in% assignment$id[assignment$half==expected_half]))
  stopifnot(nrow(obj$rows)>0L, !anyNA(obj$rows$half), all(obj$rows$half==expected_half),
    identical(obj$assignment_md5, unname(tools::md5sum(file.path(root,'data/derived/heldout/county_assignment.rds')))))
}

aggregate_oracle <- function(v) {
  if(any(!is.finite(v)) || any(v<0)) stop('Oracle requires finite leave-one-simulation-out variances for every eligible county.')
  # One threshold over ALL current eligible subgroup counties; never recompute it.
  v<=2*mean(v)
}

aggregate_summary <- function(d) {
  groups<-c('subgroup','scenario','model','diagnostic','cutoff','rule_type','missing_handling')
  base<-d[diagnostic=='unscreened',.(subgroup,scenario,model,sim,base_error=error,base_error_z=error_z)]
  d<-merge(d,base,by=c('subgroup','scenario','model','sim'),all.x=TRUE)
  d[, {
    ans<-list(n_sims=.N,n_error=sum(is.finite(error)),n_interval=sum(is.finite(coverage)),n_oracle_unavailable=if('oracle_missing_variance_count' %in% names(.SD))sum(oracle_missing_variance_count>0,na.rm=TRUE) else 0L)
    avg<-function(v) {v<-v[is.finite(v)]; c(mean=if(length(v)) mean(v) else NA_real_,se=if(length(v)>1) sd(v)/sqrt(length(v)) else NA_real_)}
    for(nm in c('bias','mae','coverage','width','n_retained','excluded_share','standardized_bias','standardized_mae')) {
      v<-switch(nm,bias=error,mae=abs(error),standardized_bias=error_z,standardized_mae=abs(error_z),get(nm))
      a<-avg(v);ans[[nm]]<-a[1];ans[[paste0(nm,'_mc_se')]]<-a[2]
    }
    for(nm in c('rmse','standardized_rmse')) {
      v<-if(nm=='rmse') error else error_z;v<-v[is.finite(v)];r<-sqrt(mean(v^2))
      ans[[nm]]<-r;ans[[paste0(nm,'_mc_se')]]<-if(length(v)>1 && r>0) sd(v^2)/(2*r*sqrt(length(v))) else NA_real_
    }
    ok<-is.finite(error)&is.finite(base_error);a<-error[ok]^2;b<-base_error[ok]^2
    ratio<-sqrt(mean(a)/mean(b)); influence<-(a-ratio^2*b)/(2*ratio*mean(b))
    ans$paired_n_sims<-length(a);ans$rmse_ratio<-ratio
    ans$rmse_ratio_mc_se<-if(length(a)>1 && ratio>0) sd(influence)/sqrt(length(a)) else NA_real_
    ans$baseline_rmse_paired<-sqrt(mean(b))
    relative_ok<-is.finite(error)&is.finite(observed_aggregate)
    e<-error[relative_ok];y<-observed_aggregate[relative_ok]
    mean_y<-mean(y);relative_rmse<-sqrt(mean(e^2))/mean_y
    relative_influence<-(e^2-mean(e^2))/(2*sqrt(mean(e^2))*mean_y)-sqrt(mean(e^2))*(y-mean_y)/mean_y^2
    ans$mean_observed_aggregate<-mean_y
    ans$mean_observed_aggregate_mc_se<-if(length(y)>1)sd(y)/sqrt(length(y)) else NA_real_
    ans$relative_rmse<-if(is.finite(mean_y)&&mean_y>0)relative_rmse else NA_real_
    ans$relative_rmse_mc_se<-if(length(e)>1 && is.finite(mean_y)&&mean_y>0)sd(relative_influence)/sqrt(length(e)) else NA_real_
    ans$relative_n_sims<-length(e)
    ans$max_excluded_share<-if(any(is.finite(excluded_share)))max(excluded_share,na.rm=TRUE) else NA_real_
    ans$oracle_main_visible<-if(.BY$diagnostic=='oracle') {if(any(is.finite(excluded_share)))max(excluded_share,na.rm=TRUE)<=.4 else FALSE} else NA
    ok_z<-is.finite(error_z)&is.finite(base_error_z);az<-error_z[ok_z]^2;bz<-base_error_z[ok_z]^2
    ratio_z<-sqrt(mean(az)/mean(bz));influence_z<-(az-ratio_z^2*bz)/(2*ratio_z*mean(bz))
    ans$standardized_paired_n_sims<-length(az);ans$standardized_rmse_ratio<-ratio_z
    ans$standardized_rmse_ratio_mc_se<-if(length(az)>1 && ratio_z>0)sd(influence_z)/sqrt(length(az)) else NA_real_
    ans
  },by=groups]
}

aggregate_manuscript_summary <- function(x) {
  x[, {
    ok<-is.finite(new_error)&is.finite(manuscript_complete_error)
    e<-new_error[ok];f<-manuscript_complete_error[ok]
    ra<-sqrt(mean(e^2));rb<-sqrt(mean(f^2))
    influence<-(e^2-mean(e^2))/(2*ra)-(f^2-mean(f^2))/(2*rb)
    se<-sd(influence)/sqrt(length(e))
    list(aggregate_rmse=sqrt(mean(new_error^2)),aggregate_rmse_mc_se=sd(new_error^2)/(2*sqrt(mean(new_error^2))*sqrt(.N)),
      manuscript_half_A_available_rmse=sqrt(mean(manuscript_available_error^2)),max_abs_replication_difference=max(abs(new_error-manuscript_available_error)),
      manuscript_half_A_complete_rmse=rb,manuscript_complete_rmse_mc_se=sd(f^2)/(2*rb*sqrt(length(f))),
      complete_window_difference=ra-rb,complete_window_difference_mc_se=se,complete_window_paired_sims=length(e),
      complete_window_agrees_within_mc_error=abs(ra-rb)<=1.96*se,
      mean_n_available=mean(n_available),mean_n_complete=mean(n_complete))
  },by=.(subgroup,scenario,model)]
}

aggregate_county_metadata <- function(root, a_ids) {
  p<-file.path(root,Sys.getenv('SAI_PANEL','data/notpublic/compact_dt.rds'))
  if(!file.exists(p)) stop('Missing noncore metadata file: ',p,' (required columns i, urb).')
  panel<-as.data.table(readRDS(p))[i %in% a_ids]
  if(!all(c('i','urb','n','year','start_year') %in% names(panel))) stop('Missing metadata columns in ',p,': ',paste(setdiff(c('i','urb','n','year','start_year'),names(panel)),collapse=', '))
  meta<-panel[, {
    u<-unique(as.character(urb));if(length(u)!=1L || is.na(u)) stop('Missing or changing urb for county ',i[1])
    pop<-n[year<start_year & is.finite(n)]
    list(noncore=u=='6.NonCore',pre_population=if(length(pop)) mean(pop) else NA_real_)
  },by=.(id=i)]
  if(!all(a_ids %in% meta$id)) stop('Missing noncore indicators for A counties: ',paste(setdiff(a_ids,meta$id),collapse=', '))
  meta
}

aggregate_predictive_A <- function(root, model, scenario, obj, cache_dir, expected_half='A', input_path=file.path(root,'data/derived/heldout/half_A_available.rds'), main_dir=file.path(root,'data/simulation/wocov'), file_prefix='new_imputation') {
  aggregate_guard(obj,root,expected_half)
  source(file.path(root,'R/screening/s3_models.R'),local=TRUE);codes<-s3_codes
  sc<-match(scenario,c('MAR','Latent level','Latent trend'));mi<-s3_seed_index[[model]]
  main_path<-s3_main_path(root,model,sc,file_prefix=='wcov_imputation')
  signature<-list(assignment=obj$assignment_md5,predictive_version=1L,
    main_path=main_path,bytes=file.info(main_path)$size,mtime=as.numeric(file.info(main_path)$mtime),input=as.numeric(file.info(input_path)$mtime))
  cache<-file.path(cache_dir,sprintf('%s_%s.rds',codes[[model]],sc))
  if(file.exists(cache)) {o<-readRDS(cache);if(identical(o$signature,signature)) {stopifnot(all(o$rows$half==expected_half));return(o$rows)}}
  message('Reading saved draws, selecting ',expected_half,' immediately: ',model,' / ',scenario)
  # Saved fits are serialized whole panels; only A is retained, before any calculation.
  a_ids<-unique(obj$rows$id);raw<-readRDS(main_path)
  draws<-grep('^y0_post_[0-9]+$',names(raw),value=TRUE)
  draws<-draws[order(as.integer(sub('y0_post_','',draws)))]
  if(length(draws)!=200L) stop('Expected the existing 200 matching predictive draws in ',main_path)
  required<-c('r','i','year','start_year_sim','y','y0_hat');stopifnot(all(required %in% names(raw)))
  a<-as.data.table(raw)[i %in% a_ids & is.finite(start_year_sim) & year>=start_year_sim & year<=start_year_sim+4,c(required,draws),with=FALSE]
  rm(raw);gc();stopifnot(all(a$i %in% a_ids),!anyDuplicated(a[,.(r,i,year)]))
  pred<-a[, {
    usable<-is.finite(y)&is.finite(y0_hat)&y0_hat>0
    lam<-as.matrix(.SD[,draws,with=FALSE]);bad<-!is.finite(lam)|lam<0
    seed<-50000000L+sc*1000000L+mi*100000L+r[1]*1000L+i[1]
    lambdas<-as.numeric(lam);lambdas[as.numeric(bad)==1]<-0
    predictive<-matrix(with_backdate_seed(seed,suppressWarnings(rpois(length(lambdas),lambdas))),nrow=nrow(lam));predictive[bad]<-NA_real_
    list(k_check=sum(usable),y_check=as.numeric(sum(y[usable])),mu_check=as.numeric(sum(y0_hat[usable])),
      draws=list(if(any(usable) && !any(bad[usable,,drop=FALSE])) colSums(predictive[usable,,drop=FALSE]) else rep(NA_real_,200)),
      legacy_y=if(sum(usable)==5 && setequal(year[usable]-start_year_sim[usable],0:4)) as.numeric(sum(y[usable])) else NA_real_,
      legacy_mu=if(sum(usable)==5 && setequal(year[usable]-start_year_sim[usable],0:4)) as.numeric(sum(y0_hat[usable])) else NA_real_)
  },by=.(sim=r,id=i)]
  pred[,half:=expected_half];saveRDS(list(signature=signature,rows=pred),cache);rm(a);gc();pred
}

run_heldout_aggregate <- function(root) {
  library(data.table)
  source(file.path(root,'R/screening/heldout_analysis.R'));source(file.path(root,'R/backdate.r'));source(file.path(root,'R/screening/heldout_exploration.R'))
  obj<-readRDS(file.path(root,'data/derived/heldout/half_A_available.rds'));aggregate_guard(obj,root)
  z<-derive_exploration_A(obj$rows)$rows
  meta<-aggregate_county_metadata(root,unique(z$id));z<-merge(z,meta,by='id',all.x=TRUE);stopifnot(!anyNA(z$noncore))
  out<-file.path(root,'_output/analysis/output/3.1a');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  cache_dir<-file.path(root,'.quarto/heldout-aggregate');dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  originals<-list.files(out,'[.]csv$',full.names=TRUE);originals<-originals[!grepl('/aggregate_',originals)]
  before<-tools::md5sum(c(originals,file.path(root,'analysis/S6_screening/s3_decision.yml')))
  rules<-fread(file.path(out,'rule_grid.csv'))
  rules<-unique(rules[,.(diagnostic,cutoff,rule_type)])
  expected<-exploration_rules()[,.(diagnostic,cutoff,rule_type)];stopifnot(fsetequal(rules,expected));rules[,rule_id:=seq_len(.N)]
  groups<-unique(z[,.(scenario,model)]);ledgers<-list();randoms<-list();comparisons<-list();count_checks<-list();draw_checks<-list()
  for(g in seq_len(nrow(groups))) {
    key<-groups[g];v<-copy(z[model==key$model & scenario==key$scenario]);pred<-aggregate_predictive_A(root,key$model,key$scenario,obj,cache_dir); random_acc<-list()
    v<-merge(v,pred[,half:=NULL],by=c('sim','id'),all.x=TRUE)
    stopifnot(all(v$k[is.na(v$k_check)]==0));v[is.na(k_check),k_check:=0L]
    stopifnot(all(v$half=='A'),all(v$k==v$k_check),all(abs(v$sum_y-v$y_check)<1e-8,na.rm=TRUE),all(abs(v$sum_mu-v$mu_check)<1e-8,na.rm=TRUE))
    valid<-v$k>=1 & vapply(v$draws,function(a)all(is.finite(a)),logical(1))
    pi<-t(vapply(v$draws[valid],function(a)quantile(a,c(.025,.975),names=FALSE,type=7),numeric(2)))
    stopifnot(all(abs(pi[,1]-v$pi_lower[valid])<1e-8),all(abs(pi[,2]-v$pi_upper[valid])<1e-8))
    draw_checks[[g]]<-data.table(scenario=key$scenario,model=key$model,county_intervals_checked=sum(valid),exact_predictive_interval_match=TRUE)
    v[,`:=`(signed_error=sum_y-sum_mu,county_standardized_error=(sum_y-sum_mu)/sqrt(sum_mu))]
    # Variance excludes the current replication in both numerator and denominator.
    totals<-v[k>=1 & is.finite(county_standardized_error),.(n=.N,s1=sum(county_standardized_error),s2=sum(county_standardized_error^2)),by=id]
    v<-merge(v,totals,by='id',all.x=TRUE)
    v[,oracle_variance:=fifelse(n>=3,pmax(0,(s2-county_standardized_error^2-(s1-county_standardized_error)^2/(n-1))/(n-2)),NA_real_)]
    for(simulation in sort(unique(v$sim))) {
      if(simulation %% 50L==0L) message(key$model,' / ',key$scenario,': simulation ',simulation)
      full<-v[sim==simulation];screen<-keep_matrix_A(full,rules)
      # Screens use the original full A pool, including the existing top-percent denominator.
      eligible<-full$k>=1 & is.finite(full$signed_error)
      cnt<-copy(screen$keys);cnt[,`:=`(scenario=key$scenario,model=key$model,sim=simulation,n_retained=colSums(screen$keep[eligible,,drop=FALSE]))]
      count_checks[[length(count_checks)+1L]]<-cnt
      for(subgroup_name in c('all_A','noncore')) {
        take<-eligible & if(subgroup_name=='noncore') full$noncore else TRUE
        u<-full[take];if(!nrow(u)) stop('Empty subgroup ',subgroup_name)
        km<-screen$keep[take,,drop=FALSE]
        oracle_missing<-sum(!is.finite(u$oracle_variance))
        oracle<-if(oracle_missing) rep(FALSE,nrow(u)) else aggregate_oracle(u$oracle_variance)
        km<-cbind(unscreened=TRUE,km,oracle=oracle)
        keys<-rbindlist(list(data.table(diagnostic='unscreened',cutoff=NA_real_,rule_type='baseline',missing_handling='none'),screen$keys[,.(diagnostic,cutoff,rule_type,missing_handling)],data.table(diagnostic='oracle',cutoff=2,rule_type='infeasible_LOSO_standardized_variance_single_pass',missing_handling='none')))
        nr<-colSums(km);ymean<-as.numeric(crossprod(km,u$sum_y))/nr;mumean<-as.numeric(crossprod(km,u$sum_mu))/nr
        error<-ymean-mumean;error[nr==0]<-NA_real_
        error_z<-error/sqrt(as.numeric(crossprod(km,u$sum_mu))/nr^2)
        mat<-do.call(rbind,u$draws);good<-apply(mat,1,function(a)all(is.finite(a)))
        mat[!is.finite(mat)]<-0
        summed<-crossprod(km,mat)/nr
        good_rule<-as.numeric(crossprod(km,!good))==0 & nr>0
        intervals<-matrix(NA_real_,ncol(km),2)
        intervals[good_rule,]<-t(apply(summed[good_rule,,drop=FALSE],1,quantile,probs=c(.025,.975),names=FALSE,type=7))
        a<-copy(keys);a[,`:=`(subgroup=subgroup_name,scenario=key$scenario,model=key$model,sim=simulation,n_eligible=nrow(u),n_retained=nr,excluded_share=1-nr/nrow(u),
          error=error,error_z=error_z,observed_aggregate=ymean,predicted_aggregate=mumean,lower=intervals[,1],upper=intervals[,2],coverage=as.numeric(ymean>=intervals[,1]&ymean<=intervals[,2]),width=intervals[,2]-intervals[,1],n_invalid_predictive_counties=as.numeric(crossprod(km,!good)))]
        a[,oracle_missing_variance_count:=ifelse(diagnostic=='oracle',oracle_missing,NA_integer_)]
        if(oracle_missing) a[diagnostic=='oracle',`:=`(n_retained=NA_real_,excluded_share=NA_real_)]
        ledgers[[length(ledgers)+1L]]<-a
        # One seeded permutation per random draw; nested sets provide paired cutoff curves.
        seed<-72000000L+g*10000L+simulation*10L+match(subgroup_name,c('all_A','noncore'))
        random_errors<-with_backdate_seed(seed,replicate(100,{ord<-sample.int(nrow(u));cumsum(u$signed_error[ord])/seq_len(nrow(u))}))
        sq<-matrix(NA_real_,nrow(keys),100);ok<-nr>0;sq[ok,]<-random_errors[nr[ok],,drop=FALSE]^2
        if(is.null(random_acc[[subgroup_name]])) random_acc[[subgroup_name]]<-list(sq=matrix(0,nrow(keys),100),base=matrix(0,nrow(keys),100),n=matrix(0,nrow(keys),100),keys=copy(keys))
        ok_sq<-is.finite(sq);finite_sq<-sq;finite_sq[!ok_sq]<-0
        random_acc[[subgroup_name]]$sq<-random_acc[[subgroup_name]]$sq+finite_sq
        random_acc[[subgroup_name]]$base<-random_acc[[subgroup_name]]$base+ok_sq*error[1]^2
        random_acc[[subgroup_name]]$n<-random_acc[[subgroup_name]]$n+ok_sq
        legacy<-is.finite(u$legacy_y)&is.finite(u$legacy_mu)
        comparisons[[length(comparisons)+1L]]<-data.table(subgroup=subgroup_name,scenario=key$scenario,model=key$model,sim=simulation,new_error=error[1],
          manuscript_available_error=mean(u$y_check-u$mu_check),manuscript_complete_error=if(any(legacy))mean(u$legacy_y[legacy]-u$legacy_mu[legacy]) else NA_real_,n_available=nrow(u),n_complete=sum(legacy))
      }
    }
    for(subgroup_name in names(random_acc)) {
      acc<-random_acc[[subgroup_name]];rr<-copy(acc$keys)
      rr[,`:=`(subgroup=subgroup_name,scenario=key$scenario,model=key$model)]
      rr<-rr[rep(seq_len(.N),each=100)]
      rr[,`:=`(draw=rep(1:100,times=nrow(acc$keys)),rmse_ratio=as.vector(t(sqrt(acc$sq/acc$base))),n_sims=as.vector(t(acc$n)))]
      randoms[[length(randoms)+1L]]<-rr
    }
    message('Aggregate group complete ',g,'/',nrow(groups),': ',key$model,' / ',key$scenario)
  }
  d<-rbindlist(ledgers);summary<-aggregate_summary(d)
  random_draws<-rbindlist(randoms);group_cols<-c('subgroup','scenario','model','diagnostic','cutoff','rule_type','missing_handling')
  bands<-random_draws[,.(random_lo=quantile(rmse_ratio,.025,na.rm=TRUE),random_hi=quantile(rmse_ratio,.975,na.rm=TRUE),random_mean=mean(rmse_ratio,na.rm=TRUE)),by=group_cols]
  summary<-merge(summary,bands,by=group_cols,all.x=TRUE)
  # Keep the scale-normalized error beside the absolute RMSE ratio in full tables.
  columns<-names(summary);columns<-setdiff(columns,c('relative_rmse','relative_rmse_mc_se'))
  where<-match('rmse_ratio_mc_se',columns)
  setcolorder(summary,append(columns,c('relative_rmse','relative_rmse_mc_se'),after=where))
  counts<-rbindlist(count_checks)
  old<-fread(file.path(out,'retention_per_simulation.csv'))[target=='err_z',.(scenario,model,sim,diagnostic,cutoff,rule_type,missing_handling,n_retained)]
  keys<-c('scenario','model','sim','diagnostic','cutoff','rule_type','missing_handling')
  count_comparison<-merge(counts,old,by=keys,suffixes=c('_aggregate','_existing'),all=TRUE)
  stopifnot(!anyNA(count_comparison$n_retained_aggregate),!anyNA(count_comparison$n_retained_existing),all(count_comparison$n_retained_aggregate==count_comparison$n_retained_existing))
  means<-counts[,.(n_retained=mean(n_retained)),by=setdiff(keys,'sim')]
  old_summary<-fread(file.path(out,'retention_summary.csv'))[target=='err_z',c(setdiff(keys,'sim'),'n_retained'),with=FALSE]
  check<-merge(means,old_summary,by=setdiff(keys,'sim'),suffixes=c('_new','_old'));stopifnot(nrow(check)==nrow(means),all(abs(check$n_retained_new-check$n_retained_old)<1e-8))
  comparison_per_sim<-rbindlist(comparisons)
  comparison<-aggregate_manuscript_summary(comparison_per_sim)
  stopifnot(all(comparison$max_abs_replication_difference<1e-8))
  strata<-merge(unique(z[,.(id,stratum,mean_pre_expected_count)]),meta,by='id')[,.(n_A_counties=.N,population_min=min(pre_population,na.rm=TRUE),population_max=max(pre_population,na.rm=TRUE),expected_count_min=if(any(is.finite(mean_pre_expected_count)))min(mean_pre_expected_count,na.rm=TRUE) else NA_real_,expected_count_max=if(any(is.finite(mean_pre_expected_count)))max(mean_pre_expected_count,na.rm=TRUE) else NA_real_),by=stratum]
  setorder(strata,stratum)
  floor_checks<-d[diagnostic=='p_two' & cutoff %in% c(.001,.0025,.005),.(identical_counts=uniqueN(n_retained)==1,identical_errors=uniqueN(error)==1),by=.(subgroup,scenario,model,sim,missing_handling)]
  shortlist<-summary[(diagnostic=='p_two' & cutoff %in% c(.01,.05,.10)) | (diagnostic=='ho_standardized' & cutoff %in% c(2,2.5)) | (diagnostic=='pre_rmse' & cutoff==10) | (diagnostic=='oracle' & oracle_main_visible)]
  minima<-summary[!diagnostic %in% c('unscreened','oracle') & is.finite(rmse_ratio)][order(rmse_ratio,cutoff),.SD[1],by=.(subgroup,scenario,model,diagnostic,missing_handling)]
  minima[,`:=`(pointwise_significantly_below_one=rmse_ratio+1.96*rmse_ratio_mc_se<1,selection_adjusted=FALSE)]
  exports<-list(aggregate_per_simulation=d,aggregate_summary=summary,aggregate_random_by_draw=random_draws,aggregate_random_bands=bands,aggregate_scorecard=shortlist,aggregate_curve_minima=minima,
    aggregate_manuscript_comparison=comparison,aggregate_manuscript_per_simulation=comparison_per_sim,aggregate_retained_count_checks=count_comparison,aggregate_predictive_checks=rbindlist(draw_checks),aggregate_size_strata=strata,aggregate_resolution_checks=floor_checks)
  for(nm in names(exports)) fwrite(exports[[nm]],file.path(out,paste0(nm,'.csv')))
  stopifnot(identical(before,tools::md5sum(names(before))))
  aggregate_guard(obj,root)
  print(minima[pointwise_significantly_below_one==TRUE,.(subgroup,scenario,model,diagnostic,cutoff,missing_handling,rmse_ratio,rmse_ratio_mc_se)])
  invisible(exports)
}

if(sys.nframe()==0L) run_heldout_aggregate(normalizePath(getwd()))

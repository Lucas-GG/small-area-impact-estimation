# Half-A revision: targets and evaluation; no models are estimated here.
post_window_summary <- function(y, mu, event, latent_draws = NULL, seed = 2032031L) {
  usable <- is.finite(y) & is.finite(mu) & mu > 0 & event %in% 0:4
  k <- sum(usable)
  summarize <- function(take, predictive = NULL) {
    if (!any(take)) return(c(raw=NA_real_,z=NA_real_,coverage=NA_real_,lower=NA_real_,upper=NA_real_,variance=NA_real_))
    sy <- sum(y[take]); sm <- sum(mu[take]); raw <- abs(sy-sm)
    out <- c(raw=raw,z=if(sm>0) raw/sqrt(sm) else NA_real_,coverage=NA_real_,lower=NA_real_,upper=NA_real_,variance=NA_real_)
    if (!is.null(predictive) && all(is.finite(latent_draws[take,,drop=FALSE]) & latent_draws[take,,drop=FALSE]>=0)) {
      total <- colSums(predictive[take,,drop=FALSE]); if(!all(is.finite(total))) return(out)
      interval <- quantile(total,c(.025,.975),names=FALSE,type=7)
      out[c('coverage','lower','upper','variance')] <- c(as.numeric(sy>=interval[1]&sy<=interval[2]),interval,var(colSums(latent_draws[take,,drop=FALSE])))
    }
    out
  }
  predictive <- NULL
  if(!is.null(latent_draws)) {
    stopifnot(nrow(latent_draws)==length(y),ncol(latent_draws)>=2)
    # Independent Poisson observation noise conditional on matching latent draw indices.
    lam <- as.numeric(latent_draws); bad <- !is.finite(lam)|lam<0; lam[bad]<-0
    predictive <- matrix(with_backdate_seed(seed,suppressWarnings(rpois(length(lam),lam))),nrow=length(y));predictive[bad]<-NA_real_
  }
  primary <- summarize(usable,predictive)
  h0 <- summarize(usable & event==0,predictive)
  old_usable <- is.finite(y)&is.finite(mu)&mu>=0&event %in% 0:4
  complete <- sum(old_usable)==5L && setequal(event[old_usable],0:4)
  k5 <- summarize(if(complete) old_usable else rep(FALSE,length(y)),predictive)
  data.table(k=as.integer(k),err_raw=primary['raw'],err_z=primary['z'],coverage=primary['coverage'],
    pi_lower=primary['lower'],pi_upper=primary['upper'],mu_post_var=primary['variance'],
    sum_y=if(k) as.numeric(sum(y[usable])) else NA_real_,sum_mu=if(k) as.numeric(sum(mu[usable])) else NA_real_,
    n_valid_predictive_draws=if(k && !is.null(predictive)) as.integer(sum(colSums(!is.finite(predictive[usable,,drop=FALSE]))==0L)) else 0L,
    h0_raw=h0['raw'],h0_z=h0['z'],h0_coverage=h0['coverage'],
    k5_raw=k5['raw'],k5_z=k5['z'],k5_coverage=k5['coverage'])
}

prepare_exploration_A <- function(root, split_seed=2032026L) {
  library(data.table)
  source(file.path(root,'R/backdate.r'))
  assignment_path<-file.path(root,'data/derived/heldout/county_assignment.rds')
  if(!file.exists(assignment_path)) stop('The fixed assignment must already exist; this revision never creates a new split.')
  assignment<-as.data.table(readRDS(assignment_path)); stopifnot(all(assignment$split_seed==split_seed))
  hash_before<-unname(tools::md5sum(assignment_path)); a_ids<-assignment[half=='A',id]
  out<-file.path(root,'data/derived/heldout');cache<-file.path(out,'post_A');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  checkpoint_paths<-list.files(file.path(root,'data/fits/backdate'),'^summary_.*[.]rds$',full.names=TRUE)
  source(file.path(root,'R/screening/s3_models.R'),local=TRUE);
  if(length(checkpoint_paths)!=3L*length(s3_models)) stop('Expected one saved no-covariate checkpoint per S3 model and scenario.')
  codes<-s3_codes
  scenarios<-c('MAR','Latent level','Latent trend'); pieces<-list(); manifests<-list()
  for(p in checkpoint_paths) {
    # The RDS format requires deserialization; immediately select A before any computation.
    cp<-as.data.table(readRDS(p))[id %in% a_ids];stopifnot(nrow(cp)>0,uniqueN(cp$model)==1,uniqueN(cp$scenario)==1)
    model<-unique(cp$model);scenario<-unique(cp$scenario); sc<-match(scenario,scenarios)
    main_path<-s3_main_path(root,model,sc)
    signature<-list(path=main_path,bytes=file.info(main_path)$size,mtime=as.numeric(file.info(main_path)$mtime),assignment=hash_before,version=1L)
    cached_path<-file.path(cache,basename(p));cached<-if(file.exists(cached_path)) readRDS(cached_path) else NULL
    if(!is.null(cached)&&identical(cached$signature,signature)) post<-cached$rows else {
      message('Reading saved main fit for half A: ',model,' / ',scenario)
      raw<-readRDS(main_path)
      needed<-c('r','i','year','start_year_sim','y','y0_hat')
      draws<-grep('^y0_post_[0-9]+$',names(raw),value=TRUE)
      if(!length(draws)) stop('STOP: no saved post-period draws in ',main_path,'; available columns: ',paste(names(raw),collapse=', '))
      stopifnot(all(needed %in% names(raw)))
      main<-as.data.table(raw)[i %in% a_ids & is.finite(start_year_sim) & year>=start_year_sim & year<=start_year_sim+4,
        c(needed,draws),with=FALSE]
      rm(raw);gc();stopifnot(all(main$i %in% a_ids),!anyDuplicated(main[,.(r,i,year)]))
      main[,event:=year-start_year_sim]
      post<-main[, {
        key_seed<-50000000L+sc*1000000L+s3_seed_index[[model]]*100000L+r[1]*1000L+i[1]
        post_window_summary(y,y0_hat,event,as.matrix(.SD[,draws,with=FALSE]),key_seed)
      },by=.(sim=r,id=i,x=start_year_sim)]
      saveRDS(list(signature=signature,rows=post),cached_path);rm(main);gc()
    }
    # Saved old complete-window totals remain the authority for k5 sensitivity.
    cp[,`:=`(old_realized_raw=abs(y_post-mu_post_mean),
      old_realized_standardized=fifelse(mu_post_mean>0,abs(y_post-mu_post_mean)/sqrt(mu_post_mean),NA_real_))]
    cp<-merge(cp,post,by=c('sim','id','x'),all.x=TRUE)
    cp[is.na(k),k:=0L]
    stopifnot(identical(is.finite(cp$k5_raw),is.finite(cp$old_realized_raw)),
      identical(is.finite(cp$k5_z),is.finite(cp$old_realized_standardized)),
      all(abs(cp$k5_raw-cp$old_realized_raw)<1e-8,na.rm=TRUE),
      all(abs(cp$k5_z-cp$old_realized_standardized)<1e-8,na.rm=TRUE))
    if(!'tail_convention' %in% names(cp)) cp[,tail_convention:=NA_character_]
    for(v in c('ppp_ho_upper','ppp_ho_lower')) if(!v %in% names(cp)) cp[,(v):=NA_real_]
    legacy<-is.na(cp$tail_convention)|!grepl('^mid-p',cp$tail_convention)
    both<-is.finite(cp$ppp_ho_upper)&is.finite(cp$ppp_ho_lower);recover<-legacy&both
    eq<-cp$ppp_ho_upper[recover]+cp$ppp_ho_lower[recover]-1;stopifnot(all(eq>=-1e-12 & eq<=1+1e-12))
    cp[recover,`:=`(ppp_ho_upper=ppp_ho_upper-.5*pmax(0,eq),ppp_ho_lower=ppp_ho_lower-.5*pmax(0,eq),tail_convention='mid-p recovered from inclusive pair')]
    cp[legacy & !both,`:=`(ppp_ho_upper=NA_real_,ppp_ho_lower=NA_real_)]
    cp[,intersect(c('ppp_ho','ppp_ho_two_sided'),names(cp)):=NULL]
    cp<-merge(cp,assignment[half=='A'],by='id',all.x=TRUE);stopifnot(all(cp$half=='A'))
    pieces[[length(pieces)+1L]]<-cp
    manifests[[length(manifests)+1L]]<-data.table(file=basename(p),main_file=main_path,bytes=file.info(p)$size,modified=as.character(file.info(p)$mtime))
  }
  x<-rbindlist(pieces,fill=TRUE);stopifnot(all(x$half=='A'),!anyDuplicated(x[,.(scenario,model,sim,id)]),identical(hash_before,unname(tools::md5sum(assignment_path))))
  # Preserve the old A input and all B files. The revised document loads only this new A file.
  saveRDS(list(rows=x,manifest=rbindlist(manifests),prepared=Sys.time(),expected_groups=3L*length(s3_models),assignment_md5=hash_before),file.path(out,'half_A_available.rds'))
  counts<-x[,.(rows=.N,old_complete_rows=sum(is.finite(old_realized_raw)),new_available_rows=sum(k>=1),
    usable_counties=uniqueN(id[k>=1]),old_usable_counties=uniqueN(id[is.finite(old_realized_raw)])),by=.(scenario,model)]
  fwrite(counts,file.path(out,'available_window_counts.csv'));print(counts)
  invisible(counts)
}

exploration_targets<-c('err_raw','err_z','h0_raw','h0_z','k5_raw','k5_z','average_raw','average_z','average_k5_raw','average_k5_z')
coverage_for_target<-function(t) if(startsWith(t,'h0')) 'h0_coverage' else if(grepl('k5',t,fixed=TRUE)) 'k5_coverage' else 'coverage'
derive_exploration_A<-function(rows,min_k_average=2L,expected_half='A') {
  stopifnot(expected_half %in% c('A','B'),all(rows$half==expected_half),min_k_average %in% 1:5)
  z<-copy(rows)
  z[,`:=`(ho_raw=y_ho-mu_ho_mean,
    ho_standardized=fifelse(mu_ho_mean>0,(y_ho-mu_ho_mean)/sqrt(mu_ho_mean),NA_real_),
    ho_uncertainty=fifelse(mu_ho_mean+mu_ho_var>0,(y_ho-mu_ho_mean)/sqrt(mu_ho_mean+mu_ho_var),NA_real_),
    p_two=pmin(1,2*pmin(ppp_ho_upper,ppp_ho_lower)))]
  averages<-z[, {
    eligible<-k>=min_k_average
    list(average_raw=mean_finite(err_raw[eligible]),average_z=mean_finite(err_z[eligible]),
      n_average_raw=sum(is.finite(err_raw[eligible])),n_average_z=sum(is.finite(err_z[eligible])),
      average_k5_raw=mean_finite(k5_raw),average_k5_z=mean_finite(k5_z),
      n_average_k5_raw=sum(is.finite(k5_raw)),n_average_k5_z=sum(is.finite(k5_z)))
  },by=.(scenario,model,id)]
  z<-merge(z,averages,by=c('scenario','model','id'),all.x=TRUE)
  list(rows=z,averages=averages)

}
mc_summary<-function(x,values,groups) {
  # Input has one value per simulation/group, after any weighted k aggregation.
  x[, {
    ans<-list()
    for(v in values) {
      a<-get(v);a<-a[is.finite(a)];n<-length(a);m<-if(n) mean(a) else NA_real_;se<-if(n>1) sd(a)/sqrt(n) else NA_real_
      ans[[v]]<-m;ans[[paste0(v,'_mc_se')]]<-se;ans[[paste0(v,'_n_sims')]]<-n
      ans[[paste0(v,'_display')]]<-if(n) sprintf('%.3f (%.3f)',m,se) else 'NA'
    };ans
  },by=groups]
}
association_stats<-function(v,group_cols) {
  v[, {
    observed<-error[is.finite(error)];boundary<-if(length(observed)) quantile(observed,.9,names=FALSE) else NA_real_
    ok<-is.finite(error)&is.finite(score);s<-score[ok];y<-error[ok];top<-y>=boundary
    usable<-length(y)>=5 && any(top) && any(!top)
    list(n_target=length(observed),n=length(y),n_top=sum(top),usable=usable,
      spearman=if(usable && sd(s)>0 && sd(y)>0) cor(s,y,method='spearman') else NA_real_,
      auc=if(usable) auc_rank(s,top) else NA_real_)
  },by=group_cols]
}
exploration_associations<-function(z) {
  out<-list();strata<-list()
  for(d in names(diagnostic_labels)) for(t in exploration_targets) {
    v<-z[if(startsWith(t,'average')) rep(TRUE,.N) else k>=1,.(scenario,model,sim,k,id,score=if(d=='p_two') 1-get(d) else abs(get(d)),error=get(t))]
    byk<-association_stats(v,c('scenario','model','sim','k'))
    byk[,`:=`(diagnostic=d,target=t)];strata[[length(strata)+1L]]<-byk
    weighted<-byk[,.(n=sum(n[usable]),spearman=if(any(is.finite(spearman))) weighted.mean(spearman,n,na.rm=TRUE) else NA_real_,
      auc=if(any(is.finite(auc))) weighted.mean(auc,n,na.rm=TRUE) else NA_real_),by=.(scenario,model,sim)]
    weighted[,evaluation:='k_stratified']
    pooled<-association_stats(v,c('scenario','model','sim'))[,.(scenario,model,sim,n,spearman,auc)]
    pooled[,evaluation:='pooled']
    combined<-rbindlist(list(weighted,pooled));combined[,`:=`(diagnostic=d,target=t)]
    out[[length(out)+1L]]<-combined
  }
  list(per_sim=rbindlist(out),per_k=rbindlist(strata))
}
exploration_rules<-function() {
  r<-rbindlist(list(data.table(diagnostic='p_two',cutoff=c(.001,.0025,.005,.01,.02,.03,.05,.075,.10,.15,.20),rule_type='p_value'),
    CJ(diagnostic=c('ho_standardized','ho_uncertainty'),cutoff=seq(1,4,.25))[,rule_type:='absolute'],
    data.table(diagnostic='pre_rmse',cutoff=c(2.5,5,10,15,20),rule_type='top_percent')))
  r[,rule_id:=seq_len(.N)];r
}
keep_matrix_A<-function(z,rules,expected_half='A') {
  stopifnot(expected_half %in% c('A','B'),all(z$half==expected_half))
  available<-vapply(seq_len(nrow(rules)),function(j) is.finite(z[[rules$diagnostic[j]]]),logical(nrow(z)))
  pass<-vapply(seq_len(nrow(rules)),function(j) {
    value<-z[[rules$diagnostic[j]]];cut<-rules$cutoff[j]
    if(rules$rule_type[j]=='p_value') value>=cut else if(rules$rule_type[j]=='absolute') abs(value)<=cut else {
      # Exact nearest requested number, with ID as a fixed tie-breaker.
      bad<-rep(FALSE,nrow(z));idx<-which(is.finite(value));nflag<-floor(length(idx)*cut/100+.5)
      if(nflag>0) bad[idx[order(-value[idx],z$id[idx])][seq_len(nflag)]]<-TRUE
      !bad
    }
  },logical(nrow(z)))
  exclude<-pass & available;exclude[is.na(exclude)]<-FALSE
  keep<-pass | !available;keep[is.na(keep)]<-FALSE
  list(keep=cbind(keep,exclude),available=cbind(available,available),
    keys=rbindlist(list(copy(rules)[,missing_handling:='keep'],copy(rules)[,missing_handling:='exclude'])))
}
random_reference<-function(y,seed,draws=100L) {
  n<-length(y);if(!n || mean(y)<=0) return(NULL)
  samples<-with_backdate_seed(seed,replicate(draws,y[sample.int(n,n,replace=FALSE)]))
  if(n==1L) samples<-matrix(samples,nrow=1L)
  sums<-apply(samples,2,cumsum);if(n==1L) sums<-matrix(sums,nrow=1L)
  list(sum=sum(y),n=n,baseline=mean(y),prefix=sums)
}
random_reduction<-function(ref,n_excluded) {
  if(is.null(ref)||n_excluded>=ref$n) return(c(mean=NA_real_,lo=NA_real_,hi=NA_real_))
  prefix<-if(n_excluded) ref$prefix[n_excluded,] else rep(0,ncol(ref$prefix))
  reductions<-1-(ref$sum-prefix)/(ref$n-n_excluded)/ref$baseline
  c(mean=mean(reductions),lo=quantile(reductions,.025,names=FALSE),hi=quantile(reductions,.975,names=FALSE))
}
retention_group_A<-function(z,rules,common_keep=NULL,seed=2032032L,targets=exploration_targets,random_baseline=TRUE,expected_half='A') {
  screens<-keep_matrix_A(z,rules,expected_half);retention<-if(is.null(common_keep)) screens$keep else common_keep
  keys<-screens$keys;rows<-list()
  for(t in targets) {
    eligible<-is.finite(z[[t]])
    y<-z[[t]][eligible];keep<-retention[eligible,,drop=FALSE];n<-length(y)
    nr<-colSums(keep);excluded<-n-nr
    base<-if(n) mean(y) else NA_real_;retained<-if(n) colSums(keep*y)/nr else rep(NA_real_,ncol(keep))
    cov<-z[[coverage_for_target(t)]][eligible];good_cov<-is.finite(cov)
    cov_n<-if(any(good_cov)) colSums(keep[good_cov,,drop=FALSE]) else rep(0,ncol(keep))
    cov_ret<-if(any(good_cov)) colSums(keep[good_cov,,drop=FALSE]*cov[good_cov])/cov_n else rep(NA_real_,ncol(keep))
    baseline_cov<-mean_finite(cov)
    ref<-if(random_baseline) random_reference(y,seed+match(t,exploration_targets)) else NULL
    randoms<-t(vapply(excluded,function(e) random_reduction(ref,e),numeric(3)))
    a<-copy(keys)
    a[,`:=`(target=t,n_total=nrow(z),n_eligible=n,n_retained=nr,missing_share=NA_real_,
      excluded_share=if(n) excluded/n else NA_real_,retained_share=if(n) nr/n else NA_real_,
      baseline_error=base,retained_error=retained,reduction=if(is.finite(base)&&base>0) 1-retained/base else NA_real_,
      baseline_coverage=baseline_cov,retained_coverage=cov_ret,n_retained_coverage=cov_n,
      coverage_missing_share=1-cov_n/nr,
      random_mean=randoms[,1],random_lo=randoms[,2],random_hi=randoms[,3])]
    # Per-rule diagnostic missing share: denominator all pseudo-treated rows, not just targets.
    a[,missing_share:=colMeans(!screens$available)]
    rows[[length(rows)+1L]]<-a
  }
  rbindlist(rows)
}
exploration_retention<-function(z,rules) {
  stopifnot(all(z$half=='A'))
  groups<-unique(z[,.(scenario,model,sim)]);rows<-list();flags<-list()
  for(g in seq_len(nrow(groups))) {
    key<-groups[g];v<-z[scenario==key$scenario & model==key$model & sim==key$sim]
    seed<-2032032L+match(key$scenario,c('MAR','Latent level','Latent trend'))*100000L+key$sim*100L
    a<-retention_group_A(v,rules,seed=seed);a[,`:=`(scenario=key$scenario,model=key$model,sim=key$sim)];rows[[g]]<-a
    screen<-keep_matrix_A(v,rules)
    flags[[g]]<-rbindlist(lapply(sort(unique(v$stratum)),function(s) {
      sub<-v$stratum==s;a<-copy(screen$keys)
      a[,`:=`(stratum=s,n_counties=sum(sub),flag_rate=colMeans(!screen$keep[sub,,drop=FALSE]),
        missing_share=colMeans(!screen$available[sub,,drop=FALSE]),scenario=key$scenario,model=key$model,sim=key$sim)];a
    }))
  }
  list(per_sim=rbindlist(rows),size=rbindlist(flags))
}
common_sample_A<-function(z,rules) {
  required<-c('TWFE','ETWFE','GB-Ens','RW-IID');  # = s3_models (R/screening/s3_models.R)
  keys<-unique(z[,.(scenario,sim)]);rows<-list()
  for(g in seq_len(nrow(keys))) {
    key<-keys[g];v<-z[scenario==key$scenario & sim==key$sim];ids<-sort(unique(v$id))
    screen_by_model<-lapply(required,function(m) {
      d<-v[model==m];d<-d[match(ids,id)];stopifnot(!anyNA(d$id));keep_matrix_A(d,rules)$keep
    })
    joint<-Reduce(`&`,screen_by_model)
    # Pair target availability too, so each model is compared on identical county sets.
    for(m in required) {
      d<-v[model==m][match(ids,id)]
      for(t in exploration_targets) {
        valid<-v[,.(ok=.N==length(required) && all(is.finite(get(t)))),by=id]
        joint_target<-joint & valid$ok[match(ids,valid$id)]
        # Compute all targets once, then select this target's result; separate target availability is enforced.
        a<-retention_group_A(d,rules,common_keep=joint_target,seed=2032033L+key$sim,targets=t,random_baseline=FALSE)
        a[,`:=`(scenario=key$scenario,model=m,sim=key$sim)];rows[[length(rows)+1L]]<-a
      }
    }
  };rbindlist(rows)
}
rank_retention_A<-function(summary,unscreened) {
  ranked<-copy(summary)[target %in% c('err_raw','err_z','h0_raw','h0_z')]
  ranked[,retained_rank:=frank(retained_error,ties.method='min',na.last='keep'),
    by=.(scenario,diagnostic,cutoff,rule_type,missing_handling,target)]
  ref<-copy(unscreened)[target %in% c('err_raw','err_z','h0_raw','h0_z')]
  ref[,unscreened_rank:=frank(baseline_error,ties.method='min',na.last='keep'),by=.(scenario,target)]
  ranked<-merge(ranked,ref[,.(scenario,model,target,unscreened_rank)],by=c('scenario','model','target'),all.x=TRUE)
  ranked[,ranking_changed:=retained_rank!=unscreened_rank];ranked
}

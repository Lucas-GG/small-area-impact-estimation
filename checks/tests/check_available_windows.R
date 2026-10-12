# Hand-built acceptance checks; no production models and no B data.
library(data.table)
source('R/backdate_helpers.R')
source('R/heldout_analysis.R')
source('R/heldout_exploration.R')
a<-post_window_summary(c(3,0,NA,NA,NA),c(2,1,2,2,2),0:4)
b<-post_window_summary(5,1,0)
none<-post_window_summary(NA_real_,2,0)
stopifnot(none$k==0,typeof(a$sum_y)==typeof(none$sum_y))
stopifnot(a$k==2,a$err_raw==0,a$err_z==0,b$k==1,b$err_raw==4,b$err_z==4)
# Complete sensitivity allows zero yearly means, matching the old saved sums.
c<-post_window_summary(c(1,2,3,4,5),c(0,2,3,4,5),0:4)
stopifnot(c$k==4,c$k5_raw==1,abs(c$k5_z-1/sqrt(14))<1e-12)
stopifnot(inherits(try(derive_exploration_A(data.table(half='B')),silent=TRUE),'try-error'))
# Targets use matching draw indices across years, and noisy count intervals are reproducible.
draws<-matrix(c(2,1,3,1,1,1,2,2),nrow=2)
c<-post_window_summary(c(3,0),c(2,1),0:1,draws,100)
stopifnot(isTRUE(all.equal(c,post_window_summary(c(3,0),c(2,1),0:1,draws,100))),
 c$mu_post_var==var(colSums(draws)),is.finite(c$coverage))
# Both policies have the same finite-target baseline but different missing retention.
z<-data.table(half='A',id=1:10,pre_rmse=c(1:9,NA),p_two=c(seq(.01,.09,.01),NA),
 ho_standardized=c(1:9,NA),ho_uncertainty=c(1:9,NA),k=2,
 err_raw=1:10,err_z=1:10,coverage=1,h0_raw=1:10,h0_z=1:10,h0_coverage=1,
 k5_raw=NA_real_,k5_z=NA_real_,k5_coverage=NA_real_,average_raw=1:10,average_z=1:10)
rules<-data.table(diagnostic='ho_standardized',cutoff=5,rule_type='absolute',rule_id=1L)
ret<-retention_group_A(z,rules,targets='err_raw')
stopifnot(ret[missing_handling=='keep',n_retained]==6,ret[missing_handling=='exclude',n_retained]==5,
 all(ret$baseline_error==5.5),all(ret$missing_share==.1))
v<-data.table(scenario='x',model='x',sim=1,k=2,id=1:10,score=1:10,error=1:10)
s<-association_stats(v,c('scenario','model','sim','k'))
stopifnot(s$n==10,s$n_top==1,s$auc==1,abs(s$spearman-1)<1e-12)
cat('Cancellation, k=1, k5, no-B assertion, predictive coverage, missing policies, and stratified statistics passed.\n')

# Common screening pairs finite target availability across all six models.
z[,`:=`(scenario='MAR',sim=1L,average_k5_raw=NA_real_,average_k5_z=NA_real_)]
models<-c('TWFE','ETWFE','RW','RW-IID','RW-ST','GB-Ens')
all_models<-rbindlist(lapply(models,function(m) copy(z)[,model:=m]))
all_models[model=='TWFE' & id==1,err_raw:=NA_real_]
common<-common_sample_A(all_models,rules)[target=='err_raw']
stopifnot(uniqueN(common[missing_handling=='keep',n_retained])==1,
 common[missing_handling=='keep',unique(n_retained)]==5,
 common[missing_handling=='exclude',unique(n_retained)]==4)
# MC standard error uses the finite simulation count.
m<-mc_summary(data.table(group='test',sim=1:3,value=c(1,3,NA_real_)),'value','group')
stopifnot(m$value==2,m$value_mc_se==1,m$value_n_sims==2)
cat('Six-model common-sample pairing and finite-simulation MC standard errors passed.\n')

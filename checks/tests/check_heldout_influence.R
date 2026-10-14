library(data.table)
source('R/heldout_analysis.R');source('R/heldout_aggregate.R');source('R/heldout_influence.R')
# Diagnostic-only mode works without post-period columns and leaves existing targets untouched.
z<-data.table(id=1:3,half='A',y_ho=c(3,4,NA),mu_ho_mean=1,mu_ho_var=0,ppp_ho_upper=.5,ppp_ho_lower=.5,pre_rmse=1)
y<-derive_heldout(z,diagnostics_only=TRUE);stopifnot(!any(c('realized_raw','realized_standardized') %in% names(y)),identical(y$ho_standardized,c(2,3,NA_real_)))
source('R/heldout_exploration.R');rule<-data.table(diagnostic='ho_standardized',cutoff=2,rule_type='absolute',rule_id=1)
stopifnot(identical(keep_matrix_A(y,rule)$keep[,1],c(TRUE,FALSE,TRUE)))
# Exact signed decomposition can have shares greater than 100% and negative cancellation.
e<-matrix(c(2,-1,4,0),2,2,byrow=TRUE);n<-rowSums(matrix(TRUE,2,2));E<-rowSums(e)/n
C<-colMeans(e*(E/n));stopifnot(abs(sum(C)-mean(E^2))<1e-12)
e2<-c(2,-1);E2<-mean(e2);C2<-E2*e2/2;stopifnot(all(C2/E2^2==c(2,-1)))
# A county absent in a replication contributes zero; its deletion leaves that replication unchanged.
e<-matrix(c(3,0,2,-1),2,2,byrow=TRUE);present<-matrix(c(TRUE,FALSE,TRUE,TRUE),2,2,byrow=TRUE)
N<-rowSums(present);E<-rowSums(e)/N;loo<-(rowSums(e)-e[,2])/(N-present[,2]);stopifnot(loo[1]==E[1],loo[2]==2)
# Independent re-sum of exported county contributions and original-file preservation.
o<-'_output/heldout/influence';read<-function(n)fread(file.path(o,paste0('influence_',n,'.csv')))
c<-read('county_MSE_contributions');s<-c[,.(sum_C=sum(C),MSE=MSE[1]),by=.(branch,half,scenario,model,diagnostic,cutoff)]
stopifnot(nrow(s)==288,all(abs(s$sum_C-s$MSE)<1e-10))
r<-read('reproduction_checks');stopifnot(nrow(r)==432,all(r$within_precision))
flags<-read('county_errors_flags');a<-readRDS('data/derived/heldout/county_assignment.rds')
stopifnot(all(flags[half=='A',id] %in% a$id[a$half=='A']),all(flags[half=='B',id] %in% a$id[a$half=='B']))
for(h in c('A','B'))stopifnot(all(flags[half==h,is.finite(e) & k>=1]))
l<-read('input_paths');stopifnot(!any(l$access=='read' & grepl('/data/simulation/',l$path)))
old<-readRDS(file.path(o,'influence_original_hashes.rds'));stopifnot(identical(old,tools::md5sum(names(old))))
# Confirm k=0 deletion tables are the full baseline, not a re-screened arm.
b<-read('full_sample');t<-read('topk_out')[k_removed==0]
k<-c('branch','half','scenario','model','diagnostic','cutoff');j<-merge(b,t,by=k,suffixes=c('_full','_top0'))
stopifnot(all(abs(j$rmse_full-j$rmse_top0)<1e-12),all(abs(j$rmse_ratio_full-j$rmse_ratio_top0)<1e-12))
cat('Diagnostic-only reuse, strict/missing-kept flags, signed/negative decomposition, absent-county deletion, 288 MSE identities, 432 reproduction comparisons, half/hash guards and unchanged original outputs agree.\n')

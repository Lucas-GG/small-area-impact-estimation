library(data.table)
source('R/heldout_aggregate.R')
source('R/heldout_exploration.R')
source('R/backdate_helpers.R')
source('R/heldout_twfe_population.R')

# Eligibility follows saved finite targets, including a zero-history county with a target.
z<-data.table(scenario='MAR',sim=1,id=1:2,model='TWFE',k=c(1,0),sum_y=c(5,NA),sum_mu=c(4,NA),positive_pre=FALSE)
p<-data.table(sim=1,id=1:2,k_check=c(1,0),y_check=c(5,NA),mu_check=c(4,NA))
e<-twfe_common_membership(z,list(MAR=p))
stopifnot(identical(e$id,1L))
bad<-copy(p);bad[2,`:=`(k_check=1,y_check=2,mu_check=1)]
stopifnot(inherits(try(twfe_common_membership(z,list(MAR=bad)),silent=TRUE),'try-error'))

# Preserve the original percentile reference, then evaluate its E subset.
v<-data.table(half='A',id=1:4,pre_rmse=c(1,2,100,200))
rule<-data.table(diagnostic='pre_rmse',cutoff=50,rule_type='top_percent',rule_id=1)
full<-keep_matrix_A(v,rule)$keep[1:2,1]
restricted<-keep_matrix_A(v[id<=2],rule)$keep[,1]
stopifnot(sum(full)==2,sum(restricted)==1,all(full==keep_matrix_A(v,rule)$keep[1:2,1]),sum(full)<=sum(keep_matrix_A(v,rule)$keep[,1]))

# Aggregation keeps signs, respects retained draw indices, and marks empty sets unavailable.
u<-data.table(half='A',id=1:2,sum_y=c(15,8),sum_mu=c(10,10),signed_error=c(5,-2))
u[,draws:=list(as.numeric(1:200)/10,as.numeric(200:1)/10)]
keys<-data.table(diagnostic=c('unscreened','one','empty'),cutoff=NA_real_,rule_type='baseline',missing_handling='none')
k<-cbind(c(TRUE,TRUE),c(TRUE,FALSE),c(FALSE,FALSE))
a<-twfe_common_group(u,keys,k)
stopifnot(a$error[1]==1.5,a$error[2]==5,a$n_retained[1]==2,a$n_retained[2]==1,!is.finite(a$error[3]))
expected<-quantile((u$draws[[1]]+u$draws[[2]])/2,c(.025,.975),names=FALSE,type=7)
stopifnot(all(abs(c(a$lower[1],a$upper[1])-expected)<1e-12))
bad_u<-copy(u);bad_u[1,half:='B']
stopifnot(inherits(try(twfe_common_group(bad_u,keys,k),silent=TRUE),'try-error'))

# Random unscreened aggregate errors equal their paired baseline for every seeded draw.
acc<-twfe_common_random_accumulate(NULL,u,keys,a$n_retained,a$error[1],123L)
acc2<-twfe_common_random_accumulate(NULL,u,keys,a$n_retained,a$error[1],123L)
stopifnot(identical(acc,acc2),all(abs(sqrt(acc$sq[1,]/acc$base[1,])-1)<1e-12))

cat('Saved-target eligibility, ambiguity rejection, preserved-screen E-subset, signed aggregation, matching predictive indices, B-row rejection and seeded-random checks passed.\n')

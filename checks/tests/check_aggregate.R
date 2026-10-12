# Hand-built acceptance checks for signed aggregation, paired uncertainty and isolation.
library(data.table)
source('R/heldout_aggregate.R')
source('R/heldout_exploration.R')
source('R/backdate_helpers.R')

# County absolute errors do not define the aggregate error: signs cancel.
y<-c(12,8);mu<-c(10,10)
stopifnot(mean(y-mu)==0,mean(abs(y-mu))==2)
stopifnot(aggregate_oracle(c(1,1,100))==c(TRUE,TRUE,FALSE))
# A single pass retains variance 5; iterating would remove it on the next pass.
stopifnot(identical(aggregate_oracle(c(1,1,5,12)),c(TRUE,TRUE,TRUE,FALSE)))

# Perfectly paired proportional errors have a constant RMSE ratio and zero paired SE.
b<-data.table(subgroup='all_A',scenario='MAR',model='test',sim=1:4,
  diagnostic='unscreened',cutoff=NA_real_,rule_type='baseline',missing_handling='none',
  error=c(1,-2,3,-4),error_z=c(1,-2,3,-4),observed_aggregate=10,coverage=1,width=2,n_retained=2,excluded_share=0)
r<-copy(b);r[,`:=`(diagnostic='p_two',cutoff=.01,rule_type='p_value',missing_handling='keep',error=error*.5,error_z=error_z*.5,n_retained=1,excluded_share=.5)]
s<-aggregate_summary(rbindlist(list(b,r)))
stopifnot(abs(s[diagnostic=='p_two',rmse_ratio]-.5)<1e-12,
  abs(s[diagnostic=='p_two',rmse_ratio_mc_se])<1e-12,
  abs(s[diagnostic=='unscreened',rmse]-sqrt(7.5))<1e-12,
  abs(s[diagnostic=='unscreened',bias]+.5)<1e-12,
  abs(s[diagnostic=='unscreened',relative_rmse]-sqrt(7.5)/10)<1e-12,
  abs(s[diagnostic=='p_two',relative_rmse]-sqrt(7.5)*.5/10)<1e-12)

# The visibility cap applies to every available replication, not only its average.
oracle<-copy(b);oracle[,`:=`(diagnostic='oracle',cutoff=2,rule_type='infeasible_LOSO_standardized_variance_single_pass',excluded_share=c(.1,.2,.3,.41))]
gated<-aggregate_summary(rbindlist(list(b,oracle)))
stopifnot(!gated[diagnostic=='oracle',oracle_main_visible])
oracle[4,excluded_share:=.4]
gated<-aggregate_summary(rbindlist(list(b,oracle)))
stopifnot(gated[diagnostic=='oracle',oracle_main_visible])

# Missing oracle replications use a matched baseline, without inventing a variance.
r2<-copy(r);r2[1,error:=NA_real_];r2[,oracle_missing_variance_count:=c(1,0,0,0)]
s2<-aggregate_summary(rbindlist(list(b,r2),fill=TRUE))
stopifnot(s2[diagnostic=='p_two',paired_n_sims]==3,
  s2[diagnostic=='p_two',n_oracle_unavailable]==1,
  abs(s2[diagnostic=='p_two',rmse_ratio]-.5)<1e-12)

# A matched manuscript window reproduces the aggregate exactly.
x<-data.table(subgroup='all_A',scenario='MAR',model='test',new_error=c(1,-2,3,-4),manuscript_available_error=c(1,-2,3,-4),manuscript_complete_error=c(1,-2,3,-4),n_available=4,n_complete=4)
c<-aggregate_manuscript_summary(x)
stopifnot(c$complete_window_agrees_within_mc_error,c$max_abs_replication_difference==0,c$complete_window_difference==0,c$complete_window_difference_mc_se==0)

# A single B row fails before any analysis, even when the assignment hash is valid.
obj<-readRDS('data/derived/heldout/half_A_available.rds')
aggregate_guard(obj,getwd())
bad<-list(rows=copy(obj$rows[1:2]),assignment_md5=obj$assignment_md5)
bad$rows[1,half:='B']
stopifnot(inherits(try(aggregate_guard(bad,getwd()),silent=TRUE),'try-error'))
bad$rows[,half:='A'];bad$assignment_md5<-'wrong'
stopifnot(inherits(try(aggregate_guard(bad,getwd()),silent=TRUE),'try-error'))

# The shared seed applies to the whole yearly matrix, before usable-year selection.
lam<-matrix(c(2,3,4,5,6,7),nrow=3)
pred<-matrix(with_backdate_seed(77,rpois(length(lam),as.numeric(lam))),nrow=3)
county<-post_window_summary(c(3,4,5),c(2,3,4),0:2,lam,77)
q<-quantile(colSums(pred),c(.025,.975),names=FALSE,type=7)
stopifnot(all(abs(c(county$pi_lower,county$pi_upper)-q)<1e-12))

cat('Signed-error cancellation, paired-ratio SE, variance oracle, A-only/hash rejection and predictive-seed checks passed.\n')

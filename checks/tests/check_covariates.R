# Configuration and fixed-rule checks only: never fit a production or synthetic model.
library(data.table)
source('R/backdate_helpers.R')
source('R/heldout_analysis.R')
# Stubs let us inspect each specification without loading native fit engines.
impute_fe <- impute_inla <- impute_xgb <- function(...) stop('No fits in this test')
fit_pm <- function(...) "RW"
fit_me <- function(...) "RW-IID"
fit_sptiv <- function(...) "RW-ST"
fit_pm_c <- function(...) "RW covariates"
fit_me_c <- function(...) "RW-IID covariates"
fit_sptiv_c <- function(...) "RW-ST covariates"
s<-backdate_specs(TRUE); old<-backdate_specs()
stopifnot(identical(s$RW$args$fit_fun,fit_pm_c),
 identical(s[['RW-IID']]$args$fit_fun,fit_me_c),
 identical(s[['RW-ST']]$args$fit_fun,fit_sptiv_c),
 identical(s[['GB-Ens']]$args$predictors,c('uemp','pov','mhinc','aian')),
 identical(old[['GB-Ens']]$args,list()))
stopifnot(identical(deparse(s$TWFE$args$fm),deparse(y0~1+uemp+pov+mhinc+aian|i+year)),
 identical(deparse(s$ETWFE$args$fm),deparse(y0~uemp+pov+mhinc+aian|
 start_year[uemp+pov+mhinc+aian]+year[uemp+pov+mhinc+aian])))
p<-data.table(i=1:4,urb='6.NonCore',y=1:4,n=100,year=2000:2003,
 start_year=Inf,y0=1:4,y0_25=1:4,uemp=1:4,pov=2:5,mhinc=3:6,AIAN=4:7)
before<-copy(p);z<-backdate_covariate_panel(p)
stopifnot(identical(p,before),identical(z$y,p$y),identical(z$n,p$n),
 !any(c('y0','y0_25') %in% names(z)))
expected<-copy(before);setnames(expected,tolower(names(expected)))
expected[, `:=`(uemp=scale(uemp),pov=scale(pov),mhinc=scale(mhinc),aian=scale(aian))]
for(v in c('uemp','pov','mhinc','aian'))
 stopifnot(isTRUE(all.equal(z[[v]],expected[[v]])))
z<-data.table(sim=1,scenario='MAR',model=c('TWFE','RW','ETWFE'),id=1,x=2008,
 ho_standardized=c(1,3,NA_real_),p_two=c(.05,.01,NA_real_))
f<-fixed_rule_flags(z,list(diagnostic='ho_standardized',cutoff=2))
stopifnot(identical(f$keep,c(TRUE,FALSE,FALSE)),identical(f$diagnostic_available,c(TRUE,TRUE,FALSE)))
f<-fixed_rule_flags(z,list(diagnostic='p_two',cutoff=.05))
stopifnot(identical(f$keep,c(TRUE,FALSE,FALSE)))
cat('Covariate scaling, specifications, default compatibility and model-specific fixed rules passed.\n')

# Mock the GB-Ens plumbing to check that the leakage audit includes covariate predictors.
predictors<-c('uemp','pov','mhinc','aian')
add_y0_history<-function(dt) {dt[,history:=shift(y0)];dt}
make_y0_context<-function(dt) data.table(context=dt$y0)
materialize<-identity
get_y0_features<-function(dt) c('history','context')
get_hazard<-function(dt,engine,predictors) {
 stopifnot(identical(predictors,c('uemp','pov','mhinc','aian')));rep(.5,nrow(dt))
}
pre_xgb<-function(dt,predictors) {
 stopifnot(identical(predictors,c('uemp','pov','mhinc','aian')))
 list(mask=c(TRUE,FALSE,FALSE,FALSE),data=dt[,c('y0',predictors),with=FALSE])
}
z<-backdate_covariate_panel(p);z[,y0:=c(1,2,NA_real_,NA_real_)]
stopifnot(isTRUE(check_backdate_features(z,203L,predictors=predictors)))
cat('GB-Ens leakage audit propagates covariate predictors.\n')

library(data.table)
source('R/heldout_analysis.R');source('R/heldout_aggregate.R');source('R/heldout_estimate.R')
o<-'_output/heldout/with_covariates/confirmation'
x<-readRDS('data/derived/heldout/with_covariates/half_A.rds');aggregate_guard(x,normalizePath('.'),'A')
stopifnot(!file.exists('data/derived/heldout/with_covariates/half_B.rds'),!file.exists(file.path(o,'run_B_authorized')),
  !eval(formals(run_half_estimate)$with_covariates),eval(formals(aggregate_predictive_A)$file_prefix)=='new_imputation')
read<-function(n)fread(file.path(o,paste0('half_A_',n,'.csv')))
a<-read('availability');stopifnot(nrow(a)==18,all(a$simulations==200))
c<-read('population_checks');stopifnot(all(c$identical_ids & c$identical_observations & c$subset_of_every_model))
f<-read('flag_accounting');stopifnot(all(f$accounted_for))
s<-read('aggregate_summary');e<-read('estimates');stopifnot(nrow(e)==8,all(abs(e$absolute_ratio_difference)<1e-12),all(abs(e$relative_ratio_difference)<1e-12),all(abs(e$coverage_change_difference)<1e-12))
t<-read('twfe_own_E_checks')
for(v in c('rmse','rmse_mc_se','relative_rmse','relative_rmse_mc_se','rmse_ratio','rmse_ratio_mc_se','coverage','coverage_mc_se','width','n_retained'))stopifnot(all(abs(t[[paste0(v,'_own')]]-t[[paste0(v,'_E')]])<1e-10,na.rm=TRUE))
l<-read('input_paths');stopifnot(!any(grepl('half_B|/post_B/|/estimate_B/',l$path)),!any(l$access=='read' & grepl('/data/simulation/wocov/',l$path)))
d<-read('aggregate_per_simulation');eq<-d[population=='common',.(same_n=uniqueN(n_retained)==1,same_y=all(is.na(observed_aggregate))||max(observed_aggregate)-min(observed_aggregate)<1e-10),by=.(subgroup,scenario,sim,diagnostic,cutoff)]
stopifnot(all(eq$same_n & eq$same_y))
cat('Covariate A availability, fixed-rule summaries, TWFE own/E equality, common populations, unique/shared flags, shared defaults and no-covariate-B input checks agree.\n')

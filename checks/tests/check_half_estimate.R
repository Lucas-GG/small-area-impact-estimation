library(data.table)
source('R/heldout_analysis.R');source('R/heldout_estimate.R')
# Decision validation uses throw-away fixtures; the fixed decision is never edited.
d<-read_decision(normalizePath('.'));stopifnot(is.list(d$estimation_plan),!('confirmation_criteria' %in% names(d)))
fixture<-tempfile('estimate-validation-');dir.create(fixture);dir.create(file.path(fixture,'analysis'))
check_invalid<-function(x,pattern){yaml::write_yaml(x,file.path(fixture,'analysis/3.1a.4_heldout_decision.yml'));err<-tryCatch(read_decision(fixture),error=function(e)conditionMessage(e));stopifnot(is.character(err),grepl(pattern,err))}
x<-d;x$estimation_plan<-NULL;check_invalid(x,'estimation_plan')
x<-d;x$estimation_plan$summaries<-'';check_invalid(x,'estimation_plan')
x<-d;x$rule$exclude_if<-'>= 2';check_invalid(x,'strict')
x<-d;x$sensitivity_rules[[1]]$cutoff<-.05;check_invalid(x,'Sensitivity')
x<-d;x$half_a_reference<-list(bad='text');check_invalid(x,'half_a_reference')
# Identically dependent models must retain the single-model SE, not divide by sqrt(5).
models<-c('TWFE','ETWFE','GB-Ens','RW','RW-IID','RW-ST')
z<-CJ(model=models,sim=1:4,diagnostic=c('unscreened','ho_standardized'))
z[,`:=`(subgroup='noncore',population='own',scenario='MAR',cutoff=ifelse(diagnostic=='unscreened',NA,2),missing_handling=ifelse(diagnostic=='unscreened','none','keep'),
 error=ifelse(diagnostic=='unscreened',c(2,2,4,4)[sim],c(1,3,2,4)[sim]),observed_aggregate=ifelse(diagnostic=='unscreened',3,2),coverage=ifelse(diagnostic=='unscreened',c(1,0,1,0)[sim],c(1,1,0,1)[sim]))]
s<-estimate_gains(z);a<-s[scenario=='MAR' & group=='Five non-TWFE models'];b<-s[scenario=='MAR' & group=='TWFE']
r<-sqrt(30/40);influence<-(c(1,3,2,4)^2-r^2*c(2,2,4,4)^2)/(2*r*mean(c(2,2,4,4)^2))
stopifnot(abs(a$absolute_ratio-r)<1e-12,abs(a$relative_ratio-1.5*r)<1e-12,a$coverage_change==.25,
 abs(a$absolute_ratio_mc_se-sd(influence)/2)<1e-12,abs(a$absolute_ratio_mc_se-b$absolute_ratio_mc_se)<1e-12)
# Audit only A exports; no B execution is invoked by this check.
o<-'_output/heldout/confirmation';r<-fread(file.path(o,'half_A_reproduction_check.csv'));stopifnot(all(r$within_numerical_precision))
l<-fread(file.path(o,'half_A_input_paths.csv'));stopifnot(!any(grepl('half_B|/post_B/|/estimate_B/',l$path)),!any(l$access=='read' & grepl('/data/simulation/',l$path)))
stopifnot(!file.exists(file.path(o,'run_B_authorized')))
cat('Revised decision validation, strict threshold, fixed sensitivities, dependent-model summary SEs, A reproduction and no-B input audit agree.\n')

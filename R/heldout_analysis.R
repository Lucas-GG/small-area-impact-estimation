# Shared evaluation functions operate only on the half already loaded by the caller.
library(data.table)
diagnostic_labels <- c(pre_rmse='In-sample RMSE',pre_rms_pearson='In-sample RMS Pearson',
  ho_raw='Held-out raw error',ho_standardized='Held-out standardized error',
  ho_uncertainty='Held-out error with model uncertainty',p_two='Two-sided predictive mid-p')
load_heldout_half <- function(root, half, input_dir = file.path(root,"data/derived/heldout")) {
  p <- file.path(input_dir,paste0('half_',half,'.rds'))
  if (!file.exists(p)) stop('Run R/heldout_prepare.R first: missing ',p)
  wanted_half <- half
  obj <- readRDS(p) # This file physically contains only the requested half.
  stopifnot(all(obj$rows$half==half))
  assignment <- readRDS(file.path(root,'data/derived/heldout/county_assignment.rds'))
  stopifnot(identical(obj$assignment_md5,unname(tools::md5sum(file.path(root,'data/derived/heldout/county_assignment.rds')))),
    all(obj$rows$id %in% assignment[half==wanted_half,id]))
  # Apply the saved assignment at load time, in addition to physical file separation.
  obj$rows <- obj$rows[id %in% assignment[half==wanted_half,id]]
  obj
}
mean_finite <- function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
derive_heldout <- function(z) {
  z <- copy(z)
  z[, `:=`(ho_raw=y_ho-mu_ho_mean,
    ho_standardized=fifelse(mu_ho_mean>0,(y_ho-mu_ho_mean)/sqrt(mu_ho_mean),NA_real_),
    ho_uncertainty=fifelse(mu_ho_mean+mu_ho_var>0,(y_ho-mu_ho_mean)/sqrt(mu_ho_mean+mu_ho_var),NA_real_),
    p_two=pmin(1,2*pmin(ppp_ho_upper,ppp_ho_lower)),
    realized_raw=abs(y_post-mu_post_mean),
    realized_standardized=fifelse(mu_post_mean>0,abs(y_post-mu_post_mean)/sqrt(mu_post_mean),NA_real_))]
  # Average absolute errors, separately by scenario/model/county, within this county half.
  averages <- z[,.(average_raw=mean_finite(realized_raw),average_standardized=mean_finite(realized_standardized),
    n_average_raw=sum(is.finite(realized_raw)),n_average_standardized=sum(is.finite(realized_standardized)),
    n_pseudo_treated=.N),by=.(scenario,model,id)]
  list(rows=merge(z,averages,by=c('scenario','model','id'),all.x=TRUE),averages=averages)
}
auc_rank <- function(score, positive) {
  np <- sum(positive); nn <- sum(!positive)
  if(np==0 || nn==0) return(NA_real_)
  (sum(rank(score,ties.method='average')[positive])-np*(np+1)/2)/(np*nn)
}
evaluate_diagnostics <- function(z, diagnostics=names(diagnostic_labels)) {
  targets <- c('realized_raw','realized_standardized','average_raw','average_standardized')
  rbindlist(lapply(diagnostics,function(d) rbindlist(lapply(targets,function(t) {
    v <- z[,.(sim,scenario,model,id,diagnostic=d,target=t,
      score=if(d=='p_two') 1-get(d) else abs(get(d)),error=get(t))]
    v[, {
      good <- is.finite(score)&is.finite(error); s<-score[good]; y<-error[good]
      # Define the target class before restricting to available diagnostic scores.
      observed <- error[is.finite(error)]
      boundary <- if(length(observed)) quantile(observed,.9,names=FALSE,type=7) else NA_real_
      top <- y >= boundary
      list(n_total=.N,n_valid=length(y),spearman=if(length(y)>2 && sd(s)>0 && sd(y)>0) cor(s,y,method='spearman') else NA_real_,
        auc=if(length(y)) auc_rank(s,top) else NA_real_,n_top=sum(top))
    },by=.(sim,scenario,model,diagnostic,target)]
  }))))
}
cutoff_performance <- function(z, rules) {
  rbindlist(lapply(seq_len(nrow(rules)),function(j) {
    d<-rules$diagnostic[j]; threshold<-rules$cutoff[j]
    v<-z[,.(sim,scenario,model,value=get(d),realized_raw,realized_standardized,average_raw,average_standardized)]
    rbindlist(lapply(c('realized_raw','realized_standardized','average_raw','average_standardized'),function(t) {
      v[, {
        good<-is.finite(value)&is.finite(get(t)); y<-get(t)[good]; val<-value[good]
        excluded<-if(d=='p_two') val<threshold else abs(val)>threshold
        base<-mean_finite(y); retained<-mean_finite(y[!excluded])
        list(diagnostic=d,cutoff=threshold,target=t,n_total=.N,n_eligible=length(y),
          unavailable_share=1-length(y)/.N,excluded_share=if(length(y)) mean(excluded) else NA_real_,
          baseline_error=base,retained_error=retained,reduction=if(is.finite(base)&&base>0) 1-retained/base else NA_real_)
      },by=.(sim,scenario,model)]
    }))
  }))
}
summarize_metrics <- function(x, columns, groups) {
  x[, c(list(simulations=.N),lapply(.SD,mean_finite)),by=groups,.SDcols=columns]
}
read_decision <- function(root) {
  p<-file.path(root,'analysis/3.1a.4_heldout_decision.yml')
  if(!file.exists(p)) stop('Fixed-rule estimate blocked: create analysis/3.1a.4_heldout_decision.yml by hand after fixing the half-A choice.')
  d<-yaml::read_yaml(p,readLines.warn=FALSE)
  if(!all(c('diagnostic','cutoff','target_definition','missing_diagnostic_handling','date','rationale') %in% names(d))) stop('Decision needs diagnostic, cutoff, target_definition, missing_diagnostic_handling, date and rationale.')
  if(any(vapply(d[c('diagnostic','cutoff','target_definition','missing_diagnostic_handling','date','rationale')],is.null,logical(1))))
    stop('The decision template is blank; record all six required fields by hand before the fixed-rule estimate.')
  stopifnot(length(d$diagnostic)==1,d$diagnostic %in% names(diagnostic_labels),
    is.numeric(d$cutoff),length(d$cutoff)==1,is.finite(d$cutoff),d$cutoff>0,
    length(d$rationale)==1,nzchar(trimws(d$rationale)),length(d$date)==1,
    !is.na(as.Date(as.character(d$date))))
  if(d$diagnostic=='p_two' && d$cutoff>=1) stop('p_two cutoff must be below 1')
  stopifnot(is.character(d$target_definition),length(d$target_definition)==1,nzchar(trimws(d$target_definition)),
    length(d$missing_diagnostic_handling)==1,d$missing_diagnostic_handling %in% c('keep','exclude'))
  required<-c('rule','primary_subgroup','primary_population','secondary_populations','sensitivity_rules','estimation_plan','half_a_reference')
  missing<-required[!required %in% names(d)]
  if(length(missing))stop('Decision missing required fields: ',paste(missing,collapse=', '))
  text<-function(x)is.character(x)&&length(x)==1L&&!is.na(x)&&nzchar(trimws(x))
  if(!is.list(d$rule)||!all(c('score','exclude_if','applied') %in% names(d$rule))||!all(vapply(d$rule[c('score','exclude_if','applied')],text,logical(1))))stop('Malformed rule: score, exclude_if and applied must be nonempty strings.')
  if(d$diagnostic!='ho_standardized'||d$cutoff!=2||gsub(' ','',d$rule$score)!='abs(ho_standardized)'||gsub(' ','',d$rule$exclude_if)!='>2')stop('Decision must record the strict abs(ho_standardized) > 2 rule.')
  if(!identical(d$primary_subgroup,'noncore')||!identical(d$primary_population,'own')||d$missing_diagnostic_handling!='keep')stop('Malformed primary population/subgroup or missing policy.')
  if(!setequal(unlist(d$secondary_populations),c('all_half','twfe_eligible_E','common_screened'))||length(unlist(d$secondary_populations))!=3)stop('Malformed secondary_populations.')
  if(!is.list(d$sensitivity_rules)||length(d$sensitivity_rules)!=2)stop('Two pre-specified sensitivity_rules are required.')
  sens<-d$sensitivity_rules
  valid_sens<-vapply(sens,function(x)is.list(x)&&text(x$diagnostic)&&is.numeric(x$cutoff)&&length(x$cutoff)==1&&is.finite(x$cutoff),logical(1))
  if(!all(valid_sens)||!setequal(vapply(sens,`[[`,character(1),'diagnostic'),c('p_two','pre_rmse')))stop('Malformed sensitivity_rules.')
  for(x in sens)if((x$diagnostic=='p_two'&&x$cutoff!=.10)||(x$diagnostic=='pre_rmse'&&(x$cutoff!=10||!identical(x$rule_type,'top_percent'))))stop('Sensitivity rules must be p_two < 0.10 and pre_rmse top 10%.')
  plan<-c('primary_estimates','summaries','secondary','comparison_with_a')
  if(!is.list(d$estimation_plan)||!all(plan %in% names(d$estimation_plan))||!all(vapply(d$estimation_plan[plan],text,logical(1))))stop('Malformed estimation_plan: require nonempty primary_estimates, summaries, secondary, comparison_with_a.')
  if(!is.list(d$half_a_reference)||!length(d$half_a_reference)||!is.numeric(unlist(d$half_a_reference))||any(!is.finite(unlist(d$half_a_reference))))stop('Malformed half_a_reference: numeric reference values required.')
  d$md5<-unname(tools::md5sum(p)); d
}

# Fixed-rule flags are model-specific; unavailable diagnostics cannot qualify a county.
fixed_rule_flags <- function(rows, decision) {
  z <- copy(rows)
  value <- z[[decision$diagnostic]]
  stopifnot(!is.null(value))
  z[, diagnostic_available := is.finite(value)]
  z[, keep := diagnostic_available & if(decision$diagnostic=='p_two')
    value >= decision$cutoff else abs(value) <= decision$cutoff]
  z[, .(sim, scenario, model, id, x, diagnostic_available, keep)]
}

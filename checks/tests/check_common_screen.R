library(data.table)
source('R/heldout_aggregate.R')
root<-normalizePath('.')
obj<-readRDS('data/derived/heldout/half_A_available.rds');aggregate_guard(obj,root)
out<-'_output/analysis/output/3.1a/common_screen'
read<-function(n)fread(file.path(out,paste0('common_screen_',n,'.csv')))
d<-read('per_simulation');membership<-read('membership');ck<-read('subset_checks')
stopifnot(nrow(d)==25200,nrow(ck)==nrow(d),all(membership$half=='A'),all(membership$id %in% obj$rows$id),
  !anyDuplicated(membership[,.(scenario,sim,diagnostic,cutoff,missing_handling,id)]),
  all(ck$n_common<=ck$n_own),all(ck$n_own==ck$n_retained),all(ck$subset_of_own))
groups<-c('scenario','sim','diagnostic','cutoff','missing_handling')
ids<-membership[order(id),.(retained_id_key=paste(id,collapse=','),n_members=.N),by=groups]
check_ids<-merge(ck,ids,by=groups,suffixes=c('_reported','_membership'),all.x=TRUE)
check_ids[is.na(n_members),`:=`(n_members=0L,retained_id_key_membership='')]
stopifnot(all(check_ids$n_common==check_ids$n_members),all(check_ids$retained_id_key_reported==check_ids$retained_id_key_membership))
e<-read('equal_population_checks');stopifnot(nrow(e)==4200,all(e$n_models==6),all(e$identical_ids & e$same_count & e$same_observed))
c<-read('contributions_per_simulation')
accounting<-c[,.(valid=sum(n_unique_flags)+n_shared_flags[1]==n_excluded[1],bounded=all(n_unique_flags<=n_model_flags)),by=groups]
stopifnot(all(accounting$valid & accounting$bounded))
b<-read('unscreened_check');stopifnot(nrow(b)==18)
metrics<-sub('_common$','',names(b)[endsWith(names(b),'_common')])
for(nm in metrics) if(is.numeric(b[[paste0(nm,'_common')]])) {
  a<-b[[paste0(nm,'_common')]];z<-b[[paste0(nm,'_10_1')]]
  stopifnot(identical(is.na(a),is.na(z)),all(abs(a-z)<1e-10,na.rm=TRUE))
}
old<-readRDS('.quarto/common-screen/original-md5.rds');stopifnot(identical(old,tools::md5sum(names(old))))
aggregate_guard(obj,root)
cat('Common membership IDs, subset bounds, own-mask reproduction, identical observed aggregates, unique/shared accounting, unscreened reproduction, A-only/hash and original-export preservation checks passed.\n')

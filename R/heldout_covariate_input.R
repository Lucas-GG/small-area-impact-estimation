# Prepare only the requested half; never fit models or construct the other half.
prepare_covariate_selected_half <- function(root,half,output,read_input=readRDS) {
  stopifnot(half %in% c('A','B'))
  assignment_path<-file.path(root,'data/derived/heldout/county_assignment.rds')
  wanted_half<-half
  assignment<-as.data.table(read_input(assignment_path));ids<-assignment[half==wanted_half,id]
  paths<-list.files(file.path(root,'data/fits/backdate_wcov'),'^summary_.*[.]rds$',full.names=TRUE)
  if(length(paths)!=18L)stop('Need all 18 saved covariate checkpoints; found ',length(paths),'. No fit is run.')
  pieces<-lapply(paths,function(p) {
    # Serialized checkpoint is read once; select the requested IDs before all calculations.
    z<-as.data.table(read_input(p));z<-z[id %in% ids]
    stopifnot(nrow(z)>0L);z
  })
  x<-rbindlist(pieces,fill=TRUE);stopifnot(!anyDuplicated(x[,.(sim,model,scenario,id)]))
  if(!'tail_convention' %in% names(x))x[,tail_convention:=NA_character_]
  if(!'ppp_ho_upper' %in% names(x))x[,ppp_ho_upper:=NA_real_]
  if(!'ppp_ho_lower' %in% names(x))x[,ppp_ho_lower:=NA_real_]
  legacy<-is.na(x$tail_convention)|!grepl('^mid-p',x$tail_convention)
  both<-is.finite(x$ppp_ho_upper)&is.finite(x$ppp_ho_lower);recover<-legacy & both
  equality<-x$ppp_ho_upper[recover]+x$ppp_ho_lower[recover]-1
  stopifnot(all(equality>= -1e-12 & equality<=1+1e-12))
  x[recover,`:=`(ppp_ho_upper=ppp_ho_upper-.5*pmax(0,equality),ppp_ho_lower=ppp_ho_lower-.5*pmax(0,equality),tail_convention='mid-p recovered from inclusive pair')]
  x[legacy & !both,`:=`(ppp_ho_upper=NA_real_,ppp_ho_lower=NA_real_,tail_convention='unavailable: legacy pair absent')]
  x[,intersect(c('ppp_ho','ppp_ho_two_sided'),names(x)):=NULL]
  x<-merge(x,assignment[id %in% ids],by='id',all.x=TRUE);stopifnot(all(x$half==half))
  manifest<-data.table(file=paths,bytes=file.info(paths)$size,modified=as.character(file.info(paths)$mtime))
  dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(rows=x,manifest=manifest,prepared=Sys.time(),expected_groups=18L,assignment_md5=unname(tools::md5sum(assignment_path))),output)
  invisible(output)
}

# Gate before any revision performance computation. Audit every observed horizon.
library(data.table)
revision_gap_rows <- function(x) {
  stopifnot(all(c('dgp','nc','model','r','i','s','y0_hat')%in%names(x)))
  x[,{
    years<-as.integer(s)
    present<-unique(years[years>=0 & years<=4 & is.finite(y0_hat)])
    n_through<-vapply(years,function(h)sum(present<=h),integer(1))
    missing_prefix<-vapply(years,function(h)length(setdiff(0:h,present)),integer(1))
    leading<-if(length(present))min(present) else NA_integer_
    list(horizon=years,n_years_through_horizon=n_through,expected_years=years+1L,
      missing_prefix=missing_prefix,missing_leading=if(is.na(leading))years+1L else pmin(leading,years+1L),
      original_cumulative_rows=seq_len(.N),out_of_order=years!=sort(years),
      duplicate_year=duplicated(years)|duplicated(years,fromLast=TRUE),
      balanced_0_4=setequal(present,0:4)&&.N==5L)
  },by=.(dgp,nc,model,r,i)]
}
run_revision_gap_audit <- function(root=getwd(),negligible_share=.001) {
  source(file.path(root,'R/placebo_revision.R'));ctx<-revision_context(root)
  before<-revision_protect(root);all<-list();statuses<-list()
  cache_dir<-file.path(root,'.quarto/placebo-revision-audit');dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  for(spec_label in c('no_covariates','with_covariates')) {
    path<-file.path(root,'data/simulation',if(spec_label=='with_covariates')'sim_fits_wcov.rds' else 'sim_fits.rds')
    message('Gap audit: ',spec_label)
    raw<-as.data.table(ctx$read(path,'observed horizon cumulative-year gap audit'))
    x<-raw[,.(dgp,nc,model,r,i,s,y0_hat)];rm(raw);gc()
    z<-revision_gap_rows(x)
    z[,`:=`(specification=spec_label,short_cumulative=n_years_through_horizon<expected_years,
      internal_missing=pmax(0L,missing_prefix-missing_leading))]
    summary<-z[,.(county_simulation_rows=.N,short_cumulative=sum(short_cumulative),
      short_cumulative_share=mean(short_cumulative),missing_leading_rows=sum(missing_leading>0),
      internal_gap_rows=sum(internal_missing>0),internal_gap_share=mean(internal_missing>0),
      out_of_order_rows=sum(out_of_order),duplicate_year_rows=sum(duplicate_year),
      original_cumulative_count_mismatch=sum(original_cumulative_rows!=expected_years)),
      by=.(specification,model,dgp,nc,horizon)]
    fwrite(summary,file.path(ctx$out,paste0('revision_',spec_label,'_gap_audit.csv')))
    fwrite(z[short_cumulative|out_of_order|duplicate_year],file.path(ctx$out,paste0('revision_',spec_label,'_gap_rows.csv')))
    signature<-list(path=normalizePath(path),bytes=file.info(path)$size,mtime=as.character(file.info(path)$mtime))
    saveRDS(list(signature=signature,rows=z),file.path(cache_dir,paste0(spec_label,'_gap_audit.rds')))
    material<-any(summary$short_cumulative_share>negligible_share) || any(summary$out_of_order_rows>0|summary$duplicate_year_rows>0)
    statuses[[spec_label]]<-data.table(specification=spec_label,negligible_threshold=negligible_share,
      max_cell_gap_share=max(summary$short_cumulative_share),overall_gap_share=sum(summary$short_cumulative)/sum(summary$county_simulation_rows),
      gap_gate=if(material)'STOP — material gaps or ordering issue' else 'PROCEED — gaps negligible')
    all[[spec_label]]<-summary
    rm(x,z);gc();message('Gap audit complete: ',spec_label)
  }
  summary<-rbindlist(all);status<-rbindlist(statuses)
  fwrite(summary,file.path(ctx$out,'revision_gap_audit_all.csv'));fwrite(status,file.path(ctx$out,'revision_gap_gate.csv'))
  stopifnot(identical(before,revision_protect(root)))
  print(status)
  invisible(list(summary=summary,status=status,proceed=all(grepl('^PROCEED',status$gap_gate))))
}

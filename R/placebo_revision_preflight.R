library(data.table)
run_placebo_revision_preflight <- function(root=getwd()) {
  source(file.path(root,'R/placebo_revision.R'));ctx<-revision_context(root)
  paths<-file.path(root,'_output/placebo_revision/phase0',c('no_covariates','with_covariates'),'revision_diagnostic_join.csv')
  a<-ctx$csv(paths[1],'no-covariate diagnostics comparison');b<-ctx$csv(paths[2],'covariate diagnostics comparison')
  z<-merge(a[,.(model,scenario,sim,id,z_no=ho_standardized)],b[,.(model,scenario,sim,id,z_cov=ho_standardized)],by=c('model','scenario','sim','id'))
  comparison<-z[,.(n_pairs=.N,n_finite=sum(is.finite(z_no)&is.finite(z_cov)),
    correlation=cor(z_no,z_cov,use='complete.obs'),mean_absolute_difference=mean(abs(z_no-z_cov),na.rm=TRUE)),by=.(model,scenario)]
  stopifnot(all(comparison$mean_absolute_difference>0))
  fwrite(comparison,file.path(ctx$out,'revision_diagnostic_specification_comparison.csv'))
  pdf<-file.path(root,'docs/NIH Manuscript for HSOR_V3.pdf');ctx$record(pdf,'printed manuscript / supplement spot-check')
  extractor<-file.path(path.expand('~'),'.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/bin/pdftotext')
  text_path<-file.path(ctx$out,'revision_manuscript_extracted.txt')
  stopifnot(file.exists(extractor),system2(extractor,c('-layout',shQuote(pdf),shQuote(text_path)))==0L)
  ctx$record(text_path,'printed appendix rows extracted from actual manuscript PDF')
  pages<-strsplit(paste(readLines(text_path,warn=FALSE),collapse='\n'),'\f')[[1]]
  lines<-strsplit(paste(pages[42:47],collapse='\n'),'\n')[[1]]
  county_group<-'NonCore';scenario<-NA_character_;rows<-list()
  pattern<-'^\\s*(TWFE|ETWFE|RW-IID|RW-ST|GB-Ens|RW)\\s+(0|4)\\s+\\([0-9.]+\\)\\s*\\([0-9.]+\\)\\s+([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)'
  for(line_index in seq_along(lines)) {
    line<-lines[line_index]
    if(grepl('^Table A 1',trimws(line))){county_group<-'Core';scenario<-NA_character_}
    if(trimws(line)%in%c('MAR','Latent level','Latent trend'))scenario<-trimws(line)
    m<-regmatches(line,regexec(pattern,line,perl=TRUE))[[1]]
    if(length(m)) {
      means<-regmatches(gsub('−','-',lines[line_index-1]),gregexpr('-?[0-9]+\\.[0-9]+',gsub('−','-',lines[line_index-1])))[[1]]
      sds<-regmatches(line,gregexpr('\\([0-9.]+\\)',line))[[1]]
      stopifnot(length(means)==4,length(sds)==4)
      rows[[length(rows)+1]]<-data.table(nc=county_group,dgp=scenario,model=m[2],s=as.integer(m[3]),
        se_diff=as.numeric(m[4]),crps=as.numeric(m[5]),cvg=as.numeric(m[6]),intv=as.numeric(m[7]),
        diff=as.numeric(means[1]),abs_diff=as.numeric(means[2]),rel_diff=as.numeric(means[3]),abs_rel_diff=as.numeric(means[4]),
        diff_sd=as.numeric(gsub('[()]','',sds[1])),abs_diff_sd=as.numeric(gsub('[()]','',sds[2])),
        rel_diff_sd=as.numeric(gsub('[()]','',sds[3])),abs_rel_diff_sd=as.numeric(gsub('[()]','',sds[4])))
    }
  }
  printed<-rbindlist(rows);stopifnot(nrow(printed)==72,!anyNA(printed$dgp))
  ref<-ctx$csv(file.path(root,'_output/placebo_revision/phase0/no_covariates/revision_S3_historical_performance.csv'),'Phase 0 S3 vs actual printed manuscript')
  ids<-c('nc','dgp','model','s');metrics<-setdiff(names(printed),ids)
  check<-merge(melt(printed,id.vars=ids,measure.vars=metrics,value.name='printed'),
    melt(ref,id.vars=ids,measure.vars=metrics,value.name='reproduced'),by=c(ids,'variable'))
  check[,`:=`(printed_display=sprintf('%.3f',printed),reproduced_display=sprintf('%.3f',reproduced))]
  check[,`:=`(display_matches=printed_display==reproduced_display,difference_at_display=round(reproduced,3)-printed)]
  rounded_mse<-ref[,.(nc,dgp,model,s,sqrt_rounded_MSE_display=sprintf('%.3f',sqrt(round(sq_diff,3))))]
  check<-merge(check,rounded_mse,by=ids)
  check[,RMSE_matches_sqrt_rounded_MSE:=fifelse(variable=='se_diff',printed_display==sqrt_rounded_MSE_display,NA)]
  fwrite(check,file.path(ctx$out,'revision_manuscript_table_spotcheck.csv'))
  availability<-data.table(specification=rep(c('no_covariates','with_covariates'),each=2),sample=rep(c('S1_full','S3_historical'),2),
    printed_table_available=c(FALSE,TRUE,FALSE,FALSE),
    note=c('Full-sample results appear as figures A2/A3, not printed numerical tables.',
      'Printed A1/A2: RMSE differences from sqrt of rounded MSE; other measures checked directly.',
      'No full-sample covariate numerical table in supplied manuscript.',
      'Covariate results appear as figure A4, not a printed numerical table.'))
  fwrite(availability,file.path(ctx$out,'revision_manuscript_table_availability.csv'))
  print(check[,.(checked=.N,display_differences=sum(!display_matches)),by=variable])
  invisible(comparison)
}

# Presentation layer: read saved half-A results; never recalculate screening statistics.
find_presentation_root<-function() {
  p<-normalizePath(getwd())
  while(!file.exists(file.path(p,'_quarto.yml'))) {q<-dirname(p);if(q==p)stop('Repository root not found');p<-q};p
}
load_presentation<-function(root,refresh_bands=TRUE) {
  .libPaths(c(file.path(root,'_output/R-library'),.libPaths()))
  library(data.table);library(ggplot2)
  source(file.path(root,'R/screening/heldout_analysis.R'))
  obj<-readRDS(file.path(root,'data/derived/heldout/half_A_available.rds'))
  stopifnot(all(obj$rows$half=='A'),sum(obj$rows$half!='A')==0,
    identical(obj$assignment_md5,unname(tools::md5sum(file.path(root,'data/derived/heldout/county_assignment.rds')))))
  old<-file.path(root,'_output/heldout/exploration/available_window')
  out<-file.path(root,'_output/analysis/output/3.1a');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  csvs<-list.files(old,'[.]csv$',full.names=TRUE)
  stopifnot(length(csvs)>=22L,file.info(file.path(old,'row_counts.csv'))$mtime>=file.info(file.path(root,'data/derived/heldout/half_A_available.rds'))$mtime)
  # Copies are byte-identical; hard links avoid duplicating the large replication ledgers.
  for(p in csvs) {
    dest<-file.path(out,basename(p))
    if(!file.exists(dest)) {if(!file.link(p,dest)) stopifnot(file.copy(p,dest))}
    else if(file.info(dest)$size!=file.info(p)$size || file.info(dest)$mtime!=file.info(p)$mtime) stopifnot(file.copy(p,dest,overwrite=TRUE))
  }
  source(file.path(root,'R/screening/heldout_random_band.R'))
  if(refresh_bands) refresh_presentation_bands(root,obj,out)
  read<-function(n) fread(file.path(out,paste0(n,'.csv')))
  list(obj=obj,out=out,counts=read('row_counts'),kd=read('k_distribution'),assoc=read('associations_summary'),
    calibration=read('calibration_summary'),retained=read('retention_summary'),ranked=read('model_rankings'),
    size=read('flag_rates_by_size'),common=read('common_sample_comparison'),grid=read('random_band_plot_grid'),
    bands=read('random_band_width_comparison'),read=read)
}
model_order<-c('TWFE','ETWFE','GB-Ens','RW-IID')
scenario_order<-c('MAR','Latent level','Latent trend')
pretty_diagnostics<-c(pre_rmse='Pre-period RMSE',pre_rms_pearson='Pre-period Pearson',ho_raw='Held-out raw',
 ho_standardized='Held-out standardized',ho_uncertainty='Held-out uncertainty',p_two='Predictive mid-p')
diag_colors<-c(pre_rmse='#8E44AD',pre_rms_pearson='#B57BC5',ho_raw='#E69F00',
 ho_standardized='#C34A36',ho_uncertainty='#007A5E',p_two='#0072B2')
model_colors<-setNames(c('#3D405B','#0072B2','#E69F00','#6B8E23','#AA4499','#CC6677'),model_order)
# This is the only shortlist definition; it presents candidates, not a selected rule.
candidate_rules<-data.table::data.table(diagnostic=c(rep('p_two',3),rep('ho_standardized',2),'pre_rmse'),
 cutoff=c(.01,.05,.10,2,2.5,10),rule_type=c(rep('p_value',3),rep('absolute',2),'top_percent'),
 rule=c('Mid-p < 0.01','Mid-p < 0.05','Mid-p < 0.10','Held-out |z| > 2','Held-out |z| > 2.5','Pre-RMSE top 10%'),
 mark=c('p .01','p .05','p .10','|z| 2','|z| 2.5','top 10%'))
shortlist_rows<-function(x) merge(x,candidate_rules,by=c('diagnostic','cutoff','rule_type'))
ordered_plot<-function(x) {x<-data.table::copy(x);x[,model:=factor(model,levels=model_order)];if('scenario'%in%names(x))x[,scenario:=factor(scenario,levels=scenario_order)];x}
pct<-function(x) sprintf('%.2f%%',100*x)
num<-function(x) sprintf('%.2f',x)
fmt_int<-function(x) format(x,big.mark=',',scientific=FALSE,trim=TRUE)
presentation_theme<-function() ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
 legend.position='bottom',strip.text=ggplot2::element_text(face='bold'),plot.caption=ggplot2::element_text(hjust=0))
auc_plot<-function(p,target_name='err_z') {
  v<-ordered_plot(p$assoc[target==target_name & evaluation=='k_stratified'])
  v[,Diagnostic:=factor(pretty_diagnostics[diagnostic],levels=rev(unname(pretty_diagnostics)))]
  ggplot(v,aes(model,Diagnostic,fill=auc))+geom_tile(color='white')+geom_text(aes(label=sprintf('%.2f',auc)),size=3)+
   facet_wrap(~scenario,nrow=1)+scale_fill_gradient2(low='#A33A3A',mid='white',high='#167C67',midpoint=.5,limits=c(0,1),name='AUC')+
   labs(x=NULL,y=NULL)+presentation_theme()+theme(axis.text.x=element_text(angle=45,hjust=1))
}
retention_plot<-function(p,target_name='err_z',policy='keep') {
  v<-ordered_plot(p$retained[target==target_name & missing_handling==policy]);marks<-shortlist_rows(v)
  g<-ordered_plot(p$grid[target==target_name]);limits<-v[,.(xmax=max(excluded_share,na.rm=TRUE)),by=.(scenario,model)]
  g<-merge(g,limits,by=c('scenario','model'));g<-g[excluded_share<=xmax+.01]
  ggplot(v,aes(excluded_share,reduction))+
   geom_ribbon(data=g,aes(x=excluded_share,ymin=random_lo_corrected,ymax=random_hi_corrected),inherit.aes=FALSE,fill='grey65',alpha=.4)+
   geom_hline(yintercept=0,color='grey45',linetype='dashed')+
   geom_line(aes(color=diagnostic,group=diagnostic),linewidth=.6)+
   geom_point(data=marks,aes(color=diagnostic),size=1.6)+
   ggrepel::geom_text_repel(data=marks,aes(label=mark,color=diagnostic),size=2.2,max.overlaps=Inf,
      min.segment.length=0,box.padding=.12,seed=2032026,show.legend=FALSE)+
   facet_grid(model~scenario)+scale_color_manual(values=diag_colors,labels=pretty_diagnostics,name=NULL)+
   scale_x_continuous(labels=scales::label_percent(accuracy=1))+scale_y_continuous(labels=scales::label_percent(accuracy=1))+
   labs(x='Excluded share',y='Reduction in mean error')+presentation_theme()
}
scorecard_table<-function(scorecard,model_name) {
  v<-scorecard[model==model_name,.(rule,excluded_share,reduction,worst_reduction,retained_coverage)]
  v<-v[match(candidate_rules$rule,rule)]
  gt::gt(as.data.frame(v)) |>
    gt::cols_label(rule='Candidate rule',excluded_share='Excluded',reduction='Mean reduction',worst_reduction='Worst reduction',retained_coverage='95% coverage') |>
    gt::fmt_percent(columns=c(excluded_share,reduction,worst_reduction,retained_coverage),decimals=2) |>
    gt::data_color(columns=c(reduction,worst_reduction),fn=function(x) scales::col_numeric(c('#F2B9B9','#FFFFFF','#AFDEC5'),domain=c(-.15,.15))(pmax(-.15,pmin(.15,x))),autocolor_text=FALSE) |>
    gt::data_color(columns=retained_coverage,fn=function(x) scales::col_numeric(c('#E5F0E8','#F5D8BC','#E8B5B5'),domain=c(0,.10))(pmin(.1,abs(x-.95))),autocolor_text=FALSE) |>
    gt::tab_options(table.font.size=12,data_row.padding=4,table.width=gt::pct(100))
}
scorecard_style_emitted<-FALSE
print_scorecards<-function(scorecard) {
  cat('::: {.panel-tabset}\n\n')
  for(m in model_order) {
    cat('### ',m,' {.unnumbered}\n\n',sep='')
    raw<-gt::as_raw_html(scorecard_table(scorecard,m),inline_css=FALSE)
    stylesheet<-regmatches(raw,regexpr('(?s)<style.*?</style>',raw,perl=TRUE))
    id_attribute<-regmatches(raw,regexpr('id="[^"]+"',raw))
    id<-sub('^id="|"$','',id_attribute);id<-sub('"$','',id)
    if(!scorecard_style_emitted) {
      cat(gsub(paste0('#',id),'.candidate-scorecard',stylesheet,fixed=TRUE),'\n')
      scorecard_style_emitted<<-TRUE
    }
    raw<-sub('(?s)<style.*?</style>','',raw,perl=TRUE)
    cat('<div class="candidate-scorecard">',raw,'</div>\n\n')
  }
  cat(':::\n\n')
}
readable_names<-function(n) {
  labels<-c(scenario='Scenario',model='Model',diagnostic='Diagnostic',cutoff='Cutoff',target='Target',evaluation='Evaluation',
   missing_handling='Missing diagnostic',rule_type='Rule type',sim='Simulation',id='County',k='Available years',stratum='Size stratum',
   spearman='Spearman correlation',auc='AUC',old_n_eligible='Old eligible rows',new_n_eligible='New eligible rows',
   old_baseline_error='Old baseline error',new_baseline_error='New baseline error')
  v<-tools::toTitleCase(gsub('_',' ',n));hit<-n %in% names(labels);v[hit]<-labels[n[hit]]
  v<-gsub('Mc Se','MC SE',v,fixed=TRUE);v<-gsub('N Sims','Finite simulations',v,fixed=TRUE);v
}
full_dt<-function(p,name) {
  x<-p$read(name)
  # All rows of report-level result tables are included, without target/cutoff filters.
  numeric_columns<-names(x)[vapply(x,is.numeric,logical(1))]
  widget<-DT::datatable(as.data.frame(x),rownames=FALSE,colnames=readable_names(names(x)),filter='top',extensions='Buttons',
    options=list(pageLength=25,lengthMenu=c(25,50,100),scrollX=TRUE,dom='Bfrtip',
      buttons=list(list(extend='csv',text='Download filtered CSV',filename=name,
        exportOptions=list(modifier=list(search='applied',page='all')))),deferRender=TRUE))
  if(length(numeric_columns)) widget<-DT::formatRound(widget,columns=numeric_columns,digits=2)
  widget
}

# CSVs already live beside the rendered HTML; no source-tree mirror is needed.
mirror_presentation_csv<-function(root,out) {
  expected<-file.path(root,'_output/analysis/output/3.1a')
  stopifnot(normalizePath(out)==normalizePath(expected))
  invisible(out)
}

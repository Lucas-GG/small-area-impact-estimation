library(data.table)
revision_export_paper_outputs <- function(root=getwd()) {
  source(file.path(root,'R/placebo_revision.R'));ctx<-revision_context(root)
  tables<-file.path(ctx$out,'tables');figures<-file.path(ctx$out,'figures')
  dir.create(tables,showWarnings=FALSE);dir.create(figures,showWarnings=FALSE)
  index<-list()
  sample_titles<-c(S1_full='Full sample (primary)',S3_historical='Historical restriction — retained counties',
    S4_own_z2='Own |z| > 2 screen — retained counties',S5_common_z2='Common |z| > 2 screen — retained counties',
    B1_balanced_full='Balanced sensitivity — full sample',B3_balanced_historical='Balanced sensitivity — historical retained counties',
    B4_balanced_own_z2='Balanced sensitivity — own-screen retained counties',B5_balanced_common_z2='Balanced sensitivity — common-screen retained counties')
  add<-function(path,place,description,specification='both',sample='not applicable',subgroup='both',metric='multiple') {
    index[[length(index)+1]]<<-data.table(filename=substring(path,nchar(ctx$out)+2L),specification,sample,subgroup,metric,intended_place=place,description)
  }
  for(spec_label in c('no_covariates','with_covariates')) {
    performance_path<-file.path(ctx$out,paste0('revision_',spec_label,'_performance.csv'))
    if(file.exists(performance_path)) {
      perf<-ctx$csv(performance_path,'paper tables and figures')
      for(sample_label in unique(perf$sample))for(subgroup_label in unique(perf$nc))for(metric_label in unique(perf$metric)) {
        z<-perf[sample==sample_label & nc==subgroup_label & metric==metric_label & s %in% c(0,4)]
        stem<-paste('revision',spec_label,sample_label,subgroup_label,metric_label,sep='_')
        table_path<-file.path(tables,paste0(stem,'_table.csv'));fwrite(z,table_path)
        place<-if(sample_label=='S1_full'&&subgroup_label=='NonCore')'main text' else if(sample_label=='S1_full')'supplement' else if(sample_label%in%c('S3_historical','B3_balanced_historical'))'appendix / response letter' else 'appendix'
        add(table_path,place,'Performance estimates and Monte Carlo SEs',spec_label,sample_label,subgroup_label,metric_label)
        if(metric_label %in% c('RMSE','crps','cvg','intv')) {
          z[,horizon:=paste('Year',s)]
          z[,dgp:=factor(dgp,levels=c('MAR','Latent level','Latent trend'))]
          z[,model:=factor(model,levels=c('TWFE','ETWFE','RW','RW-IID','RW-ST','GB-Ens'))]
          plot<-ggplot2::ggplot(z,ggplot2::aes(model,estimate))+ggplot2::geom_point()+
            ggplot2::geom_errorbar(ggplot2::aes(ymin=estimate-1.96*mc_se,ymax=estimate+1.96*mc_se),width=.2)+
            ggplot2::facet_grid(horizon~dgp,scales='free_y')+ggplot2::theme_minimal(base_size=11)+
            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,hjust=1))+
            ggplot2::labs(x=NULL,y=metric_label,title=paste(gsub('_',' ',spec_label),sample_titles[[sample_label]],subgroup_label,sep=' | '),
              caption='Error bars: estimate ± 1.96 Monte Carlo SE. S3–S5 describe retained counties.')
          figure_path<-file.path(figures,paste0(stem,'.png'))
          ggplot2::ggsave(figure_path,plot,width=11,height=6,dpi=200)
          add(figure_path,place,'Performance comparison by scenario and horizon',spec_label,sample_label,subgroup_label,metric_label)
        }
      }
    }
    for(name in c('retention','gain_vs_TWFE','upper_tail_pre_rmse','upper_tail_heldout_abs_z','poisson_diagnostics','heldout_calibration','null_benchmarks')) {
      path<-file.path(ctx$out,paste0('revision_',spec_label,'_',name,'.csv'))
      if(file.exists(path))add(path,if(name=='gain_vs_TWFE')'main text / supplement' else if(name=='retention')'appendix' else 'supplement / response letter',
        name,spec_label,if(name=='retention')'S1/S3/S4/S5 and balanced counterparts' else 'S1_full / B1_balanced_full')
    }
  }
  for(name in c('diagnostic_specification_comparison','manuscript_table_spotcheck','manuscript_table_availability')) {
    path<-file.path(ctx$out,paste0('revision_',name,'.csv'))
    if(file.exists(path))add(path,'response letter / audit',name)
  }
  index<-rbindlist(index)
  other<-setdiff(list.files(ctx$out,'[.]csv$',full.names=TRUE),file.path(ctx$out,index$filename))
  other<-other[basename(other)!='revision_output_index.csv']
  if(length(other))index<-rbind(index,data.table(filename=basename(other),specification='both',sample='audit',subgroup='both',metric='multiple',
    intended_place='response letter / audit',description='Supporting audit or input ledger'))
  fwrite(index,file.path(ctx$out,'revision_output_index.csv'))
  invisible(index)
}

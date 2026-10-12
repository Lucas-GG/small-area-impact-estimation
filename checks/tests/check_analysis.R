# Synthetic evaluation checks; no production model fit and no B data are read.
source('R/heldout_analysis.R')
stopifnot(auc_rank(c(1,2,3,4),c(FALSE,FALSE,TRUE,TRUE))==1,
          auc_rank(rep(1,4),c(FALSE,FALSE,TRUE,TRUE))==.5)
z<-data.table(sim=1L,scenario='test',model='test',id=1:10,
 pre_rmse=c(1:9,NA),realized_raw=1:10,realized_standardized=1:10,
 average_raw=1:10,average_standardized=1:10)
x<-evaluate_diagnostics(z,'pre_rmse')
# The unavailable largest-target county still defines the top-decile boundary.
stopifnot(all(x$n_valid==9),all(x$n_top==0),all(is.na(x$auc)))
z[,pre_rmse:=1:10]
x<-evaluate_diagnostics(z,'pre_rmse')
stopifnot(all(abs(x$spearman-1)<1e-12),all(x$auc==1),all(x$n_top==1))
c<-cutoff_performance(z,data.table(diagnostic='pre_rmse',cutoff=8))
stopifnot(all(c$excluded_share==.2),all(c$retained_error==4.5))
# A missing decision must fail before confirmation data can be opened.
root<-tempfile();dir.create(root)
stopifnot(inherits(try(read_decision(root),silent=TRUE),'try-error'))
cat('AUC ties, within-half target classes, thresholds, and decision gate passed.\n')

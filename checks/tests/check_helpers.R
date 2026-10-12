# Focused synthetic checks: no production fit or simulation is rerun.
library(data.table)
source('R/backdate_helpers.R')
source('R/aux.r')
source('R/fe.r')
root <- find_backdate_root()
a <- c('1'=2008, '2'=Inf)
b <- backdate_assignment(a, 3)
stopifnot(b[['1']]==2005, is.infinite(b[['2']]))
panel <- CJ(i=1:2, year=2000:2012)
panel[, `:=`(start_year=Inf, y=2L, n=100, urb='6.NonCore')]
z <- backdate_panel(panel,a,3)
stopifnot(all(is.na(z[i==1 & year>=2005]$y0)),
          all(z[i==1 & year<2005]$y0==2), all(z[i==2]$y0==2),
          all(is.infinite(z$start_year)))
e <- backdate_eligibility(panel,a,3,1,'MAR',204)
stopifnot(nrow(e)==1,e$n_pre==5,e$n_ho==3,e$n_post==5,e$y_ho==6)
z[, `:=`(y0_hat=2, y0_post_1=1, y0_post_2=3)]
s <- summarize_backdate(z,e,'TWFE',2,1000204)
stopifnot(s$mu_ho_mean==6,s$mu_ho_var==18,s$n_valid_draws==2,s$heldout_status=='ok')
stopifnot(identical(s$ppp_ho,summarize_backdate(z,e,'TWFE',2,1000204)$ppp_ho))
# Mid-p tails use the exact same count draws and split ties equally.
counts <- with_backdate_seed(1000205, rpois(2, c(3,9)))
stopifnot(s$ppp_ho_upper == mean(counts > 6)+.5*mean(counts==6),
          s$ppp_ho_lower == mean(counts < 6)+.5*mean(counts==6),
          s$ppp_ho == s$ppp_ho_upper,
          s$ppp_ho_upper+s$ppp_ho_lower == 1)
zero <- copy(z); zero[, c("y0_post_1","y0_post_2","y0_hat") := 0]
zero_e <- copy(e); zero_e[,y_ho:=0]
tied <- summarize_backdate(zero,zero_e,'TWFE',2,1000204)
stopifnot(tied$ppp_ho_upper==.5,tied$ppp_ho_lower==.5)
failed <- run_backdate_one(panel,a,3,1,'MAR',204,'TWFE',
  list(fun=function(...) stop('synthetic failure'),args=list()),2)
stopifnot(all(is.na(failed$ppp_ho_upper)),all(is.na(failed$ppp_ho_lower)))

z[i==1 & year==2005,y0_post_2:=NA_real_]
stopifnot(summarize_backdate(z,e,'TWFE',2,1000204)$heldout_status=='incomplete_draws')
# A mocked bootstrap checks the actual inference_fe body without estimating any model.
original_fit <- fit_fe
original_boot <- if(exists('cluster_boot')) cluster_boot else NULL
fit_fe <- function(.dt,fm,wts) {
  stopifnot(all(.dt$marker=='argument')); structure(list(),class='heldout_mock')
}
cluster_boot <- function(.dt) rep(1,nrow(.dt))
predict.heldout_mock <- function(object,newdata,...) rep(2,nrow(newdata))
formula.heldout_mock <- function(x,...) y0~1
dt <- data.table(i=1:2,marker='wrong_global')
passed <- data.table(i=1:2,marker='argument')
invisible(inference_fe(passed,structure(list(),class='heldout_mock'),nboot=2))
fit_fe <- original_fit
if(!is.null(original_boot)) cluster_boot <- original_boot
# Synthetic saved MAIN predictions verify joins, diagnostics and rejection of bad dates.
arr <- array(rep(a,1),dim=c(2,1,1),dimnames=list(c('1','2'),'1','1'))
main <- copy(panel)
main[, `:=`(r=1L,start_year_sim=a[as.character(i)],y0_hat=2)]
p <- tempfile(fileext='.rds');saveRDS(main,p)
m <- read_main_summary(p,panel,arr,1,'TWFE')
stopifnot(nrow(m)==1,m$y_post==10,m$mu_post_mean==10,m$pre_rmse==0,
          m$pre_rms_pearson==0)
main[i==1,start_year_sim:=2009];saveRDS(main,p)
stopifnot(inherits(try(read_main_summary(p,panel,arr,1,'TWFE'),silent=TRUE),'try-error'))
unlink(p)
cat('Masking, window aggregation, covariance-preserving variance, bootstrap argument, and main-input checks passed.\n')

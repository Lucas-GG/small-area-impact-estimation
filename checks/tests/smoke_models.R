# Optional small end-to-end test of the six shared model wrappers (synthetic data only).
# Run from the repository root; INLA's native runtime can take several minutes.
library(data.table)
invisible(lapply(list.files('R', '\\.r$', full.names=TRUE), source))
source('R/backdate_helpers.R')
panel <- CJ(i=1:8, year=1999:2014)
panel[, `:=`(start_year=2011+i%%4, urb='6.NonCore', n=100,
             y=as.integer(1+(year+i)%%3))]
W <- matrix(1,8,8); diag(W) <- 0
options(xpn.conf=list(Lset=c(2,4),graph=Matrix::Matrix(W)))
INLA::inla.setOption(num.threads=1)
a <- setNames(c(2008,2008,rep(Inf,6)),1:8)
specs <- backdate_specs()
result <- rbindlist(lapply(names(specs), function(model)
  run_backdate_one(panel,a,3,1,'MAR',204,model,specs[[model]],2,audit_features=TRUE)),fill=TRUE)
print(result[, .N, by=.(model,heldout_status)])
stopifnot(all(result$heldout_status=='ok'),all(result$n_valid_draws==2L))
cat('All six synthetic model wrappers passed.\n')

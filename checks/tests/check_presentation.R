# The two-simulation example distinguishes the order of averaging and quantiles.
source('R/heldout_random_band.R')
x<-rbind(c(-1,1),c(1,-1))
b<-simulation_averaged_band(x)
old<-rowMeans(apply(x,1,quantile,c(.025,.975)))
stopifnot(b[['mean']]==0,b[['lo']]==0,b[['hi']]==0,old[2]>old[1])
# Quantile correction preserves the mean of random reductions.
y<-rbind(c(0,.1,.2),c(.4,.3,.2))
c<-simulation_averaged_band(y)
stopifnot(abs(c[['mean']]-mean(y))<1e-12)
cat('Random-band averaging order and unchanged mean checks passed.\n')

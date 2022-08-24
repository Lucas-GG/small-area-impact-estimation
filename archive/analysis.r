

#===============================================================================
#1.-  simulate one dataset
#===============================================================================
load(file.path(path0,'data/toy.par'))

simdat0 <- with(toy.par, {sim.out.exp(
  N=N,              # total number of areas
  T=T,              # total number of periods
  n=n,              # population at risk in each county and year
  Q.scaled=Q.scaled,# (scaled) matrix with the neighborhood structure
  sigma2  =  .16,   # variance of log relative risk
  r0      =  3/1000,# baseline risk
  phi     =  .1,    # spatial dependence ('mixing') parameter
  rho     =  .8 ,   # serial correlation parameter
  t.noise =  FALSE, # additional noise
  N_t=N*.7,         # Number of treated units desired.
  T0=2)})           # The first treatment time. The rest of treatment times are equally spaced between T0 to T.

head(simdat0)

dat <- simdat0
#===============================================================================
fit.SMR0.hat <- function (dat) {
  dat      <- dat[order(dat$time,dat$id),]
  dat0     <- dat
  N        <- with(dat,length(unique(id)))
  T        <- with(dat,length(unique(time)))
  dat$SMR  <- with(dat,Y/(n*r0))

  ###0. DID
  with(dat,{DID(M=matrix(SMR,N) ,mask=matrix(O,N))})


  SMR0.hat.did <-  with(dat,{DID(M=matrix(SMR,N) ,mask=matrix(O,N))})

  ###i.	Synthetic control (SC) (Abadie et al., 2003, 2010, 2015)
  SMR0.hat.sc <-  with(dat,{adh_mp_rows( M=matrix(SMR,N) ,mask=matrix(O,N))})

  ###ii.	Staggered SC (Ben-Michael et. al. 2019) https://arxiv.org/abs/1912.03290v2
  #! replaces Elastic net (EN) (Doudchenko & Imbens, 2016)
  #dat <- simdat0
  #N        <- with(dat,length(unique(id)))
  #T        <- with(dat,length(unique(time)))
  #dat$SMR  <- with(dat,Y/(n*r0))
  dat$W    <-  1 - dat$O
  n_leads  <-  T-dat$T0[1]
  msc <- multisynth(SMR ~ W, id, time, dat,n_leads=n_leads)
  SMR0.hat.msc <-  predict(msc,att=F)
  SMR0.hat.msc <-  SMR0.hat.msc [-dim(SMR0.hat.msc )[1],-1]
  SMR0.hat.msc <-  t(apply(SMR0.hat.msc,2,na.omit))
  SMR0.hat.msc <-  merge(data.frame(id=1:N),data.frame(id=c(1:N)[msc$data$trt<Inf],SMR0.hat.msc),all=T)[,-1]
  SMR0.hat.msc <-  as.matrix(SMR0.hat.msc)

  ###iii.	Matrix completion (MC) (Athey et al., 2017) https://arxiv.org/abs/1710.10251v4
  mc <-  with(dat,{mcnnm_cv(M=matrix(SMR,N) ,mask=matrix(O,N))})
  SMR0.hat.mc <-   with(mc, L +  replicate(ncol(L),u) + t(replicate(nrow(L),v)))

  ###iv.	Disease mapping (DM) (Bauer et al., 2016)
  dat$Y0      <- dat$Y
  dat$Y0[dat$O==0]<- NA
  m0 <- inla(Y0 ~ offset(log(n*r0))
                + f(id , model="bym2", graph = gs,
                    group = time,control.group = list(model="ar1"))
               , data=dat
               , family="poisson"
               , control.predictor = list(link = 1))

  SMR0.hat.st <-  with(dat,matrix(m0$summary.fitted.values$mean/(n*r0),N))

  SMR0.hat.did.long <-  reshape(data.frame(SMR0.hat.did),direction='long',varying=1:T,v.names='SMR0.hat.did')
  SMR0.hat.sc.long  <-  reshape(data.frame(SMR0.hat.sc),direction='long',varying=1:T,v.names='SMR0.hat.sc')
  SMR0.hat.msc.long <-  reshape(data.frame(SMR0.hat.msc),direction='long',varying=1:T,v.names='SMR0.hat.msc')
  SMR0.hat.mc.long  <-  reshape(data.frame(SMR0.hat.mc),direction='long',varying=1:T,v.names='SMR0.hat.mc')
  SMR0.hat.st.long  <-  reshape(data.frame(SMR0.hat.st),direction='long',varying=1:T,v.names='SMR0.hat.st')

  rdat    <- Reduce(function(x, y) merge(x, y, by=c('id','time')),
              list(dat,
                   SMR0.hat.did.long,
                   SMR0.hat.sc.long,
                   SMR0.hat.msc.long,
                   SMR0.hat.mc.long,
                   SMR0.hat.st.long
                   ))
  row.names(rdat) <- NULL
  rdat
  }

rmse.compare <- function (dat,g) {
  with(dat,{
    SMR    <- Y/(n*r0)
    etors  <- grep('SMR0*',colnames(rdat),value = TRUE)
    rmse   <- sapply(etors, function(etor) {
        d.hat <- get(etor) - SMR
        d.hat[O==1] <- NA
        sapply(split(d.hat,g),function(d){sqrt(mean(na.omit(d)^2))})
        })
      })
   rmse
  }





#===============================================================================
load(file.path(path0,'data/case.par'))


rdat <- get.SMR0.hat(dat)

g <- case.par$map$urb[case.par$map$state%in%c('AZ','NV')]
rmse.compare(rdat,g)

#===============================================================================

plot.rmse.ctyr <- bwplot(rmse~factor(method)|factor(r0)*factor(phi)*factor(rho)*factor(N_t)*factor(T0),
         data=rmse.ctyr.dat,
         main="",
         xlab="Method", ylab="RMSE",
         panel=function(x,y,...){
           panel.bwplot(x,y, ...)
           panel.abline(h=median(y),lty=3)
         })

update(plot.rmse.ctyr,layout = c(4, 2,4))

png(file.path(path,'output/%03d.png'),width = 480*9,heigh=480*4.5,res=300)
update(plot.rmse.ctyr,layout = c(4, 2,4))
dev.off ()

#===============================================================================

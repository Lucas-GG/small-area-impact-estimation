path0 <- file.path('C:/LOCAL/ICF/R03 SAE')   
path <- file.path('C:/LOCAL/ICF/R03 SAE/STevaluation')   

#0.- load sim fucntions
setwd(path)
Rfiles <- list.files(file.path(paste0(getwd(),"/R/")), ".R")
Rfiles <- Rfiles[grepl(".R", Rfiles)]
sapply(paste0(paste0(getwd(),"/R/"), Rfiles), source)

#0.- load requiered packges
library(INLA)
#library(devtools) 
#install_github("susanathey/MCPanel")
#devtools::install_github("ebenmichael/augsynth")
library(MCPanel)
library(augsynth)

#===============================================================================
#0.- Add info to the map
#===============================================================================

load(file.path(path0,'data/case.par'))
names(case.par)

load(file.path(path0,'data/hmap'))
map <- rebuild_CRS(hmap)
map <- merge(map,tigris::fips_codes[,c('state_code','county_code','state','state_name')],by.x=c('STATEFP','COUNTYFP'),by.y=c('state_code','county_code'),all.x=T)
table(map$state)

map$urb <- 'Metro'
map$urb[map$urb13=='5.Micro']   <- 'Micro'
map$urb[map$urb13=='6.NonCore'] <- 'Noncore'


map@data$state_fips  <- map@data$STATEFP
map@data$county_fips <- map@data$COUNTYFP
map@data$id          <- map@data$COUNTYFIPS

head(map@data)
map@data <- map@data[,c('id','state_fips','county_fips','state','state_name','urb13','urb')]
head(map@data)

case.par$map <- map
save(case.par,file=file.path(path0,'data/case.par'))

#===============================================================================
#0.-  Generate toy example
#===============================================================================
library(spdep)
load(file.path(path0,'data/case.par'))

table(case.par$map$state%in%c('AZ','NV'))
plot(map)
plot(map[map$state%in%c('AZ','NV'),] ,col = "blue",add = TRUE)

smap  <- map[map$state%in%c('AZ','NV'),]
W     <- nb2mat(poly2nb(smap), style = "B")
gs    <- inla.read.graph(W)
Q = -W
diag(Q) = 0
diag(Q) = -apply(Q,1,sum)
Q[1:10,1:10]
Q.scaled <- inla.scale.model(Q, constr=list(A=matrix(1,1, dim(Q)[1]), e=0))
exp(mean(log(diag(INLA:::inla.ginv(case.par$Q.scaled)))))

n   <- case.par$n[paste(smap$id),]
N   <- dim(n)[1]
T   <- dim(n)[2]

toy.par <- list(N=N,T=T,n=n,gs=gs)
save(toy.par,file=file.path(path0,'data/toy.par'))

summary(toy.par)


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

simdat0 <- simdat0[order(simdat0$time,simdat0$id),]
#===============================================================================

rmse <- function(Y0.hat, Y, O, g) {
  N           <- dim(Y0.hat)[1]
  d.hat       <- matrix(Y,N) - Y0.hat
  d.hat[O==1] <- NA   
  sapply(split(d.hat,g),function(d){sqrt(mean(na.omit(d)^2))})
  }
  
 with(dat, {rmse(SMR0.hat.did, SMR, O, urb)}) 
detach(rdat)
 


colnames(rdat)

g <- case.par$map$urb[case.par$map$state%in%c('AZ','NV')] 
rmse.compare <- function (dat,g=g) {

  with(rdat,{
    SMR   <- Y/(n*r0)
    etors <- grep('SMR0*',colnames(rdat),value = TRUE)
    sapply(etors, function(etor) {
      d.hat <- get(etor) - SMR
      d.hat[O==1] <- NA  
      sapply(split(d.hat,g),function(d){sqrt(mean(na.omit(d)^2))})
      })
    })
  tag <- list(
    N=length(unique(id)),
    T=length(unique(time)),
    sigma2=sigma2[1], 
    phi=phi[1], 
    rho=rho[1], 
    t.noise=t.noise[1],  
    N_t=N_t[1], 
    T0=T0[1],   
  }

rmse.compare(rdat)


rmse.compare <- function (dat) {
  SMR <-  with(dat, Y/(n*r0)) 
 
  ###0. DID
  SMR0.hat.did <-  with(dat,{DID( M=matrix(SMR,N) ,mask=matrix(O,N)) })
  rmse.did <- with(dat, {rmse(SMR0.hat.did, SMR, O, urb)}) 
  
  ###i.	Synthetic control (SC) (Abadie et al., 2003, 2010, 2015)
  SMR0.hat.sc <-  with(dat,{DID( M=matrix(SMR,N) ,mask=matrix(O,N)) })
  rmse.sc <- with(dat, {rmse(SMR0.hat.sc, SMR, O, urb)}) 
  
  ###ii.	Staggered SC (Ben-Michael et. al. 2019) https://arxiv.org/abs/1912.03290v2
  #! replaces Elastic net (EN) (Doudchenko & Imbens, 2016) 
  msc <- multisynth(SMR ~ W, id, time, dat,n_leads=T-dat$T0[1] )
  SMR0.hat.msc <-  predict(msc,att=F) 
  SMR0.hat.msc <-  SMR0.hat.msc [-dim(SMR0.hat.msc )[1],-1]
  SMR0.hat.msc <-  t(apply(SMR0.hat.msc,2,na.omit))
  SMR0.hat.msc <-  merge(data.frame(id=1:N),data.frame(id=c(1:N)[msc$data$trt<Inf],SMR0.hat.msc),all=T)[,-1]
  SMR0.hat.msc <-  as.matrix(SMR0.hat.msc)
  rmse.msc <- with(dat, {rmse(SMR0.hat.msc, SMR, O, urb)})  
  
  ###iii.	Matrix completion (MC) (Athey et al., 2017) https://arxiv.org/abs/1710.10251v4
  mc <-  with(dat,{mcnnm_cv( M=matrix(SMR,N) ,mask=matrix(O,N)) })
  SMR0.hat.mc <-   with(mc, L +  replicate(ncol(L),u) + t(replicate(nrow(L),v)))
  rmse.mc <- with(dat, {rmse(SMR0.hat.mc, SMR, O, urb)}) 
  
  ###iv.	Disease mapping (DM) (Bauer et al., 2016)
  dat$Y0          <- dat$Y
  dat$Y0[dat$O==0]<- NA
  m0 <- inla(Y0 ~ offset(log(n*r0))
                + f(id , model="bym2", graph = gs,         
                    group = time,control.group = list(model="ar1"))  
               , data=dat
               , family="poisson"
               , control.predictor = list(link = 1))
   
  SMR0.hat.st <-  matrix(with(dat,m0$summary.fitted.values$mean/(n*r0)),N)
  rmse.st <- with(dat, {rmse(SMR0.hat.st, SMR, O, urb)})
   
  list(rmse.did=rmse.did,
       rmse.sc=rmse.sc,
       rmse.msc=rmse.msc,
       rmse.mc=rmse.mc ,
       rmse.st=rmse.st)
  }
  
  


get.SMR0.hat <- function (dat) {
  dat0 <- dat
  N <- length(unique(dat$id))
  T <- length(unique(dat$time))
  dat$SMR  <-  with(dat, Y/(n*r0)) 
  
  ###0. DID
  SMR0.hat.did <-  with(dat,{DID( M=matrix(SMR,N) ,mask=matrix(O,N)) })
    
  ###i.	Synthetic control (SC) (Abadie et al., 2003, 2010, 2015)
  SMR0.hat.sc <-  with(dat,{adh_mp_rows( M=matrix(SMR,N) ,mask=matrix(O,N)) })
 
  ###ii.	Staggered SC (Ben-Michael et. al. 2019) https://arxiv.org/abs/1912.03290v2
  #! replaces Elastic net (EN) (Doudchenko & Imbens, 2016) 
  dat$W    <-  1- dat$O
  msc <- multisynth(SMR ~ W, id, time, dat,n_leads=T-dat$T0[1] )
  SMR0.hat.msc <-  predict(msc,att=F) 
  SMR0.hat.msc <-  SMR0.hat.msc [-dim(SMR0.hat.msc )[1],-1]
  SMR0.hat.msc <-  t(apply(SMR0.hat.msc,2,na.omit))
  SMR0.hat.msc <-  merge(data.frame(id=1:N),data.frame(id=c(1:N)[msc$data$trt<Inf],SMR0.hat.msc),all=T)[,-1]
  SMR0.hat.msc <-  as.matrix(SMR0.hat.msc)
  
  ###iii.	Matrix completion (MC) (Athey et al., 2017) https://arxiv.org/abs/1710.10251v4
  mc <-  with(dat,{mcnnm_cv( M=matrix(SMR,N) ,mask=matrix(O,N)) })
  SMR0.hat.mc <-   with(mc, L +  replicate(ncol(L),u) + t(replicate(nrow(L),v)))
  
  ###iv.	Disease mapping (DM) (Bauer et al., 2016)
  dat$Y0          <- dat$Y
  dat$Y0[dat$O==0]<- NA
  m0 <- inla(Y0 ~ offset(log(n*r0))
                + f(id , model="bym2", graph = gs,         
                    group = time,control.group = list(model="ar1"))  
               , data=dat
               , family="poisson"
               , control.predictor = list(link = 1))
   
  SMR0.hat.st <-  matrix(with(dat,m0$summary.fitted.values$mean/(n*r0)),N)
  
  
  SMR0.hat.did.long <-  reshape(data.frame(SMR0.hat.did),direction='long',varying=1:T,v.names='SMR0.hat.did')
  SMR0.hat.sc.long  <-  reshape(data.frame(SMR0.hat.sc),direction='long',varying=1:T,v.names='SMR0.hat.sc')
  SMR0.hat.msc.long <-  reshape(data.frame(SMR0.hat.msc),direction='long',varying=1:T,v.names='SMR0.hat.msc')
  SMR0.hat.mc.long  <-  reshape(data.frame(SMR0.hat.mc),direction='long',varying=1:T,v.names='SMR0.hat.mc')      
  SMR0.hat.st.long  <-  reshape(data.frame(SMR0.hat.st),direction='long',varying=1:T,v.names='SMR0.hat.st')
  
  rdat    <- Reduce(function(x, y) merge(x, y, by=c('id','time')), 
              list(dat0,
                   SMR0.hat.did.long,
                   SMR0.hat.sc.long,
                   SMR0.hat.msc.long,
                   SMR0.hat.mc.long,
                   SMR0.hat.st.long
                   ))
  rdat
  }


#===============================================================================
load(file.path(path0,'data/case.par'))


dat <- with(case.par$map[map$state%in%c('AZ','NV'),]@data,
           data.frame(simdat0,state,urb))

rdat <- get.SMR0.hat(dat)

head(rdat,25)
           
rmse.compare(dat)


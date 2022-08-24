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



case.par$map <- map
save(case.par,file=file.path(path0,'data/case.par'))

#===============================================================================
#0.-  Generate toy example
#===============================================================================
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

n   <- case.par$n[paste(smap$COUNTYFIPS),]
N   <- dim(n)[1]
T   <- dim(n)[2]

toy.par <- list(N=N,T=T,n=n,gs=gs)
save(toy.par,file=file.path(path0,'data/toy.par'))




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
  phi     =  .3,    # spatial dependence ('mixing') parameter
  rho     =  .8 ,   # serial correlation parameter
  t.noise =  FALSE, # additional noise 
  N_t=N*.7,         # Number of treated units desired.
  T0=2)})           # The first treatment time. The rest of treatment times are equally spaced between T0 to T.

head(simdat0)

simdat0 <- simdat0[order(simdat0$time,simdat0$id),]

#===============================================================================
load(file.path(path0,'data/case.par'))


  

sdat0 <- with(case.par$map[map$state%in%c('AZ','NV'),]@data,
           data.frame(simdat0,state,urb))

rmse <- function(Y0.hat, Y, O, g) {
  N           <- dim(Y0.hat)[1]
  d.hat       <- matrix(Y,N) - Y0.hat
  d.hat[O==1] <- NA   
  sapply(split(d.hat,g),function(d){sqrt(mean(na.omit(d)^2))})
  } 

sdat0$SMR <- with(sdat0, Y/(n*r0)) 
sdat0$W <-  1- sdat0$O

###0. DID
SMR0.hat.did <-  with(sdat0,{DID( M=matrix(SMR,N) ,mask=matrix(O,N)) })
with(sdat0, {rmse(SMR0.hat.did, SMR, O, urb)}) 

###i.	Synthetic control (SC) (Abadie et al., 2003, 2010, 2015)
SMR0.hat.sc <-  with(sdat0,{DID( M=matrix(SMR,N) ,mask=matrix(O,N)) })
with(sdat0, {rmse(SMR0.hat.sc, SMR, O, urb)}) 

###ii.	Staggered SC (Ben-Michael et. al. 2019) https://arxiv.org/abs/1912.03290v2
#! replaces Elastic net (EN) (Doudchenko & Imbens, 2016) 
msc <- multisynth(SMR ~ W, id, time, sdat0,n_leads=T-sdat0$T0[1] )
SMR0.hat.msc <-  predict(msc,att=F) 
SMR0.hat.msc <-  SMR0.hat.msc [-dim(SMR0.hat.msc )[1],-1]
SMR0.hat.msc <-  t(apply(SMR0.hat.msc,2,na.omit))
SMR0.hat.msc <-  merge(data.frame(id=1:N),data.frame(id=c(1:N)[msc$data$trt<Inf],SMR0.hat.msc),all=T)[,-1]
SMR0.hat.msc <-  as.matrix(SMR0.hat.msc)
with(sdat0, {rmse(SMR0.hat.msc, SMR, O, urb)})  


###iii.	Matrix completion (MC) (Athey et al., 2017) https://arxiv.org/abs/1710.10251v4
mc <-  with(sdat0,{mcnnm_cv( M=matrix(SMR,N) ,mask=matrix(O,N)) })
SMR0.hat.mc <-   with(mc, L +  replicate(ncol(L),u) + t(replicate(nrow(L),v)))
with(sdat0, {rmse(SMR0.hat.mc, SMR, O, urb)}) 

###iv.	Disease mapping (DM) (Bauer et al., 2016) 



m0 <- inla(ifelse(sdat0$O,Y,NA) ~ offset(log(n*r0))
              + f(id , model="bym2", graph = gs,         
                  group = time,control.group = list(model="ar1"))  
             , data=sdat0
             , family="poisson"
             , control.predictor = list(link = 1))
 
SMR0.hat.st <-  matrix(with(sdat0,m0$summary.fitted.values$mean/(n*r0)),N)
with(sdat0, {rmse(SMR0.hat.st, SMR, O, urb)}) 

#ovarall, sae, county-level
# for a single point estimate (a handfull of point estiamtes) we can use the MSE  or the RMSE
# where M:mean across simulations
# but we are also intrested in the accuracy of county-levle estiamtions
# RM(MSE)  (i.e., ther root mean acroos simulation of the average square error across counties)
# this is similar to PEHE (precision in estimation of heterogeneous effects) by Hill 2010
# with county as the only categorical covariate
# what about M(RMSE); 

# in across sectional setting whre  the estiamtion variaes as a fucntion of x.

# whithin MSE is equivalent to Euclidian norm ( square root of the sum of squares) 1/sqrt(n)
# distance  
# we also have to deal with periods 


# 
# across counties mean absolute diffrece
# why? 
# we want to know how well we do with each county 
# regardless of whterh the forcast compansate 

#  empirical coverage rates of prediction intervals (PIs).


# Athey
# average root-mean-squared-error (RMSE)
# average across 



#===============================================================================
load(file.path(path0,'data/hmap'))
 

dim(coordinates(hmap))

n.area  <- dim(coordinates(hmap))[1]
n.time  <- 5
n.basis <- 4

B.x     <- splines::bs(coordinates(hmap)[,1],df=n.basis ) 
B.y     <- splines::bs(coordinates(hmap)[,2],df=n.basis )
Z       <- matrix(nrow =n.area  , ncol = n.basis ^2)
for (i in 1:n.area) {
  Z[i,] <-kronecker(B.x[i,], B.y[i,])
  }
Z[1:10,1:10]

#C <- mgcv::tensor.prod.model.matrix(list(B.x,B.y)  )
#C[1:10,1:10]
#Z[1:10,1:10]



### construct the stacked design matrix: Z matrix is the tensor-product B spline


X.mat <- matrix(NA, nrow=n.area*n.time, ncol=n.basis*n.time)
dim(X.mat )
I.tmp <- diag(1, nrow=n.time, ncol=n.time)

for (i in 1:n.area){
  X.mat[((i-1)*(n.time)+1):(i*n.time),] <- I.tmp %x% matrix(Z[i,], nrow=1)
  }
dim(X.mat) # should be (IT)*(KT)

#########--------------- type 2 interaction ---------------#########
# create vectors in the A matrix

intercept <- c(1,  rep(NA, n.basis*n.time))
idx       <- c(NA, rep(1:n.time, each=n.basis))
group     <- c(NA, rep(1:n.basis, n.time))

data.inla <- list(y=1, idx=idx, intercept=intercept, group=group)
head(data.inla)

# fit the model
formula.2 <- y  ~ -1 + intercept +
  f(idx, model = "rw2", constr=F, hyper = list(prec = list(param=c(1, 0.005))),
    replicate=group)

r.2 = inla(formula.2, data=data.inla,family="poisson", 
  control.predictor = list(A=cBind(rep(1, nrow(X.mat)), X.mat), compute=TRUE, link=1),
  control.compute=list(dic=T, cpo=TRUE), verbose=F, 
  E=E.lvec)

summary(r.2)


# extract the coefficients
b.hat.2.inla <- r.2$summary.linear.predictor[(n+2):(n+1+m), "mean"]
b.hat.mat.2.inla <- matrix(b.hat.2.inla, nrow=n.basis, ncol=n.time, byrow=F)
head(b.hat.mat.2.inla)
options(digits=3)
options(scipen=10000)
#load requiered packges
library(INLA)


#===============================================================================
###I.- Obtain parameters for simulation from case data
#===============================================================================
load('data/ldat.r')
load('data/gs')
load('data/W' )

case.par <- list()

#1.- total number of areas (1)
case.par$N <- length(unique(ldat$id))


#2.- total number of periods (there are actually 2, parallel to T0s)
#!the series for suicide is longer than for hospitlaizations
length(unique(ldat$year[!is.na(ldat$S1024)]))
length(unique(ldat$year[!is.na(ldat$Y)]))
case.par$T <- 9


#3.- population at risk in each county and year (1)
n0 <- reshape(ldat[ldat$otyp=='IP' & ldat$year >=2008,c('id','year','P1024')],
                v.names = 'P1024', idvar = 'id',
                timevar = 'year', direction = 'wide')
case.par$n <- n0[,-1]
rownames(case.par$n)<- n0$id
colnames(case.par$n)<- 2008:2016
head(case.par$n)

hist(log(case.par$n[,1]))
mean(log(case.par$n[,1]))
sd(log(case.par$n[,1]))



#4.- (scaled) matrix with the neighborhood structure (1)
#!scales Q so the geometric mean of the marginal variances is one
#!see Riabler etal 2015, https://doi.org/10.13140/RG.2.1.1899.2080, page 6
case.par$gs <- gs
Q = -INLA::inla.graph2matrix(gs)
diag(Q) = 0
diag(Q) = -apply(Q,1,sum)
case.par$Q.scaled <- INLA::inla.scale.model(Q, constr=list(A=matrix(1,1, dim(Q)[1]), e=0))
exp(mean(log(diag(INLA:::inla.ginv(case.par$Q.scaled)))))


#5.- baseline risk  (2)
#	Two levels of baseline risk (corresponding approximately with suicide to self-harm hospitalizations among youth)
with(ldat ,
  tapply(Y,list(otyp,year),sum,na.rm=T)/tapply(P1024,list(otyp,year),sum,na.rm=T))*1000

with(ldat[ldat$year>=2008,],
        tapply(Y,list(otyp),sum,na.rm=T)/tapply(P1024,list(otyp),sum,na.rm=T))*1000

with(ldat[ldat$otyp=='IP',] ,
  tapply(S1024,year,sum,na.rm=T)/tapply(P1024,year,sum,na.rm=T))*100000

case.par$r0  <- c(IP=2.9/1000,S=8.9/100000)
case.par$r0

#6.- Number of treated units desired.   (2)
N_t0 <- length(unique(ldat$id[ldat$yr1<9999]));N_t0
round(N_t0/case.par$N,1)
case.par$N_t <- c(actual=N_t0,ideal=floor(case.par$N*.5))  #actual70% and a better scenario, i.e., 50%
case.par$N_t

#7.- The first treatment time.  (2)
min.yr <- with(ldat[ldat$otyp=='IP'&ldat$yr1<9999,], min(yr1))
min.yr-c(1999,2008)
7/19*9
case.par$T0 <- c(IP=2,S=4)  #corresponding to hopitalizatins and suicide approx.
case.par$T0

#!the series for suicide is longer than for hospitlaizations
#!T0 is therefore different

#8.- (residual) variance of log relative risk (we mention only 1)
# 'residual': after taking systematic transd into account
# IP is more varible than suicide beore taking onto account trends
# residual variance is very similar though

load(file.path(path,'Ysp.ip'))
summary(Ysp.ip)
load(file.path(path,'Ysp.ip.bart'))
summary(Ysp.ip.bart)

load(file.path(path,'Ysp.s'))
summary(Ysp.s)
load(file.path(path,'Ysp.s.bart'))
summary(Ysp.s.bart)

1/6.304  ;1/6.275

case.par$sigma2 <- .16



#9.- spatial dependence ('mixing') parameter (3)
case.par$phi <- c(IP=.1,S=.3,.5)

#10.- serial correlation parameter  (3)
case.par$rho <- c(IP=.85,S=.95,.75)

summary(case.par)

case.par$urb13 <- ldat[ldat$otyp=='IP' & ldat$year ==2008,c('id','urb13')]
table(case.par$urb13$urb13)


save(case.par,file=file.path(path,'case.par'))

#===============================================================================
###II.- Prepare map for simulations
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
###III.-  Generate 'toy' example
#===============================================================================
library(spdep)
load(file.path(path,'data/case.par'))

toy.par <- case.par

#table(case.par$map$state%in%c('AZ','NV'))
#plot(map)
#plot(map[map$state%in%c('AZ','NV'),] ,col = "blue",add = TRUE)
smap  <- case.par$map[case.par$map$state%in%c('AZ','NV'),]
W     <- nb2mat(poly2nb(smap), style = "B")
Q = -W
diag(Q) = 0
diag(Q) = -apply(Q,1,sum)
#Q[1:10,1:10]
toy.par <- within(toy.par,{
  Q.scaled <- inla.scale.model(Q, constr=list(A=matrix(1,1, dim(Q)[1]), e=0))
  #exp(mean(log(diag(INLA:::inla.ginv(case.par$Q.scaled)))))
  gs    <- inla.read.graph(W)
  n   <- case.par$n[paste(smap$id),]
  N   <- dim(n)[1]
  T   <- dim(n)[2]
  T0  <- floor(T*c(2/9,4/9))
  N_t <- floor(N*(case.par$N_t/case.par$N))
  map <- smap
  })


summary(case.par)
summary(toy.par)
save(toy.par,file=file.path(path,'data/toy.par'))

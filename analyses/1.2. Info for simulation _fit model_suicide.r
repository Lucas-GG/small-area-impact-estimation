options(digits=3)
options(scipen=10000)


#===============================================================================
#.- Obtain parameters for simulation from case data
#===============================================================================
#path.h <- 'C:/LOCAL/ICF/R03 SAE/DATA'
ldat <-  readRDS('data/ldat.rds')
load('data/gs')

#.- Obtain parameters for simulation from case data

data(state.fips, package='maps')
st <- with(state.fips,tapply(abb,fips,max))
ldat$st   <- st[paste(ldat$state )]
ldat$st   <- factor(ldat$st)
ldat$otyp <- factor(ldat$otyp)
ldat$fid  <- factor(ldat$id)
ldat$fyr  <- factor(ldat$year)
ldat$rural <-  ldat$urb13%in%c('6.NonCore','5.Micro')

dim(ldat)
dat <- ldat[ldat$otyp=='IP' ,]
dim(dat)

#===============================================================================
dat$Y <- dat$S

dat$risk0 <- with(dat, ave(Y,outcome,FUN=function(x) sum(x,na.rm=T))/
                          ave(P,outcome,FUN=function(x) sum(x,na.rm=T)))
dat$E     <- with(dat, P*risk0)

tapply(dat$E,dat$outcome,summary)


#===============================================================================
options(java.parameters = "-Xmx20g")
library(bartMachine)
set_bart_machine_num_cores(6)

#-------------------------------------------------------------------------------

#standardized mortality (or outcome) ratio (SMR)
dat$smr <- dat$Y/dat$E
summary(dat$smr)
min(dat$smr[dat$smr>0],na.rm=T)
hist(dat$smr)

# log SMR
dat$lsmr <- log(dat$smr)
summary(dat$lsmr)
dat$lsmr[dat$lsmr ==-Inf] <- NA
hist(dat$lsmr)

colnames(dat)

fy <- formula(lsmr ~ E
             + yr
             + st
             + urb13
             + I(M/TOTP) #overall mortality except suicide
             + AIAN + uemp + mhinc + pui
             + opop
           )

summary(lm(fy,data=dat)) #R2~.40

yX <- model.frame(fy,data=dat[!is.na(dat$lsmr),],na.action=na.pass)
yX <- yX[!is.na(yX[, 1]),]
summary(yX)
bm <- bartMachine(yX[,-1],yX[, 1],use_missing_data =T,serialize = TRUE)
#save(bm,file=file.path(path.h,'bm',sep=''))

X <- model.frame(fy,data=dat,na.action=na.pass)[,-1]
dat$smr0 <- NULL
dat$smr0 <- predict(bm,X)

plot(dat$smr0,dat$lsmr)
summary(lm(lsmr ~ smr0,data=dat)) #R2~.64


#===============================================================================
#-------------------------------------------------------------------------------
library(lme4)

glm(Y ~ offset(log(E))  ,family=quasipoisson,data=dat))
summary(glm(Y ~ offset(log(E))+smr0  ,family=quasipoisson,data=dat))

summary(glmer(Y ~ offset(log(E)) + (1|i),family=poisson,data=dat))
summary(glmer(Y ~ offset(log(E)) + smr0 + (1|i),family=poisson,data=dat))
#varaince from .23 to .14

summary(glmer.nb(Y ~ offset(log(E)) + (1|i), data=dat))
summary(glmer.nb(Y ~ offset(log(E)) + smr0 + (1|i), data=dat))
#resiudla variance increase


#===============================================================================
library(INLA)
#.- Obtain parameters for simulation from case data

#!predict missing Y takes time
dat <- dat[!is.na(dat$Y),]


dat$yr.i         <-  as.numeric(as.factor(dat$yr))
dat$io <- 1:dim(dat)[1]

##1.- overdispersed poisson
Ypo <- inla(Y ~ offset(log(E))
              + f(io, model='iid')
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Ypo)

Ypo.bart <- inla(Y ~ offset(log(E)) + smr0
              + f(io, model='iid')
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Ypo.bart)
#! .279 to .094


##2.- hierarchical
#hirachical with overdispersion
Yh <- inla(Y ~ offset(log(E))
              + f(i , model='iid')
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Yh)

Yh.bart <- inla(Y ~ offset(log(E))  + smr0
              + f(i , model='iid')
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Yh.bart)
#! .217 to .138


Yh3 <- inla(Y ~ offset(log(E))
              + f(i , model='iid')
              + f(io, model='iid')
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Yh3)

Yh3.bart <- inla(Y ~ offset(log(E)) + smr0
              + f(i , model='iid')
              + f(io, model='iid')
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Yh3.bart)

#! smr0 exaplains most of the variation across counties .23 to .14
#!and the residual variation .03 to .017

Ypo$dic$dic
Yh$dic$dic
Yh3$dic$dic
Yh3.bart$dic$dic

##3.- spatio-temporal
Ysp.s <- inla(Y ~ offset(log(E))
              + f(i , model="bym2", graph = gs,
                  group = yr.i,control.group = list(model="ar1"))
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Ysp.ip)
save(Ysp.s,file=file.path(path.h,'Ysp.s',sep=''))
Ysp.s$dic$dic
Ysp.s$waic$waic

Ysp.s.bart <- inla(Y ~ offset(log(E)) + smr0
              + f(i , model="bym2", graph = gs,
                  group = yr.i,control.group = list(model="ar1"))
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Ysp.s.bart)
save(Ysp.s.bart,file=file.path(path.h,'Ysp.s.bart',sep=''))
Ysp.s$dic$dic
Ysp.s$waic$waic


# This feature (is documented in
# Martins  et al (2013) https://doi.org/10.1016/j.csda.2013.04.014 '4.6 Kronecker feature'
# e.g., spatially correlated innovations in a AR(1) model
# y ~ f(location, model = "besag", group = time, control.group = list(model = "ar1"))
# (they missed the graph argument requiere for besag!)
# similar use in Blangiardo et al.
# also Bakka et al. 2018 https://doi.org/10.1002/wics.1443


#===============================================================================
##4.- spatio-temporal +  adults
#===============================================================================
# equivalent to separate models
# adults and  youth have their own RE with ST structure
# adult RE loads on youth

dat$i.y     <-  dat$i
dat$i.a_s   <-  dat$i
dat$i.a_r   <-  dat$i


dat$i.y     [dat$opop=='adult']    <-   NA
dat$i.a_s   [dat$opop=='youth']    <-   NA  #shared
dat$i.a_r   [dat$opop=='youth']    <-   NA  #not shared



Yspya.s.bart <- inla(Y ~ offset(log(E)) + smr0
              + f(i.y , model="bym2", graph = gs,
                  group = yr.i,control.group = list(model="ar1"))
              + f(i.a_s , copy="i.y" , fixed = FALSE)
              + f(i.a_r , model="bym2", graph = gs,
                  group = yr.i,control.group = list(model="ar1"))
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )

summary(Yspya.s.bart)
save(Yspya.s.bart,file='data/inla/Yspya.s.bart')

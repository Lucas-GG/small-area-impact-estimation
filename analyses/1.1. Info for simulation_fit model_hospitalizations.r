options(digits=3)
options(scipen=10000)

 
#===============================================================================
#.- Obtain parameters for simulation from case data
#===============================================================================
path.h <- 'C:/LOCAL/ICF/R03 SAE/DATA'
load(file.path(path.h,'ldat.r'))
load(file.path(path.h,'gs'))
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
dat <- ldat[ldat$otyp=='IP',] 
dim(dat)



#===============================================================================
options(java.parameters = "-Xmx20g")  
library(bartMachine)
set_bart_machine_num_cores(4)

#-------------------------------------------------------------------------------

dat$smr <- dat$Y/dat$E
summary(dat$smr) 
min(dat$smr[dat$smr>0],na.rm=T)
dat$lsmr <- log(dat$smr)
summary(dat$lsmr)
dat$lsmr[dat$lsmr ==-Inf] <- NA
hist(dat$lsmr)
hist((dat$lsmr)^(1/3))

fy <- formula(lsmr ~ E
             + yr #+ fyr 
             + st  #+ region + division 
             + urb13 #+ rural
             + I(OTH/P1024)   +  I(OTHA/P25)
             + AIAN + uemp + mhinc + pui)

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


#-------------------------------------------------------------------------------
library(lme4)

qpm.ip <- glm(Y ~ offset(log(E))  ,family=quasipoisson,
            data=dat)
summary(qpm.ip)

qpm.ip.smr0 <- glm(Y ~ offset(log(E)) + smr0,family=quasipoisson,
            data=dat)
summary(qpm.ip.smr0) #dispersion is just 22%


hm.ip <-  glmer(Y ~ offset(log(E)) + (1|i),family=poisson,
            data=dat)
summary(hm.ip)

hm.ip.smr0 <-  glmer(Y ~ offset(log(E)) + smr0 + (1|i),family=poisson,
            data=dat)
summary(hm.ip.smr0) #varaince is just 14% (from .58 to .17)

hm.ed <-  glmer(Y ~ offset(log(E)) + (1|i),family=poisson,
            data=ldat[ldat$otyp=='ED',])
summary(hm.ed)
            

npm.ip <-  glmer.nb(Y ~ offset(log(E)) + (1|i) ,
            data=dat)
summary(npm.ip)
getME(npm.ip, "glmer.nb.theta")

npm.ip.smr0 <-  glmer.nb(Y ~ offset(log(E))+smr0 + (1|i) ,
            data=dat)
summary(npm.ip.smr0) #varaince  (from .57 to .09)
getME(npm.ip.smr0, "glmer.nb.theta")  #50% more residual


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
#! .61 to .11


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
# similar form .265 to .08


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

#! smr0 exaplains most of the variation across counties .26 to .05
#!and some of the residual variation .06 to .037

Ypo$dic$dic
Yh$dic$dic
Yh3$dic$dic
Yh3.bart$dic$dic

##3.- spatio-temporal
Ysp.ip <- inla(Y ~ offset(log(E))
              + f(i , model="bym2", graph = gs,
                  group = yr.i,control.group = list(model="ar1"))  
             ,data=dat,
             ,family="poisson" 
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Ysp.ip)
save(Ysp.ip,file=file.path(path.h,'Ysp.ip',sep=''))
Ysp.ip$dic$dic
Ysp.ip$waic$waic

Ysp.ip.bart <- inla(Y ~ offset(log(E)) + smr0
              + f(i , model="bym2", graph = gs,         
                  group = yr.i,control.group = list(model="ar1"))  
             ,data=dat,
             ,family="poisson" 
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Ysp.ip.bart)
save(Ysp.ip.bart,file=file.path(path.h,'Ysp.ip.bart',sep=''))
Ysp.ip$dic$dic
Ysp.ip$waic$waic


# This feature (is documented in 
# Martins  et al (2013) https://doi.org/10.1016/j.csda.2013.04.014 '4.6 Kronecker feature'
# e.g., spatially correlated innovations in a AR(1) model
# y ~ f(location, model = "besag", group = time, control.group = list(model = "ar1"))
# (they missed the graph argument requiere for besag!)
# similar use in Blangiardo et al. 
# also Bakka et al. 2018 https://doi.org/10.1002/wics.1443






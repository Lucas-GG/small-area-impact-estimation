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
ldat$otyp <- factor(ldat$otyp) # emergency vist (ED) or hospitalization (IP)
ldat$fid  <- factor(ldat$id)
ldat$fyr  <- factor(ldat$year)
ldat$rural <-  ldat$urb13%in%c('6.NonCore','5.Micro')

dim(ldat)
dat <- ldat[ldat$otyp=='IP',]#&ldat$opop=='adult',] #hospitalization (IP) only
dat <- droplevels(dat)
dim(dat)

table(ldat$st[!is.na(ldat$Y)])
table(dat$st[!is.na(dat$Y)])

#===============================================================================
options(java.parameters = "-Xmx20g")
library(bartMachine)
set_bart_machine_num_cores(6)

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


fy <- formula(lsmr ~ E
             + yr #+ fyr
             + st  #+ region + division
             + urb13 #+ rural
             + I(OTH/P) #hospitalizations not MH
             + AIAN + uemp + mhinc + pui
             + opop
           )

summary(lm(fy,data=dat)) #R2~.40
yX <- model.frame(fy,data=dat[!is.na(dat$lsmr),],na.action=na.pass)


yX <- yX[!is.na(yX[, 1]),]
summary(yX)
bm <- bartMachine(yX[,-1],yX[, 1],use_missing_data =T,serialize = TRUE)
X <- model.frame(fy,data=dat,na.action=na.pass)[,-1]


dat$smr0 <- NULL
dat$smr0 <- predict(bm,X)

plot(dat$smr0,dat$lsmr)
summary(lm(lsmr ~ smr0,data=dat)) #R2~.64

#===============================================================================
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

npm.ip.smr0 <-  glmer.nb(Y ~ offset(log(E)) + smr0 + (1|i) ,
            data=dat)
summary(npm.ip.smr0) #varaince  (from .57 to .09)
getME(npm.ip.smr0, "glmer.nb.theta")  #50% more residual


#===============================================================================
library(INLA)
#.- Obtain parameters for simulation from case data

#!predict missing Y takes time
dat <- dat[!is.na(dat$Y),]
#sdat <- dat[dat$st=='AZ'&dat$opop=='youth'&!is.na(dat$Y),]
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
dat <- dat[sample(nrow(dat),nrow(dat)),]
##3.- spatio-temporal
Ysp0.ip <- inla(Y ~ offset(log(E))
              + f(i , model='bym', graph = gs)
             ,data=dat ,
             ,family="poisson")
summary(Ysp0.ip)



Ysp.ip <- inla(Y ~ offset(log(E))
              + f(i , model="bym", graph = gs,
                  group = yr.i,control.group = list(model="ar1"))
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Ysp.ip)



#save(Ysp.ip,file=file.path(path.h,'Ysp.ip',sep=''))
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

table(dat$opop)

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
#===============================================================================



#===============================================================================
#===============================================================================
##4.- spatio-temporal +  adults
# the youth  RE loads on adult RE,
#the adutl resiudal RE is not spatial


dat$i.y     <-  dat$i
dat$i.a_s   <-  dat$i
dat$i.a_r   <-  dat$i#1:nrow(dat) #residual unestructured

dat$i.y     [dat$opop=='adult']    <-   NA
dat$i.a_s   [dat$opop=='youth']    <-   NA  #shared
dat$i.a_r   [dat$opop=='youth']    <-   NA  #not shared


Yspal.ip.bart <- inla(Y ~ offset(log(E)) + smr0
              + f(i.y , model='bym2', graph = gs,
                  group = yr.i,control.group = list(model="ar1"))
              + f(i.a_s , copy='i.y' , fixed = FALSE)
              + f(i.a_r , model='iid')
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )

summary(Yspal.ip.bart)
save(Yspal.ip.bart,file='data/inla/Yspal.ip.bart')
#===============================================================================
#===============================================================================
# equivalent to separate models
# adults and  youth have their own RE with ST structure
# adult RE loads on youth

dat$i.a     <-  dat$i
dat$i.y_s   <-  dat$i
dat$i.y_r   <-  dat$i


dat$i.a     [dat$opop!='adult']    <-   NA
dat$i.y_s   [dat$opop!='youth']    <-   NA  #shared
dat$i.y_r   [dat$opop!='youth']    <-   NA  #not shared

table(is.na(dat$i.a),is.na(dat$i.y_s))
table(is.na(dat$i.y_s),is.na(dat$i.y_s))



Yspa.ip.bart <- inla(Y ~ offset(log(E)) + smr0
              + f(i.a , model="bym2", graph = gs,
                  group = yr.i,control.group = list(model="ar1"))
              + f(i.y_s , copy="i.a" , fixed = FALSE)
              + f(i.y_r , model="bym2", graph = gs,
                  group = yr.i,control.group = list(model="ar1"))
             ,data=dat,
             ,family="poisson"
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )

summary(Yspa.ip.bart)
save(Yspa.ip.bart,file='data/inla/Yspa.ip.bart')
load('data/inla/Yspa.ip.bart')
summary(Yspa.ip.bart)

#===============================================================================
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



Yspya.ip.bart <- inla(Y ~ offset(log(E)) + smr0
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

summary(Yspya.ip.bart)
save(Yspya.ip.bart,file='data/inla/Yspya.ip.bart')
load('data/inla/Yspya.ip.bart')
summary(Yspa.ip.bart)
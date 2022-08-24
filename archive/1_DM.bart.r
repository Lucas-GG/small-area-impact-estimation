options(digits=3)
options(scipen=10000)

source('C:/LOCAL/ECOLOGICAL/HCUP/R/DM3/aux fun.r')
library(spdep)
library(Rgraphviz)
library(INLA)
inla.setOption(scale.model.default = TRUE)

path.h <- 'C:/LOCAL/ECOLOGICAL/HCUP/DATA/'
#===============================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#data

load(file=paste(path.h,'ldat.r',sep=''))
load(file=paste(path.h,'gs',sep=''))

length(unique(ldat$id))

#fit a model for Y1.
#so it can borrow from other treated counties (even large ones)
#and over time
#look at models to detect local diffreces
#what is it with the single red


tapply(ldat$E,ldat$urb13,mean)

colnames(ldat$state)

#gls$yr1 <- ave(ifelse(gls$anytr==1,gls$year,9999),gls$id,FUN=min)
ldat$evertr <- as.numeric(ldat$year>=ldat$yr1)
ldat$Y0     <- ifelse(ldat$evertr==0,ldat$Y,NA)

with(ldat,tapply(evertr,list(id,year),mean))[1:20,]
with(ldat,tapply(evertr.s,list(id,year),mean))[1:20,]


data(state.fips, package='maps')
st <- with(state.fips,tapply(abb,fips,max))
ldat$st   <- st[paste(ldat$state )]
ldat$st   <- factor(ldat$st)
ldat$otyp <- factor(ldat$otyp)
ldat$fid  <- factor(ldat$id)
ldat$fyr  <- factor(ldat$year)
ldat$rural <-  ldat$urb13%in%c('6.NonCore','5.Micro')

dim(ldat)
ldat[1:20,c('id','year','otyp')]
length(unique(ldat$id))*length(unique(ldat$year))*2

#===============================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options(java.parameters = "-Xmx20g")  
library(bartMachine)
set_bart_machine_num_cores(4)



#===============================================================================
#1.- modelling propensity
#===============================================================================
table(ldat$year)
round(tapply(ldat$id,ldat$yr1 , function(x) length(unique(x)))/516*100,1)
summary(ldat$urb13)
colnames(ldat)


summary(ldat$E)

f <- formula( ~ yr1 + id 
             + st    + region + division + urb13 + rural
             + I(Y/E)  + otyp + E 
             + I(S1024/P1024) + I(S25/P25)
             + I(OTH/P1024)   + I(OTHA/P25)
             + AIAN  + uemp + mhinc + pui
             + year)
mdat <- model.frame(f,data=ldat,na.action=na.pass)
head(mdat)
id <- unique(mdat$id)

#!only background??!can it inlcude lag outcome?

P.list <- lapply(c(2006:2016), function(year){
  #year <- 2009
  bdat <- mdat[mdat$year<year,]
  bdat[bdat$year>2006,c('I(Y/E)','I(S1024/P1024)','I(S25/P25)')] <- NA #only pre-2006 values of the outcome     
  #bdat[bdat$year==year,c('I(Y/E)','I(S1024/P1024)','I(S25/P25)')] <- NA #only lag values    
  head(bdat)
  wide <- reshape(bdat,  idvar = c('yr1','st','id','year'), timevar = 'otyp', direction = 'wide',
  v.names =c('I(Y/E)','E', 'I(OTH/P1024)','I(OTHA/P25)'))
  wdat <- reshape(wide ,  idvar = c('yr1','st','id'), timevar = 'year', direction = 'wide')
  head(wdat)
  dim(wdat)
  X    <- wdat[,!colnames(wdat)%in%c('yr1','id')]
  head(X)
  y    <-  factor((wdat$yr1)!=year,labels =c('1','0')) 
  table(y)
  bm.p <- bartMachine(X,y,use_missing_data =T,serialize = TRUE)
  P   <- cbind(id=wdat$id,p=predict(bm.p,X))
  colnames(P)[2] <- paste('e',year,sep='.')
  result <- list(bm.p,merge(data.frame(id),P,all.x=T))
  return(result)
  })
names(P.list) <- paste(c(2006:2016))

P <- Reduce(function(x,y) {merge(x,y)}, lapply(P.list, function(x) x[[2]])  ) #vectors of propensities
summary(P)
dim(P)
PM <- lapply(P.list, function(x) x[[1]])

save(P,file=paste(path.h,'P.r',sep=''))
save(PM,file=paste(path.h,'PM.r',sep=''))

investigate_var_importance(P.list[['2010']][[1]], num_replicates_for_avg = 20)
vs <- var_selection_by_permute(P.list[['2010']][[1]], bottom_margin = 10, num_permute_samples = 10)
vs$important_vars_local_names
vs$important_vars_global_max_names
vs$important_vars_global_se_names

pd_plot(P.list[['2010']][[1]], j = 'I(S25/P25).2003')


#===============================================================================
#2.- Modelling response surface
#===============================================================================
ldat <- merge(ldat,P,all.x=T,sort = F)
head(ldat)
dim(ldat)

ldat$smr <- ldat$Y0/ldat$E
summary(ldat$smr) 
min(ldat$smr[ldat$smr>0],na.rm=T)
ldat$lsmr <- log(ldat$smr)
summary(ldat$lsmr)
ldat$lsmr[ldat$lsmr ==-Inf] <- NA
hist(ldat$$lsmr)
hist((ldat$$lsmr)^(1/3))

fy <- formula(lsmr ~ E 
             + otyp 
             + fyr + yr
             + st   + region + division + urb13 + rural
             #+ I(Y/E) + I(S1024/P1024) + I(S25/P25)
             + I(OTH/P1024)   +  I(OTHA/P25)
             + AIAN + uemp + mhinc + pui
             + e.2006 + e.2007 + e.2008 + e.2009 + e.2010
             + e.2011 + e.2012 + e.2013 + e.2014 + e.2015 + e.2016
             )


#we can not used lag Y, becouse it will be missing after the start for all years
#we are using y pre 2006 trhough e
#we are including all propsenistied (even future propensities
#these are transformations of (contemporary vaues of) X



             
summary(lm(fy,data=ldat))
yX <- model.frame(fy,data=ldat[!is.na(ldat$lsmr),],na.action=na.pass)
yX <- yX[!is.na(yX[, 1]),]
summary(yX)
y  <- yX[, 1]
X  <- yX[,-1]
summary(yX)
bm <- bartMachine(X,y,use_missing_data =T,serialize = TRUE)
save(bm,file=paste(path.h,'bm.r',sep=''))

X <- model.frame(update(fy, Y~.),data=ldat,na.action=na.pass)[,-1]
dim(X)
ldat$smr0 <- NULL
ldat$smr0 <- predict(bm,X)

plot(ldat$smr0,ldat$lsmr)
 

save(ldat,file=paste(path.h,'ldatb.r',sep=''))


#===============================================================================
#3.- random part
#===============================================================================
load(file=paste(path.h,'ldatb.r',sep=''))

summary(lm(lsmr~ 1,data=ldat))
summary(lm(lsmr~ smr0,data=ldat))

m0 <- glm(Y0 ~ offset(log(E)),data=ldat, family='quasipoisson')
m1 <- glm(Y0 ~ offset(log(E)) + smr0,data=ldat, family='quasipoisson')
anova(m0,m1,test ='LRT')
names(m1)

m2 <- glm(Y0 ~ offset(log(E)) -1 + factor(i) + factor(year),data=ldat, family='quasipoisson')
anova(m1,m2,test ='LRT')

m3 <- glm(Y0 ~ offset(log(E)) + scale(AIAN) + scale(OTH/P1024)  + scale(OTHA/P25) + scale(uemp) + scale(mhinc) + scale(pui),data=ldat, family='quasipoisson')
anova(m1,m3,test ='LRT')

length(m1$y)
length(m2$y)
length(m3$y)

#-------------------------------------------------------------------------------

dat <- ldat[!is.na(ldat$Y),] 


dat$i.ip         <-  dat$i 
dat$yr.i         <-  as.numeric(as.factor(dat$yr))
dat$i.ed         <-  dat$i
dat$i.ed.res     <-  dat$i 
dat$i.ip     [dat$outcome!='ED']  <-   NA
dat$i.ed     [dat$outcome!='IP']  <-   NA
dat$i.ed.res [dat$outcome!='IP']  <-   NA

dat$io <- 1:dim(dat)[1]

##1.- overdispersed poisson
Ypo <- inla(Y0 ~ offset(log(E)) + smr0
              + f(io, model='iid')     
             ,data=dat,
             ,family="poisson" 
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Ypo)

##2.- hierarchical
#hirachical with overdispersion
Yh3 <- inla(Y0 ~ offset(log(E)) + smr0
              + f(i , model='iid')
              + f(io, model='iid')     
             ,data=dat,
             ,family="poisson" 
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Yh3)

##3.- spatio-temporal

Yspo <- inla(Y0 ~ offset(log(E)) + smr0
              + f(i , model="bym2", graph = gs,        #Leroux besagproper2
                  group = yr.i,control.group = list(model="ar1"))  
              + f(io, model='iid')   #extra randomness
             ,data=dat,
             ,family="poisson" 
             ,control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config = TRUE)
             ,control.predictor=list(compute=TRUE, link=1)
             )
summary(Yspo)




Ypo$dic$dic
Yh3$dic$dic
Yspo$dic$dic

Ypo$waic$waic
Yh3$waic$waic
Yspo$waic$waic


Ypo$cpu.used['Total']
Yh3$cpu.used['Total']
Yspo$cpu.used['Total']

save(Ypo, file=paste(path.h,'Ynn', sep='')) 
save(Yh3,file=paste(path.h,'Ynsp',sep='')) 
save(Yspo,   file=paste(path.h,'Y',   sep='')) 


 

   
#===============================================================================
#4.- Get posterior draws of Y for inference
#===============================================================================
load(file=paste(path.h,'Y',   sep='')) 
m  <- Yspo #spatio-temporal poisson

  
dat$eY0  <- NULL
dat      <- cbind(dat,eY0 = m$summary.fitted.values) 
dat$pps  <- mapply(function(Y,Y0){1-inla.pmarginal(Y, Y0)}, dat$Y, m$marginals.fitted.values)

K <- 10^4
Y00 <- t(sapply(m$marginals.fitted.values,function(x) inla.rmarginal(K,x)))
Y0 <- apply(Y00,2,function(x) sapply(x,function(l) rpois(1,l)))
cor(m$summary.fitted.values$mean,apply(Y0,1,mean))
cor(m$summary.fitted.values$'0.025quant',apply(Y0,1,quantile, .025))

dat[,paste('Y0',1:K,sep='.')] <- NULL
dat <- cbind(dat,Y0=Y0)
dim(dat)
colnames(dat) 
save(dat,file=paste(path.h,'datp.r',sep=''))






#===============================================================================
#5.- CATT (conditional average treament effect among treated)
#===============================================================================
load(file=paste(path.h,'datp.r',sep=''))

library(partykit)
#treated after start , small
pdat <- dat[dat$W==1&dat$urb13%in%c('5.Micro','6.NonCore'),]    
y <- tapply(pdat$delta,pdat$gyr,mean)

yci <- sapply(split(pdat,pdat$gyr),function(D) {
 X <- apply((D$Y-D[,paste('Y0',1:K,sep='.')])/D$E,2,mean) #E(Y) for each simulation
 c(mean(X), quantile(X,c(0.025,.975)))
 })
yl <- yci[2,]
yu <- yci[3,]
x <- -7:10
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(1,1))
plot(c(-7,10), c(-1.6,1.6), type="n",
  ylab='Stadardized hospital use ratio', 
  xlab='',axes=F)
axis(1)
axis(2)
boxplot(pdat$delta~pdat$gyr,add=T,at=x,axes=F,border='gray',outline=F,medcol=0,staplewex=0,whisklty=1)
abline(v=0,col='red',lty=3)
abline(h=0,col='red',lty=3)
polygon(c(x,rev(x)),c(yu,rev(yl)),col = adjustcolor('steelblue', alpha.f = .3),border=NA)    
lines(x,y,col = "black" , lwd = 2, cex = 4/5)
  text(c(-3,3),-1.5,c('Before GLS','After GLS'))  
  mtext("Year",side=1,line=2.75,at=0)


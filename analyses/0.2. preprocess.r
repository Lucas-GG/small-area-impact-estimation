
#===============================================================================
#2.- wide to long
#===============================================================================
dat <- readRDS('data/dat.rds')
head(dat)
colnames(dat)
# Y:youth; A: adult
# ed: ED
# SH: suicide
# OTH: non-suicide death

#ldat <- reshape(dat,
#  idvar = c('year','id','statefips'),
#  times = c('IP','ED'),
#  timevar = 'otyp',
#  varying =  list(c('SHY','SHYed'),c('OTHERY','OTHERYed'),c('OTHERA','OTHERAed')),
#  v.names =  c('Y','OTH','OTHA'),
#  direction = 'long')

colnames(dat)

dat$TOTP <- with(dat, rowSums(cbind(P009,P1024,P25)))

ldat <- reshape(dat,
  idvar = c('year','id','statefips'),
  times = c('IPY','IPA','EDY','EDA'),
  timevar = 'outcome',
  varying =  list(
    c('SHY','SHA','SHYed','SHAed'),
    c('S1024','S25','S1024','S25'),
    c('P1024','P25','P1024','P25'),
    c('OTHERY','OTHERA','OTHERYed','OTHERAed')
    ),
  v.names =  c('Y','S','P','OTH'),
  direction = 'long')

ldat$opop <-  ifelse(ldat$outcome=='IPY'|ldat$outcome=='EDY','youth','adult')
ldat$otyp <-  ifelse(ldat$outcome=='IPY'|ldat$outcome=='IPA','IP','ED')

head(ldat)
table(ldat$year)

# expected count
ldat$risk0 <- ave(ldat$Y,ldat$outcome,FUN=function(x) sum(x,na.rm=T))/
  ave(ifelse(is.na(ldat$Y),0,ldat$P),ldat$outcome,FUN=function(x) sum(x,na.rm=T))
ldat$E <- ldat$P*ldat$risk0

tapply(ldat$E,ldat$outcome,summary)

plot(ldat$E,ldat$Y)
tapply(ldat$Y/ldat$E,ldat$otyp,mean,na.rm=T)
tapply(ldat$Y/ldat$E,ldat$otyp,sd,na.rm=T)
tapply(ldat$Y/ldat$E,ldat$otyp,IQR,na.rm=T)


summary(ldat$year)
ldat$yr  <-  (ldat$year-2012)/9
ldat$yr2  <-  ldat$yr^2

#ldat[,c('yr.1','yr.2')]  <-  poly(ldat$yr,2)
#ldat <- ldat[!is.na(ldat$Y),]

summary(ldat$Y)

#===============================================================================
#sort data so it agrees with GS
library(spdep)
load('data/hmap')
nb <- poly2nb(hmap)

length(unique(attributes(nb)$region.id))
length(unique(ldat$id))

i <- 1:length(attributes(nb)$region.id);
names(i) <-  attributes(nb)$region.id
cbind(i[paste(ldat$id)],as.numeric(as.factor(ldat$id)))[1:20,]

ldat$i <-  i[paste(ldat$id)]
ldat$t  <- ldat$yr
ldat <- ldat[order(ldat$i,ldat$t),]

head(ldat)
summary(ldat$i)
length(unique(ldat$i))
max(ldat$i)

saveRDS(ldat,'data/ldat.rds')
#save(ldat,file=paste(path.h,'ldat.r',sep=''))

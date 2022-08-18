path <- file.path('C:/LOCAL/ICF/R03 SAE/')   
source(file.path(path,'R/aux fun.r'))
library(spdep)
library(Rgraphviz)
library(INLA)
#===============================================================================
data <- read.csv(file.path(path,'DATA/hcup_plus.csv'))
data$state
summary(data$state)

summary(data)
head(data)

data[is.na(data$P1024),c('id','year')]
data[data$id=='46113',c('id','year','P1024')]
data$P1024[data$id=='46113'&data$year==2016] <- 3902 #renamed Oglala 46102
data$P25[data$id=='46113'&data$year==2016]   <- 7335
dim(data)


length(unique(paste(data$id,data$year)))
data <- data [ !ave(!is.na(data$id),data$id,data$year,FUN=cumsum)>1,]
length(unique(paste(data$id,data$year)))

with(data,tapply(!is.na(SHYed),list(data$state,data$year),sum))
with(data,tapply(!is.na(SHY),list(data$state,year),sum))

#c(4,21,25,31,32)
#c(4,21,   ,26,31*,32,35,41,46,53)

#4,21,  ,32)
#AZ,KY,MA,NE,NV

#SD (SID)
#!AZ (SID and SEED)
#OR (SID)
#MI (SID)
#WA (SID)
#!NE (SID and SEED)
#!NV (SID and SEED)
#NM (SID)
#!KY (SID and SEED)
#!MA (SEED)

save(data,file=file.path(path,'dat'))
#===============================================================================
#Proximity matrix

#!whose elements are proximities or neighborhoods of the areas, i.e., a matrix 
#!with elements in [0,1], 
#!zeros on the diagonal
#!and rows adding up to 1.

rmaps <-  'C:/LOCAL/TRANSFER/MAPS/DATA/R.MAPS/'
load(paste(rmaps,'county.mapp',sep=''))
load(paste(rmaps,'state.map',sep=''))

plot(county.mapp)
summary(county.mapp)
length(unique(data$id))
table(county.mapp$COUNTYFIPS%in%unique(data$id))

plot(county.mapp[county.mapp$COUNTYFIPS%in%unique(data$id),] )

hmap <- county.mapp[county.mapp$COUNTYFIPS%in%unique(data$id),] 
plot(state.map, border = "gray")
plot(hmap, col = "blue",add = TRUE)


nb <- poly2nb(hmap)
nb2 <- correct.island(nb,hmap)
head(nb)
attributes(nb)$region.id


plot(state.map, border = "gray")
coords <- coordinates(hmap)
plot.nb(nb2, coords = coords,add = TRUE, cex = 0.5,col='red')
plot.nb(nb, coords = coords,add = TRUE, cex = 0.5)
title("polygon generated queen weights")


W <- nb2mat(nb2, style = "B")
gs <- inla.read.graph(W)
plot(gs)
head(gs)

save(gs,file=file.path(path,'gs'))
save(W,file=file.path(path,'W'))
save(hmap,file=file.path(path,'hmap'))


#===============================================================================
#-------------------------------------------------------------------------------
#wide to long
#-------------------------------------------------------------------------------

#head(ldat,50)

ldat <- reshape(data, idvar = c('year','id','statefips'),
                times = c('IP','ED',),
                timevar = 'otyp',       
                varying =  list(c('SHY','SHYed'),c('OTHERY','OTHERYed'),c('OTHERA','OTHERAed')),
                v.names =  c('Y','OTH','OTHA'),
                direction = 'long')

table(ldat$year)
ldat$E <- ldat$P1024* ave(ldat$Y,ldat$otyp,FUN=function(x) sum(x,na.rm=T))/ave(ifelse(is.na(ldat$Y),0,ldat$P1024),ldat$otyp,FUN=function(x) sum(x,na.rm=T)) 
tapply(ldat$E,ldat$otyp,summary)

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

i <- 1:length(attributes(nb)$region.id)
names(i) <-  attributes(nb)$region.id 
cbind(i[paste(ldat$id)],as.numeric(as.factor(ldat$id)))[1:20,]

ldat$i <-  i[paste(ldat$id)]
ldat$t  <- ldat$yr
ldat <- ldat[order(ldat$i,ldat$t),]

head(ldat)
summary(ldat$i)
length(unique(ldat$i))
max(ldat$i)


write.csv(ldat, paste(path.h,'ldat',sep=''))
save(ldat,file=paste(path.h,'ldat.r',sep=''))
save(gs,file=paste(path.h,'gs',sep=''))
save(W,file=paste(path.h,'W',sep=''))
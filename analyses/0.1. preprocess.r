#path <- file.path('C:/LOCAL/ICF/R03 SAE/')
#source(file.path(path,'R/aux fun.r'))
source('R/aux fun.r')
library(spdep)
library(Rgraphviz)
#library(INLA)

#===============================================================================
#1.- Proximity matrix
#===============================================================================
data <- readRDS('data/dat.rds')
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

library(INLA)
W <- nb2mat(nb2, style = "B")
gs <- inla.read.graph(W)
plot(gs)
head(gs)


save(gs, file = file.path(path, 'gs'))
save(W, file = file.path(path, 'W'))
save(hmap, file = file.path(path, 'hmap'))

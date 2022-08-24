#path <- file.path('C:/LOCAL/ICF/R03 SAE/')
#source(file.path(path,'R/aux fun.r'))

data <- read.csv('data/hcup_plus.csv')
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

saveRDS(data,'data/dat.rds')
#save(data,file=file.path(path,'dat'))
#===============================================================================

library(INLA)
load('data/inla/Yspa.ip.bart')
summary(Yspa.ip.bart)


library(inlabru)
library(sp)
#Matérn Gaussian field
#Stochastic Partial Diferential Equation
#setwd('C:/LOCAL/ST_evaluation_case')
load('data/hmap')

#The function takes the convex hull of the region and extends that to create the region for the “inner” mesh,
#with small triangles.
#Try changing the parameters of inla.mesh.2d to create different meshes, you can:
# decide how much the region should be extended by changing the parameter offset
# increase/decrease the density of the mesh by changing max.edge and cutoff

coords <- coordinates(hmap)
#mesh <- inla.mesh.2d(coords, max.edge = c(1, 5), cutoff = 1)
hull <- inla.nonconvex.hull(
  points = coords,
  convex = 2, concave = 3.5
)
mesh <- inla.mesh.2d(
  boundary = hull, max.edge = c(1, 6), # km inside and outside
  cutoff = 1, offset = c(1, 3)
) # cutoff is min edge

plot(mesh)
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")
#A    <- inla.spde.make.A (mesh, loc = coords)
#spde <- inla.spde2.matern(mesh, alpha = 2)
spde <- inla.spde2.pcmatern(mesh,
  prior.range = c(1, 0.5),# the probability of a range exceeding 500 km is 0.5
  prior.sigma = c(1, 0.5)) # the probability of a standard deviation exceeding 1 is 0.5

#===============================================================================
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

dat$yr.i         <-  as.numeric(as.factor(dat$yr))

library(splines)
#load('data/hmap')
#coords <- coordinates(hmap)
n.basis <- K^2
n.time  <- length(unique(dat$yr))
n.area  <- length(unique(dat$id))

B1 <- bs(coords[, 1], df = K)
B2 <- bs(coords[, 2], df = K)


N <- nrow(B1)
Z <- matrix(NA,N,K^2)
k <- 1
for(i in 1:K){
  for(j in 1:K){
    Z[,k]  <- B1[,i]*B2[,j]
    k <- k+1
  }
}
dim(Z)
X.mat <- matrix(NA, nrow=n.area*n.time, ncol=n.basis*n.time)
I <- diag(1, nrow=n.time, ncol=n.time)
for (i in 1:n.area){
  X.mat[((i-1)*(n.time)+1):(i*n.time),] <- I %x% matrix(Z[i,], nrow=1)
}
dim(X.mat) # should be (IT)*(KT)



#########--------------- type 2 interaction ---------------#########
# create vectors in the A matrix
intercept <- c(1,  rep(NA, n.basis*n.time))
idx       <- c(NA, rep(1:n.time, each=n.basis)) # this is the common id
group     <- c(NA, rep(1:n.basis, n.time)) #!!!!!

data.inla <- list(y=ldat$y, idx=idx, intercept=intercept, group=group,E=ldat$E)
summary(data.inla)

# fit the model
formula.2 <- y ~ -1 + intercept +
  f(idx, model = "rw2", constr=F, hyper = list(prec = list(param=c(1, 0.005))), replicate=group)

# f(idx , model="bym2", graph = gs,
#    group = yr.i,control.group = list(model="ar1"))

r.2 = inla(formula.2, data=data.inla,family="poisson",
  control.predictor =list(A=cbind(1, X.mat),compute=TRUE, link=1),
  control.compute=list(dic=T, cpo=TRUE),verbose=F, E=E)

# extract the coefficients
b.hat.2.inla <- r.2$summary.linear.predictor[(n+2):(n+1+m), "mean"]
b.hat.mat.2.inla <- matrix(b.hat.2.inla, nrow=n.basis, ncol=n.time, byrow=F)
head(b.hat.mat.2.inla)













#===============================================================================
setwd('C:/Users/Lucas/R/win-library/4.1/INLA')
data(Oral)
g = system.file("demodata/germany.graph", package="INLA")
## use data Oral to estimate a spatial field in order to simulate a
## 'realistic' dataset.
formula = Y ~ f(region, model="bym", graph=g)
result = inla(formula, data = Oral, family = "poisson", E = E)
x = result$summary.random$region$mean
n = length(x)/2
## simulate two new datasets. 'a' is the scaling between the
## log.rel.risk:
a = 2
xx = x[1:n]+1
x = c(0 + a*xx, 1 + xx/a)
E = c(Oral$E, Oral$E)
N = 2*n
y = rpois(N, lambda = E*exp(x))


## model='besag2' defines a model with length N = 2*graph->n, the
## first half is weighted with 'a' the other half is weighted with
## 1/a. here there is no unstructed terms.
idx = 1:N
mu = as.factor(rep(1:2, each=n))
formula = y ~ -1 + mu + f(idx, model="besag2", graph=g, scale.model=TRUE)
r = inla(formula, family = "poisson", data = data.frame(E,y,idx,mu)[sample(length(mu),100),], E=E, verbose=TRUE,
control.compute=list(return.marginals.predictor=TRUE))
summary(r)




plot(m1, asp = 1, main = "")

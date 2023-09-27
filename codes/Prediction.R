rm(list = ls())

library(PrevMap)
library(sf)

# Extract parameters from geostatistical model
fit.MCML.1 <- readRDS("./outputs/fit.MCML.1.rds")
print(summary(fit.MCML.1,log.cov.pars=FALSE))
#           Estimate StdErr
# sigma^2 7.9635e+00 0.1177
# phi     1.7982e+05 0.0000
# tau^2   7.0857e-01 1.4707

#import jittered coords
baseline.jits <- read.csv("./outputs/baseline.jits.csv")


#Par estimates: 
beta   <- -1.5021    # intercept
sigma2 <- 7.9635e+00 # variance of the Gaussian process
phi    <- 1.7982e+05 # scale of the spatial correlation
tau2   <- 7.0857e-01 # variance of the nugget effect
#coords<- jitter2d(data[,c("long","lat")],max=100)
coords <- baseline.jits[,c("long","lat")]
object <- fit.MCML.1 
n <- length(object$y)

# covariance matrix at observed locations
U<- as.matrix(dist(coords))
Sigma.obs <- sigma2*exp(-U/phi) 

diag(Sigma.obs) <- diag(Sigma.obs)+tau2
Sigma.obs.inv <- solve(Sigma.obs)


# covariance matrix of the grid locations
#library(sf)
grid.pred <- read.csv("./outputs/grid.pred.csv")
head(grid.pred)
grid.pred<-grid.pred[,2:3]

#U.grid.pred <- as.matrix(dist(st_coordinates(grid.pred))) 
U.grid.pred <- as.matrix(dist(grid.pred)) 
Sigma.pred <- sigma2*exp(-U.grid.pred/phi) 

#Cross correlation 
#U.grid.obs <- as.matrix(pdist(st_coordinates(grid.pred), coords))
U.grid.obs <- as.matrix(pdist(grid.pred, coords))
C <- sigma2*exp(-U.grid.obs/phi) # cross-correlation


#S.samples
object$mu <- object$D %*% beta
mcmc.1 <- control.mcmc.MCML(n.sim=50000,burnin=5000,thin=45,h=(1.65)/(nrow(D)^(1/6)))

S.samples <- Laplace.sampling(object$mu, Sigma.obs, object$y, 
                              object$units.m, 
                              control.mcmc=mcmc.1 , 
                              plot.correlogram = T, 
                              messages = T)

# conditional mean
n.sim <- 1000 # number of simulations
S.obs <- S.samples$samples
#S.obs <-S.sim.res$samples

# conditional mean at prediction locations
cond.mean <- sapply(1:n.sim, function(i) beta+C%*%Sigma.obs.inv%*%(S.obs[i,]-object$mu))



# conditional covariance
Sigma.cond <- t(chol(Sigma.pred-C%*%Sigma.obs.inv%*%t(C)))

# S.pred.samples <- sapply(1:n.sim, function(i) cond.mean[,i]+
#                            Sigma.cond%*%rnorm(nrow(st_coordinates(grid.pred))))
S.pred.samples <- sapply(1:n.sim, function(i) cond.mean[,i]+
                           Sigma.cond%*%rnorm(nrow(grid.pred)))

# prevalence samples on the grid
#prev.samples <- exp(beta_new+S.pred.samples)/(1+exp(beta.new+S.pred.samples))
prev.samples <- exp(S.pred.samples)/(1+exp(S.pred.samples))


#r<- rasterFromXYZ(cbind(st_coordinates(grid.pred), apply(prev.samples,1,mean)))
r<- rasterFromXYZ(cbind(grid.pred, apply(prev.samples,1,mean)))
plot(r)

library(sf)
map<-st_read(".\\raw_data\\map.shp")
map<- st_transform(map,3857)

plot(st_geometry(map),add=T)
writeRaster(r, ".\\figures\\Prediction_from_algorithm.tif")

r

#Problem: 
#Expected prevalence at far away places from survey locations
1/(1+exp(-beta)) #0.1821125 

#predicted prevalence at a far away place is much higher than 1/(1+exp(-beta))
points(2357592,-511841,col="red",pch=16)
extract(r,data.frame(x=2357592,y=-511841)) #0.5227119




#Prediction with PrevMap
?spatial.pred.binomial.MCML

pred.MCML <- spatial.pred.binomial.MCML(fit.MCML.1 ,
                                        #grid.pred = st_coordinates(grid.pred),
                                        grid.pred = grid.pred,
                                        control.mcmc=mcmc.1,
                                        scale.predictions = "prevalence")



#r.prevmap<- rasterFromXYZ(cbind(st_coordinates(grid.pred), pred.MCML$prevalence$predictions ))
r.prevmap<- rasterFromXYZ(cbind(grid.pred, pred.MCML$prevalence$predictions ))
plot(r.prevmap)
plot(st_geometry(map),add=T)
writeRaster(r.prevmap, ".\\figures\\Prediction_from_PrevMap.tif")

r.prevmap

par(mfrow=c(1,3))
plot(r,main="Prediction algorithm")
plot(st_geometry(map),add=T)

plot(r.prevmap,main="PrevMap")
plot(st_geometry(map),add=T)

plot(r-r.prevmap,main="Difference between\nPrediction algorithm and PrevMap")
plot(st_geometry(map),add=T)







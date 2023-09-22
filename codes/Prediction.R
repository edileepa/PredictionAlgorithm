
# Extract parameters from geostatistical model
fit.MCML.1 <- readRDS(".\\outputs\\fit.MCML.1.rds")
print(summary(fit.MCML.1,log.cov.pars=FALSE))
#           Estimate StdErr
# sigma^2 7.9635e+00 0.1177
# phi     1.7982e+05 0.0000
# tau^2   7.0857e-01 1.4707

#import jittered coords
baseline.jits <- read.csv(".\\outputs\\baseline.jits.csv")


#Par estimates: 
beta   <- -1.5021    # intercept
sigma2 <- 7.9635e+00 # variance of the Gaussian process
phi    <- 1.7982e+05 # scale of the spatial correlation
tau2   <- 7.0857e-01 # variance of the nugget effect
#coords<- jitter2d(data[,c("long","lat")],max=100)
coords <- baseline.jits
object <- fit.MCML.1 
n <- length(object$y)

# covariance matrix at observed locations
U<- as.matrix(dist(coords))
Sigma.obs <- sigma2*exp(-U/phi) 

diag(Sigma.obs) <- diag(Sigma.obs)+tau2
Sigma.obs.inv <- solve(Sigma.obs.inv)


# covariance matrix of the grid locations
U.grid.pred <- as.matrix(dist(st_coordinates(grid.pred))) 
Sigma.pred <- sigma2*exp(-U.grid.pred/phi) 

#Cross correlation 
U.grid.obs <- as.matrix(pdist(st_coordinates(grid.pred), coords))
C <- sigma2*exp(-U.grid.obs/phi) # cross-correlation


#S.samples
object$mu <- object$D %*% beta


S.samples <- Laplace.sampling(object$mu, Sigma.obs, object$y, 
                              object$units.m, 
                              control.mcmc=mcmc.1 , 
                              plot.correlogram = T, 
                              messages = T)

# conditional mean
n.sim <- 1000 # number of simulations
S.obs <- S.samples$samples
S.obs <-S.sim.res$samples

cond.mean <- sapply(1:n.sim, function(i) C%*%Sigma.obs.inv%*%(S.obs[i,]))

# conditional covariance
Sigma.cond <- t(chol(Sigma.pred-C%*%Sigma.obs.inv%*%t(C)))


C <- sigma2 * exp(- U.pred.coords / phi)
A <- C %*% Sigma.inv

mu.cond <- mu.pred + A %*% (t(S.sim) - as.numeric(object$mu))
sd.cond <- sqrt(sigma2 - rowSums(A * C))


S.pred.samples <- sapply(1:n.sim, function(i) cond.mean[,i]+
                           Sigma.cond%*%rnorm(nrow(st_coordinates(grid.pred))))

# prevalence samples on the grid
#prev.samples <- exp(beta_new+S.pred.samples)/(1+exp(beta.new+S.pred.samples))
prev.samples <- exp(S.pred.samples)/(1+exp(S.pred.samples))


r<- rasterFromXYZ(cbind(st_coordinates(grid.pred), apply(prev.samples,1,mean)))
plot(r)
plot(st_geometry(map),add=T)


r

#Problem: 
#Expected prevalence at far away places from survey locations
1/(1+exp(-beta)) #0.1821125 

#predicted prevalence at a far away place is much higher than 1/(1+exp(-beta))
points(2357592,-511841,col="red",pch=16)
extract(r,data.frame(x=2357592,y=-511841)) #0.5021257













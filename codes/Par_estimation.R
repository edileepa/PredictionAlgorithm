#no duplicated coords
table(duplicated(data[,3:4]))


#parameter estimation - GLM
glm.fit <- glm(cbind(y,n-y) ~ 1,data=data,family = binomial())
print(summary(glm.fit))

#assess spatial coorelation
library(geoR)

geodata<- as.geodata(obj=cbind(data[,3:4],glm.fit$residuals ), coords.col = 1:2, data.col = 3)
plot(geodata)

v <- variog(geodata)
plot(v)

v1<-variofit(v,ini.cov.pars=c(2,200000))
lines(v1)
v1

?likfit
v2 <- likfit(geodata,ini.cov.pars=c(2,200000))
v2
lines(v2,col="red")



library(PrevMap)
#round 1
mcmc.0 <- control.mcmc.MCML(n.sim=10000,burnin=2000,thin=8,h=(1.65)/(nrow(D)^(1/6)))
par.0 <- c(coef(glm.fit),c(3.5,200000,0.5))     #c(beta,sigma2,phi,tau)
start.0 <- c(200000,0.5/3.5 )     #starting values of phi 

cat("\npar0.0",par.0)
cat("\nstart.0",start.0,"\n\n")

fit.MCML.0<-binomial.logistic.MCML(y ~ 1,
                                 units.m=~n, 
                                 par0= par.0,
                                 coords= ~ long + lat,
                                 data=data,
                                 control.mcmc=mcmc.0,
                                 kappa=0.5,
                                 plot.correlogram = FALSE,
                                 #fixed.rel.nugget=10e-07,
                                 start.cov.pars=start.0)
print(summary(fit.MCML.0,log.cov.pars=FALSE))
saveRDS(fit.MCML.0,".\\outputs\\fit.MCML.0.rds")

mcmc.1 <- control.mcmc.MCML(n.sim=50000,burnin=5000,thin=45,h=(1.65)/(nrow(D)^(1/6)))
par.1 <- coef(fit.MCML.0)     #c(beta,sigma2,phi,tau)
start.1 <- c(coef(fit.MCML.0)[3], coef(fit.MCML.0)[4]/coef(fit.MCML.0)[2] )     #starting values of phi 

cat("\npar0.1",par.1)
cat("\nstart.1",start.1,"\n\n")

fit.MCML.1<-binomial.logistic.MCML(y ~ 1,
                                   units.m=~n, 
                                   par0= par.1,
                                   coords= ~ long + lat,
                                   data=data,
                                   control.mcmc=mcmc.1,
                                   kappa=0.5,
                                   plot.correlogram = FALSE,
                                   #fixed.rel.nugget=10e-07,
                                   start.cov.pars=start.1)
print(summary(fit.MCML.1,log.cov.pars=FALSE))
saveRDS(fit.MCML.1,".\\outputs\\fit.MCML.1.rds")

# Covariance parameters Matern function (kappa=0.5) 
#           Estimate StdErr
# sigma^2 7.9635e+00 0.1177
# phi     1.7982e+05 0.0000
# tau^2   7.0857e-01 1.4707




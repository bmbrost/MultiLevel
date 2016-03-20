rm(list=ls())

library(splines)
library(lmer)
library(nlme)

T <- 100  # number of observations
time <- cumsum(rgamma(T,shape=1.1,scale=4.5))  # time covariate
hr <- ((time)-24*floor(time/24))
day <- ceiling(time/24)

y <- sin(12-hr) + sin(day)
plot(time,y,type="l")


###
### Confirm MCMC algorithm estimates parameters well
###

X <- cbind(1,day,hr)
beta <- c(0.1,0.25,0.15)

int <- 10  # interval between knots
knots <- seq(0,max(time),by=int)
Z <- bs(time,knots=knots,degree=3,intercept=FALSE)  # cubic spline
sigma.alpha <- 1
alpha <- rnorm(ncol(Z),0,sigma.alpha)

sigma <- 1
y <- rnorm(T,X%*%beta+Z%*%alpha,sigma)
plot(time,y,type="l")

par(mfrow=c(3,1))
plot(X%*%beta,type="l")
plot(Z%*%alpha,type="l")
plot(time,y,type="l")

Id <- factor(rep(1,length(y)))
fit <- lme(y~day+hr,random=list(Id=pdIdent(~Z-1)))
y.hat.lme <- fit$fitted[,2]
lines(time,y.hat.lme,col=2)
beta.hat.lme <- as.matrix(fixef(fit))
alpha.hat.lme <- t(as.matrix(ranef(fit)))
plot(alpha,alpha.hat.lme);abline(a=0,b=1)

source('~/Documents/git/GLMM/normal.glmm.mcmc.R', chdir = TRUE)
start <- list(beta=beta,alpha=rep(0,ncol(Z)),sigma=sigma,sigma.alpha=sigma.alpha)
# hist(sqrt(1/rgamma(1000,1,,2)))
priors <- list(sigma.beta=2,r.sigma=2,q.sigma=1,r.sigma.alpha=2,q.sigma.alpha=1)
out1 <- normal.glmm.mcmc(y,X,Z,priors=priors,start=start,sigma.alpha=NULL,n.mcmc=1000)
out1$DIC

matplot(out1$beta,type="l",lty=1);abline(h=beta,col=1:3)
matplot(out1$alpha[,44],type="l",lty=1);abline(h=alpha[44],col=1:3)
matplot(out1$sigma,type="l");abline(h=sigma)
matplot(out1$sigma.alpha,type="l");abline(h=sigma.alpha)

alpha.hat <- apply(out1$alpha,2,mean)
plot(alpha.hat,alpha.hat.lme);abline(a=0,b=1)
plot(alpha.hat,alpha);abline(a=0,b=1)

alpha.quant <- t(apply(out1$alpha,2,quantile,c(0.025,0.975)))
sum(alpha>alpha.quant[,1]&alpha<alpha.quant[,2])

y.hat <- apply(out1$y.hat,1,mean)
plot(time,y,type="l")
lines(time,y.hat,col=3)
lines(time,y.hat.lme,col=2)

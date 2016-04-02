###
### Simulate grouped binomial distributed data and fit model using 
### binomial.varying.coef.mcmc.R
###

rm(list=ls())

library(mvtnorm)

logit <- function(x){
	log(x/(1-x))
}

expit <- function(x){
	exp(x)/(1+exp(x))
}


N <- 200  # number of trials per group
I <- 50  # number of events per group
J <- 10  # number of groups
N <- rep(N,each=I*J)

g <- rep(1:J,each=I)  # grouping variable

##########################################################
### Population-level process model parameters
##########################################################

mu.beta <- matrix(c(-0.5,1),,1)  # mean of betas
qX <- nrow(mu.beta)
rho <- -0.15  # correlation between betas
Sigma <- diag(qX)*0.75  # variance-covariance of betas
Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho

##########################################################
### Simulate group-level process model parameters
##########################################################

beta <- t(rmvnorm(J,mu.beta,Sigma))  # betas for each group
plot(t(beta))


##########################################################
### Simulate count data
##########################################################

X <- cbind(1,rnorm(I*J,0,sd=1))  # design matrix
beta.tmp <- t(beta[,g])
p <- expit(rowSums(X*beta.tmp))  # binomial probability
hist(p,breaks=100)
z <- rbinom(J*I,N,p)  # successes per event
hist(z)

##########################################################
### Fit model
##########################################################

source('~/Documents/git/Multilevel/nested/binomial/binomial.varying.coef.mcmc.R')
start <- list(beta=beta,mu.beta=mu.beta,Sigma=Sigma)
priors <- list(sigma.beta=5,S0=diag(qX),nu=qX+1)
tune <- list(beta=rep(1,J))
# tune <- list(beta=out1$tune$beta)
out1 <- binomial.varying.coef.mcmc(z,N,X,g,priors,start,tune,adapt=TRUE,1000)
out1$tune

# Examine estimates for beta_j
g.idx <- 10  # group idx for plotting beta_j
matplot(out1$beta[,,g.idx],type="l",lty=1);abline(h=beta[,g.idx],col=1:qX,lty=2)

# Examine estimates for mu.beta
matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2)

# Examine estimates for Lambda
matplot(cbind(out1$Sigma[1,1,],out1$Sigma[1,2,]),type="l")
abline(h=c(Sigma[1,1],Sigma[1,2]),lty=2,col=1:qX)





##########################################################
### Simulate data for Poisson GLM
##########################################################

H <- 10  # number of events
N <- 1000  # number of trials per event
X <- cbind(1,rnorm(H))  # design matrix
qX <- ncol(X)
beta <- c(-3,0.15)  # coefficients
p <- expit(X%*%beta)  # probability
z <- rbinom(H,N,p)  # observed successes
sum(z)

##########################################################
### Fit model
##########################################################

source('~/Documents/git/GLM/binomial/binomial.glm.mcmc.R')
priors <- list(sigma.beta=5)
tune <- list(beta=0.25)
start <- list(beta=coef(glm(cbind(z,N-z) ~ 0+X, family=binomial("logit"))))
out1 <- binomial.glm.mcmc(z,N=N,X,priors,start,tune,n.mcmc=50000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)
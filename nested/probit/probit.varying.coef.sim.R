###
### Simulate binary data and fit model using probit.glmm.mcmc.R
###

rm(list=ls())
library(mvtnorm)

I <- 100  # number of observations per group
J <- 10  # number of groups
g <- rep(1:J,each=I)  # grouping variable

mu.beta <- c(0,-4)  # mean of betas
k <- length(mu.beta)
rho <- -0.15  # correlation between betas
Lambda <- diag(k)*1  # variance-covariance of betas
Lambda[1,2] <- Lambda[2,1] <- Lambda[1,1]*Lambda[2,2]*rho

beta <- rmvnorm(J,mu.beta,Lambda)  # betas for each group
plot(beta)

# Define covariates
X <- cbind(1,rnorm(I*J))
qX <- ncol(X)

# Simulate data
beta.tmp <- beta[g,]
# p <- pnorm(X%*%beta.tmp)
p <- pnorm(rowSums(X*beta.tmp))
hist(p);summary(p)

y <- rbinom(I*J,1,p)
table(y)
plot(p,y)

# Fit model
source('~/Documents/git/Multilevel/nested/probit/probit.varying.coef.mcmc.R', chdir = TRUE)
start <- list(beta=beta,mu.beta=mu.beta,Lambda=Lambda)
priors <- list(sigma.beta=5,S0=diag(qX),nu=qX+1)
out1 <- probit.varying.coef.mcmc(y,X,g,priors,start,5000)

# Examine estimates for beta_j
g.idx <- 9  # group idx for plotting beta_j
matplot(out1$beta[,,g.idx],type="l",lty=1);abline(h=beta[g.idx,],col=1:qX,lty=2)

# Examine estimates for mu.beta
matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2)

# Examine estimates for Lambda
matplot(cbind(out1$Lambda[1,1,],out1$Lambda[1,2,]),type="l")
abline(h=c(Lambda[1,1],Lambda[1,2]),lty=2,col=1:qX)

# Examine estimates for v
boxplot(pnorm(out1$v),col=8,outline=FALSE)
points(y,col=3)


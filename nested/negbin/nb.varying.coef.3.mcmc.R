#
#
# Bayesian negative binomial generalized linear mixed model
#
# Function name: nb.varying.coef.3.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 11 April 2016
#
# Model statement:
#	y_ijk ~ NB(lambda_ijk,theta_i)
#	log(y_ijk) = x_ijk%*%alpha_ij
# 	alpha_ij ~ N(beta_j,Sigma_alpha)
# 	beta_i ~ N(mu_beta,Sigma_beta)
#	mu_beta ~ N(0,sigma.beta^2*I)
#	Sigma_alpha ~ Wish(S_0,nu)
#	Sigma_beta ~ Wish(S_0,nu)
#	theta ~ Gamma(a,b)
#	note: E[y_ijk]=lambda_ijk and Var[y_ijk]=lambda_ijk+lambda_ijk^2/theta
#	note: E[theta]=a/b and Var[theta]=a/(b^2)
#
# Reference:
#
# Required R packages: mvtnorm
#
# Inputs:
#
# y - vector of length n containing count observations corresponding to each row in
#	the design matrix X. Note that value y[1] corresponds to X[1,], y[2] corresponds
#	to X[2,], etc.
# X - design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# g1 - variable that defines groups of observations in y
# g2 - variable that defines subgroups (within groups) of observations in y
# priors - list of priors containing the following elements:
#	1. sigma.beta - standard deviation of normal prior on mu.beta
#	2. S0 - scale matrix for the inverse-Wishart prior on Sigma
#	3. nu - degrees of freedom for the IW prior on Sigma
#	4. a - shape parameter of gamma prior for alpha 
#	5. b - rate parameter of gamma prior for alpha
# start - list of starting values containing the following elements:
#	1. alpha - list of matrices containing starting values for subgroup-level coefficients
#	2. beta - matrix of starting values for group-level coefficients
#	3. mu.beta - matrix of starting values for population-level coefficients
#	4. Sigma_alpha - list of variance-covariance matrices for alphas
#	5. Sigma_beta - variance-covariance matrix for betas
#	6. theta - vector of over-dispersion parameter for negative binomial distribution
# tune - list of tuning parameters containing the following elements:
#	1. alpha - tuning parameter for Metropolis-Hastings update on alpha
#	2. theta - tuning parameter for Metropolis-Hastings update on theta
# adapt - switch to enable adapative tuning (TRUE/FALSE)
# Ta - interval at which to save subgroup-level samples (alpha,Sigma.alpha)
# n.mcmc - number of desired MCMC iterations
#
#

nb.varying.coef.3.mcmc <- function(y,X,g1,g2,priors,start,tune,adapt=TRUE,Ta=50,
	n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	library(mvtnorm)

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	###
	###  Setup variable
	###
# browser()
	g1 <- as.numeric(g1)
	g2 <- as.numeric(g2)

	N <- length(unique(g1))  # number of groups
	n.j <- tapply(g2,g1,function(x) length(unique(x)))
	g1.idx <- sapply(sort(unique(g1)),function(x) which(g1==x),simplify=FALSE)
		# indexes of observations in y by group
	n <- length(y)  # total number of observations
	qX <- ncol(X)
	
	###
	###  Starting values, and priors
	###

	beta <- matrix(start$beta,qX,N)
	alpha <- start$alpha
	theta <- start$theta
	
	lambda <- numeric(n)  # intensity of Poisson process
	for(i in 1:n){
		lambda[i] <- exp(t(X[i,])%*%alpha[[g1[i]]][,g2[i]])  
	}

	Sigma.alpha <- start$Sigma.alpha
	Sigma.alpha.inv <- lapply(Sigma.alpha,solve)

	mu.beta <- matrix(start$mu.beta,qX,1)
	Sigma.beta <- start$Sigma.beta
	Sigma.beta.inv <- solve(Sigma.beta)
	
	Sigma.mu.beta <- priors$sigma.mu.beta^2*diag(qX)
	Sigma.mu.beta.inv <- solve(Sigma.mu.beta)
	S0 <- priors$S0
	nu <- priors$nu

	###
	### Create receptacles for output
	###
  
	beta.save <- array(0,dim=c(n.mcmc,qX,N))
	alpha.save <- sapply(1:N,function(x) 
		array(0,dim=c(n.mcmc/Ta,qX,n.j[x])),simplify=FALSE)
	theta.save <- matrix(0,n.mcmc,N)
	mu.beta.save <- matrix(0,n.mcmc,qX)
	Sigma.beta.save <- array(0,dim=c(qX,qX,n.mcmc))
	Sigma.alpha.save <- array(0,dim=c(qX,qX,n.mcmc/Ta))
	Sigma.alpha.save <- lapply(1:N,function(x) Sigma.alpha.save)

	keep <- list(alpha=lapply(1:N,function(x) rep(0,n.j[x])),theta=rep(0,N))
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning

	###
	###  Begin MCMC loop
	###

	for(k in 1:n.mcmc){
		if(k%%1000==0) cat(k,"");flush.console()	

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			# browser()
			keep.tmp$alpha <- lapply(keep.tmp$alpha,function(x) x/Tb)
			tune$alpha <- sapply(1:N,function(x) 
				get.tune(tune$alpha[[x]],keep.tmp$alpha[[x]],k),simplify=FALSE)
			keep.tmp$alpha <- lapply(keep.tmp$alpha,function(x) x*0)			
			keep.tmp$theta <- keep.tmp$theta/Tb
			tune$theta <- get.tune(tune$theta,keep.tmp$theta,k)
			keep.tmp$theta <- keep.tmp$theta*0

	   	} 	

		for(i in 1:N){  # loop through groups
			for (j in 1:n.j[i]){  # loop through subgroups within groups

				idx <- which(g1==i&g2==j)	
				
				###
				### Update alpha_ij
				### 
				
				y.tmp <- y[idx]
				alpha.tmp <- alpha[[i]][,j]
				alpha.star <- c(rmvnorm(1,alpha.tmp,diag(qX)*(tune$alpha[[i]][j])^2))
				lambda.star <- exp(X[idx,]%*%alpha.star)  # intensity of Poisson process
		  		mh.0 <- sum(dnbinom(y.tmp,size=theta[i],mu=lambda[idx],log=TRUE))+
					sum(dmvnorm(alpha.tmp,beta[,i],Sigma.alpha[[i]],log=TRUE))
				mh.star <- sum(dnbinom(y.tmp,size=theta[i],mu=lambda.star,log=TRUE))+
					sum(dmvnorm(alpha.star,beta[,i],Sigma.alpha[[i]],log=TRUE))
				if(exp(mh.star-mh.0)>runif(1)){
					alpha[[i]][,j] <- alpha.star
					lambda[idx] <- lambda.star
					keep$alpha[[i]][j] <- keep$alpha[[i]][j]+1
					keep.tmp$alpha[[i]][j] <- keep.tmp$alpha[[i]][j]+1
				}
			}  # end loop through subgroups

			###
			### Update theta (over-dispersion parameter)
			### 
# browser()
			idx <- g1.idx[[i]]
			theta.star <- rnorm(1,theta[i],tune$theta[i])
			if(theta.star>0){
				lambda.tmp <- lambda[idx]
				y.tmp <- y[idx]
		  		mh.0.theta <- sum(dnbinom(y.tmp,size=theta[i],mu=lambda.tmp,log=TRUE))+
					dgamma(theta[i],shape=priors$a,rate=priors$b,log=TRUE)
				mh.star.theta <- sum(dnbinom(y.tmp,size=theta.star,mu=lambda.tmp,log=T))+
					dgamma(theta.star,shape=priors$a,rate=priors$b,log=TRUE)
				if(exp(mh.star.theta-mh.0.theta)>runif(1)){
					theta[i] <- theta.star
					keep$theta[i] <- keep$theta[i]+1
					keep.tmp$theta[i] <- keep.tmp$theta[i]+1
				}		
			}

			###
			### Update Sigma.alpha_i
			###	
	
			# browser()		
		  	Sn <- S0+crossprod(t(alpha[[i]])-matrix(beta[,i],n.j[i],qX,byrow=TRUE))
			Sigma.alpha[[i]] <- solve(rWishart(1,nu+n.j[i],solve(Sn))[,,1])
			Sigma.alpha.inv[[i]] <- solve(Sigma.alpha[[i]])
		
			###
			### Update beta_i
			### 

			alpha.sum <- apply(alpha[[i]],1,sum)
			A.inv <- solve(n.j[i]*Sigma.alpha.inv[[i]]+Sigma.beta.inv)
			b <- Sigma.alpha.inv[[i]]%*%alpha.sum+Sigma.beta.inv%*%mu.beta		    
		    beta[,i] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

			###
			### Save subgroup-level samples
			###
			
			if(k%%Ta==0){
				k.tmp <- k/Ta
				alpha.save[[i]][k.tmp,,] <- alpha[[i]]
				Sigma.alpha.save[[i]][,,k.tmp] <- Sigma.alpha[[i]]
			} 
		}  # end loop through groups
					
		###
		### Update Sigma.beta
		###	

		# browser()		
	  	Sn <- S0+crossprod(t(beta)-matrix(mu.beta,N,qX,byrow=TRUE))
		Sigma.beta <- solve(rWishart(1,nu+N,solve(Sn))[,,1])
		Sigma.beta.inv <- solve(Sigma.beta)
				
	  	###
	  	### Sample mu_beta
	  	###

		beta.mean <- apply(beta,1,sum)
		A.inv <- solve(N*Sigma.beta.inv+Sigma.mu.beta.inv)
		b <- Sigma.beta.inv%*%beta.mean
	    mu.beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

		###
		###  Save population and group-level samples 
	    ###

	  	beta.save[k,,] <- beta
		mu.beta.save[k,] <- mu.beta
		Sigma.beta.save[,,k] <- Sigma.beta
		theta.save[k,] <- theta
	}

	keep$alpha <- lapply(keep$alpha,function(x) round(x/n.mcmc,2))
	keep$theta <- round(keep$theta/n.mcmc,2)
	print("alpha acceptance rate:") 
	cat("\n")
	print(keep$alpha)
	cat("\ntheta acceptance rate:",keep$theta) 
	
	###
	### Write output
	###

	# Create list of ending values
# browser()
	end <- list(beta=beta.save[k,,],mu.beta=matrix(mu.beta.save[k,],qX,1),
		theta=theta.save[k,],Sigma.beta=Sigma.beta.save[,,k],
		alpha=lapply(alpha.save,function(x) x[k/Ta,,]),
		Sigma.alpha=lapply(Sigma.alpha.save,function(x) x[,,k/Ta]))

	
	list(alpha=alpha.save,beta=beta.save,mu.beta=mu.beta.save,theta=theta.save,
		Sigma.beta=Sigma.beta.save,Sigma.alpha=Sigma.alpha.save,keep=keep,end=end,
		y=y,X=X,g1=g1,g2=g2,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}
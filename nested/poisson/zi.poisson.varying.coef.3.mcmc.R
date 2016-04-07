#
#
# Bayesian binomial generalized linear mixed model
#
# Function name: binomial.varying.coef.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 02 April 2016
#
# Model statement:
#	z_ij ~ Binomial(N_ij,p_ij)
#	logit(p_ij) = x_ij%*%beta_j
# 	beta_j ~ N(mu_beta,Sigma)
#	mu_beta ~ N(0,sigma.beta^2*I)
#	Sigma ~ Wish(S_0,nu)
#
# Reference:
#
# Required R packages: mvtnorm
#
# Inputs:
#
# z - vector of length n containing the number of successes during event i in group j.
#	Order of elements in z match order of rows in design matrix X, i.e., 
# 	z[1] corresponds to X[1,], z[2] corresponds to X[2,], etc.
# N - vector of length n containing the number of trials during event i in group j.
# X - design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# g - variable that defines groups of observations in z
# priors - list of priors containing the following elements:
#	1. sigma.beta - Standard deviation of normal prior on mu.beta
#	2. S0 - Scale matrix for the inverse-Wishart prior on Sigma
#	3. nu - Degrees of freedom for the IW prior on Sigma
# start - list of starting values containing the following elements:
#	1. beta - Vector of starting values for coefficients
#	2. mu.beta - Vector of starting values for mean of betas
#	3. Sigma - Variance-covariance matrix for betas
# tune - list of tuning parameters containing the following elements:
#	1. beta - Tuning parameter for Metropolis-Hasting update on beta
# adapt - Switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - Number of desired MCMC iterations
#
#

zi.poisson.varying.coef.3.mcmc <- function(y,X,g1,g2,priors,start,tune,adapt=TRUE,n.mcmc=1000){

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

	I <- length(unique(g1))  # number of groups
	J <- tapply(g2,g1,function(x) length(unique(x)))
	# g1.idx <- sapply(sort(unique(g1)),function(x) which(g1==x),simplify=FALSE)
		# indexes of observations in y by group
	n <- length(y)  # total number of observations
	qX <- ncol(X)
	y0 <- which(y==0)  # zero valued observations by group
	z <- rep(1,n)  # latent mixture component indicator
	
	###
	###  Starting values, and priors
	###

	beta <- matrix(start$beta,qX,I)
	alpha <- start$alpha
	p <- start$p
	
	lambda <- numeric(n)  # intensity of Poisson process
	for(i in 1:n){
		lambda[i] <- exp(t(X[i,])%*%alpha[[g1[i]]][,g2[i]])  
	}

# browser()

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
  
	beta.save <- array(0,dim=c(n.mcmc,qX,I))
	alpha.save <- sapply(1:I,function(x) array(0,dim=c(n.mcmc,qX,J[x])),simplify=FALSE)
	p.save <- sapply(1:I,function(x) matrix(0,n.mcmc,J[x]),simplify=FALSE)
	z.save <- matrix(0,n.mcmc,n)
	mu.beta.save <- matrix(0,n.mcmc,qX)
	Sigma.beta.save <- array(0,dim=c(qX,qX,n.mcmc))
	Sigma.alpha.save <- lapply(1:I,function(x) Sigma.beta.save)

	keep <- list(alpha=lapply(1:I,function(x) rep(0,J[x])))
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
			tune$alpha <- sapply(1:I,function(x) 
				get.tune(tune$alpha[[x]],keep.tmp$alpha[[x]],k),simplify=FALSE)
			keep.tmp$alpha <- lapply(keep.tmp$alpha,function(x) x*0)			
	   	} 	

		for(i in 1:I){  # loop through groups

			for (j in 1:J[i]){  # loop through subgroups within groups

				idx <- which(g1==i & g2==j)
					
				###
				### Update z_ijk
				###

				y0 <- idx[which(y[idx]==0)]
				n.y0 <- length(y0)
				p.tmp <- p[[i]][j]*exp(-lambda[y0])
				p.tmp <- p.tmp/(p.tmp+1-p[[i]][j])  
				z[y0] <- rbinom(n.y0,1,p.tmp)
		
				###
				### Update p_ij
				###
		
				n.k <- length(idx)
				z.tmp <- z[idx]	
				sum.z <- sum(z.tmp)
				p[[i]][j] <- rbeta(1,sum.z+priors$a,n.k-sum.z+priors$b)

				###
				### Update alpha_ij
				### 
				
				z1 <- which(z.tmp==1)
				y.tmp <- y[idx[z1]]
				alpha.tmp <- alpha[[i]][,j]
				alpha.star <- c(rmvnorm(1,alpha.tmp,diag(qX)*tune$alpha[[i]][j]))
				lambda.star <- exp(X[idx,]%*%alpha.star)  # intensity of Poisson process
		  		mh.0 <- sum(dpois(y.tmp,lambda[idx[z1]],log=TRUE))+
					sum(dmvnorm(alpha.tmp,beta[,i],Sigma.alpha[[i]],log=TRUE))
				mh.star <- sum(dpois(y.tmp,lambda.star[z1],log=TRUE))+
					sum(dmvnorm(alpha.star,beta[,i],Sigma.alpha[[i]],log=TRUE))
				if(exp(mh.star-mh.0)>runif(1)){
					alpha[[i]][,j] <- alpha.star
					lambda[idx] <- lambda.star
					keep$alpha[[i]][j] <- keep$alpha[[i]][j]+1
					keep.tmp$alpha[[i]][j] <- keep.tmp$alpha[[i]][j]+1
				}
			}  # end loop through subgroups
# browser()
			###
			### Update Sigma.alpha_i
			###	
	
			# browser()		
		  	Sn <- S0+crossprod(t(alpha[[i]])-matrix(beta[,i],J[i],qX,byrow=TRUE))
			Sigma.alpha[[i]] <- solve(rWishart(1,nu+J[i],solve(Sn))[,,1])
			Sigma.alpha.inv[[i]] <- solve(Sigma.alpha[[i]])
		
			###
			### Update beta_i
			### 

			alpha.sum <- apply(alpha[[i]],1,sum)
			A.inv <- solve(J[i]*Sigma.alpha.inv[[i]]+Sigma.beta.inv)
			b <- Sigma.alpha.inv[[i]]%*%alpha.sum+Sigma.beta.inv%*%mu.beta		    
		    beta[,i] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

			###
			### Save subgroup-level samples
			###
			
			alpha.save[[i]][k,,] <- alpha[[i]]
			Sigma.alpha.save[[i]][,,k] <- Sigma.alpha[[i]]
			p.save[[i]][k,] <- p[[i]]
		}  # end loop through groups
					
		###
		### Update Sigma.beta
		###	

		# browser()		
	  	Sn <- S0+crossprod(t(beta)-matrix(mu.beta,I,qX,byrow=TRUE))
		Sigma.beta <- solve(rWishart(1,nu+I,solve(Sn))[,,1])
		Sigma.beta.inv <- solve(Sigma.beta)
				
	  	###
	  	### Sample mu_beta
	  	###

		beta.mean <- apply(beta,1,sum)
		A.inv <- solve(I*Sigma.beta.inv+Sigma.mu.beta.inv)
		b <- Sigma.beta.inv%*%beta.mean
	    mu.beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

		###
		###  Save population and group-level samples 
	    ###

	  	beta.save[k,,] <- beta
		z.save[k,] <- z
		mu.beta.save[k,] <- mu.beta
		Sigma.beta.save[,,k] <- Sigma.beta
	}

	keep$alpha <- lapply(keep$alpha,function(x) round(x/n.mcmc,2))
	print("alpha acceptance rate:") 
	cat("\n")
	print(keep$alpha)
	
	###
	### Write output
	###

	list(alpha=alpha.save,beta=beta.save,mu.beta=mu.beta.save,p=p.save,z=z.save,
		Sigma.beta=Sigma.beta.save,Sigma.alpha=Sigma.alpha.save,keep=keep,
		y=y,X=X,g1=g1,g2=g2,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}
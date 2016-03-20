probit.varying.coef.mcmc <- function(y,X,g,priors,start,n.mcmc){

	###
	### Brian M. Brost (08 MAR 2016)
	### Mixed effects regression for binary data using probit link
	###

	###
	### Model statement:
	### y_ij=0,u_ij<=0
	### y_ij=1,u_ij>0
	### v_ij ~ N(x_j*beta_j,1)
	### beta_j ~ N(mu.beta,Lambda)
	###	mu.beta ~ N(0,sigma.beta^2*I)
	###	Lambda ~ Wish(S0,nu)
	###	
	
	###
	###  Libraries and Subroutines
	###
	
	truncnormsamp <- function(mu,sig2,low,high,nsamp){
	  flow=pnorm(low,mu,sqrt(sig2)) 
	  fhigh=pnorm(high,mu,sqrt(sig2)) 
	  u=runif(nsamp) 
	  tmp=flow+u*(fhigh-flow)
	  x=qnorm(tmp,mu,sqrt(sig2))
	  x
	}
	
	###
	###  Preliminary Variables
	###

# browser()
	J <- length(unique(g))  # number of groups
	g <- as.numeric(g)
	g.idx <- sapply(sort(unique(g)),function(x) which(g==x),simplify=FALSE)
		# indexes of observations in y by group
	n.j <- unlist(lapply(g.idx,length))  # number of observations per group
	n <- length(y)  # total number of observations
	qX <- ncol(X)
	y1 <- (y==1)
	y0 <- (y==0)
	y1.sum <- sum(y1)
	y0.sum <- sum(y0)
	v <- numeric(n)

	###
	### Starting values and priors
 	###
	
	beta <- matrix(start$beta,J,qX)
	mu.beta <- start$mu.beta
	Lambda <- start$Lambda
	Lambda.inv <- solve(Lambda)

	Sigma.beta <- priors$sigma.beta^2*diag(qX)
	Sigma.beta.inv <- solve(Sigma.beta)
	S0 <- priors$S0
	nu <- priors$nu

	###
	### Create receptacles for output
	###
  
	beta.save <- array(0,dim=c(n.mcmc,qX,J))
	mu.beta.save <- matrix(0,n.mcmc,qX)
	Lambda.save <- array(0,dim=c(qX,qX,n.mcmc))
	v.save <- matrix(0,n.mcmc,n)

	###
	###  Gibbs Loop
	###
	
	for(k in 1:n.mcmc){
		if(k%%1000==0) cat(k," ");flush.console()
# browser()	
		###
		### Sample v (auxilliary variable for probit regression)
	  	###

		linpred <- rowSums(X*beta[g,])
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	
	  	###
	  	### Sample beta_j
	  	###

		# Using g
		# A.inv <- lapply(1:J,function(x) solve(t(X[g==x,])%*%X[g==x,]+Lambda.inv))
		# b <- t(sapply(1:J,function(x) t(v[g==x])%*%X[g==x,]))+
			# matrix(t(mu.beta)%*%Lambda.inv,J,2,byrow=TRUE)
	  	# beta <- t(sapply(1:J,function(x)
	  		# A.inv[[x]]%*%b[x,]+t(chol(A.inv[[x]]))%*%matrix(rnorm(qX),qX,1)))
		
		# Using g.idx
		A.inv <- lapply(g.idx,function(x) solve(t(X[x,])%*%X[x,]+Lambda.inv))		
		b <- lapply(g.idx,function(x) t(v[x])%*%X[x,]+t(mu.beta)%*%Lambda.inv)
		beta <- t(sapply(1:J,function(x)
	  		A.inv[[x]]%*%t(b[[x]])+t(chol(A.inv[[x]]))%*%matrix(rnorm(qX),qX,1)))

	  	###
	  	### Sample mu_beta
	  	###

		beta.mean <- apply(beta,2,sum)
		A.inv <- solve(J*Lambda.inv+Sigma.beta.inv)
		b <- Lambda.inv%*%beta.mean
	    mu.beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
	  	
	  	###
	  	### Sample Lambda
	  	###
	  
	  	Sn <- S0+crossprod(beta-matrix(mu.beta,J,qX,byrow=TRUE))
		Lambda <- solve(rWishart(1,nu+J,solve(Sn))[,,1])
		Lambda.inv <- solve(Lambda)
	
	  	###
	  	### Save Samples 
	  	###

	  	beta.save[k,,] <- t(beta)
	  	v.save[k,] <- v
		mu.beta.save[k,] <- mu.beta
		Lambda.save[,,k] <- Lambda
		
	}
	
	###
	###  Write output 
	###
	
	list(beta=beta.save,v=v.save,mu.beta=mu.beta.save,Lambda=Lambda.save,
		y=y,X=X,g=g,start=start,priors=priors,n.mcmc=n.mcmc)
}

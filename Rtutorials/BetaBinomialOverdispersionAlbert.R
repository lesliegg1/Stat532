
###The Beta-Binomial is a Beta-mixture of Binomials.  
###  It can be helpful when individual study units vary greatly in the probability 
###  of a success (like the cities here).  The units are allowed to have unequal 
###   probabilities that follow a Beta dist'n.   The distribution of eta is 
###   describing a distribution from which the observed proportions could
###    have been drawn.  This is actually an example of a hierarchical model.
### See Beta-Binomial example in Chapter 5 in BDA3.  We will apply some of the 
### Chapter 5 results to this example.  See BetaBinomialPriorInvestigation.R for
### more info on the priors.  The prior used here for K is not immediately 
### intuitive, but provides a good setting to talk more about priors on variance
### parameters in hierarchical setting and the issues that can arise.
### NOTE to me:  We used this as example for multivariate rejection sampling
### and then transitioned into hieararchical models with it (given ties to Ch. 5 BDA3)

#We can assume the counts come from a Beta-Binomial model with mean eta 
# and precision K  (See Appendix page 584 in Gelman et al 
##    for parameterization)
##     n=nj,  beta=K(1-eta), alpha=K*eta
##    K=(alpha+beta) ##prior total number of "trials"
##    eta = alpha/(alpha + beta) ##proportion success out of prior "trials"
  
 ##This function evaluates the likelihood for eta and K 
  BB.loglik.fun <- function(etaK.vec, n.vec, y.vec) {
     ll.out <- sum(lchoose(n.vec, y.vec) + 
         lbeta((etaK.vec[2]*etaK.vec[1] + y.vec),(etaK.vec[2]*(1-etaK.vec[1])+n.vec-y.vec))- 
         lbeta(etaK.vec[2]*etaK.vec[1], etaK.vec[2]*(1-etaK.vec[1])))
  	 return(ll.out)
  	}
  	
  
  #### Applied problem  (Tsutakawa et al. 1985) - Estimate rate of death from
  ## stomach cancer for at risk males between ages 45-64 for the largest cities in 
  ## Missouri.  Here are the mortality rates for 20 of the cities
  ##  (number at risk (nj), number of cancer deaths (yj))
  
  nj.vec <- c(1083,855,3461,657,1208,1025, 527, 1668, 583, 582, 917, 857, 680, 917, 53637,
               874, 395, 581, 588, 383)
  yj.vec <- c(0,0,2,0,1,1,0,2,1,3,0,1,1,1,54,0,0,1,3,0)
  
   
  ## We need to assign priors for eta and K - let's use the vague prior 
  ##  porportional to  g(eta,K) propto (1/(eta*(1-eta)))*(1/(1+K)^2)
  
   BB.prior.fun <- function(etaK.vec) {
   	 (1/(etaK.vec[1]*(1-etaK.vec[1])))*(1/((1+etaK.vec[2])^2))
   	}
   	
   ## We could combine them into one function to find the posterior
   
   BB.logpost1<- function(etaK.vec, n.vec, y.vec) {
  	 ll <- sum(lbeta(etaK.vec[2]*etaK.vec[1] + y.vec, etaK.vec[2]*(1-etaK.vec[1])+n.vec-y.vec)- 
                 lbeta(etaK.vec[2]*etaK.vec[1], etaK.vec[2]*(1-etaK.vec[1])))
  	 lprior <- -log(etaK.vec[1]) - log(1-etaK.vec[1]) - 2*log(1+etaK.vec[2])              
  	 lpost.out <- ll + lprior
  	 return(lpost.out)
  	}

    eta.vec <- seq(0.0001,0.003,length=200)
   	K.vec <- seq(1,15000, length=200)
   	etaK.grid <- expand.grid(eta.vec, K.vec)
   	lp.grid <- apply(etaK.grid, 1, BB.logpost1, n.vec=nj.vec, y.vec=yj.vec)
   	lp.mat <- matrix(lp.grid, nrow=200, ncol=200)
   	p.mat <- exp(lp.mat - max(lp.mat)) 
   	contour(eta.vec, K.vec, p.mat, ylab="K", xlab=expression(eta))
   
   #### Now, let's transform the parameters to the real line as they do
   ### for the Beta-Binomial hierarchical example in BDA3 Chapter 5 
   ## theta1 <- logit(eta)
   ## theta2 <- log(K)
        	
   	BB.logpost <- function(theta.vec, n.vec, y.vec) {
   	  eta <- exp(theta.vec[1])/(1+exp(theta.vec[1]))
   	  K <- exp(theta.vec[2])
   	  ll <- sum(lbeta(K*eta + y.vec, K*(1-eta) + n.vec - y.vec) - 
                lbeta(K*eta, K*(1-eta)))
   	  #log.prior <- -log(eta) - log(1-eta) - 2*log(1+K)  
   	  #log.Jacob <- (theta.vec[1]+theta.vec[2]) - 2*log(1+exp(theta.vec[1]))
   	  log.rest <- theta.vec[2] - 2*log(1+exp(theta.vec[2])) #combined prior and Jacob
   	  #trans.log.post <- ll + log.prior + log.Jacob
   	  trans.log.post <- ll + log.rest
   	  return(trans.log.post)
   	}
   	
   	   	
   	theta1.vec <- seq(-8,-5,length=200)
   	theta2.vec <- seq(5,17, length=200)
   	theta.grid <- expand.grid(theta1.vec, theta2.vec)
   	logpost.grid <- apply(theta.grid, 1, BB.logpost, n.vec=nj.vec, y.vec=yj.vec)
   	logpost.mat <- matrix(logpost.grid, nrow=200, ncol=200)
   	post.mat <- exp(logpost.mat - max(logpost.mat)) 
   	contour(theta1.vec, theta2.vec, post.mat, ylab="log(K)", 
              xlab=expression(logit(eta)))
   	
   	
    contour(theta1.vec, theta2.vec, logpost.mat, ylab="log(K)", 
      xlab=expression(logit(eta)), levels=seq(-600,-500,by=5))  #not as nice
   	

   	### From contour plot, see that mode is near (-7,6) -can use this as a starting
   	## point to find mode using the Nelder-Mead Method   	
   	optim.out <- optim(c(-7,6), BB.logpost, n.vec=nj.vec, y.vec=yj.vec, 
                control=list(fnscale=-100), method="Nelder-Mead", hessian=TRUE)
   	
   	optim.out
   	optim.out$par   	
   	Var <- solve(-optim.out$hessian)
   	Var
   	
   	## Using this information we could use a multi-variate normal approximation for 
    ##  the posterior distribution,  though we would rather sample from it directly.  
    ##  We can also use this information to come
   	##  up with a proposal distribution for a rejection sampling algorithm.  If 
   	##  we use the normal it will probably be too skinny in the tails and may 
   	##  not bound our posterior, so let's try a t-distribution on very few d.f. 
   	##  to get fatter tails.
    ##  There is a nice multivariate t function in the LearnBayes package
    ##   install.packages("LearnBayes") 
       library(LearnBayes)
       help(dmt)
   	
   	#### NOW let's use REJECTION SAMPLING to obtain samples from the joint
   	#### posterior distribution of logit(eta) and log(K)
   	
   	### STEP 1:  Find a proposal density of a simple form that, when multiplied
   	# by an appropriate constant covers the posterior density of interest 
   
  	#We essentially want to maximize log(g(theta|y)) - log(p(theta))
  	#  Let's use optim() again to do this, by first writing a function to maximize
  	
  	post.prop.diff <- function(theta, n.vec, y.vec, t.params) {
  		post.part <- BB.logpost(theta, n.vec=n.vec, y.vec=y.vec)
  		proposal.part <- dmt(theta, mean=t.params$m, S=t.params$var, 
                             df=t.params$df, log=TRUE)
  		d <- post.part - proposal.part
  		return(d)
  		} 
	
	t.params.set <- list(m=c(-6.8, 7.6), var=2*Var, df=4)
	
	d.out <- optim(c(-7, 7), post.prop.diff, n.vec=nj.vec, y.vec=yj.vec, 
                   t.params=t.params.set, control=list(fnscale=-10))
	
	d.out$par   #-6.8899, 12.46 So, max value of d occurs at (-6.88, 12.46)
	            # -6.886769  7.507773  (for some starting values)
	            
	post.prop.diff(c(-5, 7), n.vec=nj.vec, y.vec=yj.vec, t.params=t.params.set)
	d.max <- post.prop.diff(c(-6.8899, 12.46), n.vec=nj.vec, y.vec=yj.vec, 
             t.params=t.params.set) #-569.3335
	dmax <- post.prop.diff(c(-6.8868, 7.5077),n.vec=nj.vec, y.vec=yj.vec, 
            t.params=t.params.set) #-570.07
	
  ##### NOW, finally we obtain a sample using rejection sampling...we will do 
  ###   it in vector form, rather than a loop
	
	#1. Take a lot of draws from the proposal distribution (a multivariate-t)
	  n.draws <- 10000  ##will not end up with 10000 draws b/c using rejection
	  prop.draws <- rmt(n.draws, mean=t.params.set$m, S=t.params.set$var, 
                        df=t.params.set$df)

      #2. Evaluate the posterior and the proposal distribution at all the values
	  log.post <- apply(prop.draws, 1, BB.logpost, n.vec=nj.vec, y.vec=yj.vec)
	  log.g.theta <- dmt(prop.draws, mean=t.params.set$m, S=t.params.set$var, 
                         df=t.params.set$df, log=TRUE)

      #3. Calculate the ratio
	  accept.probs <- exp(log.post - log.g.theta - d.max)

      #4. Pick off the draws that we accept
	  theta.draws <- prop.draws[runif(n.draws) <= accept.probs,]
	  
	  #What percent did we accept?
	   length(theta.draws[,1])/n.draws #only 28% -- inefficient?, but easy to find?
	  
	    contour(theta1.vec, theta2.vec, post.mat, ylab="log(K)", 
                   xlab=expression(logit(eta)))
	    points(theta.draws[,1], theta.draws[,2], pch=1, col=4, cex=0.75)

        ##Inference about the mean, eta?
           eta.draws <- exp(theta.draws[,1])/(1+exp(theta.draws[,1]))
           hist(eta.draws, nclass=100, col="lightgray", freq=FALSE, 
                       xlab=expression(eta))
           median(eta.draws) #0.001081160
           abline(v=yj.vec/nj.vec, col=3) # Add the observed proportions
           quantile(eta.draws, c(.025, .975)) #0.0006410274 0.0020984802 
                            #may need more draws to better capture the right-hand
                            #  tail of the distn

  #### NEW MATERIAL to bring in BDA3 Chapter 5 stuff.
  #####So, using rejection sampling we obtained draws from p(logit(eta), log(K)|y).
  #### This follows the steps on the bottom of Page 112 to sample from the
  #### joint posterior distribution of the parameters (pi's) AND the hyperparameters
  ### For step (1), they describe using grid approximation and we just used rejection
  ### sampling instead.
  
  #Step 2.(a) Transform the lth draw of (logit(eta), log(K)) to (alpha,beta) to
  # get a draw from their marginal posterior distribution
   theta.to.ab.fun <- function(theta.vec) {
   	   eta <- exp(theta.vec[1])/(1+exp(theta.vec[1]))
   	   K <- exp(theta.vec[2])
   	   a <- eta*K
   	   b <- (1-eta)*K
   	   c(a,b)
    }
   alphabeta_draws <- t(apply(theta.draws, 1, theta.to.ab.fun))
   head(alphabeta_draws)
   
  #Step 2 (b) For each j=1,....,J, sample the pi's from their conditional posterior
  # distribution  pi.j|alpha,beta,y ~ Beta(alpha+y.j, beta+n.j-y.j)
  
    get.pis.fun <- function(ab, yvec, mvec) {
    	n <- length(yvec)
    	rbeta(n, ab[1]+yvec, ab[2]+mvec-yvec)
    }
  
   pis.draws <- t(apply(alphabeta_draws, 1, get.pis.fun, yvec=yj.vec, mvec=nj.vec))
   
   ###So, could use these draws to get the pi's you wanted to do the posterior
   ### predictive checks for the Problem 1 on Assignment 8.
   n.draws <- dim(pis.draws)[1]
   
   par(mfrow=c(5,2))
    for (j in 1:10){
    	hist(pis.draws[,j], nclass=100, col=gray(0.8), main="", freq=FALSE, xlim=c(0,0.008))
    	abline(v=yj.vec[j]/nj.vec[j], lwd=2, col="blue") #observed proportion   
     }
    
   par(mfrow=c(5,2))
    for (j in 11:20){
    	hist(pis.draws[,j], nclass=100, col=gray(0.8), main="", freq=FALSE, xlim=c(0,0.008))
    	abline(v=yj.vec[j]/nj.vec[j], lwd=2, col="blue") #observed proportion   
     }


     
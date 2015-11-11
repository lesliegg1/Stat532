
### Albert put prior directly on eta and K
### Beta(0,0) and (1/((K+1)^2)) 
### p(eta,K) propto (eta*(1-eta))^{(-1)}*(1/((K+1)^2))

### First, good to realize that K=(alpha + beta) is 
###  a precision parameter (though not equal to THE precision
###  of a beta distributed random variable).  As alpha and beta
###  increase, the variance decreases.
###  Exercise 5.9 in BDA 3 says that the variance is approximately
###  1/(alpha+beta), which is consistent with prec approx. alpha+beta
###  This approximation makes sense if thinking about alpha/beta=1,
###  and easiest to show via algebra if start with precision instead
###  of variance (see my notes from Fall 15) 

## Question:  WHY did they choose (1/((K+1)^2))
##   1. Makes sense to put a 1 in denominator to cover improper case where
##       beta=0 and alpha=0  (like including 0 prior trials)
##   2. Compare 1/(1+K) to 1/(1+K)^2 curves below.  Squared version has a 
##       little more mass on larger values of precision (smaller variances)

##Let's just explore some shapes for precisions

p1_K <- function(K) {1/K}  # Equivalent to log(K) propto 1
 #To keep from blowing up at zero, we will add 1 to K
 # to instead approximate variance with 1/(alpha+beta+1) = 1/(K+1)
p2_K <- function(K) {1/(K+1)} #From transform p(var) \propto 1/var
p3_K <- function(K) {1/((K+1)^2)} #From Albert and from transform p(var) \propto 1
p4_K <- function(K) {1/(2*(K+1)^(3/2))} #From transform p(SD) \propto 1

curve(p1_K(x), from=0.001, to=10, n=2000, lwd=2, col=2, ylim=c(0,8))
 curve(p2_K(x), col=3, lwd=2, add=TRUE, n=2000)
 curve(p3_K(x), col=4, lwd=3, add=TRUE, n=2000)
 curve(p4_K(x), col=5, lwd=2, add=TRUE, n=2000)
 legend(5,6, legend=c("1/k", "1/(k+1)", "1/((k+1)^2)", "1/(2*(K+1)^(3/2))"), col=c(2,3,4,5), 
         lwd=c(2,2,2,2), bty="n")
   #Note behavior of the functions near 0
   #Do not want prior to blow up at zero -- leads to infinite variance


#### Now, in BDA3 Chapter 5, they have a discussion of what priors on K and eta
#### lead to improper posterior distributions.  Let's look at the posterior
#### contours for some of these.  We will go ahead and use the Rat Tumor Data
#### to look at posterior contours under different prior distributions
#### I essentially implement the suggestions given in Section 5.3 for the Rat
#### tumor study 

#Load data set
#setwd("/Users/meganhiggs/Documents/STAT 532 - Bayesian Fall2015/Computation")
setwd("~/Documents/Stat532/Rtutorials/")
rattumor.data <- read.table("rats.asc", header=TRUE)
rattumor.data

#### 1. Implement using the prior obtained from assuming a uniform prior on the
####    SD = (1/sqrt(K)) = (1/sqrt(alpha+beta))
####    Working through the transformation gives 
####      p(alpha,beta) propto (alpha + beta)^(-5/2)
####    or p(logit(eta), log(K)) propto alpha*beta*(alpha+beta)^(-5/2)
####    See Exercise 5.9 in BDA3 for more details 
####    This prior should lead to PROPER Posterior
    
	# Write a function to get the log of the unnormalized posterior, so
	# we can look at contour lines (make sure not going to infinity) 
		l.unpost.BDA3 <- function(phi,yvec,mvec) {
    		n <- length(yvec)
    		b <- exp(phi[2])/(1+exp(phi[1]))
    		a <- b*exp(phi[1])
    		llik <- n*(lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    		         sum(lgamma(a + yvec) + lgamma(b + mvec - yvec) - 
    		             lgamma(a + b + mvec))
    		lprior <- -5/2*log(a+b) 
    		lJacob <- log((1+ b^2)) - log(a*b*(a+b)) 
    		out <- llik + lprior + lJacob
    		return(out)
		}


	#Calculate l.unpost over a grid
 		phi.1 <- seq(-2.3,-1.3,length=100) #logit(eta) = log(alpha/beta) 
 		phi.2 <- seq(1,5,length=100) #log(K) = log(alpha+beta)  
 		phi.grid <- expand.grid(phi.1,phi.2)
 		l.unpost <- apply(phi.grid, 1, l.unpost.BDA3, yvec=rattumor.data[,1], 
 						mvec=rattumor.data[,2])
 		unpost <- exp(l.unpost - max(l.unpost))
 		unpost.mat <- matrix(unpost, nrow=length(phi.1), ncol=length(phi.2))

##Plot contours of unnormalized posterior on the transformed scale
 		contour(phi.1, phi.2, unpost.mat, levels=seq(0.05, 0.95, 0.1),
           ylab=expression(log(alpha + beta)),xlab=expression(log(alpha/beta)))
    

#### 2. Implement using the prior obtained using the 1/(K+1)^2 from Albert
#### and the stomach cancer example we did
####    p(eta, K) = (eta*(1-eta))^(-1)*(1/(K+1)^2)
####    This prior should lead to PROPER Posterior
    
	# Write a function to get the log of the unnormalized posterior, so
	# we can look at contour lines (make sure not going to infinity) 
		l.unpost.Albert <- function(phi,yvec,mvec) {
    		n <- length(yvec)
    		b <- exp(phi[2])/(1+exp(phi[1]))
    		a <- b*exp(phi[1])
    		llik <- n*(lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    		         sum(lgamma(a + yvec) + lgamma(b + mvec - yvec) - 
    		             lgamma(a + b + mvec))
    		lprior <-  -log(a/(a+b)) - log(1-(a/(a+b))) - 2*log(1+a+b) 
    		lJacob <- log(b) - 2*log(a) - 2*log(1+(b/a)) + log(a) - log(1+(b/a))
    		out <-  llik + lprior + lJacob
    		return(out)
		}

		l.unpostAlbert <- apply(phi.grid, 1, l.unpost.Albert, yvec=rattumor.data[,1], 
 						mvec=rattumor.data[,2])
 		unpostAlbert <- exp(l.unpostAlbert - max(l.unpostAlbert))
 		unpost.mat.Albert <- matrix(unpostAlbert, nrow=length(phi.1), ncol=length(phi.2))

 		contour(phi.1, phi.2, unpost.mat.Albert, levels=seq(0.05, 0.95, 0.1),
           ylab=expression(log(alpha + beta)),xlab=expression(log(alpha/beta)))
           
        ## Plot contour from first prior over the top?
          contour(phi.1, phi.2, unpost.mat, levels=seq(0.05, 0.95, 0.1), col="blue", add=TRUE)
        ## Results are pretty close to the prior used in BDA3    


#### 3. Implement using the prior obtained using 
####    p(logit(eta), log(K)) propto 1
###    This is uniform on APPROXIMATELY the log(SD) which can typically
###    lead to improper posterior distributions because Equation 5.8 in BDA3,
###    which is given as the llik part of the code in the l.unpost functions
###    has an infinite integral as the population SD approaches 0.
####   See Problem 5.9(a) in BDA3: 

		l.unpost.fun3 <- function(phi,yvec,mvec) {
    		n <- length(yvec)
    		b <- exp(phi[2])/(1+exp(phi[1]))
    		a <- b*exp(phi[1])
    		llik <- n*(lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    		         sum(lgamma(a + yvec) + lgamma(b + mvec - yvec) - 
    		             lgamma(a + b + mvec))
    		lprior <-  log(1) 
    		lJacob <- log((1+ b^2)) - log(a*b*(a+b)) ##same Jacob as BDA3
    		out <-  llik + lprior + lJacob
    		return(out)
		}
		
	    l.unpost3 <- apply(phi.grid, 1, l.unpost.fun3, yvec=rattumor.data[,1], 
 						mvec=rattumor.data[,2])
 		unpost3 <- exp(l.unpost3 - max(l.unpost3))
 		unpost.mat.3 <- matrix(unpost3, nrow=length(phi.1), ncol=length(phi.2))


 		contour(phi.1, phi.2, unpost.mat.3, levels=seq(0.05, 0.95, 0.1),
           ylab=expression(log(alpha + beta)),xlab=expression(log(alpha/beta)),col="green")
        contour(phi.1, phi.2, unpost.mat, levels=seq(0.05, 0.95, 0.1), col="blue", add=TRUE)
        contour(phi.1, phi.2, unpost.mat.Albert, levels=seq(0.05, 0.95, 0.1), col=1,
                add=TRUE)
           #Albert is in between the uniform and the BDA3 prior
           
        #Contour looks okay for this problem BUT will be improper as (alpha+beta) -> infty
        # (prec goes to infinity and variance goes to 0).  In this limit p(y|alpha,beta) has
        # infinite integral if (alpha/beta = c) AND
        # the integral of log(alpha+beta) also goes to infinity - so posterior goes to 
        # infinity.  So, if there was very little variance in the pi's in the data
        # there could be problems. (See notes on yellow paper)
        #FUTURE: Simulate some artificial data with very little pop. variance
        #         to illustrate the problem with actual data
        

#### 4. Implement using the prior obtained using p(alpha, beta) propto 1
#### With this prior, we are essentially checking if p(y|alpha,beta) has
#### finite integral  (Allison's question from class on Friday)
####  
       l.unpost.fun4 <- function(phi,yvec,mvec) {  
			n <- length(yvec)
    		b <- exp(phi[2])/(1+exp(phi[1]))
    		a <- b*exp(phi[1])
    		llik <- n*(lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    		         sum(lgamma(a + yvec) + lgamma(b + mvec - yvec) - 
    		             lgamma(a + b + mvec))
    		lprior <-  log(1) 
    		lJacob <- log(1)
    		out <-  llik + lprior + lJacob
    		return(out)
		}
		
	    l.unpost4 <- apply(phi.grid, 1, l.unpost.fun4, yvec=rattumor.data[,1], 
 						mvec=rattumor.data[,2])
 		unpost4 <- exp(l.unpost4 - max(l.unpost4))
 		unpost.mat.4 <- matrix(unpost4, nrow=length(phi.1), ncol=length(phi.2))
 		
 		cont.levels <- c(seq(0.0001,.0009,0.0001), seq(0.01, 0.99, 0.1))
 		contour(phi.1, phi.2, unpost.mat.4, levels=cont.levels,
           ylab=expression(log(alpha + beta)),xlab=expression(log(alpha/beta)), col="orange")
        contour(phi.1, phi.2, unpost.mat.3, levels=cont.levels, col="green", add=TRUE)
        contour(phi.1, phi.2, unpost.mat, levels=cont.levels, col="blue", add=TRUE)
        contour(phi.1, phi.2, unpost.mat.Albert, levels=cont.levels, col=1,
                add=TRUE)

		### Plot the function  p(y|alpha,beta)
		### CHECK THIS -- did in a hurry
		 lik.ab <- function(ab,yvec,mvec) {  
			n <- length(yvec)
			a <- ab[1]
    		b <- ab[2]	
    		llik <- n*(lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    		         sum(lgamma(a + yvec) + lgamma(b + mvec - yvec) - 
    		             lgamma(a + b + mvec))
    		return(llik)
		   }
		
		a <- seq(1,20,length=400)  
 		b <- seq(1,60,length=600)   
 		ab.grid <- expand.grid(a,b)
 		llik.ab.rat <- apply(ab.grid, 1, lik.ab, yvec=rattumor.data[,1], mvec=rattumor.data[,2])
 		lik.ab.rat <- exp(llik.ab.rat - max(llik.ab.rat))
 		lik.ab.rat.mat <- matrix(lik.ab.rat, nrow=length(a), ncol=length(b))
   
		contour(a, b, lik.ab.rat.mat, ylab=expression(beta), xlab=expression(alpha),
		          levels=seq(0.0001, 0.95, 0.05))
		  abline(0,1,col="blue") 
		  
		#If alpha=beta, then posteriors going out to ....
		        
		#Try for fake data simulated from same prob, m=25, n=20
		 #set.seed(1342)
		 #y.fake <- rbinom(20,size=25,prob=0.1)
		 y.fake <- rep(2,10)
         llik.ab.fake <- apply(ab.grid, 1, lik.ab, yvec=y.fake, mvec=rep(25,10))
         lik.ab.fake <- exp(llik.ab.fake - max(llik.ab.fake))
 		 lik.ab.fake.mat <- matrix(lik.ab.fake, nrow=length(a), ncol=length(b))
		 contour(a, b, lik.ab.fake.mat, ylab=expression(alpha), xlab=expression(beta),
		          levels=seq(0.0001, 0.8, 0.01))
		  abline(0,1,col="blue") 
		  
	 
##### *******************************************  #####
#### 5. Implement using the prior  p(eta, K) propto 1
####   This prior should also lead to an IMPROPER Posterior   
####   
    
	# Write a function to get the log of the unnormalized posterior, so
	# we can look at contour lines (make sure not going to infinity) 
		l.unpost.fun5 <- function(phi,yvec,mvec) {  
			n <- length(yvec)
    		b <- exp(phi[2])/(1+exp(phi[1]))
    		a <- b*exp(phi[1])
    		llik <- n*(lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    		         sum(lgamma(a + yvec) + lgamma(b + mvec - yvec) - 
    		             lgamma(a + b + mvec))
    		lprior <-  log(1) 
    		lJacob <- log(b) - 2*log(a) - 2*log(1+(b/a)) + log(a) - log(1+(b/a))
    		out <-  llik + lprior + lJacob
    		return(out)
		}
		
	    l.unpost5 <- apply(phi.grid, 1, l.unpost.fun5, yvec=rattumor.data[,1], 
 						mvec=rattumor.data[,2])
 		unpost5 <- exp(l.unpost5 - max(l.unpost5))
 		unpost.mat.5 <- matrix(unpost5, nrow=length(phi.1), ncol=length(phi.2))
 		
 		contour(phi.1, phi.2, unpost.mat.5, levels=seq(0.01, 0.99, 0.01),
           ylab=expression(log(alpha + beta)),xlab=expression(log(alpha/beta)))
 		contour(phi.1, phi.2, unpost.mat.4, col="magenta", levels=seq(0.01, 0.99, 0.01))
 		    #Really close to upost.mat.5, but they are not exactly the same
        contour(phi.1, phi.2, unpost.mat.3, levels=seq(0.05, 0.95, 0.1), col="orange", add=TRUE)
        contour(phi.1, phi.2, unpost.mat, levels=seq(0.05, 0.95, 0.1), col="blue", add=TRUE)
        contour(phi.1, phi.2, unpost.mat.Albert, levels=seq(0.05, 0.95, 0.1), col="green",
                add=TRUE) 
        

##### *******************************************  #####
#### 6. What about (logit(eta), log(K)) ~ Unif[-10^10, 10^10]x[-10^10, 10^10]
####
#### Technically IS a proper prior distribution, but if used, most of the posterior
#### mass will be in the range of alpha and beta near infinity, which leads to
#### same problems in terms of posterior being sensitive to supposedly
#### "non-informative" prior
####  TO-DO -- Add in the future


###### ***************************************************#######
#### Now compare results for center of posterio from prior used in BDA 3 to the 
####  the PLUG IN method of estimating alpha and beta
### (a) from first 70 experiments (see page 102 in BDA3 for description)
### (b) from all 71 experiments (the empirical Bayes analysis)

##(a). Plot the point obtained using method 2 of plugging in
## point estimates from the first 70 experiments.  More detail 
## about this is given in Chapter 5 BDA3.  Just finding point
## estimates for alpha and beta via MoM using the first 70
## experiments and then transforming 
	historical.avg <- mean(rattumor.data[1:70,1]/rattumor.data[1:70,2])
	historical.var <- var(rattumor.data[1:70,1]/rattumor.data[1:70,2])
	c(historical.avg, historical.var)
	#Using result on page 583 in BDA3 for method of moments
	hist.K <- ((historical.avg*(1-historical.avg)) - historical.var)/historical.var
	hist.eta <- historical.avg
	#can go all the way to alpha and beta from these
	hist.alpha <- hist.eta*hist.K
	hist.beta <- (1-hist.eta)*hist.K
	c(hist.alpha, hist.beta)  #1.356149, #8.615058 
	
	##ADD the point to the contour plot of the posterior
	 contour(phi.1, phi.2, unpost.mat, levels=seq(0.05, 0.95, 0.1),
           ylab=expression(log(alpha + beta)),xlab=expression(log(alpha/beta)))

	 points(log(hist.eta), log(hist.K), col="green", cex=3, pch=16)  

 
 ###(b) Compare to the empirical Bayes point estimate using all 71 experiments    
     all.avg <- mean(rattumor.data[,1]/rattumor.data[,2])
     all.var <- var(rattumor.data[,1]/rattumor.data[,2])
     all.K <- ((all.avg*(1-all.avg)) - all.var)/all.var
     all.eta <- all.avg
     
     ##ADD the point to the contour plot of the posterior
     points(log(all.eta), log(all.K), col="lightblue", cex=3, pch=17)
    



## EXTRA    #####
################################################
### Graphically investigating the var(theta) approx 1/(alpha + beta) claim
####
###See algebra in notes and gone over in class
fun1 <- function(a,c=1){c+(c/a)}
fun2 <- function(a,c=1){(1/c)+(1/a)}
curve(fun1(x,c=1), from=0, to=10, col=3, lwd=4, ylim=c(0,10))
curve(fun2(x,c=1), from=0, to=10, col=4, lwd=2, add=TRUE)

#### How does variance change as alpha + beta changes?
curve(dbeta(x,0.5,0.5), from=0, to=1, col=4, lwd=5, ylim=c(0,3))
curve(dbeta(x,1,1), add=TRUE, col=5, lwd=3)
curve(dbeta(x,2,2), add=TRUE, col=3, lwd=3)
curve(dbeta(x,4,4), add=TRUE, col=2, lwd=2)
curve(dbeta(x,10,10), add=TRUE, col=1, lwd=1)


### is variance of the beta approximately 1/(a+b)?
### Just some code I was playing around with -- not cleaned up
var.beta <- function(a,b) { (a*b)/(((a+b)^2)*(a+b+1))}
var.beta.approx <- function(a,b) {1/(a+b)} #From Gelman Problem 5.8
var.beta.approx2 <- function(a,b) {1/(a+b+1)^2} #better approx of variance and happens to be
# prior used by Albert?
sd.beta.approx <- function(a,b) {(a+b)^(-0.5)}

curve(var.beta(x,b=0.01), from=0, to=1, col=2, ylim=c(0,5))
curve(var.beta.approx(x,b=0.01), add=TRUE, col=3)
curve(var.beta.approx2(x,b=0.01), add=TRUE, col=4)

curve(var.beta.approx(x,b=0.01), from=0, to=1, col=3)
curve(var.beta(x,b=0.01), add=TRUE, col=2)
curve(var.beta.approx2(x,b=0.01), add=TRUE, col=4)

a.vec <- seq(0.01, 5, length=100)
b.vec <- seq(0.01, 5, length=100)
var.beta.out <- outer(a.vec, b.vec, var.beta)
var.beta.approx.out <- outer(a.vec, b.vec, var.beta.approx)
var.beta.approx2.out <- outer(a.vec, b.vec, var.beta.approx2)
hist((var.beta.out/var.beta.approx2.out), nclass=100, col=gray(0.5))
hist((var.beta.approx2.out/var.beta.out), nclass=100, col=gray(0.5))
hist(var.beta.out - var.beta.approx2.out, nclass=100, col=gray(0.5))

        
        
##################################################################
### Compare INVERSE GAMMA on sigma2 to GAMMA on precision     ####
###
 curve(dgamma(x,0.1,0.1), from=0, to=10, col=2, lwd=2, xlab=expression(tau))
 curve(dgamma(x,1000,1000),from=0, to=10, col=2, lwd=2, xlab=expression(tau))
        
 a <- 0.1 #0.1
 b <- 0.1  #0.1
 tau_draws <- rgamma(1000, a, b)
 summary(tau_draws)
 hist(tau_draws, col=gray(0.8), nclass=100, xlim=c(0,20), 
             freq=FALSE, xlab=expression(tau), main="Tau draws from Gamma")
        curve(dgamma(x,a,b), from=0, to=20, col=2, lwd=2,  add=TRUE)
        curve(1/x, from=0, to=20, col=3, lwd=2, add=TRUE)
        curve((1/(1+x)), from=0, to=20, col=4, lwd=2, add=TRUE)
        curve((1/((1+x)^2)), from=0, to=20, col=5, lwd=3, add=TRUE)
        var(tau_draws)
        
  sig2_draws <- 1/tau_draws
  var(sig2_draws)  #huge
        hist(sig2_draws, col=gray(0.8), nclass=100, xlim=c(0,10000), 
             freq=FALSE, xlab=expression(sigma^2), main="Sigma2 draws from Inv-Gamma")
        
  sum(tau_draws<0.0001)  ##causing plotting issues b/c they are so huge
   sig2_draws2 <- 1/tau_draws[tau_draws>0.0001]  #approx.
   hist(sig2_draws2, col=gray(0.8), nclass=100, xlim=c(0,500), 
             freq=FALSE, xlab=expression(sigma^2), main="Sigma2 draws from Inv-Gamma")
        curve(1/x, from=0, to=10000, col=3, lwd=2, add=TRUE)
        curve((1/((1+x)^2)), from=0, to=70, col=5, lwd=3, add=TRUE)
        
        
    min(tau_draws)
    1/min(tau_draws)
        
      

###########################################
### Play around with hyperparameters for binomial logit normal #####
##############################################

## Assuming variance of mu is related to precision/variance on the logit(pi)'s
sig2 <- 1
tau <- 1/sig2
k0 <- 1
theta.0 <- 0
curve(dnorm(x,mean=theta.0, sd=sqrt(sig2/k0)), from=-5, to=5, lwd=2, col=2)
  curve(dnorm(x,mean=theta.0, sd=sqrt(sig2/10)), lwd=2, col=2, add=TRUE)
  curve(dnorm(x,mean=theta.0, sd=sqrt(sig2/0.1)), lwd=2, col=2, add=TRUE)
  curve(dnorm(x,mean=theta.0, sd=sqrt(10/k0)), lwd=2, col=3, add=TRUE)
  curve(dnorm(x,mean=theta.0, sd=sqrt(100/k0)), lwd=2, col=3, add=TRUE)

curve(dgamma(x,0.1,0.1), from=0, to=10, col=2, lwd=2, xlab=expression(tau))
curve(dgamma(x,1000,1000),from=0, to=10, col=2, lwd=2, xlab=expression(tau))

a <- 1  #0.1
b <- 1  #0.1
curve(dgamma(x,a,b), from=0, to=20, col=2, lwd=2)
tau_draws <- rgamma(1000, a, rate=b)
hist(tau_draws, col=gray(0.8), nclass=100, xlim=c(0,20), 
     freq=FALSE, xlab=expression(tau), main="Tau draws from Gamma")
curve(dgamma(x,a,rate=b), from=0, to=20, col=2, lwd=2,  add=TRUE)

sig2_draws <- 1/tau_draws
hist(sig2_draws, col=gray(0.8), nclass=100, 
     freq=FALSE, xlab=expression(sigma^2), main="Sigma2 draws from Inv-Gamma")
d.invgamma <- function(x,a,b) {
   (b^(a)/gamma(a))*x^(-(a+1))*exp(-b/x)
}
curve(d.invgamma(x,a,b), col=5, lwd=2, from=0, to=100, add=TRUE)





#### IGNORE THIS --- started something that I didn't finish  ####
###########################################
### Write function evaluate p(y|pi) and p(y|phi) ###
### for different values of pi (or phi)  #####
##############################################

log.lik.fun.pi <- function(pie,y,m) {
                sum(dbinom(y, size=m, prob=pie, log=TRUE))
               }

log.lik.fun.pi(0.5,y.data, m=20)

log.lik.fun.phi <- function(phi,y,m) {
                pie <- expit(phi)
                sum(dbinom(y, size=m, prob=pie, log=TRUE))
               }

log.lik.fun.phi(logit(0.5),y.data, m=20)  #do they match?
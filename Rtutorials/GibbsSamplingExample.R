

#### Very simple Gibbs Sampling algorithm to get draws for mu and sigma2 from a 
#### Normal data model Yi ~ iid N(mu, sig2)
####  with the following prior distributions:
###  mu ~ Normal(theta, tau2)
###  sig2 ~ Inv-Gamma(a,b)
###  This is the semi-conjugate prior distribution described in Section 3.4
### It is NOT a conjugate family for the normal likelihood and the posterior
###  does not follow any standard parametric form.  In section 3.4, Gelman etal.
###  show how one can get useful results by using the conditional posterior of 
###  of mu|sig2,y (i.e. the complete conditional distribution!) and then computing
###  p(sig2|y) on a discrete grid (READ section 3.4) so that we can sample from 
###  p(sig2|y) and then p(mu|sig2,y).
### NOW, let's use a Gibbs sampling algorithm to solve the problem instead.
###  This is convenient because we already can take samples from the normal
### distribution of p(mu|sig2,y) and the complete cond. for sig2 (p(sig2|mu,y)) 
### also has an easy form to sample from, since it is just a Inv-Gam distribution!
### Therefore, the Gibbs sampling algorithm gets us around having to use the 
### grid approximation to sample from the marginal posterior of sig2.

set.seed(1543)
data <- rnorm(21, mean=30, sd=3)
ybar <- mean(data)
s2 <- var(data)
n <- length(data)

#### First, find the complete conditional distribution for mu: (Equation 3.11)
## mu|sig2,y,s2 ~ Normal(mu.cc.mean, mu.cc.var)
 # mu.cc.mn <- ((theta/tau2) + ((n*ybar)/sig2))/((1/tau2) + (n/sig2))       
 # mu.cc.var <- 1/((1/tau2) + (n/sig2))
 #Note these complete conditional params are functions of the unknown sigma2!
  
### Second, find the complete conditional distribution for sig2:
###    sig2|mu,y ~ InvGamma(sig.cc.a, sig.cc.b)  or  1/sig2|mu,y ~ Gamma(sig.cc.a, sig.cc.b)
 # a and b are hyperparameters for prior on sigma2 (see above)
 #  sig.cc.a <- a + n/2
 #  sig.cc.b <- ((1/2)*sum((data-mu)^2))+b
   
   
#Now write the Gibbs sampling algorithm:

 #set hyperparameters
 # Using here what many use as "default" priors.  I am not condoning or suggesting they
 # are the most appropriate prior parameters to use.
  theta <- 0
  tau2 <- 100
   a <- 0.001
   b <- 0.001

  nsim <- 10000
  mu.vec <- numeric(nsim)
  sig2.vec <- numeric(nsim)
  mu.vec[1] <- 20  #don't really need this if start by drawing mu
  sig2.vec[1] <- 1  
  
  for (t in 2:nsim) {
  	
  	  sig2.cur <- sig2.vec[t-1]
  	  
  	  #Draw a value of mu given the current value of sig2
  	  mu.cc.mn <- ((theta/tau2) + ((n*ybar)/sig2.cur))/((1/tau2) + (n/sig2.cur))
  	  mu.cc.var <- 1/((1/tau2) + (n/sig2.cur))
  	  mu.cur <- rnorm(1, mean=mu.cc.mn, sd=sqrt(mu.cc.var))
  	  mu.vec[t] <- mu.cur
  	
  	  #Now draw a value of sigma2 given the new value of mu
  	  sig.cc.a <- a + (n/2)
  	  sig.cc.b <- ((1/2)*sum((data-mu.cur)^2))+b
  	  phi <- rgamma(1, sig.cc.a, sig.cc.b)
  	  sig2.vec[t] <- 1/phi 	  
  	}

  burn.in <- 1000
  mu.post <- mu.vec[(burn.in+1):nsim]
  sig2.post <- sig2.vec[(burn.in+1):nsim]

  dev.new()
  plot(mu.post, sig2.post, type="n", main="Draws from Joint Posterior Distribution",
        xlab=expression(mu),  ylab=expression(sigma^2))
     points(mu.post, sig2.post, pch=1, col=4)
     
  dev.new()
  hist(mu.post, col=gray(0.85), nclass=100, freq=F, main="Posterior for mu")
    lines(density(mu.post), col=2)
    abline(v=30, lwd=3, col="darkgreen") #true mean
    abline(v=ybar, lwd=3, col=2)  #sample mean
    
  dev.new()
  hist(sig2.post, col=gray(0.85), nclass=100, freq=F, main="Posterior for sigma2")
    lines(density(sig2.post), col=2)
    abline(v=9, lwd=3, col="darkgreen") #true sig2
    abline(v=s2, lwd=3, col=2)  #sample mean
    
   
  ### Assessing convergence - an intro #####
  ## To get Rhat - we would need to run at least 2 parallel chains ###


    #Make Sample Path plots
     par(mfrow=c(2,1))
     plot(1:nsim, mu.vec, type="n", main="Sample Path History Plot", 
     ylab=expression(mu),  xlab="iteration")
       lines(1:nsim, mu.vec, col="red")
     plot(1:nsim, sig2.vec, type="n", main="", ylab=expression(sigma^2),
            xlab="iteration")
       lines(1:nsim, sig2.vec, col="blue")


    #Make Sample Path plots
     par(mfrow=c(2,1))
     plot(2000:2500, mu.vec[2000:2500], type="n", main="Sample Path History Plot", 
           ylab=expression(mu), xlab="iteration")
       lines(2000:2500, mu.vec[2000:2500], col="red")
     plot(2000:2500, sig2.vec[2000:2500], type="n", main="", ylab=expression(sigma^2),
            xlab="iteration")
       lines(2000:2500, sig2.vec[2000:2500], col="blue")

   ### Calculate the sample lag 1 autocorrelation - Both are small!
     acf(mu.post, plot=FALSE)[1]
     acf(sig2.post, plot=FALSE)[1]


   ### Plot white noise just to remind us what it looks like
     par(mfrow=c(2,1))
     plot(2000:2500, mu.vec[2000:2500], type="n", main="", ylab=expression(mu),
            xlab="iteration")
       lines(2000:2500, mu.vec[2000:2500], col="red")
     plot(2000:2500, mu.vec[2000:2500], type="n", main="", ylab=expression(mu),
            xlab="iteration")
       lines(2000:2500, rnorm(501, mean(mu.post), sd=sd(mu.post)), col=1)


  ###Use the CODA package
    #library(coda)
    #help(package="coda") #for list of functions
    mu.sig2.post <- mcmc(cbind(mu.post, sig2.post))
    summary(mu.sig2.post)
    plot(mu.sig2.post) 
    effectiveSize(mu.sig2.post)
    #gelman.diag(mu.sig2.post)  #Can't do b/c need at least 2 chains
    #gelman.plot(mu.sig2.post)
    

    



  
 

#### R-code for INTRODUCTION to Posterior Predictive Distributions in 
### one-parameter normal model.

#### Posterior predictive distribution for Normal model with known variance 
####  (sig2=1) and n=1 where y=7
###  Obtain computationally and then check analytically:

 norm.lik.fun <- function(mu, sig2=1, y) {
 	out <- dnorm(y, mean=mu, sd=sqrt(sig2))
 	return(out)
 	}
 	
 #Use normal prior distribution with mu.0=5, tau2.0=5
 
  ###Computationally find posterior distribution and draw samples from it
  ## Just doing this for practice! (really no need b/c know posterior and can
  ## easily draw from it)
  
   mu.vec <- seq(3,10, length=1000) ##grid over mu
   prior.vec <- apply(cbind(mu.vec), 1, function(mu) {dnorm(mu,mean=5,sd=sqrt(5))})
   lik.vec <- apply(cbind(mu.vec), 1, norm.lik.fun, y=7)
   lik.x.prior.vec <- lik.vec*prior.vec   #multiply them together
   interval.width <- mu.vec[2]-mu.vec[1]
   norm.constant <- sum(lik.x.prior.vec*interval.width)
   posterior.vec <- lik.x.prior.vec/norm.constant
   
  curve(dnorm(x, mean=5, sd=sqrt(5)), from=0, to=10, xlab=expression(mu), ylab="Density", 
           lty=1, lwd=3, col=2, ylim=c(0,0.75))  #prior
     lines(mu.vec, posterior.vec, lty=1, col="orange", lwd=3)  #add to previously made plot
     abline(v=7) #the one piece of data
     curve(dnorm(x,mean=(20/3), sd=sqrt(5/6)), lty=2, col="magenta", lwd=3, add=TRUE) 
         #analytical posterior distribution  It matches!  :)
     legend(1,0.7, legend=c("prior", "comp. post", "analyt. post"), lwd=c(3,3,3), lty=c(1,1,2),
              col=c(2,"orange","magenta"), bty="n")

 
 #Obtain draws from the posterior predictive distribution: We really want to 
 #  repeat these two steps over and over again: 
 #  (1) draw a value of mu from the posterior distribution
 #  (2) draw a value of y from the p(y|mu) using the value of mu draw in (1)
 # 
 #In other words, for each draw from the posterior distribution, obtain a draw from the 
 # posterior pred. distribution.  We could do this in a loop, but it is easier 
 # to just do it all at once in vector/matrix form
 
   ##RUN ONE LINE AT A TIME TO SEE WHAT's GOING ON
   # Draw samples of mu from the posterior distribution 
   n.draws <- 2000
   post.samps <-  sample(mu.vec, size=2000, prob=posterior.vec, replace=TRUE) 
   hist(post.samps, freq=FALSE, col="lightgray", nclass=75, add=TRUE)  #Look reasonable?
   
   #Now, for each draw of mu from the posterior obtained above,
   #draw from the normal distribution with those means!
   
   post.pred.samps <- apply(cbind(post.samps), 1, function(mu) rnorm(1,mean=mu, sd=1))
   hist(post.pred.samps, freq=FALSE, col="yellow", nclass=75, add=TRUE)  #Look reasonable?
   curve(dnorm(x,mean=(20/3), sd=sqrt(1 + (5/6))), lty=1, col="green", lwd=4, add=TRUE) 
       #analytical posterior predictive distribution  It matches!  :)
 
   plot(post.pred.samps, post.samps, ylab=expression(mu), xlab="y pred", xlim=c(0,11), ylim=c(0,10))
     hist(post.pred.samps, freq=FALSE, col="yellow", nclass=75, add=TRUE) 
 
     hist(post.pred.samps, freq=FALSE, col="yellow", nclass=75)  

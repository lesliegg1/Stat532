
### Based on example 3.7 in Carlin & Louis and data from Bliss(1935)
library(mvtnorm)

### Data: Number of adult flour beetles killed after 5 hours of exposure
##   to various levels of gaseous carbon disulphide (CS2)

dose <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
killed <- c(6, 13, 18, 28, 52, 53, 61, 60)
exposed <- c(59,60,62,56,63,59,62,60)

### Some exploratory plots 

prop.k <- killed/exposed
emp.logits <- log(prop.k/(1-prop.k))

#Plot data and check linearity assumptions
par(mfrow=c(2,1))
plot(dose, prop.k, pch=16, col=3, cex=1.2)
plot(dose, emp.logits, pch=16, col=3, cex=1.2) #note missing one because p.hat=1


#Fit model quickly in glm()
options(show.signif.stars=FALSE)

yes.no <- cbind(killed, (exposed-killed))
glm.out <- glm(yes.no~dose, family=binomial)
summary(glm.out)  #Check for overdispersion (Resid Dev=11.232 vs 6 d.f.)

dev.new()
plot(dose, emp.logits, pch=16, col=3, cex=1.2) 
abline(glm.out) #model seems reasonable except for a bit of overdispersion
                # 1-pchisq(11.232,6) [1] 0.08146544  -dispersion approx 11.232/6
                #

#We could use quasilikelihood to account for the overdispersion
yes.no <- cbind(killed, (exposed-killed))
quasi.glm.out <- glm(yes.no~dose, family=quasibinomial)
summary(quasi.glm.out)  #1.671141= estimated dispersion param from Pearson resids

  #Do conclusions change?

##### Let's fit the generalized logit model suggested by Prentice (1976) ##
## kill ~ Binomial(n.exposed, pi(kill|dose))
## pi(kill|dose) = [exp(x)/(1+exp(x))]^m1
## Let x = (dose-mu)/sigma    sigma>0, m1>0
### Do you know how to do this easily in standard software? ###
### This changes the number of parameters from 2 to 3 - the extra
###  is needed to model the dispersion (a variance parameter)


### Specify prior distributions
# m1 ~ gamma(a0, b0)  a0=.25, b0=0.25
# mu ~ Normal(c0, d02)
# sig2 ~ InvGamma(eo, fo)

### Is the posterior available in closed form?
### What are the three complete conditional distributions? Are
###  they available in closed form?  NO

### Let's try to understand the model better by plotting logits of pi^1/m1
##  Convince yourself that all this makes sense! (or convince me it doesn't!)
emp.pm1.logits <- function(m1) {log((prop.k^(1/m1))/(1-(prop.k^(1/m1))))} #prop.k = observed proportions
dev.new()
plot(dose, emp.logits, pch=16, col=1, cex=1, ylab="logits", type="b", ylim=c(-10,10)) #note missing one because p.hat=1
 points(dose, emp.pm1.logits(exp(-1.05)), col=2, type="b", pch=16)
 points(dose, emp.pm1.logits(0.5), col=3, type="b", pch=16)
 points(dose, emp.pm1.logits(1.5), col=4, type="b", pch=16)
 points(dose, emp.pm1.logits(2), col=5, type="b", pch=16)
 legend(1.80, -4, legend=c("m1=exp(-1.05)", "m1=0.5", "m1=1 (observed logits)", 
           "m1=1.5", "m1=2.0"), lwd=c(1,1,1,1,1), pch=c(16,16,16,16,16), 
            col=c(2,3,1,4,5), bty="n")

 lm(emp.pm1.logits(exp(-1.05))[1:7]~dose[1:7])$coef
    #  (Intercept)   dose[1:7] 
    # -98.75    54.51125
    # This implies that the slope (1/sigma) should be around 55 if m1 is around exp(-1.05)=0.35
    #   so that sig=(1/54.5)=0.0183 and log(sig)=log(0.0183) = -4.00, 
    # The intercept (-mu/sig) should be around -100 -->mu=100/(54.5) = 0.545
  

## But, the above still doesn't tell what value for m1 is reasonable from the data
## Let's go back to just plotting the observed proportions
pi.dose.fun <- function(w.vec, mu, sigma, m1){
	((exp((w.vec-mu)/sigma))/(1 +exp((w.vec-mu)/sigma)))^m1
}

plot(dose, prop.k, pch=16, type="b", col=1, cex=1.2, ylim=c(0,1))
 points(dose, pi.dose.fun(dose, mu=1.77 ,sig=0.0183, m1=1), col=4, type="b", pch=18) #Logistic regression
 #points(dose, pi.dose.fun(dose, mu=1.54 ,sig=exp(-4.87), m1=exp(2.18)), col=5, type="b", pch=18) #Old #posterior mode
 points(dose, pi.dose.fun(dose, mu=1.835, sig=(1/54.5), m1=0.35), col=2, type="b", pch=16) #Linear regression logit scale (above)
 points(dose, pi.dose.fun(dose, mu=1.812, sig=exp(-4.01), m1=exp(-1.06)), col=3, type="b", pch=16) #mode of posterior
 #points(dose, pi.dose.fun(dose, mu=1.55, sig=54.5, m1=0.35), col=4, type="b", pch=16) #
 #points(dose, pi.dose.fun(dose, mu=1.55, sig=54.5, m1=0.35), col=5, type="b", pch=16)
 legend(1.79, 0.3, legend=c("Observed proportions",  "From Logistic Regression", "Old Posterior Mode","Crude linear regression check", "Mode of posterior"), 
     lwd=c(1,1,1,1,1),pch=c(16,18,16,16,16), col=c(1,4,5,2,3), bty="n")
  
 #mode of posterior using optim 1.812353 -4.010371 -1.056079

##Write function to compute log likelihood for theta=(mu, sig2, m1) given the data
llik.fun <- function(theta.vec, dose.vec, y.vec, n.vec) {
       mu <- theta.vec[1]
       sig <- theta.vec[2]
       m1 <- theta.vec[3]
       x.vec <- (dose.vec - mu)/sig
	   llik.vec <- m1*y.vec*(x.vec-log(1+exp(x.vec))) + (n.vec - y.vec)*log(1-((exp(x.vec)/(1+exp(x.vec)))^m1))
       out <- sum(llik.vec)
       return(out)
	}

 #test function - see if it matches for mu=0, sig=1, m1=1 which is just standard logistic regression 
    llik.fun(c(0,1,1), dose.vec=dose, y.vec=killed, n.vec=exposed)
    sum(log(dbinom(killed, exposed, prob=(exp(dose)/(1+exp(dose)))))) - sum(lchoose(exposed, killed))
   #they match!
   
   llik.fun(c(1.2, 0.03, 0.35), dose.vec=dose, y.vec=killed, n.vec=exposed)

###Transformed the prior parameters to the real line
## I worked out the prior distributions for the transformed parameters
##  on paper (you should do the same to check my work!). Transformed parameters:
##  phi1 = mu
##  phi2 = log(sigma) or (0.5)*log(sigma^2) ##CHECK THIS
##  phi3 = log(m1)
##  Priors are rather vague, but proper (See Carlin & Gelfand (1991b) for more info on specification)

### This first function takes the original params as inputs
l.prior.fun <- function(theta.vec, a0=0.25, b0=0.25, c0=2, d0=10, e0=2.000004, f0=0.001) {
      phi1 <- theta.vec[1]
      phi2 <- log(theta.vec[2])
      phi3 <- log(theta.vec[3])
      log.p.phi1 <- log(dnorm(phi1, mean=c0, sd=d0)) 
      log.p.phi2 <- (-2*e0*phi2) - (f0*(exp(-2*phi2)))
      log.p.phi3 <- (phi3*a0) - (b0*exp(phi3))
      log.p <- log.p.phi1 + log.p.phi2 + log.p.phi3 #assuming priors independent
      return(log.p)
    }

### This second function takes the transformed params as inputs
l.prior.fun2 <- function(phi.vec, a0=0.25, b0=0.25, c0=2, d0=10, e0=2.000004, f0=0.001) {
      phi1 <- phi.vec[1]
      phi2 <- phi.vec[2]
      phi3 <- phi.vec[3]
      log.p.phi1 <- log(dnorm(phi1, mean=c0, sd=d0)) 
      log.p.phi2 <- -2*e0*phi2 - (f0*(exp(-2*phi2)))
      log.p.phi3 <- phi3*a0 - (b0*exp(phi3))
      log.p <- log.p.phi1 + log.p.phi2 + log.p.phi3 #assuming priors independent
      return(log.p)
    }


   l.prior.fun(c(0.5, 0.03, 0.3))
   l.prior.fun2(c(0.5, log(0.03), log(0.3)))


l.unpost.fun <- function(theta.vec, dose.vec, y.vec, n.vec) {
     llik <- llik.fun(theta.vec, dose.vec=dose.vec, y.vec=y.vec, n.vec=n.vec)
     lp <- l.prior.fun(theta.vec)
     lout <- llik + lp
     return(lout)
    }

l.unpost.fun2 <- function(phi.vec, dose.vec, y.vec, n.vec) {
     theta.vec <- c(phi.vec[1], exp(phi.vec[2]), exp(phi.vec[3]))
     llik <- llik.fun(theta.vec, dose.vec=dose.vec, y.vec=y.vec, n.vec=n.vec)
     lp <- l.prior.fun2(phi.vec)
     lout <- llik + lp
     return(lout)
    }

 l.unpost.fun(c(1.8, 1, 1), dose.vec=dose, y.vec=killed, n.vec=exposed)
 l.unpost.fun2(c(1.8,log(1), log(1)), dose.vec=dose, y.vec=killed, n.vec=exposed)
   

### Plot the priors: (should also do this on transformed scale)
  library(pscl) #For the densiigamma() function with same parameterization as Gelman et al.
    par(mfrow=c(1,3))
    curve(dnorm(x,2,10), xlim=c(-30,30), main=expression(mu) )
    curve(densigamma(x,2.000004,0.001), xlim=c(0.001,.01), main=expression(sigma^2)) #prior mean of .001 and sd of .5 for sig2
      ### should plot transformed version too
    curve(dgamma(x,0.25,0.25), xlim=c(0,4), main=expression(m[1]))



### Now, let's proceed since the priors seem at least reasonable - they are vague, but
##  cover the appropriate range of plausible values for the problem. The complete 
##  conditional distributions are not available in closed form, so we cannot just 
##  proceed using a Gibbs sampling algorithm.

 ##Let's first assume we have no information about the covariance between phi1,phi2, and phi3
 ## and use a Multivariate Normal proposal distribution with a diag. var-cov matrix
 ## We could just play around with values to see what seems to be working in terms of a reasonable
 ## acceptance rate, but it's often easier to use optim() to get a good guess at the 
 ## mode and curvature of the posterior, thus giving us a good starting point

   l.unpost.fun(c(1.8, exp(-4), exp(-1)), dose.vec=dose, y.vec=killed, n.vec=exposed)
               
   optim.out <- optim(c(1.77, log(0.03), log(0.35)), l.unpost.fun2, dose.vec=dose, 
                    y.vec=killed, n.vec=exposed, control=list(fnscale=-100), 
                    method="Nelder-Mead", hessian=TRUE)
   	
   optim.out$par  #MODE [1]  1.812353 -4.010371 -1.056079	
   VarCov <- solve(-optim.out$hessian)
   VarCov
             [,1]         [,2]         [,3]
[1,]  0.0001073514 -0.001418677 -0.003012326
[2,] -0.0014186772  0.034213984  0.051382611
[3,] -0.0030123256  0.051382611  0.099574752


   #### Program the Multivariate M-H algorithm
   ## We are actually using a symmetric proposal distribution, so this
   ## could be programmed simply as a multivariate Metropolis, but I'm
   ## including the jumping distribution part into the M-H ratio so that
   ## you can see how it is done - You can check that the parts from the
   ## proposal distribution are equal and cancel out.

  ###Let's run three separate chains (which I make explicit using a lot
  ### of code below.

  ### FIRST USE A DIAGONAL VAR-COV MATRIX for the PROPOSAL DISTN
   nsim <- 10000  #number of iteration
   phi.mat1 <- matrix(NA, nrow=nsim, ncol=3) #chain 1
   phi.mat2 <- matrix(NA, nrow=nsim, ncol=3) #chain 2
   phi.mat3 <- matrix(NA, nrow=nsim, ncol=3) #chain 3
   jump.vec1 <- numeric(nsim-1) #keep track of when we jump (accept candidates)
   jump.vec2 <- numeric(nsim-1) 
   jump.vec3 <- numeric(nsim-1)
   phi.mat1[1,] <- c(1.8, -4, -1)  #1.812353 -4.010371 -1.056079
   phi.mat2[1,] <- c(2, -3, -1)
   phi.mat3[1,] <- c(1, -2,  -0.5)
   ind.varcov <- diag(diag(VarCov))

   for (i in 2:nsim) {
     phi.cur1 <- phi.mat1[i-1,]
     phi.cur2 <- phi.mat2[i-1,]
     phi.cur3 <- phi.mat3[i-1,]
     phi.cand1 <- rmvnorm(1, mean=phi.cur1, sigma=2*ind.varcov)
     phi.cand2 <- rmvnorm(1, mean=phi.cur2, sigma=2*ind.varcov)
     phi.cand3 <- rmvnorm(1, mean=phi.cur3, sigma=2*ind.varcov)

     log.r.num1 <- l.unpost.fun2(phi.cand1, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cur1, mean=phi.cand1, sigma=2*ind.varcov, log=TRUE)
     log.r.num2 <- l.unpost.fun2(phi.cand2, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cur2, mean=phi.cand2, sigma=2*ind.varcov, log=TRUE)
     log.r.num3 <- l.unpost.fun2(phi.cand3, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cur3, mean=phi.cand3, sigma=2*ind.varcov, log=TRUE)

     log.r.denom1 <- l.unpost.fun2(phi.cur1, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cand1, mean=phi.cur1, sigma=2*ind.varcov, log=TRUE)
     log.r.denom2 <- l.unpost.fun2(phi.cur2, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cand2, mean=phi.cur2, sigma=2*ind.varcov, log=TRUE)
     log.r.denom3 <- l.unpost.fun2(phi.cur3, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cand3, mean=phi.cur3, sigma=2*ind.varcov, log=TRUE)
 
     log.r1 <- log.r.num1 - log.r.denom1
     log.r2 <- log.r.num2 - log.r.denom2
     log.r3 <- log.r.num3 - log.r.denom3

     p.accept1 <- min(1, exp(log.r1))
     p.accept2 <- min(1, exp(log.r2))
     p.accept3 <- min(1, exp(log.r3))

     u.vec <- runif(3)
     ifelse(u.vec[1] <= p.accept1, phi.mat1[i,]<- phi.cand1, phi.mat1[i,] <- phi.cur1)
     ifelse(u.vec[2] <= p.accept2, phi.mat2[i,]<- phi.cand2, phi.mat2[i,] <- phi.cur2)
     ifelse(u.vec[3] <= p.accept3, phi.mat3[i,]<- phi.cand3, phi.mat3[i,] <- phi.cur3)

     jump.vec1[i-1] <- ifelse(u.vec[1] <= p.accept1, 1, 0)
     jump.vec2[i-1] <- ifelse(u.vec[2] <= p.accept2, 1, 0) 
     jump.vec3[i-1] <- ifelse(u.vec[3] <= p.accept3, 1, 0)
    }  

   ###Clearly I could have done a better job programming my multiple parallel chains in fewer
   ##   lines, but wanted to make it explicit

   ##Look at chains
     par(mfrow=c(3,1))
     plot(seq(1:nsim), phi.mat1[,1], type="l", ylab=expression(phi[1]), col=4, ylim=c(1,3))
        lines(seq(1:nsim), phi.mat2[,1], col=2)
        lines(seq(1:nsim), phi.mat3[,1], col=3)
     plot(seq(1:nsim), phi.mat1[,2], type="l", ylab=expression(phi[2]), col=4, ylim=c(-5,-2))
        lines(seq(1:nsim), phi.mat2[,2], col=2)
        lines(seq(1:nsim), phi.mat3[,2], col=3)
     plot(seq(1:nsim), phi.mat1[,3], type="l", ylab=expression(phi[3]), col=4, ylim=c(-2,2))
        lines(seq(1:nsim), phi.mat2[,3], col=2)
        lines(seq(1:nsim), phi.mat3[,3], col=3)

   ##Acceptance rates?  For multivariate - around 10-15% is okay
     mean(jump.vec1)   
     mean(jump.vec2)
     mean(jump.vec3)
     
     

   ##What do you think of the convergence???
   ### Think about correlation among parameters - use plots and calculations.
      pairs(~ phi.mat1[,1] + phi.mat1[,2] + phi.mat1[,3])
      var(phi.mat1)  #We could use these to create a better chain (though!
      var(phi.mat2) 
      var(phi.mat3) 
      #We can use the above estimated variance-covariance matrices to
      # help improve efficiency of the chain or use original optimization

       
    #### ******************************************************** ####
    #### *******************************************************  #####
   ### NOw, redo including off-diagonals of the variance covariance matrix
   nsim <- 1000  #number of iteration
   phi.mat1 <- matrix(NA, nrow=nsim, ncol=3) #chain 1
   phi.mat2 <- matrix(NA, nrow=nsim, ncol=3) #chain 2
   phi.mat3 <- matrix(NA, nrow=nsim, ncol=3) #chain 3
   jump.vec1 <- numeric(nsim-1) #keep track of when we jump (accept candidates)
   jump.vec2 <- numeric(nsim-1) 
   jump.vec3 <- numeric(nsim-1)
   phi.mat1[1,] <- c(1.8, -4, -1)  #1.543270 -4.874497  2.179276
   phi.mat2[1,] <- c(1.5, -3.8, -1.2)
   phi.mat3[1,] <- c(1.8, -4.5,  -0.8)
  
   for (i in 2:nsim) {
     phi.cur1 <- phi.mat1[i-1,]
     phi.cur2 <- phi.mat2[i-1,]
     phi.cur3 <- phi.mat3[i-1,]
     phi.cand1 <- rmvnorm(1, mean=phi.cur1, sigma=2*VarCov)
     phi.cand2 <- rmvnorm(1, mean=phi.cur2, sigma=2*VarCov)
     phi.cand3 <- rmvnorm(1, mean=phi.cur3, sigma=2*VarCov)

     log.r.num1 <- l.unpost.fun2(phi.cand1, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cur1, mean=phi.cand1, sigma=2*VarCov, log=TRUE)
     log.r.num2 <- l.unpost.fun2(phi.cand2, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cur2, mean=phi.cand2, sigma=2*VarCov, log=TRUE)
     log.r.num3 <- l.unpost.fun2(phi.cand3, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cur3, mean=phi.cand3, sigma=2*VarCov, log=TRUE)

     log.r.denom1 <- l.unpost.fun2(phi.cur1, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cand1, mean=phi.cur1, sigma=2*VarCov, log=TRUE)
     log.r.denom2 <- l.unpost.fun2(phi.cur2, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cand2, mean=phi.cur2, sigma=2*VarCov, log=TRUE)
     log.r.denom3 <- l.unpost.fun2(phi.cur3, dose.vec=dose, y.vec=killed, n.vec=exposed) + 
                     dmvnorm(phi.cand3, mean=phi.cur3, sigma=2*VarCov, log=TRUE)
 
     log.r1 <- log.r.num1 - log.r.denom1
     log.r2 <- log.r.num2 - log.r.denom2
     log.r3 <- log.r.num3 - log.r.denom3

     p.accept1 <- min(1, exp(log.r1))
     p.accept2 <- min(1, exp(log.r2))
     p.accept3 <- min(1, exp(log.r3))

     u.vec <- runif(3)
     ifelse(u.vec[1] <= p.accept1, phi.mat1[i,]<- phi.cand1, phi.mat1[i,] <- phi.cur1)
     ifelse(u.vec[2] <= p.accept2, phi.mat2[i,]<- phi.cand2, phi.mat2[i,] <- phi.cur2)
     ifelse(u.vec[3] <= p.accept3, phi.mat3[i,]<- phi.cand3, phi.mat3[i,] <- phi.cur3)

     jump.vec1[i-1] <- ifelse(u.vec[1] <= p.accept1, 1, 0)
     jump.vec2[i-1] <- ifelse(u.vec[2] <= p.accept2, 1, 0) 
     jump.vec3[i-1] <- ifelse(u.vec[3] <= p.accept3, 1, 0)
    }  

   ##Look at chains
     par(mfrow=c(3,1))
     plot(seq(1:nsim), phi.mat1[,1], type="l", ylab=expression(phi[1]), col=4, ylim=c(1,2))
        lines(seq(1:nsim), phi.mat2[,1], col=2)
        lines(seq(1:nsim), phi.mat3[,1], col=3)
     plot(seq(1:nsim), phi.mat1[,2], type="l", ylab=expression(phi[2]), col=4, ylim=c(-6,-3))
        lines(seq(1:nsim), phi.mat2[,2], col=2)
        lines(seq(1:nsim), phi.mat3[,2], col=3)
     plot(seq(1:nsim), phi.mat1[,3], type="l", ylab=expression(phi[3]), col=4, ylim=c(0,5))
        lines(seq(1:nsim), phi.mat2[,3], col=2)
        lines(seq(1:nsim), phi.mat3[,3], col=3)

   ##Acceptance rates?
     mean(jump.vec1)  
     mean(jump.vec2)
     mean(jump.vec3)





   ## Use coda to help 

    library(coda)
    help(package="coda") #for list of functions
    phi.post1 <- mcmc(phi.mat1)
    phi.post2 <- mcmc(phi.mat2)
    phi.post3 <- mcmc(phi.mat3)
    phi.mcmc <- mcmc.list(list(phi.post1, phi.post2, phi.post3))
    summary(phi.mcmc)
    plot(phi.mcmc) 
   
    #Look at some other convergence diagnostics

    effectiveSize(phi.mcmc)
    gelman.diag(phi.mcmc)  #Can't do b/c need at least 2 chains
    gelman.plot(phi.mcmc)
    geweke.diag(phi.mcmc)
    heidel.diag(phi.mcmc)
    raftery.diag(phi.mcmc)
      


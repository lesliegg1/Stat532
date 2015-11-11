
#### Schools example in Chapter 5 ####
####  Normal Hierarchical Model   ####
## TO DO:  Make into an .Rnw file and bring ###
##  in ggplot2 for better graphics          ##

  ### The data we are going to use are estimates of effects (y_j) obtained from 
  ###  separate linear regression models, along with the standard errors of 
  ###  those effects (i.e. we will not be considering the the estimates as
  ###  constants, but rather assuming we know the error associated with them
  ###  The standard errors are considered known (sigma_j) - this is akin to the 
  ###  Chapter 5 results for normal hierarchical models where the sample means
  ###  are modeled rather than the individual observations using theory 
  ###  of sufficient statistics                                                        
  
 ### The data
  y.j <- c(28, 8, -3, 7, -1, 1, 18, 12)
  sig.j <- c(15, 10, 16, 11, 9, 11, 10, 18)
  data.mat <- cbind(y.j, sig.j)

 ### Model
 #  y.j ~ N(theta.j, sig.j^2)
 #  theta.j ~ N(mu, tau2)
 #  mu propto 1
 #  tau propto 1  (NOT log(tau) propto 1) OR  tau ~ folded t


  ### Plot the "raw data" #####
  ## 
   interval.fun <- function(y.sd, m) {y.sd[1] + c(-1,1)*y.sd[2]*m}
   est.pm.1SD <- apply(data.mat, 1, interval.fun, m=1)  #2x8
   est.pm.2SD <- apply(data.mat, 1, interval.fun, m=2)  #
   est.pm.3SD <- apply(data.mat, 1, interval.fun, m=3)
   
   
   plot(seq(min(est.pm.3SD[1,])-0.5, max(est.pm.3SD[2,])+0.5,length=8), seq(0,9,length=8), type="n", xlab=" ", ylab=" ",yaxt="n") 
      mtext(c("A","B","C","D","E","F","G","H"), side=2, at=1:8, line=1,las=2)
      
      #Separate estimates
      points(y.j, 1:8, pch="|", cex=1.2, col="purple")
      segments(est.pm.3SD[1,], 1:8, est.pm.3SD[2,], 1:8, lwd=1, col=1)
      segments(est.pm.2SD[1,], 1:8, est.pm.2SD[2,], 1:8, lwd=3, col="purple")
      segments(est.pm.1SD[1,], 1:8, est.pm.1SD[2,], 1:8, lwd=5, col=4)

      ##Add lines for completely pooled estimate
      abline(v=7.7, lwd=2, col="orange")
      abline(v=(7.7 + c(-1,1)*4.1*2), lwd=2, lty=2, col="orange")

     # points(post.mns, 1:n.p, pch="|", cex=1, col=4)  #For later


  ##Discussion:  Probability statement implied by separate estimates = 
  ##  probability that the true effect in A is more than 28 is 0.5
  ## For the pooled model, probability that the true effect in A
  ##  is less than 7.7 is about 0.5.  Are either of these reasonable?


 #### Let's find the approximate posterior distribution  ########################
 #### We will first use the results that were derived in Chapter 5 to          ##
 #### obtain draws from the joint posterior distribution using factorization   ##
 ####  p(theta, mu, tau|y)=p(theta|mu,tau,y)p(mu|tau,y)p(tau|y)                ##
 ################################################################################

 ### Overall STRATEGY:
 # (1) Obtain samples from p(tau|y) using the following result: 
 #       p(tau|y) = p(mu,tau|y)/p(mu|tau,y) #propto [p(y|mu,tau)p(tau)]/p(mu|tau,y)  
 #       Sample using grid approximation over tau
 # (2) Simulate from p(mu|tau,y) using normal dist'n and values of tau from (1)
 # (3) Simulate from p(theta|mu,tau,y) using normal distribution and value of 
 #       tau from (1) and mu from (2)

## PRIOR for tau  ################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
##Before we begin, we must specify a prior distribution for tau 

## Think about tau within the context of the problem:
## Here, tau=PRIOR SD
## The effects are measured as point on the SAT testt (200 - 800) with 
## an avg 500, gives realistic upper limit on tau of about 100 so that 
## avg +/- 3SD = (200, 800).
## We will compare the p(tau) \propto 1  (as used first in text)
## and folded-t as suggested in Gelman et al. (2006)
## Could also compare p(tau2) = InverseGamma(1,1) and InverseGamma(0.001,0.001)
##   (See Figure 5.9)

 # install.packages("LearnBayes") #make sure installed before sourcing "folded_t_functions.R"
 setwd("~/Documents/Stat532/Rtutorials/8 schools")
 source("folded_t_functions.R")  

  ##Plot the density of the folded t
 	curve(d.tfold(x, df=2), xlim=c(0,5), ylim=c(0,0.85))  #plot for diff dfs
 	 for (i in 1:20) { curve(d.tfold(x, df=i), xlim=c(0,5), col=i, add=TRUE)} 

   #function to go with pg 520 in Gelman (2006) - Scaled 
   #plot it for diff s - NOTE: not standardized! 
    curve(un.tfold(x, df=1, s=1), xlim=c(0,50), ylim=c(0,1))  
    for (i in 1:20) { curve(un.tfold(x, df=1, s=i), xlim=c(0,50), col=i, add=TRUE)} 
	
   #function to obtain random draws from a scaled, non-central folded t
     hist(r.tfold(10000, df=4, mn=0, s=1), nclass=100, col=gray(.5), freq=F, xlim=c(0,5), ylim=c(0,0.85))
   	 curve(d.tfold(x, df=4), xlim=c(0,5), ylim=c(0,0.85), add=TRUE, 
               col=2, lwd=2)

     curve(un.tfold(x, df=1, s=5), xlim=c(0,50), ylim=c(0,1))  #Use for prior
     
      

          
   ## Function to obtain value proportional to log(p(tau|y)) to
   ## use to sample from the marginal posterior of tau
   ## We will use Equation 5.21 on page 117 and 
   ## a folded-t prior on 1 df and s=5

   lp.tauGy.lik <- function(tau, y.j, sig.j) {  #log scale
         V.mu.inv <- sum(1/(sig.j^2 + tau^2))
         mu.hat <- (sum((1/(sig.j^2 + tau^2))*y.j))/V.mu.inv
         part1 <- -0.5*log(V.mu.inv)
         part2 <- sum(-0.5*log(sig.j^2 + tau^2) - 
                       0.5*((y.j - mu.hat)^2)/(sig.j^2 + tau^2))
         out <- part1 + part2
         return(out)
       }
       
   l.prior.t1 <- function(tau) { log(un.tfold(tau,df=1,s=5)) }          

          
     ## Function to combine lik and prior pieces and exponentiate
     lp.tauGy.t1 <- function(tau, y.j, sig.j) {lp.tauGy.lik(tau, y.j=y.j, sig.j=sig.j) + l.prior.t1(tau)}
     p.tauGy.t1 <- function(tau) { exp(lp.tauGy.lik(tau, y.j=y.j, sig.j=sig.j) + l.prior.t1(tau)) }

 
    ## Function to obtain value proportional to log(p(tau|y)) = marg post of tau
    ## Using improper uniform prior on tau
    l.prior.unif <- function(tau) {log(1)}

     ## Simple function to exponentiate results of above function
     lp.tauGy.unif <- function(tau,y.j, sig.j) { lp.tauGy.lik(tau, y.j=y.j, sig.j=sig.j) + l.prior.unif(tau)}
     p.tauGy.unif <- function(tau) { exp(lp.tauGy.lik(tau, y.j=y.j, sig.j=sig.j) + l.prior.unif(tau))}
       
       
   ### Let's look at the prior vs. marginal posterior for tau       
   ### Plot p(tau|y) and p(tau) for the folded t and 
   ### plot p(tau|y) for the uniform
      grid.tau <- seq(0.001,50, length=5000) #play around with upper limit 
                                              #to see how things change
      width <- grid.tau[2] - grid.tau[1]
      grid.tau.g.y<- apply(cbind(grid.tau), 1, p.tauGy.t1)  #for folded-t prior
      norm.grid.tau <- grid.tau.g.y/sum(grid.tau.g.y*width)

      grid.tau.g.y2<- apply(cbind(grid.tau), 1, p.tauGy.unif) #unif prior on tau
      norm.grid.tau2 <- grid.tau.g.y2/sum(grid.tau.g.y2*width)

      #grid.tau.prior <- apply(cbind(grid.tau), 1, d.tfold, df=1, ncp=0)
      grid.tau.prior <- apply(cbind(grid.tau), 1, un.tfold, df=1, s=5)
       sum(grid.tau.prior*width)  
       norm.tau.prior <- grid.tau.prior/sum(grid.tau.prior*width) #grid version
   
      plot(grid.tau, norm.grid.tau, type="n", xlim=c(0,50), ylim=c(0,0.22), xlab="tau",
            ylab="density")
       lines(grid.tau, norm.tau.prior, col=2, lwd=2, lty=1) #folded-t prior
       lines(grid.tau, norm.grid.tau, col=3, lwd=2)  #marg. post with foldedt    
       lines(grid.tau, norm.grid.tau2, col=4, lwd=3, lty=1)
       legend(15,0.16, legend=c("folded-t prior (df=1, s=5)","p(tau|y) (p(tau) = folded-t)", 
               "p(tau|y) (p(tau)= unif)"),lwd=c(2,2,2), lty=c(1,1,1), col=c(2,3,4), bty="n")                  
          ##Compare to Figure 5.5 in text
          ## NOTE: Values near zero are most plausible, zero is mode
          ## Values of tau larger than 10 are half as likely as near tau=0
          ## Be sure to look at and understand Figures 5.6 and 5.7
     
      ###  Question:  What might plausible values of tau be?       
      ###  We already know upper limit on tau is about 100 from the context
      ###  of the problem. Also is helpful to stare at equations 5.17 and 
      ###  5.20 to help think about what values of tau are reasonable. 
      ###  What is tau the SD of?  How different could true SAT-coaching
      ###  effects at the school level really be?  Maybe 30 points?
      ###  this means that tau is probably not greater than 5? (15 +/- 3*5)
      
      ###  NOTE: For a classical analysis we could get an estimate of tau2 = MSB - MSW,
      ###  which can be negative!
 
 
    ### Now, let's think about where we can get into trouble with improper posteriors
    ### Since p(mu|tau,y) and p(theta|mu,tau,y) are both proper distributions (normals),
    ### if there is going to be a problem, it is going to be with p(tau|y)
    ### i.e. p(theta,mu,tau|y) has finite integral IFF p(tau|y) has finite integral
    
    ### What about p(mu, log(tau)) propto 1  (equivalent to p(mu,tau) propto 1/tau)?
     l.prior.logunif <- function(tau) {log(1/tau)}

     ## Simple function to exponentiate results of above function
     p.tauGy.logunif <- function(tau) { exp(lp.tauGy.lik(tau, y.j=y.j, sig.j=sig.j) + l.prior.logunif(tau))}
      grid.tauGy.logunif <- apply(cbind(grid.tau), 1, p.tauGy.logunif) #unif prior on tau
      norm.grid.logunif <- grid.tauGy.logunif/sum(grid.tauGy.logunif*width)
       
      ##ZOOM in on lower values of tau
      plot(grid.tau, norm.grid.logunif, type="n", xlim=c(0,4), ylim=c(0,0.4), xlab="tau",
            ylab="density")
       lines(grid.tau, norm.grid.logunif, col=5, lwd=3, lty=1) #unif on log(tau)
       lines(grid.tau, norm.tau.prior, col=2, lwd=2, lty=1) #folded-t prior
       lines(grid.tau, norm.grid.tau, col=3, lwd=2)  #marg. post with foldedt    
       lines(grid.tau, norm.grid.tau2, col=4, lwd=3, lty=1) #unif on tau
       legend(1.75,0.4, legend=c("folded-t prior (df=1, s=5)","p(tau|y) (p(tau) = folded-t)", 
               "p(tau|y) (p(tau)= unif)", "p(tau|y)  (p(log(tau))= unif)" ),lwd=c(3,3,3,3),
                lty=c(1,1,1,1), col=c(2,3,4,5), bty="n")
 
       ###Could show that p(tau|y) for unif on log(tau) is not proper (See Exercise 5.10)
       ### Math: Everything in 5.21 (ignoring p(tau)) approaches a nonzero constant limit
       ### as tau goes to zero (see the p(tau|y) with p(tau) propto 1 in above plot)
       ### [i.e. the marginal likelihood p(y|tau) approaches finite non-zero value as tau->0]
       ### 
       ### Therefore, the behavior of posterior near zero depends on the prior.
       ### Since p(log(tau)) = p(tau) propto 1/tau, and 1/tau is not integrable for any
       ### small interval including tau=0
       ### SO, take-home message is that if p(tau) is integrable near zero, then the
       ### posterior should be integrable near zero as well
       ### What about the other extreme?  What is happening as tau goes to infinity?
       ### The exponential term in 5.21 is less than or equal to 1, and with some algebra
       ### one can find the upper bound for p(tau|y) for tau large is 
       ### p(tau)*sqrt(1/J)*(tau^(-(J-1))) where J is the number of groups
       ### Therefore, if p(tau) propto 1, it IS integrable as long as there are more than
       ### 2 groups (J>2), but Gelman et al. worry that with small J (4 or 5) it overestimates
       ### tau, resulting in less than optimal amount of pooling. For 1 or 2 groups is
       ### improper, but actually results in tau=infinity and no pooling.
 
       ### FOR MORE ABOUT PRIORS on tau, read section 5.7 and Gelman et al. (2006)
       ### Particularly pay attention to the commonly used priors
       ### Easy to think that if they work at lowest data level, they will work
       ### for variance params at next level, but that is not the case
       ### Page 129 quote: "In a hierarchical model the data can never rule out a 
       ### group-level variance of zero, and so the prior distribution cannot put an
       ### infinite mass in this area".  Page 130 in BDA3 goes through different priors
       ### for the Schools example (TO BE ADDED IN FUTURE)
 

 ###### NOW, we will move away from just focusing on tau, to getting draws
 ###### from the joint posterior distribution p(theta,mu,tau|y)
 ###### To do this, we will use the Metropolis algorithm instead of grid
 ###### approximation to get draws of tau from p(tau|y)
 
 
   ##STEP 1 - draw tau from p(tau|y) for a particular p(tau)### 
   ## We will sample from the two different p(tau|y)'s from  ##
   ## the folded-t and uniform priors together in same loop ##  
   
    nsamp <- 5000
    tau.vec <- numeric(nsamp)  #save draws u
    tau.vec2 <- numeric(nsamp) #save draws using uniform prior

    tau.vec[1] <- 10  #starting value
    tau.vec2[1] <- 10
     
     
    for (t in 2:nsamp) {
     tau.cur <- tau.vec[t-1] ##FOLDED-t sampling
     #tau.cand <- r.tfold(1, df=2, mn=tau.cur) #jumping distn centered on current  
     tau.cand <- r.tfold(1, df=2, s=4, mn=0)     
     l.p.cand <- lp.tauGy.t1(tau.cand, y.j, sig.j)
     l.p.cur <- lp.tauGy.t1(tau.cur, y.j, sig.j)
     l.j.cand <- log(un.tfold(tau.cand, df=2, s=4))
     l.j.cur <- log(un.tfold(tau.cur, df=2, s=4))
     #l.j.cand <- log(d.tfold(tau.cand,df=2,ncp=tau.cur))  #not correct
     #l.j.cur <- log(d.tfold(tau.cur,df=2,ncp=tau.cand))   #not correct
     log.rat <- l.p.cand - l.p.cur + l.j.cur - l.j.cand
     r <- min(1, exp(log.rat))
     ifelse(runif(1)<r, tau.vec[t] <- tau.cand, tau.vec[t] <- tau.cur)

     tau.cur2 <- tau.vec2[t-1]  #UNIFORM sampling
     #tau.cand2 <- r.tfold(1, df=2, mn=tau.cur2) #jumping distn centered on curr
     tau.cand2 <- r.tfold(1, df=2, s=4, mn=0)
     l.p.cand2 <- lp.tauGy.unif(tau.cand2, y.j, sig.j)
     l.p.cur2 <- lp.tauGy.unif(tau.cur2, y.j, sig.j) 
     l.j.cand2 <- log(un.tfold(tau.cand2, df=2, s=4))
     l.j.cur2 <- log(un.tfold(tau.cur2, df=2, s=4))  
     #l.j.cand2 <- log(d.tfold(tau.cand2,df=2,ncp=tau.cur2)) #not correct
     #l.j.cur2 <- log(d.tfold(tau.cur2,df=2,ncp=tau.cand2)) #not correct
     log.rat2 <- l.p.cand2 - l.p.cur2 + l.j.cur2 - l.j.cand2 
     r2 <- min(1, exp(log.rat2)) 
     ifelse(runif(1)<r2, tau.vec2[t] <- tau.cand2, tau.vec2[t] <- tau.cur2)
    }

   #Plot histograms of posterior draws of tau from the two priors
   #Overlay curves obtained in first part of the code 
   par(mfrow=c(2,1))
    hist(tau.vec, nclass=20, col=gray(.5), freq=F, main="Folded-t prior", 
            xlab="tau|y", xlim=c(0,50), ylim=c(0,0.2))
      lines(grid.tau, norm.grid.tau, col=3, lwd=2)  #marg. post with foldedt 
    hist(tau.vec2, nclass=40, col=gray(.5), freq=F, main="Unif prior", 
             xlab="tau|y", xlim=c(0,50), ylim=c(0,0.2))
      lines(grid.tau, norm.grid.tau2, col=4, lwd=3, lty=1) #unif on tau
             

   ##Look at sample path plots to see where the draws came from
   dev.new()
   plot(1:nsamp, tau.vec2, type="n", main="Sample path plots", ylim=c(0,60))  
     lines(1:nsamp, tau.vec2, col=3)
     lines(1:nsamp, tau.vec, col=2) #looks more reasonable to me?
     legend(100,58, legend=c("folded-t prior (df=1, s=5)","Uniform prior"), 
                 lwd=c(2,2), lty=c(1,1), col=c(2,3), bty="n")

   plot(1:nsamp, log(tau.vec2), type="n", main="Sample path plots", ylim=c(-5,5))  
     lines(1:nsamp, log(tau.vec2), col=3)
     lines(1:nsamp, log(tau.vec), col=2) #looks more reasonable to me?
     legend(100,5.2, legend=c("folded-t prior (df=1, s=5)","Uniform prior"), 
                 lwd=c(2,2), lty=c(1,1), col=c(2,3), bty="n")


   


  ###STEP 2: Obtain samples from p(mu|tau,y) = Normal(mu.hat,V.mu)
  #Function to easily get a draw from mu|tau,y  
  #conditional posterior of mu in Normal (Equation 5.20 in text)
  
   r.mu.g.tau.y <- function(tau, y.j, sig.j) {
     V.mu <- 1/(sum(1/(sig.j^2 + tau^2)))
     mu.hat <- (sum((1/(sig.j^2 + tau^2))*y.j))*V.mu
     draw <- rnorm(1, mu.hat, sd=sqrt(V.mu))
     return(draw)
   }
    #r.mu.g.tau.y(1, y.j=y.j, sig.j=sig.j)

  #Actually get the draws under the two different priors for tau
   mu.vec <- apply(cbind(tau.vec), 1, r.mu.g.tau.y, y.j=y.j, sig.j=sig.j)
   mu.vec2 <- apply(cbind(tau.vec2), 1, r.mu.g.tau.y, y.j=y.j, sig.j=sig.j)

   #Take a look at results for mu|y under the two different priors
    par(mfrow=c(2,1))
    hist(mu.vec, nclass=50, col=gray(.7), freq=F, main="Folded-t prior", 
            xlab="mu|y", xlim=c(-10,30))
    hist(mu.vec2, nclass=50, col=gray(.7), freq=F, main="Improper Unif. prior", 
              xlab="mu|y", xlim=c(-10,30))

   #Look at summary measures to help compare the marginal distributions of mu|y
   # under the two priors
    summary(mu.vec)
    summary(mu.vec2)
    quantile(mu.vec, c(0.025, .975))   
    quantile(mu.vec2, c(0.025, .975))  
    mean(mu.vec>28) #0   
    mean(mu.vec2>28) #0.003 

   ###STEP 3: Obtain samples from (theta_{j}|mu,tau,y)=N(theta.j.hat, Vj)
   #Function to get a draw
   r.thetaj.g.mu.tau.y <- function(mu.tau, y.j, sig.j) {
     V.j <- 1/((1/(sig.j^2)) + (1/(mu.tau[2]^2)))
     theta.hat.j <- ((y.j/(sig.j^2)) + (mu.tau[1]/(mu.tau[2]^2)))*V.j 
     theta.draws <- rnorm(8, theta.hat.j, sqrt(V.j))
     return(theta.draws)
   }

   #Apply function to get vector of draws for each theta.j
   theta.mat <- t(apply(cbind(mu.vec,tau.vec), 1, r.thetaj.g.mu.tau.y, y.j=y.j, 
                    sig.j=sig.j))
   theta.mat2 <- t(apply(cbind(mu.vec2,tau.vec2), 1, r.thetaj.g.mu.tau.y, 
                      y.j=y.j, sig.j=sig.j))

    
   
  ###################################################################
  ### Now, let's look at school level results and further summarize # 
   
   burn.in <- 500  #Because based on Metropolis draws for tau

  #Summarize the draws, make plots in text
   
   school <- c("A","B","C","D","E","F","G","H")
   dev.new()
   par(mfrow=c(4,2))
    for (j in 1:8) {
       hist(theta.mat[-(1:burn.in),j], col=gray(.5), nclass=100, main=paste(school[j]), 
             xlim=c(-20,40), freq=FALSE)
     }

   dev.new()
   par(mfrow=c(4,2))
    for (j in 1:8) {
       hist(theta.mat2[-(1:burn.in),j], col=gray(.5), nclass=100, main=paste(school[j]), 
             xlim=c(-20,40), freq=FALSE)
     }
     
     
   #### Replot so can directly compare the same school for different priors  
   dev.new()
    par(mfrow=c(2,4)) #Compare results for schools 1 to 4
    for (j in 1:4){
      hist(theta.mat[-(1:burn.in),j], col=gray(.5), nclass=100, main=paste(school[j]), 
             xlim=c(-20,40), freq=FALSE, xlab=" ") } #Folded-t
    for (j in 1:4){
      hist(theta.mat2[-(1:burn.in),j], col=gray(.5), nclass=100, main=paste(school[j]), 
             xlim=c(-20,40), freq=FALSE, xlab=" ")} #Unif 
        
   dev.new()
    par(mfrow=c(2,4))  #Compare results for schools 5 to 8
    for (j in 5:8){
      hist(theta.mat[-(1:burn.in),j], col=gray(.5), nclass=100, main=paste(school[j]), 
             xlim=c(-20,40), freq=FALSE, xlab=" ")} #Folded-t  
    for (j in 5:8){
      hist(theta.mat2[-(1:burn.in),j], col=gray(.5), nclass=100, main=paste(school[j]), 
             xlim=c(-20,40), freq=FALSE, xlab=" ")} #Unif

    
   ### Find the posterior probability that
   ###   school A effect is as large as 28 points:
       mean(theta.mat[,1]>28) #approx 0.01       
       mean(theta.mat2[,1]>28) #approx 0.05

    ## Figure 5.8 (b) The posterior dist'n for largest effect max(theta_{j})
    ## Note - these are NOT posterior predictive distributions
      max.effect <- apply(theta.mat[-(1:burn.in),], 1, max)
      max.effect2 <- apply(theta.mat2[-(1:burn.in),], 1, max)

      dev.new()
      par(mfrow=c(1,2))
      hist(max.effect, col=gray(.5), nclass=100, main="max effect- folded-t", xlim=c(-20,70))
       abline(v=28, col="orange", lwd=2, freq=FALSE)
      hist(max.effect2, col=gray(.5), nclass=100, main="max effect- unif", xlim=c(-20,70))
        abline(v=28, col="orange", lwd=2, freq=FALSE)

      mean(max.effect > 28.4) #approx 0.016
      mean(max.effect2 > 28.4) #approx 0.076


   ### Can we estimate the Pr(theta_A > theta_C|y)?  Posterior probability that school A is
   ##   more effective than school C?
      mean(theta.mat[,1] > theta.mat[,3]) # approx 0.62
    
    
  #### Plot the separate effects vs. the theta.j's  ###
  ## TO DO: update these with varying line width caterpillar plots
  
    mu.dummy <- seq(-20,50, length=500)
    dev.new()
    par(mfrow=c(1,1))
    plot(mu.dummy, seq(0,9, length=500), type="n", xlim=c(-50,70), ylab="", 
          xlab="School `effect'", yaxt="n", main="Folded-t prior 95% intervals")
      mtext(rev(school), side=2, at=c(1,2,3,4,5,6,7,8), adj=0, las=2, line=1)
      abline(v=mean(mu.vec[-(1:burn.in)]))
      segments(-0.3, 0, 16.0, 0, lwd=3, col="orange")  #complete pooling intervals
      points(7.9, 0, pch=18, col=1, cex=1.25)
      text(-30, 0, "Pooled")
      abline(h=seq(0.5,7.5,by=1), lty=2, col=gray(.3))
      for (i in 1:8) {
      	 ci <- y.j[9-i] + c(-1,1)*2*sig.j[9-i]
      	 pi <- quantile(theta.mat[-(1:burn.in),9-i],c(0.025, 0.975))
      	 segments(ci[1], i-0.15, ci[2], i-0.15, lwd=3, col=2)
      	 points(y.j[9-i], i-0.15, pch=17, cex=1.25, col=1)
      	 segments(pi[1], i+0.15, pi[2], i+0.15, lwd=3, col=4)
      	 points(mean(theta.mat[-(1:burn.in),9-i]), i+0.15, cex=1.25, pch=16, col=1)
      	 }
     
     dev.new() 	 
     plot(mu.dummy, seq(0,9, length=500), type="n", xlim=c(-50,70), ylab="", 
        xlab="School `effect'", yaxt="n", main="Uniform on SD prior 95% intervals")
      mtext(rev(school), side=2, at=c(1,2,3,4,5,6,7,8), adj=0, las=2, line=1)
      abline(v=mean(mu.vec2[-(1:burn.in)]))
      segments(-0.3, 0, 16.0, 0, lwd=3, col="orange")  #complete pooling intervals
      points(7.9, 0, pch=18, col=1, cex=1.25)
      text(-30, 0, "Pooled")
      abline(h=seq(0.5,7.5,by=1), lty=2, col=gray(.3))
      for (i in 1:8) {
      	 ci <- y.j[9-i] + c(-1,1)*2*sig.j[9-i]
      	 pi <- quantile(theta.mat2[-(1:burn.in),9-i],c(0.025, 0.975))
      	 segments(ci[1], i-0.15, ci[2], i-0.15, lwd=3, col=1)
      	 points(y.j[9-i], i-0.15, pch=17, cex=1.25, col=1)
      	 segments(pi[1], i+0.15, pi[2], i+0.15, lwd=3, col="magenta")
      	 points(mean(theta.mat2[-(1:burn.in),9-i]), i+0.15, cex=1.25, pch=16, col=1)
      	 }
      	 
      	 
    dev.new() ##ALL On same plot	 
     plot(mu.dummy, seq(0,9, length=500), type="n", xlim=c(-50,70), ylab="", 
        xlab="School `effect'", yaxt="n", main="Compare Priors")
      mtext(rev(school), side=2, at=c(1,2,3,4,5,6,7,8), adj=0, las=2, line=1)
      abline(v=mean(mu.vec2[-(1:burn.in)]))
      segments(-0.3, 0, 16.0, 0, lwd=3, col="orange")  #complete pooling intervals
      points(7.9, 0, pch=18, col=1, cex=1.25)
      text(-30, 0, "Pooled")
      abline(h=seq(0.5,7.5,by=1), lty=2, col=gray(.3))
      for (i in 1:8) {
      	 ci <- y.j[9-i] + c(-1,1)*2*sig.j[9-i]
      	 segments(ci[1], i-0.2, ci[2], i-0.2, lwd=3, col=1)
      	 points(y.j[9-i], i-0.2, pch=17, cex=1.25, col=1)
      	 
      	 pi.unif <- quantile(theta.mat2[-(1:burn.in),9-i],c(0.025, 0.975))
      	 segments(pi.unif[1], i+0, pi.unif[2], i+0, lwd=3, col="magenta")
      	 points(mean(theta.mat2[-(1:burn.in),9-i]), i+0, cex=1.25, pch=16, col=1)
      	 
      	 pi.t <- quantile(theta.mat[-(1:burn.in),9-i],c(0.025, 0.975))
      	 segments(pi.t[1], i+0.2, pi.t[2], i+0.2, lwd=3, col=4)
      	 points(mean(theta.mat[-(1:burn.in),9-i]), i+0.2, cex=1.25, pch=16, col=1)
      	 }
      	legend(40,9.5, legend=c("Folded-t","Unif","Separate"), col=c(4,"magenta",1),
      	         lwd=c(2,2,2), lty=c(1,1,1), bty="n")
 
 

######## EXTRA - different methods of obtaining draws #############
## Since we have covered Gibbs sampling and used JAGS and Stan, ###
## Let's look at other ways we can obtain draws from the posterior#           
## In BDA3 (came out after I programmed the above stuff) provides #
## a lot of code for the Schools example.  SEE PAGES 594 - 604!!

 ### Gibbs sampling version 
 ###  complexity to the one above?
 ### The following code is on page 596 in the new text
 
 
J <- 8
y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)

theta_cc_update <- function() {
	theta_hat <- ((mu/tau^2) + (y/sigma^2))/((1/tau^2)+(1/sigma^2))
	V_theta <- 1/((1/tau^2)+(1/sigma^2))
	rnorm(J,theta_hat, sqrt(V_theta)) }

mu_cc_update <- function() {
	rnorm(1, mean(theta), tau/sqrt(J)) }

tau_cc_update <- function(){
	sqrt(sum((theta-mu)^2)/rchisq(1,J-1)) }

chains <- 4
iter <- 2000
sims <- array(NA, c(iter,chains,J+2))
dimnames(sims) <- list(NULL, NULL, c(paste("theta[", 1:8, "]", sep=""), "mu", "tau"))
for (m in 1:chains){
	mu <- rnorm(1, mean(y), sd(y))
	tau <- runif(1, 0, sd(y))
	for (t in 1:iter){
		theta <- theta_cc_update()
		mu <- mu_cc_update()
		tau <- tau_cc_update()
		sims[t,m,] <- c(theta,mu,tau)
	}
}


##See BDA3 for writing your own HMC code, or 
## we can use Stan

##### Use Stan    ######
library(rstan)

J <- 8
y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)
schools.fit <- stan(file="8schools.stan", data=c("J", "y","sigma"), iter=1000, chains=4)
print(schools.fit)
plot(schools.fit) 

schools.fit1 <- stan(fit=schools.fit, data=c("J", "y","sigma"), iter=2000, chains=4)

schools.sim <- extract(schools.fit1, permuted=TRUE) #combines chains and permutes values, after warmup
names(schools.sim) 

dim(schools.sim$mu)  
dim(schools.sim$theta)

hist(schools.sim$tau, nclass=100, col=gray(0.8), freq=FALSE)
hist(schools.sim$mu, nclass=100, col=gray(0.8), freq=FALSE)
mean(schools.sim$theta[,1] > schools.sim$theta[,3])

     
      
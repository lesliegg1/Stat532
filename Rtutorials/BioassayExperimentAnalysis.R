
### This code follows Section 3.7 in Gelman et al.: Analysis of Bioassay Expt.

### Description of the problem (from pg 74):
## In the development of drugs and other chemical compounds, acute toxicity
## test or bioassay experiments are commonly performed on animals. Such
## experiments proceed by administering various does levels of the compound
## to batches of animals.  The animals' responses are typically characterized by
## a dichotomous outcome: for example, alive or dead, tumor or no tumor.  An
## experiment of this kind gives rise to data of the form
## (xi, ni, yi) where xi=ith of k dose levels
##                    ni=number of animals given the ith dose
##                    yi=number of animals out of ni with positive outcome


## Data from Racine et al. (1986) (Table 3.1 in Gelman et al)

log.dose <- c(-0.86, -0.30, -0.05, 0.73)
n.dose <- c(5,5,5,5)
y.dose <- c(0,1,3,5)

data.frame(log.dose,n.dose,y.dose)


#1. Think about exchangeability of the 5 animals within each dose group:
#  Independent with equal probabilities of death? Not too bad of assump? If yes, 
#  then Binomial dist'n useful:
# Treat outcomes of the 4 different doses as independent after accounting for 
#   dose by allowing different probs for each dose? Not too bad of assump?
# We will proceed with a "yes" for both of these and assume yi|pi ~ Bin(ni,pi)


 #1.  Specify a prior distribution for the pi's:
 #    Gelman suggests the simplest way to proceed is assuming the pi's are exchangeable in their
 #    prior dist'n and using some non-informative density: p(p1,p2,p3,p4) prop.to 1
 #    What posterior distributions would the pi's then have?

 #   What is the flaw of the above prior distribution?  We KNOW the dose for each group
 #   and DO NOT expect the probabilities of death to be exchangeable since they should
 #   depend on dose!  So, let's model the changes in pi we expect due to dose using a 
 #   logit(pi) = alpha + beta*dosei  (Why wouldn't we just use pi = a + b*dosei?)
 #   Now we just have a logistic regression model!
# NOTE: dose is already log transformed


   inv.logit <- function(x) { exp(x)/(1 + exp(x)) }

   LR.loglik.fun <- function(ab.vec,x.vec,n.vec,y.vec) {
       ll.out <- sum(y.vec*log(inv.logit(ab.vec[1] + ab.vec[2]*x.vec)) + 
                    (n.vec-y.vec)*log(1-inv.logit(ab.vec[1] + ab.vec[2]*x.vec)))
       return(ll.out)
     }


    LR.lik.fun <- function(ab.vec,x.vec,n.vec,y.vec) {
       out <- prod( ((inv.logit(ab.vec[1] + ab.vec[2]*x.vec))^y.vec)*
                        (1-inv.logit(ab.vec[1] + ab.vec[2]*x.vec))^(n.vec-y.vec) )
       return(out)       
     }


 
 ## Now we need to think about a prior distribution for (alpha, beta) instead of p.  To proceed with
 ##  a simple analysis we will use p(alpha,beta) propto 1  (independent and locally uniform priors)


     #b/c of the prior we used, we don't need to define another function for lik.x.prior


  ## Now, we want to compute the joint posterior distribution over a grid of points. It's nice to get
  ##  a crude idea of where alpha and beta will lie.  To do this you can perform ordinary regression on
  ## the proportions, or go ahead and run logistic regression to find the MLEs.

    ### Let's go ahead and perform the standard non-Bayesian analysis:
        yes.no <- cbind(y.dose, n.dose-y.dose)
        glm.out <- glm(yes.no ~ log.dose, family=binomial)
        summary(glm.out)
            #           Estimate Std. Error z value Pr(>|z|)
          #(Intercept)   0.8466     1.0191   0.831    0.406
          #log.dose      7.7488     4.8728   1.590    0.112


   ### Define the grid over which to compute the posterior density
         a.grid <- seq(-5,8, length=200)
         b.grid <- seq(-10,35, length=200)
         ab.grid <- expand.grid(a.grid,b.grid)
         lik.x.p.vec1 <- apply(ab.grid, 1, LR.loglik.fun, x.vec=log.dose, n.vec=n.dose, y.vec=y.dose)
         lik.x.p.vec <- exp(lik.x.p.vec1 - max(lik.x.p.vec1))
         lik.x.p.mat <- matrix(lik.x.p.vec, nrow=length(a.grid), ncol=length(b.grid))
         norm.const <- sum((a.grid[2]-a.grid[1])*(b.grid[2]-b.grid[1])*lik.x.p.vec)
         post.vec <- lik.x.p.vec/norm.const        
         post.mat <- lik.x.p.mat/norm.const

         contour(a.grid, b.grid, post.mat, xlab=expression(alpha), ylab=expression(beta),
                     levels=seq(0.05, 0.95, by=0.05)*max(post.mat), drawlabels=FALSE)
         abline(v=0.8466, col=3) #MLE's
         abline(h=7.7488, col=3)


   ### Now obtain draws from the posterior distribution:

     ###Method 1: Use the same grid approximation sampling method we did for univariate problem     ## 
         grid.rows <- 1:length(post.vec)
         samp.rows <- sample(grid.rows, size=10000, prob=post.vec, replace=TRUE)
        
         post.samps <- ab.grid[samp.rows,]

         #Now let's make our posterior draws more continuous
         a.interval <- a.grid[2]-a.grid[1]
         b.interval <- b.grid[2]-b.grid[1]
         post.samps.cont <- post.samps + c(runif(1, -a.interval/2, a.interval/2), runif(1, -b.interval/2, b.interval/2))
         dim(post.samps)

         contour(a.grid, b.grid, post.mat, xlab=expression(alpha), ylab=expression(beta),
                     levels=seq(0.05, 0.95, by=0.05)*max(post.mat), drawlabels=FALSE)
            points(post.samps.cont[,1], post.samps.cont[,2], col=3, cex=0.5)



    ### Method 2:  Using the conditional posterior for mu and marginal posterior for sig2 as 
     ##   Gelman et al. describe on page 92

        #1. Compute the marginal posterior distribution of alpha by numerically summing over beta
        #   in the discrete distribution computed on grid  (post.mat)

            a.marg.post <- apply(post.mat, 1, sum)
         
        #2.a. Draw alpha from the marginal posterior obtained in 1
            a.post.samps <- sample(a.grid, 10000, a.marg.post, replace=TRUE)
           
         #  b. Draw betas from the discrete conditional distribution p(beta|alpha,y)
  
             b.post.samps <- numeric(length(a.post.samps))
             for (i in 1:length(a.post.samps)){
                 b.post.samps[i] <- sample(b.grid, 1, prob=post.mat[a.grid==a.post.samps[i],])
                 }
             
            contour(a.grid, b.grid, post.mat, xlab=expression(alpha), ylab=expression(beta),
                     levels=seq(0.05, 0.95, by=0.05)*max(post.mat), drawlabels=FALSE)
                points(post.samps.cont[,1], post.samps.cont[,2], col=3, cex=0.5) #Method 1 points
                points(a.post.samps, b.post.samps, col=4, cex=0.5)  #Method 2 points
            

       ### ISSUES to consider when using a grid approximation:
         # -Size of grid: too small misses features of posterior distribution, too large can miss smaller features
         #     falling b/t grid points
         # -Beware of overflow and underflow operations - computers can only keep track of so many digits
         #   -**Compute log of posterior density, subtract off maximum, and then exponentiate
         #       as I did above 
         #


   ###########################################
   ### Posterior distribution of the LD50
     #  A common parameter of interest in bioassay studies is LD50 (Lethal Dose 50 = Dose at which 
     #  prob. of death is 50%)
     #  LD50:  E(Yi/ni)=logit^-1(alpha + beta*xi)=0.5 -->  xi= -alpha/beta
     #  What is the distribution of LD50?
     #  Given that we have posterior draws of alpha and beta, and LD50 is just a function of
     #  alpha and beta, it is "easy peasy lemon squeezy"(quote from my kindergartener) 
     #  to get p(LD50|y)
     
      LD50.draws <- -a.post.samps/b.post.samps
      min(LD50.draws)
      hist(LD50.draws, col="lightgray", main="Posterior for LD50", xlim=c(-1,1), nclass=300, freq=FALSE) 
     
     #For this problem, if beta<=0, LD50 is not a useful concept
     
     ## Practical inferences
     #  What is the probability that the drug is harmful (beta<0)?
        mean(b.post.samps>0) #1  Pr(beta>0) > 0.9999  (On one run I did get one posterior draw < 0)
     
        mean(LD50.draws>0) #= 0.0
        quantile(LD50.draws, c(0.025, 0.975)) ##How would you have obtained a CI or SE in likelihood
                                               # analysis?
        
     #Now, can compute the posterior of LD50 conditional on beta>0
     # If we have non-positive values for beta take only those positive values to summarize
     
        LD50.post.cond <- LD50.draws[b.post.samps>0] #draws conditional on beta>0
        hist(LD50.post.cond, col="lightgray", main="Posterior for LD50 given Beta>0", xlim=c(-1,1), 
                nclass=200, freq=FALSE, xlab="LD50|beta>0 (log g/ml)") 
        mean(LD50.post.cond)
        median(LD50.post.cond)
        quantile(LD50.post.cond, c(.05,.95)) 
        
        #Need to backtransform to get back to original scale of dose
        
        LD50.post.exp <- exp(LD50.post.cond)
        hist(LD50.post.exp, col="lightgray", main="LD50|beta>0: backtransformed", xlim=c(0,3), 
             nclass=200, freq=FALSE, xlab="LD50|beta>0 (g/ml)") 
        mean(LD50.post.exp)
        median(LD50.post.exp)
        
quantile(LD50.post.exp, c(.05,.95)) 
        
       #What is the posterior probability that LD50 is less than 0.75, which is close to the
       # second dose used in the study?
         mean(LD50.post.exp<0.75)  #0.0167
        
      
     
  

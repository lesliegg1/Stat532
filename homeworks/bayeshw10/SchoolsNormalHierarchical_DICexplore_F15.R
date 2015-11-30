 ###############################################################
 ### First load data and run MCMC code to get posterior samples
 ### From code in SchoolsNormalHierarchicalExample_F15.R  #####
 ### theta.mat (folded-t) and theta.mat2 (unif) for tau   #####
 ### Updated Fall 2015                                    #####
 ##############################################################

 ### Think about residuals
 
 ### theta.mat is 5000 x 8

  resid.fun <- function(theta.j, y.j){y.j - theta.j}

  resid.out <- apply(theta.mat, 1, resid.fun, y.j=y.j) #8x5000
  dev.new()
  par(mfrow=c(3,3), mar=c(4,4,2,2))
   for (s in 1:8){
     plot(theta.mat[s,], resid.out[,s], pch=16, cex=1.2, ylim=c(-20,35), xlim=c(-15,30))
    }
 
 ###Versus just plugging in the posterior mean of theta
  theta.mns <- apply(theta.mat,2,mean)
  plot(theta.mns, (y.j-theta.mns), pch=16, cex=1.2, ylim=c(-20,35), xlim=c(-15,30),
         main="Using posterior MEAN", col=2)

 ### Think about average wt'd sum of squared resids (approx deviance)
  avg.wt.sq.resid.fun <- function(theta.j, y.j, sig.j) {
     mean(((y.j-theta.j)^2)/(sig.j^2))
    }
  avg.sq.resid.pm <- avg.wt.sq.resid.fun(theta.mns, y.j=y.j, sig.j=sig.j) #at posterior means
  avg.sq.resid.pool <- avg.wt.sq.resid.fun(rep(mean(y.j),8), y.j=y.j, sig.j=sig.j) #at overall mean
  

  ##Do for EACH posterior draw to get posterior distribution

   post.avg.sq.resid <- apply(theta.mat, 1, avg.wt.sq.resid.fun, y.j=y.j, sig.j=sig.j)
   
   par(mfrow=c(1,1))
   hist(post.avg.sq.resid, nclass=100, col=gray(0.9), main="Posterior of Average Weighted Squared Resids",
          freq=FALSE)
     abline(v=avg.sq.resid.pm, lwd=2, col=2)  #Using posterior mean as plug-in
     abline(v=avg.sq.resid.pool, lwd=2, col=4)  #using overall mean as plug-in
     legend(1.5,3.0, legend=c("posterior means", "overall mean"), col=c(2,4), 
                  lwd=c(2,2), lty=c(1,1), bty="n")


 ## CHANGED IN _F15 version
 ## Now each y.rep is only used with its generating theta, rather than moving
 ## over all theta for a y.rep like we do for the observed data.
 ## BUT, can also do for y.rep's (for these 8 schools)
     n_draws <- length(theta.mat[,1])
     y.rep.mat <- matrix(NA, nrow=n_draws, ncol=8)
       for (j in 1:8){
          y.rep.mat[,j] <- rnorm(n_draws, mean=theta.mat[,j], sd=sig.j[j])
         }
   
    #Write function to calculate wt'd MSE for different draws of (theta, yrep)  
      T.yrep.mse.fun <- function(theta_yrep, sig.j) {
    	  avg.wt.sq.resid.fun(theta.j=theta_yrep[1:8], y.j=theta_yrep[9:16], sig.j=sig.j)
             }
             
     theta_yrep_mat <- cbind(theta.mat,y.rep.mat) #5000 x 16
     choose.draws <- sample(500:5000, size=1000) #choose draws randomly after 500      
     T.yrep.theta <- apply(theta_yrep_mat[choose.draws,], 1, T.yrep.mse.fun, sig.j=sig.j)
  
    par(mfrow=c(1,2))
     plot(post.avg.sq.resid[choose.draws], T.yrep.theta, pch=16, cex=1.2, xlim=c(0,5), 
            ylim=c(0,5), xlab="Posterior MSE (y_obs)", ylab="Posterior MSE (y_preds)")
         abline(a=0, b=1)
         abline(v=c(avg.sq.resid.pm,avg.sq.resid.pool), col="orange")
      #Why would it make sense for ypreds to have more larger values for this application?
    
     #Proportion of posterior predictive draws such that the T(yrep,theta) is greater
     # than or equal to T(yobs)?
      mean(T.yrep.theta >= post.avg.sq.resid[choose.draws]) #.707
      

  ### Or histograms of differences?
    dev.new()
     hist((post.avg.sq.resid[choose.draws]-T.yrep.theta), col=gray(0.8), nclass=100, xlim=c(-3,3),
           main="T(y,theta)-T(yrep,theta)", xlab="Difference")
     abline(v=0, col=2)
     
     avg.sq.resid.pm <- avg.wt.sq.resid.fun(theta.mns, y.j=y.j, sig.j=sig.j) #at posterior means
  avg.sq.resid.pool 
  

 ### More general: Compare deviance among models:

    Dev.fun <- function(theta.j, y.j, sig.j) {
          -2*sum(log(dnorm(y.j, theta.j, sig.j)))
          }

     Dev.fun(y.j, y.j, sig.j)


     #1. Calculate D.theta.hat
        burn.in <- 500
       #a. Find a theta.hat - use posterior means      
           post.mns <- apply(theta.mat[-(1:burn.in),],2,mean)
           post.mns2 <- apply(theta.mat2[-(1:burn.in),],2,mean) 

       #b. Plug in posterior means to deviance function
           D.theta.hat <- Dev.fun(post.mns, y.j, sig.j)     #58.35
           D.theta.hat2 <- Dev.fun(post.mns2, y.j, sig.j)   #57.40 ("matches" book)

     #2. Calculate D.avg.hat
       #a. Apply the deviance function to every set of posterior draws of theta
           Dev.draws <- apply(theta.mat[-(1:burn.in),],1,Dev.fun, y.j=y.j, sig.j=sig.j)
           Dev.draws2 <- apply(theta.mat2[-(1:burn.in),],1,Dev.fun, y.j=y.j, sig.j=sig.j)
           par(mfrow=c(2,1))
           hist(Dev.draws, col=gray(.5), nclass=100, main="Deviance (folded-t)", xlim=c(50,80), 
       			xlab=expression(D(y,theta)), freq=F, ylim=c(0,0.4))
                  abline(v=D.theta.hat, col=2, lwd=2)
                  text(56,0.4, expression(D[hat(theta)]))
           hist(Dev.draws2, col=gray(.5), nclass=100, main="Deviance (uniform)", xlim=c(50,80), 
                xlab=expression(D(y,theta)),freq=F, ylim=c(0,0.4))
                 abline(v=D.theta.hat2, col=2, lwd=2)
                 text(55,0.4, expression(D[hat(theta)]))


        #b. Average the deviance over all the draws 
            D.avg.hat <- mean(Dev.draws)     #60.227
            D.avg.hat2 <- mean(Dev.draws2)   #60.2916
            
            abline(v=D.avg.hat2, col="orange", lwd=2) #add to histogram for uniform
             text(65, 0.4, expression(bar[D(theta)]))
             
           
          
      #3. Calculate the effective number of parameters (pD) using the two
      #    different methods in the book (pg 181,182): pD1 and pD2
            pD1 <- D.avg.hat - D.theta.hat       #1.9
            pD1.2 <- D.avg.hat2 - D.theta.hat2   #2.78

             #pD2 is approximately 1/2 the posterior variance of the deviance
            pD2 <- 0.5*var(Dev.draws)  #1.82
            pD2.2 <- 0.5*var(Dev.draws2) #2.48


      #4. Combine to get DIC= D.pred.avg.hat = 2*D.avg.hat - D.theta.hat
      #                                      = D.avg.hat + p.D1
           DIC.foldedt <- D.avg.hat + pD1    #62.03
           DIC.uniform <- D.avg.hat2 + pD1.2   #63.03
   

      #5. Get DIC for the no pooling (tau=infinity) model (separate means model)
             D.theta.hat.nopool <- Dev.fun(y.j, y.j, sig.j)  #54.64

			 #Drawing nsamp-burn.in draws of 8 draws from separate means and variances
             theta.mat.nopool <- t(replicate(nsamp-burn.in, rnorm(8, mean=y.j, sd=sig.j))) #4500x8
            
             Dev.draws.nopool <- apply(theta.mat.nopool,1, Dev.fun, y.j=y.j, sig.j=sig.j)
             D.avg.hat.nopool <- mean(Dev.draws.nopool)  #62.65
             pD1.nopool <- D.avg.hat.nopool - D.theta.hat.nopool  #8.00

             DIC.nopool <- D.avg.hat.nopool + pD1.nopool  #70.66       
             


      #6. Get DIC for the complete pooling (tau=0) model (single mean model)
             mu.hat.tau0 <- sum(y.j/(sig.j^2))/sum(1/(sig.j^2)) #See page 136
             sig.mu.tau0 <- sqrt(1/(sum(1/(sig.j^2))))
             D.theta.hat.pool <- Dev.fun(rep(mean(mu.vec[burn.in:5000]),8), y.j, sig.j) #59.41682

             theta.mat.pool <- t(replicate(nsamp-burn.in, rnorm(8, mean=mu.hat.tau0, sd=sig.mu.tau0)))

             Dev.draws.pool <- apply(theta.mat.pool,1, Dev.fun, y.j=y.j, sig.j=sig.j)
             D.avg.hat.pool <- mean(Dev.draws.pool)  #60.3
             pD1.pool <- D.avg.hat.pool - D.theta.hat.pool  #0.99

             DIC.pool <- D.avg.hat.pool + pD1.pool   #61.3



      ####Compare posterior deviance distributions between uniform, folded t hiearch. and nopooling and complete pooling 
          par(mfrow=c(2,2))
           hist(Dev.draws, col=gray(.5), nclass=100, main="Deviance (folded-t)", xlim=c(50,80), 
             xlab=expression(D(y,theta)), freq=F, ylim=c(0,0.4))
               abline(v=c(D.theta.hat, D.avg.hat, DIC.foldedt), col=c(2,4,5), lwd=2)
               text(56,0.4, expression(D[hat(theta)]))
           hist(Dev.draws2, col=gray(.5), nclass=100, main="Deviance (uniform)", xlim=c(50,80), 
             xlab=expression(D(y,theta)),freq=F, ylim=c(0,0.4))
               abline(v=c(D.theta.hat2, D.avg.hat2, DIC.uniform), col=c(2,4,5), lwd=2)
               text(55,0.4, expression(D[hat(theta)]))
           hist(Dev.draws.nopool, col=gray(.5), nclass=100, main="Deviance (No pooling )", xlim=c(50,80), 
             xlab=expression(D(y,theta)), freq=F, ylim=c(0,0.4))
               abline(v=c(D.theta.hat.nopool, D.avg.hat.nopool, DIC.nopool), col=c(2,4,5), lwd=2)
               text(56,0.4, expression(D[hat(theta)]))
               text(65,0.4, expression(bar(D(theta))))
               text(72,0.4, "DIC")
           hist(Dev.draws.pool, col=gray(.5), nclass=100, main="Deviance (Complete Pool)", xlim=c(50,80), 
             xlab=expression(D(y,theta)),freq=F, ylim=c(0,0.4))
               abline(v=c(D.theta.hat.pool, D.avg.hat.pool, DIC.pool), col=c(2,4,5), lwd=2)
               text(55,0.4, expression(D[hat(theta)]))


      #### Make table of results like Table 6.2
         out.table <- data.frame( 
               D.theta.hat=round(c(D.theta.hat,D.theta.hat2,D.theta.hat.nopool,D.theta.hat.pool), digits=2),
               D.avg.hat= round(c(D.avg.hat,D.avg.hat2,D.avg.hat.nopool,D.avg.hat.pool), digits=2),
               pD=round(c(pD1, pD1.2, pD1.nopool, pD1.pool),digits=2),
               DIC=round(c(DIC.foldedt, DIC.uniform, DIC.nopool, DIC.pool), digits=2),
               row.names=c("Folded-t", "Uniform", "No pooling", "Complete pool")
                       )
         out.table
         
      #### Also think of DIC as trying to get at Expected MSE using posterior predictive dist'n
      ##  Averaging over all values in y.rep.  Can we do that?
      ##  We already did it for the weighted MSE, let's just do it for unweighted
      ## avg.sq.resid = (1/n)*sum((yi.rep - E(yi.rep|y))^2)
      
        sum.sq.resid.yrep.fun <- function(y.rep, E.yrep, sig.j) {
           #mean(((y.rep-E.yrep)^2))
           sum(((y.rep-E.yrep)/sig.j)^2)
          }
          
        yrep.post.mns <- apply(y.rep.mat, 2, mean)
         #sum.sq.resid.yrep.fun(y.rep.mat[3,], E.yrep=yrep.post.mns, sig.j=sig.j) #test
        
        T2.yrep.theta <- apply(y.rep.mat, 1, sum.sq.resid.yrep.fun, E.yrep=yrep.post.mns, sig.j=sig.j)
        T2.yrep.DICcompare <- -2*(sum(-0.5*(log(2*pi*(sig.j^2)))) - 0.5*T2.yrep.theta)
        mean(T2.yrep.DICcompare)  #64.32
        
      ### Or can get this by calculating expected deviance of yreps at posterior means
      ###Approximate Eyrep[D(yrep,theta.ha)] = Eyrep[-2*log(p(yrep|thetaha))]
          Dev.pred.fun <- function(yrep,theta.hat, sig.j) {
          	  -2*sum(dnorm(yrep, theta.hat, sig.j, log=TRUE))
             }
          D.pred <- apply(y.rep.mat,1,Dev.pred.fun, theta.hat=post.mns, sig.j=sig.j)
          dev.new()
          hist(D.pred, col=gray(.5), nclass=100, main="Post Pred. Deviance at Posterior Means", 
             xlab="Post Pred Deviance", freq=F, xlim=c(40,150))
             abline(v=DIC.foldedt, col=c(5,"purple","purple"), lwd=2)
          D.pred.avg <- mean(D.pred) #64.32
                
                
         par(mfrow=c(2,1))
         hist(Dev.draws, col=gray(.5), nclass=100, main="Deviance (folded-t)", xlim=c(50,80), 
             xlab=expression(D(y,theta)), freq=F, ylim=c(0,0.4))
               abline(v=c(D.theta.hat, D.avg.hat, DIC.foldedt), col=c(2,4,5), lwd=2)
               text(56,0.4, expression(D[hat(theta)]))
         hist(T2.yrep.DICcompare, col=gray(.5), nclass=150, main="Post. Pred. MSE (DIC scale)", 
             xlab="post pred MSE (Dev. scale)", freq=F, xlim=c(50,80))
               abline(v=c(DIC.foldedt, mean(T2.yrep.DICcompare), median(T2.yrep.DICcompare)),
                col=c(5,"purple","magenta"), lwd=2)
               abline(v=D.pred.avg, col="orange", lwd=2)
    
                 
      ### Posterior predictive loss - MSPE- focus on error between yi's and y.rep's
          SSPE.fun <- function(y.rep, y, sig.j) {
          	  #mean(((y-y.rep)^2))
          	  sum(((y-y.rep)/sig.j)^2)
          }   
          T.SSPE <- apply(y.rep.mat, 1, SSPE.fun, y=y.j, sig.j=sig.j)
          T.SSPE.DevScale <- -2*(sum(-0.5*(log(2*3.1415926*(sig.j^2)))) - 0.5*T.SSPE)
          mean(T.SSPE.DevScale) #68.17

   

      #### Let's look at WAIC - Chapter 7 in Gelman et al.
       ## Watanabe-Akaike or "widely available information criterion"
       ## More fully Bayesian approach for estimating out-of-sample expectation
       ## starting with computed log pointwise posterior predictive density, then
       ## addeing correction for effective number of parameters to adjust for overfitting
       ## Advantage = averages over the posterior distribution rather than conditioning on pt. est.
       
       #Can use dev.draws and dev.draws2 from earlier
       
        l.p.y.g.theta.sepj.fun <- function(theta.j, y.j, sig.j) {
             dnorm(y.j, theta.j, sig.j, log=TRUE)
          }
        #l.p.y.g.theta.sepj.fun(theta.mat[2,], y.j, sig.j)
        
       p.y.g.theta.sepj.fun <- function(theta.j, y.j, sig.j) {
              dnorm(y.j, theta.j, sig.j)
          }
        #p.y.g.theta.sepj.fun(theta.mat[2,], y.j, sig.j)
       
       lp.y.g.theta.sepj <- apply(theta.mat, 1, l.p.y.g.theta.sepj.fun, y.j=y.j, sig.j=sig.j) #2nd part of pWAIC
       dim(lp.y.g.theta.sepj) #8 x 5000
       
       ### Calculate p.y.g.theta for each theta, then average, then take log to get computed lppd?
       p.y.g.theta.sepj <- apply(theta.mat, 1, p.y.g.theta.sepj.fun, y.j=y.j, sig.j=sig.j) #calculate for each theta
       mean.dev.sepj <- apply(p.y.g.theta.sepj, 1, mean) #average over theta
       log.mean.dev.sepj <- log(mean.dev.sepj) #take log
       est.lppd <- sum(log.mean.dev.sepj) 
       
       
       #### Calculate the pWAIC1
       pD.WAIC1 <- 2*(sum(log.mean.dev.sepj - mean(lp.y.g.theta.sepj)))
       pD.WAIC1
       
       #### Calculate the pWAIC2
       var.sepj <- apply(lp.y.g.theta.sepj, 1, var)
       pD.WAIC2 <- sum(var.sepj)  #See Equation 7.12 in text #1.009749
       pD.WAIC2
            
       ### Calculate WAIC2 on deviance scale
       WAIC2 <- -2*(est.lppd - pD.WAIC2)  #61.40
       WAIC2
       WAIC1 <- -2*(est.lppd - pD.WAIC1)
       WAIC1
       
       ### For the Uniform prior  ###
       lp.y.g.theta.sepj.2 <- apply(theta.mat2, 1, l.p.y.g.theta.sepj.fun, y.j=y.j, sig.j=sig.j) #2nd part of pWAIC
       p.y.g.theta.sepj.2 <- apply(theta.mat2, 1, p.y.g.theta.sepj.fun, y.j=y.j, sig.j=sig.j) #calculate for each theta
       mean.dev.sepj.2 <- apply(p.y.g.theta.sepj.2, 1, mean) #average over theta
       log.mean.dev.sepj.2 <- log(mean.dev.sepj.2) #take log
       est.lppd.2 <- sum(log.mean.dev.sepj.2) 
       var.sepj.2 <- apply(lp.y.g.theta.sepj.2, 1, var)
       pD.WAIC2.2 <- sum(var.sepj.2)  #See Equation 7.12 in text #1.009749
       WAIC2.2 <- -2*(est.lppd.2 - pD.WAIC2.2)



       #### Look at approximate classical estimate of tau2 - using Deviance from one mean model (red) vs.
       ####  SS of 8 mean model (full)   tau2.hat <- MSB - MSW
          ExtraDev <- (59.35-54.64)
          ExtraDf <- (8-1)
          FullDev <- 54.54
          FullDf <- 8
          tau.2.est <- (ExtraDev/ExtraDf) - (FullDev/FullDf)
      


        library(xtable)
        xtable(out.table)
           
     #### Conditional AIC - 
     #### How to count parameters for the hierachical model? 
     ### What if used hierarchical model, but only interested in the overall mean and not the means of the 
     ####  individual schools?
     
       #Calculate AIC:   AIC = -2log(p(y|mu.hat,tau.hat)) + 2p
                 
            AIC.hm.uniform <- D.theta.hat2 + 2*2     #61.548 = 2 params are mu and tau2?
            AIC.hm.foldedt <- D.theta.hat + 2*2      #62.404 = 2 params are mu and tau2
            AIC.pool <- D.theta.hat.pool + 2*1       #61.349
            AIC.nopool <- D.theta.hat.nopool + 2*8   #70.641
            
             out.table <- data.frame( 
               D.theta.hat=round(c(D.theta.hat,D.theta.hat2,D.theta.hat.nopool,D.theta.hat.pool), digits=2),
               D.avg.hat= round(c(D.avg.hat,D.avg.hat2,D.avg.hat.nopool,D.avg.hat.pool), digits=2),
               pD=round(c(pD1, pD1.2, pD1.nopool, pD1.pool),digits=2),
               DIC=round(c(DIC.foldedt, DIC.uniform, DIC.nopool, DIC.pool), digits=2),
               AIC=round(c(AIC.hm.foldedt, AIC.hm.uniform, AIC.nopool, AIC.pool), digits=2),
               WAIC=round(c(WAIC2, WAIC2.2, NA, NA), digits=2),
               row.names=c("Folded-t", "Uniform", "No pooling", "Complete pool")
                       )
         out.table
            
      #### Bayes factor between the pooled vs. no pooled model?
        ## Can't do it under the improper prior for theta
        ## Ratio unstable for large values of tau2 and undefined for tau=inf   
        ## See Section 7.4 in text  
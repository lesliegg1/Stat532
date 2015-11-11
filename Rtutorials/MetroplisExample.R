

##### Let's do an example where we know the answer to demonstrate
#####   the Metropolis and Metropolis-Hastings algorithms

##Setting:  Suppose  (mu1,mu2)|y ~ BivariateNormal((y1, y2), I2)
##  That is, we have two unknown, independent means with equal, known (sig2=1) variance
##  (y1,y2)=(10,5) (we only have one observation from each group)
##  Clearly, we have other methods for sampling from the normal distribution
##  but let's suppose we do not.


###***************************####
##    METROPOLIS                ##
###***************************####


#1. Specify a symmetric jumping/proposal distribution
#    Let's choose a standard Normal distribution centered on
#    the current value of theta and scaled by sf
#    That is,  mu.cand = mu.cur + sf*Z  where Z~N(0,1)

  draw.cand.fun <- function(cur, sd.scale) {
       mu.draw <- c(rnorm(1, cur[1], sd=sd.scale), rnorm(1, cur[2], sd=sd.scale))
       return(mu.draw)
  }
  #draw.cand.fun(c(6,2), sd.scale=1)  #check the function

#2. Write a function to evaluate the posterior distribution
#   at (mu1,mu2) #write as prod. of two norms since independent
  ## Should really do this on the log scale

  post.fun <- function(mu.vec,y.vec) {
     out <- prod((1/sqrt(2*pi))*exp(-0.5*(mu.vec-y.vec)^2))
     return(out)
    }

#2. Get a starting value for mu

    mu.start <- c(6,2)

#3. Define a few things before we start

    y.obs <- c(10,5)  #observed data --- just n=1 for each

   nsim <- 500  #number of iteration
   mu.mat <- matrix(NA, nrow=nsim, ncol=2) #empty matrix for mu values
   r.vals <- numeric(nsim-1)  #keep track of r's
   jump.vec <- numeric(nsim-1) #keep track of when we jump (accept candidates)
   mu.mat[1,] <- mu.start


plot(mu.mat[,1], mu.mat[,2], type="n", ylim=c(0,10), xlim=c(0,15), 
            main="Bivariate Sample Path Plot", xlab=expression(mu[1]), ylab=expression(mu[2]))
  mu1.vec <- seq(5,15, length=400)
  mu2.vec <- seq(1,10, length=400)
  mu1.mu2.grid <- expand.grid(mu1.vec,mu2.vec)
  out.grid <- matrix(apply(mu1.mu2.grid, 1, post.fun, y.vec=c(10,5)), nrow=400)
  contour(mu1.vec, mu2.vec, out.grid, main="Bivariate Sample Path Plot", 
               xlab=expression(mu[1]), ylab=expression(mu[2]))  
           points(mu.start[1],mu.start[2], col=2, pch=16, cex=1.2)

 for (i in 2:nsim) {
   mu.cur <- mu.mat[i-1,]
   mu.cand <- draw.cand.fun(mu.cur, sd.scale=0.25)
    points(mu.cand[1],mu.cand[2], pch=18, cex=1, col=2)
   
   post.cur <- post.fun(mu.cur, y.vec=y.obs)
   post.cand <- post.fun(mu.cand, y.vec=y.obs)
   
   r <- post.cand/post.cur

   p.accept <- min(1, r)
   u <- runif(1)
   ifelse(u <= p.accept, mu.mat[i,]<- mu.cand, mu.mat[i,] <- mu.cur)

   #Alternatively just draw from a Bernoulli with p=r
   #ind.accept <- rbinom(1,1,r)
   #ifelse(ind.accept==1, mu.mat[i,]<- mu.cand, mu.mat[i,] <- mu.cur)
    points(mu.mat[i,1],mu.mat[i,2], pch=16, cex=1.2, col=3)
    lines(c(mu.mat[i-1,1], mu.mat[i,1]), c(mu.mat[i-1,2], mu.mat[i,2]), col=4)
    
   jump.vec[i-1] <- ifelse(u <= p.accept, 1, 0)
   r.vals[i-1] <- round(r, digits=3)
   }

  ## Look at Bivariate sample path plot
  #dev.new()
  plot(mu.mat[,1], mu.mat[,2], type="n", ylim=c(0,10), xlim=c(0,15), 
            main="Bivariate Sample Path Plot", xlab=expression(mu[1]), ylab=expression(mu[2]))
     lines(mu.mat[,1], mu.mat[,2], col=3)  
     points(mu.mat[1,1], mu.mat[1,2], pch=16, col=4, cex=1.2)
     points(y.obs[1], y.obs[2], col=2, cex=1.4, pch=17)

  #What proportion of draws were accepted?
   mean(jump.vec) 


  par(mfrow=c(2,1))
  plot(seq(1:nsim), mu.mat[,1], type="l", ylab=expression(mu[1]), xlab="iteration", col=4)
  plot(seq(1:nsim), mu.mat[,2], type="l", ylab=expression(mu[2]), xlab="iteration", col=4)

  burn.in <- 100  #Drop the first 100, b/c hasn't converged yet!
  par(mfrow=c(2,1))
  plot(seq((burn.in+1),nsim, by=1), mu.mat[(burn.in+1):nsim,1], type="l", xlab="iteration",
           ylab=expression(mu[1]), col=4)
  plot(seq((burn.in+1),nsim, by=1), mu.mat[(burn.in+1):nsim,2], type="l", xlab="iteration",
           ylab=expression(mu[2]), col=4)




  #############################################################
  ##### Now, let's do it with the unnormalized log posterior distribution
  ###   AND log Metropolis ratio  
  ##############################################################

   log.unpost.fun <- function(mu.vec, y.vec) {
      out <- sum(-0.5*(mu.vec - y.vec)^2)
      return(out)
     }

   y.obs <- c(10,5)

   nsim <- 10000  #number of iteration
   mu.mat <- matrix(NA, nrow=nsim, ncol=2) #empty matrix for mu values
   r.vals <- numeric(nsim-1)  #keep track of r's
   jump.vec <- numeric(nsim-1) #keep track of when we jump (accept candidates)
   mu.mat[1,] <- mu.start

 for (i in 2:nsim) {
   mu.cur <- mu.mat[i-1,]
   mu.cand <- draw.cand.fun(mu.cur, sd.scale=(2.4/sqrt(length(mu.cur))))  
   log.post.cur <- log.unpost.fun(mu.cur, y.vec=y.obs)
   log.post.cand <- log.unpost.fun(mu.cand, y.vec=y.obs)
   log.r <- log.post.cand - log.post.cur

   p.accept <- min(1, exp(log.r))
   u <- runif(1)
   ifelse(u <= p.accept, mu.mat[i,]<- mu.cand, mu.mat[i,] <- mu.cur)

   jump.vec[i-1] <- ifelse(u <= p.accept, 1, 0)
   r.vals[i-1] <- round(r, digits=3)
   }

  ## Look at Bivariate sample path plot
  par(mfrow=c(1,1))
  plot(mu.mat[,1], mu.mat[,2], type="n", ylim=c(0,10), xlim=c(0,15), 
            main="Bivariate Sample Path Plot", xlab=expression(mu[1]), ylab=expression(mu[2]))
     lines(mu.mat[,1], mu.mat[,2], col=3)  
     points(mu.mat[1,1], mu.mat[1,2], pch=16, col=4, cex=1.2)
     points(y.obs[1], y.obs[2], col=2, cex=1.4, pch=17)

  mean(jump.vec)

  par(mfrow=c(2,1))
  plot(seq(1:nsim), mu.mat[,1], type="l", ylab=expression(mu[1]), col=4)
  plot(seq(1:nsim), mu.mat[,2], type="l", ylab=expression(mu[2]), col=4)

  burn.in <- 100  #Drop the first 100, b/c hasn't converged yet!
  par(mfrow=c(2,1))
  plot(seq((burn.in+1):nsim), mu.mat[(burn.in+1):nsim,1], type="l", ylab=expression(mu[1]), col=4)
  plot(seq((burn.in+1):nsim), mu.mat[(burn.in+1):nsim,2], type="l", ylab=expression(mu[2]), col=4)

  par(mfrow=c(2,1)) 
  hist(mu.mat[(burn.in+1):nsim,1], col="lightgray", nclass=20, xlab=expression(mu[1]), main="",
           freq=FALSE, xlim=c(0,15))
  hist(mu.mat[(burn.in+1):nsim,2], col="lightgray", nclass=20, xlab=expression(mu[2]), main="",
           freq=FALSE, xlim=c(0,15))

  ### What about the posterior distribution of mu_{1} - mu_{2} ??

    mu.diff <- mu.mat[(burn.in+1):nsim,1] - mu.mat[(burn.in+1):nsim,2]
    dev.new()
    hist(mu.diff, col="lightgray", nclass=20, xlab=expression(mu[1]-mu[2]), main="",
           freq=FALSE, xlim=c(0,12))
    abline(v=mean(mu.diff), col=2, lwd=2)
    abline(v=quantile(mu.diff, c(.025, 0.975)), col=4, lwd=2)
    quantile(mu.diff, c(.025, 0.975))  #3.214828 7.230235
    mean(mu.diff)  #4.9997
    curve(dnorm(x, mean=(y.obs[1]-y.obs[2]), sd=sqrt(2)), add=TRUE, col=5)
 
    mean(mu.diff>0)  #.999293



   

### From Section 4.3.4 in Link & Barker
### Updated Fall 2014  #####################

expit <- function(x) {exp(x)/(1+exp(x))}
logit <- function(x) { log(x/(1-x))}

# Generate data under the model  ##
set.seed(459)
n <- 10
m <- 25
phi.data <- rnorm(10, mean=0.2, sd=1)
y.data <- rbinom(10, size=m, prob=expit(phi.data))



##########################################################
## Write functions to sample from complete conditional ###
##  distributions                                      ###
##########################################################


## Complete conditional for precision of phi: phi ~ N(mu,prec=tau) where
##  p(tau) = Gamma(a.0, b.0)
draw.cc.tau <- function(mu, phi.vec, a.0, b.0, mu.0,k.0=1) {
                 n <- length(phi.vec)
                 a.gamma <- a.0 + ((n+1)/2)
                 b.gamma <- b.0 + ((1/2)*sum((phi.vec - mu)^2)) + ((k.0/2)*((mu - mu.0)^2))
                 rgamma(1, a.gamma, b.gamma)
                 }

draw.cc.tau(0,phi.data, a.0=0.001, b.0=0.001, mu.0=0)

## Complete conditional for mu: phi ~ N(mu,prec=tau), where
##  mu|tau ~ Normal(mu.0, prec=(K.0*tau)

draw.cc.mu <- function(tau, phi.vec, mu.0, K.0=1) {
                 n <- length(phi.vec)
                 cc.mean <- ((n/(n+K.0))*mean(phi.vec)) + ((K.0/(n+K.0))*mu.0)
                 cc.sd <- sqrt(1/(tau*(K.0+n)))
                 rnorm(1, cc.mean, cc.sd)
  }

draw.cc.mu(1.5,phi.data,mu.0=0, K.0=1)


## Complete conditional distribution for phi's

phi.cc <- function(phi.i, y.i, mu, tau) {
           (expit(phi.i)^(y.i))*((1-expit(phi.i))^(m-y.i))*exp((-tau/2)*((phi.i-mu)^2))
  }
phi.cc(phi.data[1], y.data[1], 0, 1)
log(phi.cc(phi.data[1], y.data[1], 0, 1))


log.phi.cc <- function(phi.i, y.i, mu, tau) {
               y.i*log(expit(phi.i)) + (m-y.i)*log(1-expit(phi.i)) - ((tau/2)*((phi.i-mu)^2))
           }
log.phi.cc(phi.data[1], y.data[1], 0, 1)

  ### Now what?  We can use Metropolis-Hastings or Rejection

 #### Let's try rejection sampling first  #####
 ## Let's use a Beta(y.i, m-y.i) distribution for logit(phi.i)=pi.i's
 ## This should make sense given form of the complete conditional
 ## ********************
 ## This function will allow us to actually control the number of samples
 ## we end up with using rejection sampling.  For each run of the function
 ## it will continue to propose values until one is accepted.

  reject.draw.cc.phi <- function(y.i, mu, tau) {
       keep <- 0
       n.propose <- 0
       while (keep==0) {
       pi.cand <- rbeta(1, y.i, (m-y.i))
       phi.cand <- logit(pi.cand)
       prob.accept <- exp((-tau/2)*((phi.cand-mu)^2))
       keep <- ifelse(runif(1)<prob.accept, 1, 0)
       n.propose <- n.propose + 1
       #print(n.reject)   
      }
      return(c(n.propose,phi.cand))    
     }

     

   ### Check it?  ####
    test.out <- t(replicate(10000, reject.draw.cc.phi(y.data[1], mu=0, tau=1)))
       #first column is the number of rejections BEFORE accepting one for each run through
       #second column is the actual values
    head(test.out)
    test.draws <- test.out[,2]
    sum(test.out[,1]-1)/sum(test.out[,1])  #prop. rejects = rejects/total proposals
    
    #find normalizing constant computationally using numerical integration
    ## We do not NEED to do this for this problem, but this is an example
    ## to refresh how we could get it for this complete conditional distribution
    ## and it allows us to easily check our Rejection and M-H algorithms
      phi.grid <- seq(-2,4,length=400)
      un.phi.dens <- apply(cbind(phi.grid), 1, phi.cc, y=y.data[1], mu=0, tau=1)
      step.size <- phi.grid[2] - phi.grid[1]
      nc <- sum(un.phi.dens*step.size)
      cc.phi.dens <- un.phi.dens/nc
      
       #Check rejection sampling
      hist(test.draws, nclass=100, freq=FALSE, col=gray(0.9)) #Rejection draws
       lines(phi.grid, cc.phi.dens, col=3, lwd=2)  #approx. density
       
       #compare to this histogram - Is this one effective?
      hist(test.draws, nclass=10, freq=FALSE, col=gray(0.9)) #Rejection draws
       lines(phi.grid, cc.phi.dens, col=3, lwd=2)  #approx. density
       
       #Compare to proposal distribution
       cand.draws <- logit(rbeta(100000, y.data[1], m-y.data[1]))
       lines(density(cand.draws),col=2,lwd=3)
       

 #### Now, let's try M-H  #####
 ## Let's use a Beta(y.i, m-y.i) distribution for logit(phi.i)=pi.i's
 ## This should make sense given form of the complete conditional
 ## 
    MH.draw.cc.phi <- function(phi.cur, y.i, mu, tau) {
     pi.cand <- rbeta(1,y.i, m-y.i)
     phi.cand <- logit(pi.cand)
     r <- exp((-tau/2)*(((phi.cand-mu)^2) - ((phi.cur-mu)^2))) #figured out analytically b/c lots of stuff cancels!
     ind.jump <- ifelse(runif(1)<r, 1, 0)
     phi.cur <- ifelse(ind.jump==1, phi.cand, phi.cur)
     return(c(phi.cur,ind.jump))
    }

   #check M-H  
     n.MH <- 10000
     MH.test.draws <- numeric(n.MH)
     MH.test.draws[1] <- logit(0.5)
     for (k in 2:n.MH) {
         MH.test.draws[k] <- MH.draw.cc.phi(MH.test.draws[k-1], y.i=y.data[1], mu=0, tau=1)[1]
       }
     hist(MH.test.draws, nclass=100, freq=FALSE, col=gray(0.9))
       lines(phi.grid, cc.phi.dens, col=3, lwd=2)  #looks pretty good!


#####################################################
### Now set up the actual Gibbs sampling algorithm ##
#####################################################

##1. Set values for hyperparameters
## More thought should be put into these in real life!
## Plot the priors and think about making them weakly
## informative. #****# SEE CODE AT END FOR PLAYING WITH PRIORS
a.0 <- 0.01
b.0 <- 0.01
K.0 <- 1
mu.0 <- 0


##2. Set up vectors and matrices to store results in
n.gibbs <- 100
mu.vec <- numeric(n.gibbs)
tau.vec <- numeric(n.gibbs)
phi.mat <- matrix(nrow=n.gibbs, ncol=length(y.data))
jump.vec <- numeric(n.gibbs-1)

##3. Set initial values
mu.vec[1] <- 2  #initial value for mu
tau.vec[1] <- 0.5  #initial value for tau
phi.mat[1,] <- rep(0,length(y.data)) #logit(0.5)=0

##4. NOW Set up the d steps for each iteration ##

for (t in 2:n.gibbs) {

 #(a) Update tau
 tau.vec[t] <- draw.cc.tau(mu.vec[t-1], phi.mat[t-1,], a.0=a.0, b.0=b.0, mu.0=mu.0)

 #(b) Update mu
 mu.vec[t] <- draw.cc.mu(tau.vec[t], phi.mat[t-1,], mu.0=mu.0, K.0=K.0)

 #(c) Update the vector of phi's

 ### Rejection sampling version  ###
  # for (j in 1:length(y.data)) {
      # phi.mat[t,j] <- reject.draw.cc.phi(y.data[j], mu=mu.vec[t], tau=tau.vec[t])[1]
    # }

 ### M-H version  ####
  for (j in 1:length(y.data)) {
    draw.and.jump <- MH.draw.cc.phi(phi.mat[t-1,j], y.i=y.data[j], mu=mu.vec[t], tau=tau.vec[t])
    phi.mat[t,j] <- draw.and.jump[1]
    jump.vec[t-1] <- draw.and.jump[2]
   }
  }

mean(jump.vec)
jump.vec[1:10]

par(mfrow=c(3,1))
plot(1:n.gibbs, mu.vec, type="n", main=expression(mu), ylab=expression(mu),  xlab="iteration")
  lines(1:n.gibbs, mu.vec, col="red")
    abline(h=0.2, lwd=2, col="darkgreen") #true mu
    abline(h=mean(mu.vec), lwd=2, col=2)  #posterior mean
    abline(h=mean(logit(y.data/m)), lwd=2, col=4) #logit of average of observed proportions

plot(1:n.gibbs, tau.vec, type="n", main=expression(tau), ylab=expression(tau),  xlab="iteration")
  lines(1:n.gibbs, tau.vec, col="red")
    abline(h=1, lwd=2, col="darkgreen") #true tau
    abline(h=mean(tau.vec), lwd=2, col=2)  #posterior mean
    abline(h=1/var(logit(y.data/m)), lwd=2, col=4)

#SD scale, instead of precision or variance
plot(1:n.gibbs, sqrt(1/tau.vec), type="n", main=expression(sqrt(1/tau)), ylab=expression(sqrt(1/tau)),  xlab="iteration")
  lines(1:n.gibbs, sqrt(1/tau.vec), col="red")
    abline(h=1, lwd=3, col="darkgreen") #true tau
    abline(h=sqrt(1/mean(tau.vec)), lwd=2, col=2)  #posterior mean
    abline(h=sd(logit(y.data/m)), lwd=2, col=4)  #observed sd in observed phi's


par(mfrow=c(5,2))
  for (j in 1:10){
    plot(1:n.gibbs, phi.mat[,j], type="l", main=bquote(phi[.(j)]), ylab=bquote(phi[.(j)]),  
          xlab="iteration", ylim=c(-3,3))
    abline(h=phi.data[j], lwd=1, col="darkgreen") #true phi
    abline(h=logit(y.data[j]/m), lwd=1, col="blue") #observed proportion   
}
par(mfrow=c(3,1))
for (j in 1:3){
  plot(1:n.gibbs, phi.mat[,j], type="l", main=bquote(phi[.(j)]), ylab=bquote(phi[.(j)]), 
       xlab="iteration", ylim=c(-3,3))
  abline(h=phi.data[j], lwd=1, col="darkgreen") #true phi
  abline(h=logit(y.data[j]/m), lwd=1, col="blue") #observed proportion   
}


burn.in <- 1:10
par(mfrow=c(1,2))
  hist(mu.vec[-burn.in], col=gray(0.85), nclass=100, freq=F, main="Posterior for mu", xlim=c(-1,2))
    abline(v=0.2, lwd=3, col="darkgreen") #true mu
    abline(v=mean(mu.vec[-burn.in]), lwd=3, col=2)  #posterior mean
    abline(v=mean(logit(y.data/m)), lwd=2, col=4) #logit of average of observed proportions
                                                  #like empirical phi.bar
    legend(0.5,4, legend=c("truth", "posterior mean", "empirical avg"), 
              col=c("darkgreen",2,4), lty=c(1,1,1), lwd=c(2,2,2), bty="n")

  hist(tau.vec[-burn.in], col=gray(0.85), nclass=100, freq=F, main="Posterior for tau")
    abline(v=1, lwd=2, col="darkgreen") #true tau
    abline(v=mean(tau.vec[-burn.in]), lwd=2, col=2)  #posterior mean
    abline(v=(1/var(logit(y.data/m))), lwd=2, col=4) #observed data summary



#dev.new()
par(mfrow=c(2,5))
  for (j in 1:10){
    hist(phi.mat[-burn.in,j], main=bquote(phi[.(j)]), ylab=bquote(phi[.(j)]),  
          xlab="iteration", xlim=c(-3,3), col=gray(0.9))
    abline(v=phi.data[j], lwd=1, col="darkgreen") #true phi
    abline(v=y.data[j]/m, lwd=1, col="blue") #observed proportion 
    abline(v=mean(phi.mat[-burn.in,j]), lwd=1, col=2) #observed proportion  
}

## Apply this function to a matrix where the number of columns=number of parameters
##  and number of rows in MCMC iterations
my.caterpillar <- function(mat) {
  n.p <- dim(mat)[2]
  post.mns <- apply(mat, 2, mean)
  q.out <- apply(mat, 2, function(x) quantile(x, c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995)))
  dev.new()
  plot(seq(min(q.out[1,])-0.5,max(q.out[7,])+0.5,length=n.p), 1:n.p, type="n", xlab=" ", ylab=" ",
        yaxt="n")
   segments(q.out[1,], 1:n.p, q.out[7,], 1:n.p, lwd=2, col=1)
   segments(q.out[2,], 1:n.p, q.out[6,], 1:n.p, lwd=3, col=1)
   segments(q.out[3,], 1:n.p, q.out[5,], 1:n.p, lwd=6, col=1)
   points(q.out[4,], 1:n.p, pch="|", cex=1.5, col=2)
   points(post.mns, 1:n.p, pch="|", cex=1, col=4)
  }


my.caterpillar(phi.mat)
 abline(v=0.2, lwd=1, col="darkgreen") #true mu
 abline(v=mean(mu.vec[-burn.in]), col="green", lwd=2) #posterior mean of mu
 points(logit(y.data/m), 1:length(y.data), pch="|", cex=2, col="orange")
 for (j in 1:10) {
   text(-3.5, j, labels=bquote(phi[.(j)]))
   }
 


###########################################
### Play around with hyperparameters  #####
##############################################

## Assuming variance of beta.0 is related to precision/variance on the logit(pi)'s
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





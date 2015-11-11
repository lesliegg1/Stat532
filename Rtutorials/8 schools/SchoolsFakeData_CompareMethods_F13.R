

#### Schools example in Chapter 5 ####
####  Normal Hierarchical Model   ####

### TO DO -- clean up code and make into .Rnw file  ###

  ### The data we are going to use are estimates of effects (y_j) obtained from 
  ###  separate linear regression models, along with the standard errors of 
  ###  those effects (i.e. we will not be considering the the estimates as
  ###  constants, but rather assuming we know the error associated with them
  ###  The standard errors are considered known (sigma_j) - this is akin to the 
  ###  Chapter 5 results for normal hierarchical models where the sample means
  ###  are modeled rather than the individual observations using theory 
  ###  of sufficient statistics    
  
#### Here we are going to generate FAKE data for individuals (effects of coaching) 
####  from the 8 schools so
####  that we can compare starting "from scratch" on this analysis and compare
####  results to those obtained using                                                     
  
 ### The actual data used in Gelman et al. - estimated effects and SE's
  y.j <- c(28, 8, -3, 7, -1, 1, 18, 12)
  sig.j <- c(15, 10, 16, 11, 9, 11, 10, 18)
  data.mat <- cbind(y.j, sig.j)
  
 n.j <- c(30, 100, 20, 100, 50, 60, 100, 50) 
 n_j <- n.j
 n.tot <- sum(n.j)  #610
 sig.y.j <- sig.j*sqrt(n.j)
 
 set.seed(2567)
 y.ij <- numeric(n.tot)
 school <- rep(NA,n.tot)
 n.index <- c(0, cumsum(n.j)) +1
 for (j in 1:8) {
   y.ij[n.index[j]:(n.index[j+1]-1)] <- rnorm(n.j[j], mean=y.j[j], sd=sig.y.j[j])
   school[n.index[j]:(n.index[j+1]-1)] <- rep(LETTERS[j], n.j[j])
 }
 school <- factor(school)
 
 ### Plot the raw fake data ####
 boxplot(y.ij ~ school, col="gray", var.width=TRUE)
   points(school, y.ij)
   points(1:8, y.j, col="red", cex=2, pch=20)
 
 
### Plot the "raw data" #####
 ## Use lm() to get averages and SE's per school
 ## Separate means, NO POOLING analysis
  options(show.signif.stars=FALSE)
  lm.out <- lm(y.ij ~ school -1)
  summary(lm.out)
 
  lm.out.1SE <- confint(lm.out, level=0.68)
  lm.out.2SE <- confint(lm.out, level=0.95)
  lm.out.3SE <- confint(lm.out, level=0.99)

  
  ### For the original data 
   interval.fun <- function(y.sd, m) {y.sd[1] + c(-1,1)*y.sd[2]*m}
   est.pm.1SD <- apply(data.mat, 1, interval.fun, m=1)  #2x8
   est.pm.2SD <- apply(data.mat, 1, interval.fun, m=2)
   est.pm.3SD <- apply(data.mat, 1, interval.fun, m=3)
   
   
  ### Compare results with data generating values
   plot(seq(min(est.pm.3SD[1,])-0.5, max(est.pm.3SD[2,])+0.5,length=8), seq(0,9,length=8), type="n", xlab=" ", ylab=" ",yaxt="n") 
      mtext(c("A","B","C","D","E","F","G","H"), side=2, at=1:8, line=1,las=2)
      points(y.j, 1:8, pch="|", cex=1.2, col="purple")
     
      segments(est.pm.3SD[1,], 1:8, est.pm.3SD[2,], 1:8, lwd=1, col=1)
      segments(est.pm.2SD[1,], 1:8, est.pm.2SD[2,], 1:8, lwd=3, col="purple")
      segments(est.pm.1SD[1,], 1:8, est.pm.1SD[2,], 1:8, lwd=5, col=4)
      
      segments(lm.out.3SE[,1], (1:8)-0.25, lm.out.3SE[,2], (1:8)-0.25, lwd=1, col=1) 
      segments(lm.out.2SE[,1], (1:8)-0.25, lm.out.2SE[,2], (1:8)-0.25, lwd=3, col="magenta")
      segments(lm.out.1SE[,1], (1:8)-0.25, lm.out.1SE[,2], (1:8)-0.25, lwd=5, col=4)
      points(coef(lm.out), (1:8)-0.25, pch="|", cex=1.2, col=3)
      
  ### Plot NO pooling vs. Complete pooling results
    ###
    lm.1mean <- lm(y.ij ~ 1)
    summary(lm.1mean)
  
  
   plot(seq(min(est.pm.3SD[1,])-0.5, max(est.pm.3SD[2,])+0.5,length=8), seq(0,9,length=8), type="n", xlab=" ", ylab=" ",yaxt="n") 
      mtext(c("A","B","C","D","E","F","G","H"), side=2, at=1:8, line=1,las=2)
      points(coef(lm.out), (1:8)-0.25, pch="|", cex=1.2, col=3)
      
      segments(lm.out.3SE[,1], (1:8)-0.25, lm.out.3SE[,2], (1:8)-0.25, lwd=1, col=1) 
      segments(lm.out.2SE[,1], (1:8)-0.25, lm.out.2SE[,2], (1:8)-0.25, lwd=3, col="magenta")
      segments(lm.out.1SE[,1], (1:8)-0.25, lm.out.1SE[,2], (1:8)-0.25, lwd=5, col=4)
     
      ##Add lines for completely pooled estimate
      abline(v=coef(lm.1mean), lwd=2, col="orange")
      abline(v=confint(lm.1mean), lwd=2, lty=2, col="orange")

     # points(post.mns, 1:n.p, pch="|", cex=1, col=4)  #For later
     
  #### Now, let's just do the random effects version using lmer()
  library(lme4)
  library(arm) #Gelman and Hill functions
  
  lmer.out <- lmer(y.ij~1 + (1|school))
  display(lmer.out)
  coef(lmer.out)
  fixef(lmer.out)
  ranef(lmer.out)
  se.ranef(lmer.out)
  
  int.lmer.1 <- cbind(coef(lmer.out)$school - 1*se.ranef(lmer.out)$school,coef(lmer.out)$school + 1*se.ranef(lmer.out)$school)
  int.lmer.2 <- cbind(coef(lmer.out)$school - 2*se.ranef(lmer.out)$school,coef(lmer.out)$school + 2*se.ranef(lmer.out)$school)
  int.lmer.3 <- cbind(coef(lmer.out)$school - 3*se.ranef(lmer.out)$school,coef(lmer.out)$school + 3*se.ranef(lmer.out)$school)
 
  segments(int.lmer.3[,1], (1:8)-0, int.lmer.3[,2], (1:8)-0, lwd=1, col=1) 
  segments(int.lmer.2[,1], (1:8)-0, int.lmer.2[,2], (1:8)-0, lwd=3, col="magenta")
  segments(int.lmer.1[,1], (1:8)-0, int.lmer.1[,2], (1:8)-0, lwd=5, col=4)
  points(coef(lmer.out)$school[,1], (1:8), pch="|", cex=1.2, col=3)
  
  
#####################################  
#### Now get results from Stan ######
#####################################

library(rstan)
J <- 8
y_A <- y.ij[school=="A"]
y_B <- y.ij[school=="B"]
y_C <- y.ij[school=="C"]
y_D <- y.ij[school=="D"]
y_E <- y.ij[school=="E"]
y_F <- y.ij[school=="F"]
y_G <- y.ij[school=="G"]
y_H <- y.ij[school=="H"]

schools.fake.fit <- stan(file="schools_fake.stan.R", data=c("y_A", "y_B","y_C","y_D","y_E","y_F","y_G","y_H","J"), iter=1000, chains=4)

#schools.fake.fit2 <- stan(file="schools_fake.stan.R", data=c("y_A", "y_B","y_C","y_D","y_E","y_F","y_G","y_H","J","n_j"), iter=1000, chains=4)
print(schools.fake.fit)
plot(schools.fake.fit) 
schools.fake.fit1 <- stan(fit=schools.fake.fit, data=c("y_A", "y_B","y_C","y_D","y_E","y_F","y_G","y_H","J"), iter=10000, chains=4)
print(schools.fake.fit1)

schools.fake.t.fit <- stan(file="schools_fake_halft.stan.R", data=c("y_A", "y_B","y_C","y_D","y_E","y_F","y_G","y_H","J"), iter=1000, chains=4)
schools.fake.t.fit1 <- stan(fit=schools.fake.t.fit, data=c("y_A", "y_B","y_C","y_D","y_E","y_F","y_G","y_H","J"), iter=10000, chains=4)
print(schools.fake.t.fit1)
traceplot(schools.fake.t.fit1, "mu")
traceplot(schools.fake.t.fit1, "tau")

schools.sim <- extract(schools.fake.fit1, permuted=TRUE) #combines chains and permutes values, after warmup
names(schools.sim)

dim(schools.sim$theta)
post.theta <- schools.sim$theta
post.int.1 <- apply(post.theta, 2, function(x) quantile(x,c(0.025,0.975)))
post.int.2 <- apply(post.theta, 2, function(x) quantile(x,c(0.16,0.84)))
post.int.3 <- apply(post.theta, 2, function(x) quantile(x,c(0.005,0.995)))
post.mean <- apply(post.theta, 2, mean)
post.median <- apply(post.theta, 2 ,median)

t.schools.sim <- extract(schools.fake.t.fit1, permuted=TRUE) #combines chains and permutes values, after warmup
t.post.theta <- t.schools.sim$theta
t.post.int.1 <- apply(t.post.theta, 2, function(x) quantile(x,c(0.025,0.975)))
t.post.int.2 <- apply(t.post.theta, 2, function(x) quantile(x,c(0.16,0.84)))
t.post.int.3 <- apply(t.post.theta, 2, function(x) quantile(x,c(0.005,0.995)))
t.post.mean <- apply(t.post.theta, 2, mean)
t.post.median <- apply(t.post.theta, 2 ,median)


#### FIRST, let's compare posteriors for the folded-t vs. the improper uniform prior as we did previously
plot(seq(min(post.int.3[1,])-5, max(post.int.3[2,])+5,length=8), seq(0.25,8.75,length=8), type="n", xlab=" ", ylab=" ",yaxt="n") 
      mtext(c("A","B","C","D","E","F","G","H"), side=2, at=1:8, line=1,las=2)     
       abline(h=seq(1.5,7.5,by=1), lty=2)
       
       abline(v=mean(schools.sim$mu), col="red")
       abline(v=mean(t.schools.sim$mu), col="purple")


### Bayesian hierarchical model - uniform prior on tau and sigma
      segments(post.int.3[1,], (1:8)+0.25, post.int.3[2,], (1:8)+0.25, lwd=1, col=1)
      segments(post.int.1[1,], (1:8)+0.25, post.int.1[2,], (1:8)+0.25, lwd=3, col="orange")
      segments(post.int.2[1,], (1:8)+0.25, post.int.2[2,], (1:8)+0.25, lwd=5, col="red") 
      points(post.mean, (1:8)+0.25, pch="|", cex=1.2, col=1)
      points(post.median, (1:8)+0.25, pch="|", cex=1.2, col=1)
      
      ### Bayesian hierarchical model - folded-t on tau and improper uniform on sigma
      segments(t.post.int.3[1,], (1:8)-0.25, t.post.int.3[2,], (1:8)-0.25, lwd=1, col=1)
      segments(t.post.int.1[1,], (1:8)-0.25, t.post.int.1[2,], (1:8)-0.25, lwd=3, col="magenta")
      segments(t.post.int.2[1,], (1:8)-0.25, t.post.int.2[2,], (1:8)-0.25, lwd=5, col="purple") 
      points(t.post.mean, (1:8)-0.25, pch="|", cex=1.2, col=1)
      points(t.post.median, (1:8)-0.25, pch="|", cex=1.2, col=1)
      
      legend(-45, 9.25, bty="n", legend=c("Improper Uniform", "folded-t(df=1,0,s=5)"), lwd=c(5,5), col=c(2,"purple"))

##### Now plot stan results compared to lm() and lmer()
plot(seq(min(post.int.3[1,])-0.5, max(post.int.3[2,])+0.5,length=8), seq(0.25,8.75,length=8), type="n", xlab=" ", ylab=" ",yaxt="n") 
      mtext(c("A","B","C","D","E","F","G","H"), side=2, at=1:8, line=1,las=2)
      
       abline(h=seq(1.5,7.5,by=1), lty=2)
      
      ##Add lines for completely pooled estimate
      abline(v=coef(lm.1mean), lwd=2, col=1)
      abline(v=mean(schools.sim$mu), col=1)
      abline(v=fixef(lmer.out), col=2)
      
      ### No pooling - lm output
      segments(lm.out.3SE[,1], (1:8)-0.25, lm.out.3SE[,2], (1:8)-0.25, lwd=1, col=1) 
      segments(lm.out.2SE[,1], (1:8)-0.25, lm.out.2SE[,2], (1:8)-0.25, lwd=3, col=gray(0.5))
      segments(lm.out.1SE[,1], (1:8)-0.25, lm.out.1SE[,2], (1:8)-0.25, lwd=5, col=gray(0.1))  
      points(coef(lm.out), (1:8)-0.25, pch="|", cex=1.2, col=1)
      
      ### Linear mixed model - "random effects" - lmer output
      segments(int.lmer.3[,1], (1:8)-0, int.lmer.3[,2], (1:8)-0, lwd=1, col=1) 
      segments(int.lmer.2[,1], (1:8)-0, int.lmer.2[,2], (1:8)-0, lwd=3, col="magenta")
      segments(int.lmer.1[,1], (1:8)-0, int.lmer.1[,2], (1:8)-0, lwd=5, col=4)
      points(coef(lmer.out)$school[,1], (1:8), pch="|", cex=1.2, col=4)
      
      ### Bayesian hierarchical model - uniform prior on tau and sigma
      segments(post.int.3[1,], (1:8)+0.25, post.int.3[2,], (1:8)+0.25, lwd=1, col=1)
      segments(post.int.1[1,], (1:8)+0.25, post.int.1[2,], (1:8)+0.25, lwd=3, col="orange")
      segments(post.int.2[1,], (1:8)+0.25, post.int.2[2,], (1:8)+0.25, lwd=5, col="red") 
      points(post.mean, (1:8)+0.25, pch="|", cex=1.2, col=1)
      points(post.median, (1:8)+0.25, pch="|", cex=1.2, col=1)
      
      legend(-55, 9.25, bty="n", legend=c("No pooling", "lmer()", "Bayesian hier."), lwd=c(5,5,5), col=c(gray(0.1), 4, 2))
      
    
  
     
      
      
      


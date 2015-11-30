### STILL NEED TO GO THROUGH THIS 


### This is meant to directly following obtaining posterior draws
### by running the code in SchoolsNormalHierarchicalExample_F15.R

##### Let's do some model checking (Chapter 6)  ####
## What aspects of reality are NOT captured by the model?
## As with any statistical model -can get misleading results by using a poor model.
## Check adequacy of fit of the model to the data AND the plausibility of the model
##  for purposes for which the model will be used (quote from Gelman- Chapt. 6)
## MODEL= sampling distribution + prior distribution + hierarchical structure 
##     + explanatory variables + etc..
## See Section 6.8 after reading the chapter

 ### Assumptions to think about
  ## (1) normality of y.j|theta.j, sig.j  (sig.j known)
  ## (2) exchangeability related to theta.j's - b/f seeing the results, is 
  ##       there any desire to include school-specific differences in the model, 
  ##        such as two schools being more similar than two other schools - 
  ##          similarity of students, teachers, etc.?
  ## (3) normality of theta.j|mu, tau
  ## (4) prior on mu (improp. uniform)
  ## (5) prior on tau (improp. unifrom, folded-t, etc.)
  

 ### Check 1:  Do the inferences about theta.j seem practically reasonable?  
 ###  yes, they do not violate our knowledge about SAT prep courses.     


 ### Check 2:  Inference about preditive values?
 ##  Let's simulate a hypothetical replication of experiments for post pred. dist'n.

   #Draw yj.rep from a normal distribution with mean theta.j and sd=sig.j
   #  This represents new data from the SAME 8 schools/expts.

       y.rep.mat <- matrix(NA, nrow=(nsamp-burn.in), ncol=8)
       for (j in 1:8){
          y.rep.mat[,j] <- rnorm((nsamp-burn.in), mean=theta.mat[-(1:burn.in),j], sd=sig.j[j])
         }
      school <- c("A","B","C","D","E","F","G","H")
      dev.new()
      par(mfrow=c(4,2))
      for (j in 1:8) {
          hist(y.rep.mat[,j], col=gray(.5), nclass=100, main=paste(school[j]), 
            xlim=c(-50,50))
          abline(v=y.j[j], col=2, lwd=2)
        }  
      max(y.rep.mat)  #74   74/8 corresponds to getting 9 or 10 more questions correct
      min(y.rep.mat)  #-56  56/8 (corresponds to missing 7 questions)
      #neither one of these is completely unrealistic
      
  ## Plot connecting point from same draws from posterior
  
    dev.new()
    par(mfrow=c(1,1))
    plot(1,1, type="n", ylim=c(-50,80), xlim=c(0.5,9.5), ylab="", xaxt="n", xlab="", main="Posterior replications")
     mtext(school, side=1, at=c(1,2,3,4,5,6,7,8), adj=0, las=1, line=0.5)
     lines(1:8, y.j, typ="b", col=2, lwd=4)
     for (s in 1:length(y.rep.mat[,1])) {
     	lines(1:8, y.rep.mat[s,], col=gray(0.7))
        }
     lines(1:8, y.j, typ="b", col=2, lwd=4)
     


 ### More formal posterior predictive checks:  Does model fit important aspects of the data?
  ##.  Example question:  Is the largest observed outcome of 28 points consistent with the model?
  ##  To address this, let's calculate max.j(y.rep.j) - If all lie below 28 points
  ##  then that would indicate that the model does not capture this important
  ##  aspect of the data - We don't want it to shrink school A in too far

   ##Max
    max.j <- apply(y.rep.mat, 1, max)
    dev.new()
    hist(max.j, col=gray(.5), nclass=100, main="Max y.rep", xlim=c(-10,90))
    abline(v=28, col=2, lwd=2)

   ##Min
     min.j <- apply(y.rep.mat, 1, min)

   ##Average
     avg.j <- apply(y.rep.mat, 1, mean)

   ##Sample SD
     s.j <- apply(y.rep.mat, 1, sd) ##Observed SD of effects


   #Plots of the comparisions
    dev.new()
    par(mfrow=c(2,2))
      hist(max.j, col=gray(.5), nclass=100, main="Max", xlim=c(-20,90), 
              xlab="y(rep)", freq=F)
        abline(v=max(y.j), col=2, lwd=2)
        text(50,0.04, paste(mean(max.j>max(y.j))))
      hist(min.j, col=gray(.5), nclass=100, main="Min", xlim=c(-70,10), 
               xlab="y(rep)", , freq=F)
        abline(v=min(y.j), col=2, lwd=2)
        text(6,0.04, paste(mean(min.j>min(y.j))))
      hist(avg.j, col=gray(.5), nclass=100, main="Sample Average", xlim=c(-40,40), 
             xlab="y(rep)", , freq=F)
        abline(v=mean(y.j), col=2, lwd=2)
         text(25,0.06, paste(mean(avg.j>mean(y.j))))
       hist(s.j, col=gray(.5), nclass=100, main="SD of observed effects", xlim=c(0,40),
              xlab="y(rep)", , freq=F)
        abline(v=sd(y.j), col=2, lwd=2)
          text(25,0.1, paste(mean(s.j>sd(y.j))))


   ### Now let's think about further assess sensitivity - Do other models also provide just 
   ###  as good of fit but different conclusions??

     ###Prior on tau:  we already looked at this some by considering the folded-t 
     ##    and uniform, but could consider different folded-t's or other distributions

     ### Normal dist'n for theta.j - could also consider the robust t-dist'n
      ##   See Section 17.4 for more details

     ### Normal likelihood?
      # No need to seriously challenge it - justification is CLT and designs 
      # of studies.  We would need access to original data and not just the 
      # estimates to formally check this
      # This is no different than for classical statistical analysis!


##### Use Stan    ######
library(rstan)
schools.fit <- stan(file="schools.stan.R", data=c("J", "y","sigma"), iter=1000, chains=4)
print(schools.fit)
plot(schools.fit) 

schools.fit1 <- stan(fit=schools.fit, data=c("J", "y","sigma"), iter=2000, chains=4)

schools.sim <- extract(schools.fit1, permuted=TRUE) #combines chains and permutes values, after warmup
names(schools.sim) 

dim(schools.sim$mu)  
dim(schools.sim$theta)

hist(schools.sim$tau)
hist(schools.sim$mu)
mean(schools.sim$theta[,1] > schools.sim$theta[,3])

### Now, can do posterior predictive simulations and graphs directly in R
###. Posterior predictive replications for NEW DATA FROM SAME SCHOOLS
n.sims <- length(schools.sim$lp__)
y.rep <- array(NA, c(n.sims,J))
for (s in 1:n.sims){
	y.rep[s,] <- rnorm(J, schools.sim$theta[s,], sigma)
}

dim(y.rep)

par(mfrow=c(5,4), mar=c(4,4,2,2))
hist(y, xlab="", main="Observed")
for (s in 1:19){
	hist(y.rep[s,], xlab="", main=paste("y.rep",s))
}  #not very satisfying, but can be for larger number of groups

#Compute numerical test statistic - difference between best and second best of 8 coaching programs

ppc.1 <- function(x) {
	x.sort <- rev(sort(x))
	return(x.sort[1]-x.sort[2])
}
t.y <- ppc.1(y)
t.yrep <- apply(y.rep, 1, ppc.1)

par(mfrow=c(1,1))
cat("T(y)=", round(t.y,1), " and T(y.rep) has mean", round(mean(t.yrep),1), "and sd", round(sd(t.yrep),1),
     "\nPr (T(y.rep) > T(y)) =", round(mean(t.yrep>t.y),2), "\n")
     
hist0 <- hist(t.yrep, xlim=range(t.y,t.yrep), nclass=100, col=gray(0.9), xlab="T(y.rep)")
  abline(v=t.y, col=2, lwd=2)
  text(t.y, 0.9*max(hist0$count), "T(y)", adj=0)
  
##### Now, let's get replicated data for NEW schools, not just new data from SAME schools
### Have to make an assumption about the variances of the estimates (the sigma2.j's).  We will just
### assume we get another 8 schools with the same sigma2.j's as the first 8 schools

theta.rep <- array(NA, c(n.sims,J))
y.rep2 <- array(NA, c(n.sims,J))
for (s in 1:n.sims){
	theta.rep[s,] <- rnorm(J, schools.sim$mu[s], schools.sim$tau[s])
	y.rep2[s,] <- rnorm(J, theta.rep[s,], sigma)
}

##Intro example for Rejection Sampling

### Suppose we want to sample from the following pdf:
##  t(Y) = 1/8 if 0 < y <= 1
##         1/4 if 1 < y <= 2
##         1/2 if 2 < y <= 3
##         1/8 if 3 < y <= 4


## Plot it
plot(-2:6, 1:9, type="n", xlim=c(-2,6), ylim=c(0,1.2), xaxp=c(-2,6, 8), ylab="t(y)", xlab="y")
abline(v=0, lty=2)
abline(h=0, lty=1)
segments(c(0,1,2,3), c(1/8,1/4,1/2,1/8), c(1,2,3,4), c(1/8,1/4,1/2,1/8), lwd=2)
points(c(0,1,2,3), c(1/8,1/4,1/2,1/8), pch=1)
points(c(1,2,3,4), c(1/8,1/4,1/2,1/8), pch=16)


### Define a g(y)=1 over [-1,4]
segments(-1,1, 4,1, lwd=2, col="darkgray")
text(4,1.05, "g(y)")

### We can easily normalize g(y), so that it is Uniform(-1,4)
segments(-1, (1/5), 4, (1/5), lwd=2, col="blue")
text(4, (1/5)+.05, "g*(y)") #normalized g(y)

### REJECTION SAMPLING ####
 ### STEP 1 - Draw samples from g*(y) = Unif(-1,4)

g.star.draws <- runif(100000, -1, 4)
hist(g.star.draws, nclass=20, col="lightgray", freq=FALSE)
abline(h=1/5, lwd=2, col=2)

 ### STEP 2 (a) - Get probabilities for accepting each draw
 ### To do this we need a function to calculate t(y) quickly

ty.fun <- function(y) {
  if(y<0) p.y <- 0
  if(y>0 & y<=1) p.y <- 1/8
  if(y>1 & y<=2) p.y <- 1/4     
  if(y>2 & y<=3) p.y <- 1/2
  if(y>3 & y<=4) p.y <- 1/8
  if(y>4) p.y <- 0
  return(p.y)
}
ty.fun(3.5)

 ### STEP 2 (b) -  Calculate t(y) for each value drawn from g*(y)
    t.y <- apply(cbind(g.star.draws), 1, ty.fun)

 ### STEP 3 - calculate the probability of acceptance for each proposed value
    prob.accept <- t.y/((1/2)*1)  ## t(y)/cg(y) where c=1/2 and g(y)=1 for all y
    
 ### STEP 4 - Use the probabilities from STEP 3 to accept or refect candidate values
 # We can do this by just drawing U ~ Unif(0,1) and if U < prob then accept, and U > prob, then reject
    U.vec <- runif(100000)
    ind.accept <- (U.vec < prob.accept)
    y.samps <- g.star.draws[ind.accept==1]
  
   #check what proportion were accepted
    mean(ind.accept) #what proportion did we accept?
    mean(ind.accept[g.star.draws <= 3 & g.star.draws > 2])
    mean(ind.accept[g.star.draws <= 1 & g.star.draws > 0])
    mean(ind.accept[g.star.draws <= 2 & g.star.draws > 1])
   
    
### STEP 5 -- Plot samples and t(y) to check
    #Plot the draws:
    hist(y.samps, nclass=50, col="lightgray", freq=FALSE)
    
    #Plot 
    segments(c(0,1,2,3), c(1/8,1/4,1/2,1/8), c(1,2,3,4), c(1/8,1/4,1/2,1/8), lwd=2, col=2)
    points(c(0,1,2,3), c(1/8,1/4,1/2,1/8), pch=1, col=2)
    points(c(1,2,3,4), c(1/8,1/4,1/2,1/8), pch=16, col=2)

## Option 2 ------------------------------------------------
## Could also start with g*(y) -- normalized version and find the c such that the ratio
## is less than or equal to 1.  So, find c such that p(y)/cg*(y) <= 1

### c= 5/2  because p(y)/g*(y) <- (1/2)/(1/5) = 5/2
### Then cg*(y) = (5/2)*g*(y) = (5/2)*(1/5)=1/2, so end up with same ratio

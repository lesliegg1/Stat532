prior <- function(x){
  ifelse(dunif(x, 0, .384)>0, 1/2, 
         ifelse(dunif(x, .384, .486)>0, 50*x-18.75, 
                ifelse(dunif(x, .484, .584)>0, -50*x+29.975,
                       ifelse(dunif(x, .584, 1)>0, 1/2, 0))))
}

x <- seq(0,1, .00001)
plot(x, prior(x), type="l", ylim=c(0, 6))

y <- 0.00001*prior(x)
lhsum <- sum(y)

###use small grid size for sampling from it

##sample with probability proportional to height of function
draws <- sample(x, 1000, prob=y)
#make histogram and overlay prior
hist(draws, nclass=100, freq=F)
lines(x, prior(x))


mean(draws)
median(draws)
#find probability of being less than 0.385, should be close to 0.2
##Matt idea
empcdf <- cumsum(draws/length(x))
pl.385 <- sum(ifelse(draws<=0.385,1,0))/1000
#prob btwn .385 and .485
pbetwn <- (sum(ifelse(draws<0.485,1,0))-pl.385)/1000

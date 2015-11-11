mu.post.fun <- function(prior.mean, prior.sd, ybar, n, data.sd=1, x.lo, x.hi, ymax=1.4) {
  post.prec <- (1/(prior.sd)^2) + (n/(data.sd^2))
  post.mean <- ((prior.mean/(prior.sd^2)) + ((ybar*n)/(data.sd^2)))/(post.prec)
  post.sd <- sqrt(1/post.prec)
  post.int <- qnorm(c(.05,.95), mean=post.mean, sd=post.sd)
  cat("posterior mean=", post.mean, "\n")
  cat("90% posterior interval = ", post.int, "\n")
  #return(c(post.mean,post.int))
  curve(dnorm(x,mean=prior.mean, sd=prior.sd), col="red", from=x.lo, to=x.hi, ylim=c(0,ymax), lwd=2,
        ylab="density", xlab=expression(mu))
  curve(dnorm(x,mean=ybar, sd=data.sd/sqrt(n)), add=TRUE, col="blue", lwd=2)
  curve(dnorm(x,mean=post.mean, sd=post.sd), add=TRUE,  col="purple", lwd=4)
  legend(x.lo, ymax, legend=c("Sampling Dist'n", "Prior","Posterior"), 
         lwd=c(2,2,4), col=c("blue","red","purple"), bty="n")
}

mu.post.fun(19, 5, ybar=5, n=12, data.sd=10, x.lo=5, x.hi=35, ymax=0.5)

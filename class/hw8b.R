post.thetaGpos <- function(x, a=10, b=90) {(.82*x + 0.1)*(x^(a-1))*(1-x)^{(b-1)}}
curve(post.thetaGpos, from=0, to=1, a=10, b=90)

grid.theta <- matrix(seq(0,1,length=1000), ncol=1000, nrow=1)
un.post <- apply(grid.theta, 2, post.thetaGpos, a=2, b=2)
norm.post <- un.post/(sum(un.post*(grid.theta[2]-grid.theta[1])))

play.w.Betaprior <- function(alpha, beta, x.low=0, x.hi=1) {
  grid.theta <- matrix(seq(0,1,length=1000), ncol=1000, nrow=1)
  un.post <- apply(grid.theta, 2, post.thetaGpos, a=alpha, b=beta)
  norm.post <- un.post/(sum(un.post*(grid.theta[2]-grid.theta[1])))
  plot(grid.theta, norm.post, col=4, lwd=2, type="l", xlim=c(x.low,x.hi))
  curve(dbeta(x,alpha,beta), add=TRUE, col="magenta")
}


##Blue is posterior, prior is magenta
play.w.Betaprior(alpha=10, beta=90, x.low=0, x.hi=0.25)
play.w.Betaprior(alpha=5, beta=100, x.low=0, x.hi=0.25)
play.w.Betaprior(alpha=1, beta=1)
play.w.Betaprior(alpha=10, beta=10)
play.w.Betaprior(alpha=2, beta=2)
play.w.Betaprior(alpha=5, beta=5)
play.w.Betaprior(alpha=0.01, beta=0.01)
play.w.Betaprior(alpha=0.1, beta=10)
play.w.Betaprior(alpha=10, beta=0.1)

play.w.Betaprior(alpha=50, beta=5, x.low=0.75, x.hi=1.0)

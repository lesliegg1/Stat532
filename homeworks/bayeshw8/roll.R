

roll <- function(x0, grad, posterior, y = 5, n = 10, mu = 0, sigma = 2, eps = 0.1, L = 100, p.t = rnorm(1, 0, 1))
{
	x <- numeric(L)

	for (i in 1:L)
	{
		x0 <- x0 + p.t * eps / 2
		p.t <- p.t - grad(x0, y, n, mu, sigma) * eps / 2
		x0 <- x0 + p.t * eps / 2
		x[i] <- x0
	}

	xm <- min(x); xM <- max(x)
	r <- xM - xm
	for (i in 1:length(x))
	{
		plot(xx <- seq(xm - r/2, xM + r/2, 0.01), -log(posterior(xx, y, n, mu, sigma)), type="l")
		points(x[i], -log(posterior(x[i], y, n, mu, sigma)), col="red", cex=2, pch=16)
		Sys.sleep(0.05)
	}
}
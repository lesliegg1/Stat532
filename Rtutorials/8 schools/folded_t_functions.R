library(LearnBayes)
 #following function does not incorporate a scale parameter
 #function to evaluate density of half-t
   d.tfold <- function(x, df) {2*dt(x,df=df) } 
 	

   #function to go with pg 520 in Gelman (2006)
   un.tfold <- function(x, df, s=1) { (1 + (1/df)*((x/s)^2))^(-0.5*(df+1)) } 
   	
   	
   #function to obtain random draws from a scaled, non-central folded t
   r.tfold <- function(n, df, s=1, mn=0) { #get random draws from folded t
   	    z <- abs(rnorm(n, mn))
   	    x <- rchisq(n, df=df)
   	    theta <- s*z*sqrt(df/x)
   	    return(theta)
   	}
   	
alpha <- seq(-5, 10, length=20)
beta <- seq(-10, 40, length=20)

draw <- expand.grid(alpha, beta)

draw$priorprobs <- 1/length(draw$Var1)



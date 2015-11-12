data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0,upper=1> y[N];
}
parameters {
  real alpha;
  real<lower=0> beta;
}
model {
  alpha ~ cauchy(0, 2);
  beta ~ cauchy(0, 2.5);
  y ~ bernoulli_logit(alpha + beta * x);
}
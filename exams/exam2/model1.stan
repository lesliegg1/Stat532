data {
  int y_1; 
  int y_2;
  int y_3;
  int y_4;
  int y_5;
}

parameters {
  real<lower=0> m;
  real<lower=0> lambda;
  real<lower=0, upper=1> pi;
  real<lower=0, upper=1> eta;
  real<lower=0> kappa;
}

model {
  pi ~ beta(eta*kappa, kappa*(1-eta)); 
  m ~ poisson(100);
  y_1 ~ binomial(m, pi);
  y_2 ~ binomial(m, pi);
  y_3 ~ binomial(m, pi);
  y_4 ~ binomial(m, pi);
  y_5 ~ binomial(m, pi);
}

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
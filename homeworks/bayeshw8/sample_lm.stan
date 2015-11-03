data {
    int<lower=0> N;
    int<lower=0> p;
    vector[N] y;
    matrix[N, p] x;
}

parameters {
    vector[p] beta;
    real<lower=0> sigma;
}

model {
    beta[1] ~ normal(100, 1);
    y ~ normal(x * beta, sigma);
}

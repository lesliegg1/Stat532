// There were a couple small mistakes in here, but they are fixed now :)
data {
    int<lower=0> N;
    int<lower=0> Ncens;
    vector<lower=0>[N-Ncens] y;
    real<lower=0> C;
}

parameters {
    real<lower=0> lambda;
}

transformed parameters{
    real<lower=0> theta;

    theta <- 1 / lambda;
}

model {
    for (i in 1:(N-Ncens))
    {
     y[i] ~ exponential(lambda);
    }

    increment_log_prob(Ncens * exponential_ccdf_log(C, lambda));
}

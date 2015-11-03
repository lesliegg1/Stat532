data {
}

parameters {
    real theta;
}

transformed parameters {
}

model {
    // increment_log_prob(normal_log(theta, 0, 1));
    // theta ~ normal(0, 1);
    increment_log_prob(- theta ^ 2 / 2);
}

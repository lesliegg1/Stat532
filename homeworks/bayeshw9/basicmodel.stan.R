data {
  int J; 
  real y_1[5]; 
  real y_2[10];
  real y_3[30];
  real y_4[30];
  real y_5[20];
  real y_6[25];
  real y_7[50];
  real y_8[10];
}

parameters {
  real mu; 
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_y;
  real alpha[J];
}

transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] <- mu + alpha[j];
}

model {
  alpha ~ normal(0, sigma_alpha);
  
  y_1 ~ normal(theta[1], sigma_y);
  y_2 ~ normal(theta[2], sigma_y);
  y_2 ~ normal(theta[3], sigma_y);
  y_3 ~ normal(theta[4], sigma_y);
  y_4 ~ normal(theta[5], sigma_y);
  y_5 ~ normal(theta[6], sigma_y);
  y_6 ~ normal(theta[7], sigma_y);
  y_7 ~ normal(theta[8], sigma_y); 
}
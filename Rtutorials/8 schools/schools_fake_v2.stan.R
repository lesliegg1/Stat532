
#THIS DOES NOT WORK using n_j's 

data {
    int J; 
    int n_j;
    real y_A[n_j[1]]; 
    real y_B[n_j[2]];
    real y_C[n_j[3]];
    real y_D[n_j[3]];
    real y_E[n_j[4]];
    real y_F[n_j[5]];
    real y_G[n_j[6]];
    real y_H[n_j[7]];
  }
  
  parameters {
    real mu; 
    real<lower=0> tau;
    real<lower=0> sigma_y;
    real eta[J];
  }
  
  transformed parameters {
    real theta[J];
    for (j in 1:J)
      theta[j] <- mu + tau * eta[j];
  }
  
  model {
  	eta ~ student_t(1,0,5);
    y_A ~ normal(theta[1], sigma_y);
    y_B ~ normal(theta[2], sigma_y);
    y_C ~ normal(theta[3], sigma_y);
    y_D ~ normal(theta[4], sigma_y);
    y_E ~ normal(theta[5], sigma_y);
    y_F ~ normal(theta[6], sigma_y);
    y_G ~ normal(theta[7], sigma_y);
    y_H ~ normal(theta[8], sigma_y); 
  }
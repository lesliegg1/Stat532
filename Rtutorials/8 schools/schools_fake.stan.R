  data {
    int J; 
    real y_A[30]; 
    real y_B[100];
    real y_C[20];
    real y_D[100];
    real y_E[50];
    real y_F[60];
    real y_G[100];
    real y_H[50];
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
   eta ~ normal(0, 1);
    
    y_A ~ normal(theta[1], sigma_y);
    y_B ~ normal(theta[2], sigma_y);
    y_C ~ normal(theta[3], sigma_y);
    y_D ~ normal(theta[4], sigma_y);
    y_E ~ normal(theta[5], sigma_y);
    y_F ~ normal(theta[6], sigma_y);
    y_G ~ normal(theta[7], sigma_y);
    y_H ~ normal(theta[8], sigma_y); 
  }
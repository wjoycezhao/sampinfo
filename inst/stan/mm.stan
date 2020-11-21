
data{
  // memory
  int deltaM_value; //  number of option categories
  int condition_value; // number of option categories
  int<lower=2> C; // number of option categories
  int<lower=0> K; // number of features; if 0 then no betas
  int<lower=1> N; // number of time points
  int<lower=1> tNo[N]; // thought No.x; used to reset decay
  matrix[C,K] X[N]; // feature matrix, for different time points
  int<lower=1, upper=C> cID[N]; // response cluster ID, 1-7 for 3-cluster solution
  int<lower=2> MC; // ID of the neutral category
}

parameters{
  real<lower=0, upper=1> deltaM [deltaM_value == 8]; // decay
  vector[C-1] alpha_raw; // base rate
  vector[K] beta_raw; // group-level mean for beta
}

transformed parameters {
  vector[C] alpha;
  vector[K] beta; // group-level mean for beta
  matrix[C,N] theta;
  matrix[C,K * (deltaM_value == 8)] X_acc;

  alpha[1:(C-1)] = alpha_raw * 2;
  alpha[C] = 0;
  beta = beta_raw * 5;

  if (K > 0 && deltaM_value == 8){
      for (n in 1:N){
            if (tNo[n] == 1){
              X_acc = to_matrix(rep_array(0,C,K));
            }else{
              X_acc *= deltaM[1];
            }
            X_acc += X[n];
            theta[,n] = X_acc * beta;
      }
      theta += rep_matrix(alpha,N);
  }else if (K > 0){
      for (n in 1:N){
            theta[,n] = X[n] * beta;
      }
      theta += rep_matrix(alpha,N);
  }else{
      theta = rep_matrix(alpha,N);
  }

  if(condition_value == 1){
      for (n in 1:N){
          if (tNo[n] == 1){
              if(cID[n] < MC){
                  theta[(MC + 1):C, n] += to_vector(rep_array(-999999, MC - 1));
              }else if (cID[n] > MC){
                  theta[1:(MC - 1), n] += to_vector(rep_array(-999999, MC - 1));
              }
          }
      }
  }
}
model{
  if(deltaM_value == 8){
    deltaM[1] ~ uniform(0,1);
  }
  alpha_raw ~ normal(0,1);
  beta_raw ~ normal(0,5);
  for (n in 1:N){
    cID[n] ~ categorical_logit(theta[,n]);
  }
}

generated quantities {
  vector[N] log_lik;
  real dev_m;
  dev_m = 0;
  for (n in 1:N){
    log_lik[n] = categorical_logit_lpmf(cID[n]|theta[,n]);
    dev_m += - 2*log_lik[n];
  }
}

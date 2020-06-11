
data{
  // memory
  int deltaM_value; // number of option categories
  int<lower=2> C; // number of option categories
  int<lower=0> K; // number of features; if 0 then no betas
  int<lower=1> N; // number of time points
  int<lower=1> tNo[N]; // thought No.x; used to reset decay
  matrix[C,K] X[N]; // feature matrix, for different time points
  int<lower=1, upper=C> cID[N]; // response cluster ID, 1-7 for 3-cluster solution
}

parameters{
  real<lower=0> deltaM [deltaM_value == 9]; // decay
  vector[C-1] alpha_raw; // base rate
  vector[K] beta; // group-level mean for beta
}

transformed parameters {
  vector[C] alpha;
  matrix[C,N] theta;
  matrix[C,K * (deltaM_value == 9)] X_acc;

  alpha[1:(C-1)] = alpha_raw/2;
  alpha[C] = 0;

  if (K > 0 && deltaM_value == 9){
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
}

model{
  if(deltaM_value==9){
    deltaM[1] ~ normal(0.5,0.8)T[0,];
  }
  beta ~ std_normal();
  alpha_raw ~ std_normal();
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

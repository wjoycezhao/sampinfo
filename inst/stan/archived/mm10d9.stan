
data{
  // memory
  int<lower=2> C; // number of option categories
  int<lower=0> K; // number of features; if 0 then no betas
  int<lower=1> N; // number of time points
  int<lower=1> S; // number of participants
  int<lower=1> Q; // number of questions
  int<lower=0> sID[N]; // participant ID variable
  int<lower=0> qID[N]; // question ID variable
  int<lower=1> tNo[N]; // thought No.x; used to reset decay
  matrix[C,K] X[N]; // feature matrix, for different time points
  int<lower=1, upper=C> cID[N]; // response cluster ID, 1-7 for 3-cluster solution
}

parameters{
  real<lower=0> deltaM; // decay
  matrix[C-1,Q] alpha_raw; // base rate
  vector[K] beta_mu; // group-level mean for beta
  vector<lower=0>[K] beta_s_sd; // sd for beta across individuals
  vector<lower=0>[K] beta_q_sd; // sd for beta across questions
  matrix[K,Q] beta_q_raw; // question-level beta
  matrix[K,S] beta_s_raw; // individual-level beta
}

transformed parameters {
  matrix[C,Q] alpha;
  matrix[C,N] theta;
  matrix[K,S] beta[Q];
  if (K > 0){
    for (q in 1:Q){
    beta[q] = rep_matrix(beta_mu,S) +
              diag_pre_multiply(beta_s_sd,beta_s_raw) +
              rep_matrix(beta_q_sd .* beta_q_raw[,q],S);
    }
  }

  alpha[1:(C-1),] = alpha_raw/2;
  alpha[C,] = to_row_vector(rep_array(0,Q));

  {
    if (K > 0){
      matrix[C,K] X_acc;
      for (n in 1:N){
        if (tNo[n] == 1){
          X_acc = to_matrix(rep_array(0,C,K));
        }else{
          X_acc *= deltaM;
        }
        X_acc += X[n];
        theta[,n] = X_acc * beta[qID[n],,sID[n]];
      }
      theta += alpha[,qID];
    }else{
      theta = alpha[,qID];
    }
  }
}

model{
  deltaM ~ normal(0.5,0.8)T[0,];
  beta_mu ~ std_normal();
  beta_q_sd ~ std_normal();
  beta_s_sd ~ std_normal();
  to_vector(beta_q_raw) ~ std_normal();
  to_vector(beta_s_raw) ~ std_normal();
  to_vector(alpha_raw) ~ std_normal();

  for (n in 1:N){
    cID[n] ~ categorical_logit(theta[,n]);
  }
}

generated quantities {
  vector[N] log_lik;
  real dev_m;
  matrix[K,Q] beta_qdev;
  matrix[K,S] beta_sdev;
  vector[K] beta_q_temp;
  vector[K] beta_s_temp;
  matrix[K,Q] beta_qmean;
  matrix[K,S] beta_smean;
  dev_m = 0;
  for (n in 1:N){
    log_lik[n] = categorical_logit_lpmf(cID[n]|theta[,n]);
    dev_m += - 2*log_lik[n];
  }
  beta_qdev = diag_pre_multiply(beta_q_sd, beta_q_raw);
  beta_sdev = diag_pre_multiply(beta_s_sd,beta_s_raw);
  if(K>0){
   for (k in 1:K) {
     beta_q_temp[k] = mean(to_array_1d(beta_qdev[k,]));
     beta_s_temp[k] = mean(to_array_1d(beta_sdev[k,]));
     }
   beta_qmean = rep_matrix(beta_mu+beta_s_temp,Q) + beta_qdev;
   beta_smean = rep_matrix(beta_mu+beta_q_temp,S) + beta_sdev;
  }
}

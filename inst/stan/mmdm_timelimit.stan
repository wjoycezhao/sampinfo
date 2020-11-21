
functions{
///
      vector[] decision_l (int N, real ar_value,
                          vector p_n, vector p_y, vector p_end){
            //variable to hold prob
          vector[N] prob_d[3];
          // to say no
          prob_d[1]=  p_end .* (ar_value * (p_n +  .5 * (1-p_n-p_y)) +
                              (1-ar_value)*(p_n .* (1-p_y) + .5*p_y .* p_n + .5*(1-p_n) .* (1-p_y)));
          // to continue
          prob_d[2] =  (1 - p_end);

          // to say yes
          prob_d[3] = p_end .* (ar_value * (p_y + .5*(1-p_n-p_y)) +
                              (1-ar_value)*(p_y .* (1-p_n) + .5*p_y .* p_n + .5*(1-p_n) .* (1-p_y)));
          return prob_d;
      }
      int sign(int x) {
          int x_sign;
          if(x<0){x_sign=-1;}
          if(x==0){x_sign = 0;}
          if(x>0){x_sign=1;}
          return x_sign;
      }
}

data{
  // memory
  int deltaM_value; // number of option categories
  int condition_value; // number of option categories
  int<lower=2> C; // number of option categories
  int<lower=0> K; // number of features; if 0 then no betas
  int<lower=1> N; // number of time points
  int<lower=1> tNo[N]; // thought No.x; used to reset decay
  matrix[C,K] X[N]; // feature matrix, for different time points
  int<lower=1, upper=C> cID[N]; // response cluster ID, 1-7 for 3-cluster solution
  int<lower=2> MC; // ID of the neutral category
  // decision
  int ar_value;
  int binary_value;
  int deltaD_value;
  int terminate[N];
  int rating[N];
  int rating_y[N];
  int rating_n[N];
  // simulations
  matrix[C,K] X_sim[C]; // feature matrix, for categoreis
  vector[7] cluster_rating_m[C];
  int max_tNo_prd;
}

parameters{
  // memory
  real<lower=0, upper=1> deltaM [deltaM_value == 8]; // decay
  vector[C-1] alpha_raw; // base rate
  vector[K] beta_raw; // group-level mean for beta
  // decision
  real<lower=0, upper=1> deltaD [deltaD_value == 8];
  real<lower=0,upper=1> p_end;
  real<lower=0> sigma_raw;
}

transformed parameters {
  // memory
  vector[C] alpha;
  vector[K] beta; // group-level mean for beta
  matrix[C,N] theta;
  matrix[C,K * (deltaM_value == 8)] X_acc;
  // decision
  vector[N] p_n;
  vector[N] p_y;
  real<lower=0> sigma;

  // memory
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
  // decision
  sigma = 35* sigma_raw;
{
  real utility_n;
  real utility_y;
  if (deltaD_value == 8){
      for (n in 1:N){
        if (tNo[n] == 1){
          utility_n = rating_n[n];
          utility_y = rating_y[n];
          }else{
          utility_n *= deltaD[1];
          utility_n += rating_n[n];
          utility_y *= deltaD[1];
          utility_y += rating_y[n];
          }
       p_n[n] = 1 - normal_cdf(utility_n, 0, sigma);
       p_y[n] = normal_cdf(utility_y, 0, sigma);
      }
  }else{
      for (n in 1:N){
       p_n[n] = 1 - normal_cdf(rating_n[n], 0, sigma);
       p_y[n] = normal_cdf(rating_y[n], 0, sigma);
      }
  }
}
}

model{
  // decision
  vector[N] decision_p[3];
  // memory
  if(deltaM_value == 8){
    deltaM[1] ~ uniform(0,1);
  }
  alpha_raw ~ normal(0,1);
  beta_raw ~ normal(0,5);
  for (n in 1:N){
    cID[n] ~ categorical_logit(theta[,n]);
  }

  // decision
  if(deltaD_value == 8){
    deltaD[1] ~ uniform(0,1);
  }
  sigma_raw ~ normal(0, 1)T[0,];
  decision_p = decision_l (N, ar_value, p_n, p_y, to_vector(rep_array(p_end,N)));
  for (n in 1:N){
      // decision
      target += log(decision_p[terminate[n]+2,n]);
  }
}

generated quantities {
  // memory
  vector[N] log_lik_m;
  real dev_m;
  // decision
  vector[N] log_lik_d;
  vector[N] decision_p[3];
  vector[N] decision_no;
  vector[N] decision_0;
  vector[N] decision_yes;
  real dev_d;
  // simulations
  int cID_prd[max_tNo_prd];
  int decision_prd;
  int terminate_prd;

  // memory WAIC
  dev_m = 0;
  for (n in 1:N){
    log_lik_m[n] = categorical_logit_lpmf(cID[n]|theta[,n]);
    dev_m += - 2*log_lik_m[n];
  }
  // decision WAIC
  dev_d = 0;
  decision_p = decision_l (N, ar_value, p_n, p_y, to_vector(rep_array(p_end,N)));
  decision_no =  decision_p[1];
  decision_0 = decision_p[2];
  decision_yes =  decision_p[3];
  for (n in 1:N){
      log_lik_d[n] = log(decision_p[terminate[n]+2,n]);
      dev_d += - 2*log_lik_d[n];
  }
  // simulations
cID_prd = rep_array(99,max_tNo_prd);
{
  matrix[C,K * (deltaM_value == 8)] X_acc_prd;
  vector[C] theta_prd;
  int rating_prd;
  int rating_n_prd;
  int rating_y_prd;
  real utility_n_prd;
  real utility_y_prd;
  real p_n_prd;
  real p_y_prd;
  vector[1] decision_p_prd[3];
  vector[3] decision_p_prd1;
  int tNo_prd;

  tNo_prd= 1;
  terminate_prd = 0;
  while (terminate_prd == 0 && tNo_prd < max_tNo_prd + 1){
    if (K > 0 && deltaM_value == 8){
      if (tNo_prd == 1){
        X_acc_prd = to_matrix(rep_array(0,C,K));
      }else{
        X_acc_prd *= deltaM[1];
        X_acc_prd += X_sim[cID_prd[tNo_prd - 1]];
      }
      theta_prd = X_acc_prd * beta;
      theta_prd += alpha;
    }else if (K > 0 && deltaM_value == 0){
      if (tNo_prd == 1){
        X_acc_prd = to_matrix(rep_array(0,C,K));
      }else{
        X_acc_prd = X_sim[cID_prd[tNo_prd - 1]];
      }
      theta_prd = X_acc_prd * beta;
      theta_prd += alpha;
    }else if (K > 0 && deltaM_value == 1){
      if (tNo_prd == 1){
        X_acc_prd = to_matrix(rep_array(0,C,K));
      }else{
        X_acc_prd += X_sim[cID_prd[tNo_prd - 1]];
      }
      theta_prd = X_acc_prd * beta;
      theta_prd += alpha;
    }else{
      theta_prd = alpha;
    }

    cID_prd[tNo_prd] = categorical_logit_rng(theta_prd);

    rating_prd = categorical_rng(cluster_rating_m[cID_prd[tNo_prd]])-4;

    if (binary_value == 1) {rating_prd = sign(rating_prd);}
    if (ar_value == 0){
      rating_n_prd = min(rating_prd,0);
      rating_y_prd = max(rating_prd,0);
    }else{
      rating_n_prd = rating_prd;
      rating_y_prd = rating_prd;
    }

    if (deltaD_value == 8){
          if (tNo_prd == 1){
            utility_n_prd = rating_n_prd;
            utility_y_prd = rating_y_prd;
          }else{
            utility_n_prd *= deltaD[1];
            utility_n_prd += rating_n_prd;
            utility_y_prd *= deltaD[1];
            utility_y_prd += rating_y_prd;
          }
    }else if (deltaD_value == 0){
            utility_n_prd = rating_n_prd;
            utility_y_prd = rating_y_prd;
    }else if (deltaD_value == 1){
          if (tNo_prd == 1){
            utility_n_prd = rating_n_prd;
            utility_y_prd = rating_y_prd;
          }else{
            utility_n_prd += rating_n_prd;
            utility_y_prd += rating_y_prd;
          }
    }
    p_n_prd = 1 - normal_cdf(utility_n_prd, 0, sigma);
    p_y_prd = normal_cdf(utility_y_prd, 0, sigma);
    decision_p_prd = decision_l (1, ar_value, to_vector(rep_array(p_n_prd,1)), to_vector(rep_array(p_y_prd,1)), to_vector(rep_array(p_end,1)));
    for(i in 1:3){decision_p_prd1[i] = decision_p_prd[i,1];}
    terminate_prd = categorical_rng(decision_p_prd1) - 2;
    tNo_prd += 1;
  }
}
}


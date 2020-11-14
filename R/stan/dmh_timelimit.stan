// with deltaM; when deltaM is constrained to a fixed number x, lpdf = log(1/(1-1))=0 see binary -> p_end /3; random p_end/10
// https://groups.google.com/forum/#!topic/stan-users/pIEJOTEtbAU
// generate prediction for support distribution

functions{
// 
      vector decision_l (int ar_value, int terminate, int N, vector p_n, vector p_y){
            //variable to hold prob
          vector[N] prob_d; 
          if(ar_value == 0){
              if (terminate == -1){
                  prob_d = p_n .* (1-p_y) + .5 * p_y .* p_n + .5 * (1-p_y) .* (1-p_n);
              }
              else if (terminate == 1){
                  prob_d = p_y .* (1-p_n) + .5 * p_y .* p_n + .5 * (1-p_y) .* (1-p_n);
              }            
          }else{
              if (terminate == -1){
                  prob_d = p_n + .5 * (1-p_y-p_n);
              }
              else if (terminate == 1){
                  prob_d = p_y + .5 * (1-p_y-p_n);
              }  
          }
         return (prob_d);
}
}


data{
  // decision
  int deltaD_value;
  int ar_value;
  int<lower=1> N;
  int S;
  int Q;
  int sID[N];
  int qID[N];
  int tNo[N];
  int rating_n[N];
  int rating_y[N];
  int terminate[N]; 
  int N_continue;
  int N_no;
  int N_yes;
  int idx_continue[N_continue];
  int idx_no[N_no];
  int idx_yes[N_yes];
}

parameters{
  // decision
  real deltaD_mu_raw [deltaD_value==8];
  real<lower=0> deltaD_s_sd [deltaD_value==8];
  real<lower=0> deltaD_q_sd [deltaD_value==8];
  vector[S*(deltaD_value==8)] deltaD_s_raw;
  vector[Q*(deltaD_value==8)] deltaD_q_raw;
  real p_end_mu_raw;
  real<lower=0> p_end_s_sd;
  real<lower=0> p_end_q_sd;
  vector[S] p_end_s_raw;
  vector[Q] p_end_q_raw;
  real sigma_mu_raw;
  real<lower=0> sigma_s_sd;
  real<lower=0> sigma_q_sd;
  vector[S] sigma_s_raw;
  vector[Q] sigma_q_raw;
}
transformed parameters {
  // decision
  vector[S] deltaD[Q*(deltaD_value==8)];
  vector[S] p_end[Q];
  vector[S] sigma[Q];
  
  if (deltaD_value == 8){
    for (q in 1:Q){
      deltaD[q] = inv_logit(deltaD_mu_raw[1] + 
                      deltaD_s_sd[1] * deltaD_s_raw + 
                      deltaD_q_sd[1] * deltaD_q_raw[q]);
    }
  }
  for (q in 1:Q){
    p_end[q] = inv_logit(p_end_s_raw +
                      p_end_q_sd * p_end_q_raw[q]);
    sigma[q] = exp(sigma_mu_raw + 
                sigma_s_sd * sigma_s_raw + 
                sigma_q_sd * sigma_q_raw[q]);
  }
}

model{
  vector[N] p_n;
  vector[N] p_y;
  real utility_n;
  real utility_y;
  // decision
  // decision
  if (deltaD_value == 8){
    deltaD_mu_raw[1] ~ normal(0, 2);
    deltaD_s_sd[1] ~ normal(0, 1)T[0,];
    deltaD_q_sd[1] ~ normal(0, 1)T[0,];
    deltaD_s_raw ~ std_normal();
    deltaD_q_raw ~ std_normal();
  }
  
  p_end_mu_raw ~ normal(-1.2,.3);
  p_end_s_sd ~ normal(0,0.3)T[0,];
  p_end_q_sd ~ normal(0,0.3)T[0,];
  p_end_s_raw ~ normal(p_end_mu_raw,p_end_s_sd);
  p_end_q_raw ~ std_normal();

  sigma_mu_raw ~ normal(0, 2);
  sigma_s_sd ~ normal(0, 2)T[0,];
  sigma_q_sd ~ normal(0, 2)T[0,];
  sigma_s_raw ~ std_normal();
  sigma_q_raw ~ std_normal();
  
  if (deltaD_value == 9){
      for (n in 1:N){
        if (tNo[n] == 1){
          utility_n = rating_n[n];
          utility_y = rating_y[n];
          }else{
          utility_n *= deltaD[qID[n],sID[n]];
          utility_n += rating_n[n];
          utility_y *= deltaD[qID[n],sID[n]];
          utility_y += rating_y[n];
          }
       p_n[n] = 1 - normal_cdf(utility_n, 0, sigma[qID[n],sID[n]]);
       p_y[n] = normal_cdf(utility_y, 0, sigma[qID[n],sID[n]]);
      }
  }else{
      utility_n = 0;
      utility_y = 0;
      for (n in 1:N){
       p_n[n] = 1 - normal_cdf(rating_n[n], 0, sigma[qID[n],sID[n]]);
       p_y[n] = normal_cdf(rating_y[n], 0, sigma[qID[n],sID[n]]);
      }    
  }

  target += sum(log(decision_l (ar_value, -1, N_no, p_n[idx_no], p_y[idx_no])));
  target += sum(log(decision_l (ar_value, 1, N_yes, p_n[idx_yes], p_y[idx_yes])));
  // change this to multiplication
  for (n in 1:N){
    if(terminate[n] == 0){
      target += log(1 - p_end[qID[n],sID[n]]);
    }else{
      target += log(p_end[qID[n],sID[n]]);
    }
  }
}
generated quantities {
  real deltaD_mu[deltaD_value==8];
  real p_end_mu;
  real sigma_mu;
  
  vector[S*(deltaD_value==8)] deltaD_smean;
  vector[Q*(deltaD_value==8)] deltaD_qmean;
  vector[Q] p_end_qmean;
  vector[S] p_end_smean;
  vector[Q] sigma_qmean;
  vector[S] sigma_smean;
  
  vector[N] p_n;
  vector[N] p_y;
  real utility_n;
  real utility_y;
  vector[N] log_lik;
  real dev_d;
  vector[N] decision_prob[3];
  
  if(deltaD_value == 8){
    for (s in 1:S){
      deltaD_smean[s] = mean(to_array_1d(deltaD[,s]));
    }
    for (q in 1:Q){
      deltaD_qmean[q] = mean(to_array_1d(deltaD[q]));
    }
      deltaD_mu[1] = mean(to_array_1d(deltaD_smean)); 
  }
  
  for (s in 1:S){
    p_end_smean[s] = mean(to_array_1d(p_end[,s]));
    sigma_smean[s] = mean(to_array_1d(sigma[,s]));
  }
  for (q in 1:Q){
    p_end_qmean[q] = mean(to_array_1d(p_end[q]));
    sigma_qmean[q] = mean(to_array_1d(sigma[q]));
  }
  
  p_end_mu = mean(p_end_smean);
  sigma_mu = exp(sigma_mu_raw + sigma_q_sd^2/2 + sigma_s_sd^2/2);    

    if (deltaD_value == 9){
      for (n in 1:N){
        if (tNo[n] == 1){
          utility_n = rating_n[n];
          utility_y = rating_y[n];
          }else{
          utility_n *= deltaD[qID[n],sID[n]];
          utility_n += rating_n[n];
          utility_y *= deltaD[qID[n],sID[n]];
          utility_y += rating_y[n];
          }
       p_n[n] = 1 - normal_cdf(utility_n, 0, sigma[qID[n],sID[n]]);
       p_y[n] = normal_cdf(utility_y, 0, sigma[qID[n],sID[n]]);
      }
  }else{
      utility_n = 0;
      utility_y = 0;
      for (n in 1:N){
       p_n[n] = 1 - normal_cdf(rating_n[n], 0, sigma[qID[n],sID[n]]);
       p_y[n] = normal_cdf(rating_y[n], 0, sigma[qID[n],sID[n]]);
      }    
  }
  decision_prob[1] = decision_l (ar_value, -1, N, p_n, p_y);
  decision_prob[3] = decision_l (ar_value, 1, N, p_n, p_y);
  dev_d = 0;
  for (n in 1:N){
    decision_prob[1,n] *= p_end[qID[n],sID[n]];
    decision_prob[2,n] = (1 - p_end[qID[n],sID[n]]);
    decision_prob[3,n] *= p_end[qID[n],sID[n]];
    log_lik[n] = log(decision_prob[terminate[n]+2, n]);
    dev_d += - 2*log_lik[n];
  }
}

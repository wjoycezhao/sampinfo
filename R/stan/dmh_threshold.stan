// with deltaM; when deltaM is constrained to a fixed number x, lpdf = log(1/(1-1))=0 see binary -> threshold /3; random threshold/10
// https://groups.google.com/forum/#!topic/stan-users/pIEJOTEtbAU
// generate prediction for support distribution

functions{
// 
      vector decision_l (int ar_value, int terminate, int N, vector p_n, vector p_y){
            //variable to hold prob
          vector[N] prob_d; 
          if(ar_value == 0){
              if (terminate == 0 ){//continue
                  prob_d =  (1-p_n) .* (1-p_y);
              }
              else if (terminate == -1){
                  prob_d = p_n .* (1-p_y) + .5 * p_y .* p_n;
              }
              else if (terminate == 1){
                  prob_d = p_y .* (1-p_n) + .5 * p_y .* p_n ;
              }            
          }else{
              if (terminate == 0 ){//continue
                  prob_d =  (1 - p_n - p_y);
              }
              else if (terminate == -1){
                  prob_d = p_n;
              }
              else if (terminate == 1){
                  prob_d = p_y;
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
  real threshold_mu_raw;
  real<lower=0> threshold_s_sd;
  real<lower=0> threshold_q_sd;
  vector[S] threshold_s_raw;
  vector[Q] threshold_q_raw;
  real sigma_mult_mu_raw;
  real<lower=0> sigma_mult_s_sd;
  real<lower=0> sigma_mult_q_sd;
  vector[S] sigma_mult_s_raw;
  vector[Q] sigma_mult_q_raw;
}
transformed parameters {
  // decision
  vector[S] deltaD[Q*(deltaD_value==8)];
  vector[S] threshold[Q];
  vector[S] sigma_mult[Q];
  vector[S] sigma[Q];
  
  if (deltaD_value == 8){
    for (q in 1:Q){
      deltaD[q] = inv_logit(deltaD_mu_raw[1] + 
                      deltaD_s_sd[1] * deltaD_s_raw + 
                      deltaD_q_sd[1] * deltaD_q_raw[q]);
    }
  }
  
  for (q in 1:Q){
    threshold[q] = exp(threshold_mu_raw + 
                       threshold_s_sd * threshold_s_raw + 
                       threshold_q_sd * threshold_q_raw[q]);
  
    sigma_mult[q] = exp(sigma_mult_s_raw +
                        sigma_mult_q_sd * sigma_mult_q_raw[q]);
    sigma[q] = sigma_mult[q] .* threshold[q];
  }

}

model{
  vector[N] p_n;
  vector[N] p_y;
  real utility_n;
  real utility_y;
  
  // decision
  if (deltaD_value == 8){
    deltaD_mu_raw[1] ~ normal(0, 2);
    deltaD_s_sd[1] ~ normal(0, 1)T[0,];
    deltaD_q_sd[1] ~ normal(0, 1)T[0,];
    deltaD_s_raw ~ std_normal();
    deltaD_q_raw ~ std_normal();
  }
  
  threshold_mu_raw ~ normal(1,2);
  threshold_s_sd ~ normal(0,2)T[0,];
  threshold_q_sd ~ normal(0,2)T[0,];
  threshold_s_raw ~ std_normal();
  threshold_q_raw ~ std_normal();

  sigma_mult_mu_raw ~ normal(-0.7, 0.7);
  sigma_mult_s_sd ~ normal(0, 0.7)T[0,];
  sigma_mult_q_sd ~ normal(0, 0.7)T[0,];
  sigma_mult_s_raw ~ normal(sigma_mult_mu_raw,sigma_mult_s_sd);//!!
  sigma_mult_q_raw ~ std_normal();


  if (deltaD_value == 8){
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
       p_n[n] = 1 - normal_cdf(utility_n, -threshold[qID[n],sID[n]], sigma[qID[n],sID[n]]);
       p_y[n] = normal_cdf(utility_y, threshold[qID[n],sID[n]], sigma[qID[n],sID[n]]);
      }
  }else{
      utility_n = 0;
      utility_y = 0;
      for (n in 1:N){
       p_n[n] = 1 - normal_cdf(rating_n[n], -threshold[qID[n],sID[n]], sigma[qID[n],sID[n]]);
       p_y[n] = normal_cdf(rating_y[n], threshold[qID[n],sID[n]], sigma[qID[n],sID[n]]);
      }    
  }
  
  target += sum(log(decision_l (ar_value, 0, N_continue, p_n[idx_continue], p_y[idx_continue])));
  target += sum(log(decision_l (ar_value, -1, N_no, p_n[idx_no], p_y[idx_no])));
  target += sum(log(decision_l (ar_value, 1, N_yes, p_n[idx_yes], p_y[idx_yes])));
  
}
generated quantities {
  real deltaD_mu[deltaD_value==8];
  real threshold_mu;
  real sigma_mult_mu;
  real sigma_mu;
  
  vector[S*(deltaD_value==8)] deltaD_smean;
  vector[Q*(deltaD_value==8)] deltaD_qmean;
  vector[Q] threshold_qmean;
  vector[S] threshold_smean;
  vector[Q] sigma_mult_qmean;
  vector[S] sigma_mult_smean;
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
  
  
  threshold_mu = exp(threshold_mu_raw + threshold_q_sd^2/2 +threshold_s_sd^2/2); 
  sigma_mult_mu = exp(sigma_mult_mu_raw + sigma_mult_q_sd^2/2 + sigma_mult_s_sd^2/2);
  sigma_mu = sigma_mult_mu * threshold_mu;


  for (s in 1:S){
    threshold_smean[s] = mean(to_array_1d(threshold[,s]));
    sigma_mult_smean[s] = mean(to_array_1d(sigma_mult[,s]));
    sigma_smean[s] = mean(to_array_1d(to_vector(threshold[,s]) .* to_vector(sigma_mult[,s])));
  }
  for (q in 1:Q){
    threshold_qmean[q] = mean(to_array_1d(threshold[q]));
    sigma_mult_qmean[q] = mean(to_array_1d(sigma_mult[q]));
    sigma_qmean[q] = mean(to_array_1d(threshold[q] .* sigma_mult[q]));
  }  
 
    if (deltaD_value == 8){
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
       p_n[n] = 1 - normal_cdf(utility_n, -threshold[qID[n],sID[n]], sigma[qID[n],sID[n]]);
       p_y[n] = normal_cdf(utility_y, threshold[qID[n],sID[n]], sigma[qID[n],sID[n]]);
      }
  }else{
      utility_n = 0;
      utility_y = 0;
      for (n in 1:N){
       p_n[n] = 1 - normal_cdf(rating_n[n], -threshold[qID[n],sID[n]], sigma[qID[n],sID[n]]);
       p_y[n] = normal_cdf(rating_y[n], threshold[qID[n],sID[n]], sigma[qID[n],sID[n]]);
      }    
  }
  decision_prob[1] = decision_l (ar_value, -1, N, p_n, p_y);
  decision_prob[2] = decision_l (ar_value, 0, N, p_n, p_y);
  decision_prob[3] = decision_l (ar_value, 1, N, p_n, p_y);
  dev_d = 0;
  for (n in 1:N){
    log_lik[n] = log(decision_prob[terminate[n]+2, n]);
    dev_d += - 2*log_lik[n];
  }
}

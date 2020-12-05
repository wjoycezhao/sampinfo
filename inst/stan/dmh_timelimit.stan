
functions{
//
      vector[] decision_l (int N, real ar_value,
                          vector p_n, vector p_y, vector p_end){
            //variable to hold prob
          vector[N] prob_d[3];
          // to say no
          prob_d[1]=  p_end .* (ar_value * (p_n +  .5 * (1-p_n-p_y)) +
                              (1-ar_value)*(p_n .* (1-p_y) + .5*p_y .* p_n + .5*(1-p_n) .* (1-p_y)));
          // to continue;
          prob_d[2] =  (1 - p_end);

          // to say yes
          prob_d[3] = p_end .* (ar_value * (p_y + .5*(1-p_n-p_y)) +
                              (1-ar_value)*(p_y .* (1-p_n) + .5*p_y .* p_n + .5*(1-p_n) .* (1-p_y)));
          return prob_d;
}
}


data{
  // decision
  int deltaD_value; //8 for flexible decay parameter
  int ar_value; //0 if absolute accumulation; 1 if relative accumulation
  int<lower=1> N; //number of time points
  int<lower=1> S; // number of participants
  int<lower=1> Q; // number of questions
  int<lower=0> sID[N]; // participant ID variable
  int<lower=0> qID[N]; // question ID variable
  int<lower=1> tNo[N]; // thought No.x; used to reset representations
  int rating_n[N]; // supports for no (accumulatied supports if deltaD_value!=8)
  int rating_y[N]; // supports for yes (accumulatied supports if deltaD_value!=8)
  int terminate[N]; // 0 for continue; 1 for yes; -1 for no
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
  vector[N] decision_p[3];
  vector[N] p_end_qs;
  // decision
  if (deltaD_value == 8){
    deltaD_mu_raw[1] ~ normal(0, 1);
    deltaD_s_sd[1] ~ normal(0, 1)T[0,];
    deltaD_q_sd[1] ~ normal(0, 1)T[0,];
    deltaD_s_raw ~ std_normal();
    deltaD_q_raw ~ std_normal();
  }

  p_end_mu_raw ~ normal(0,1);
  p_end_s_sd ~ normal(0,1)T[0,];
  p_end_q_sd ~ normal(0,1)T[0,];
  p_end_s_raw ~ normal(p_end_mu_raw,p_end_s_sd);
  p_end_q_raw ~ std_normal();

  sigma_mu_raw ~ normal(0, 1);
  sigma_s_sd ~ normal(0, 1)T[0,];
  sigma_q_sd ~ normal(0, 1)T[0,];
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

  for (n in 1:N){
    p_end_qs[n] = p_end[qID[n],sID[n]];
  }
  decision_p = decision_l (N, ar_value, p_n, p_y, p_end_qs);
  for (n in 1:N){
      // decision
      target += log(decision_p[terminate[n]+2,n]);
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
  vector[N] decision_p[3];
  vector[N] p_end_qs;

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

  dev_d = 0;
  for (n in 1:N){
    p_end_qs[n] = p_end[qID[n],sID[n]];
  }
  decision_p = decision_l (N, ar_value, p_n, p_y, p_end_qs);
  for (n in 1:N){
    log_lik[n] = log(decision_p[terminate[n]+2,n]);
    dev_d += - 2*log_lik[n];
  }
}

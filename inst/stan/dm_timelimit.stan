// with deltaM; when deltaM is constrained to a fixed number x, lpdf = log(1/(1-1))=0
// https://groups.google.com/forum/#!topic/stan-users/pIEJOTEtbAU
// generate prediction for support distribution

functions{
//
      vector[] decision_l (int N, real ar_value,
                          vector terminate, vector p_n, vector p_y, vector p_end){
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
}


data{
  // decision
  real ar_value;
  int deltaD_value;
  int<lower=1> N;
  int tNo[N];
  int terminate[N];
  int rating_y[N];
  int rating_n[N];
}
parameters{
  // decision
  real<lower=0, upper=1> deltaD [deltaD_value == 8];
  real<lower=0,upper=1> p_end;
  real<lower=0> sigma_raw;
}
transformed parameters {
  // decision
  vector[N] p_n;
  vector[N] p_y;
  real<lower=0> sigma;
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
  if(deltaD_value == 8){
    deltaD[1] ~ uniform(0,1);
  }
  sigma_raw ~ normal(0, 1)T[0,];
  decision_p = decision_l (N, ar_value, to_vector(terminate), p_n, p_y, to_vector(rep_array(p_end,N)));
  for (n in 1:N){
      // decision
      target += log(decision_p[terminate[n]+2,n]);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] decision_p[3];
  vector[N] decision_0;
  vector[N] decision_no;
  vector[N] decision_yes;
  real dev_d;
  dev_d = 0;
  decision_p = decision_l (N, ar_value, to_vector(terminate), p_n, p_y, to_vector(rep_array(p_end,N)));
  decision_no =  decision_p[1];
  decision_0 = decision_p[2];
  decision_yes =  decision_p[3];
  for (n in 1:N){
      log_lik[n] = log(decision_p[terminate[n]+2,n]);
      dev_d += - 2*log_lik[n];
  }
}

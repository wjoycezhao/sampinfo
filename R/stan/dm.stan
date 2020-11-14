// with deltaM; when deltaM is constrained to a fixed number x, lpdf = log(1/(1-1))=0
// https://groups.google.com/forum/#!topic/stan-users/pIEJOTEtbAU
// generate prediction for support distribution

functions{
// 
      vector[] decision_l (int N, real ar_value, real random_value,
                          vector terminate, vector p_n, vector p_y, vector p_end){
            //variable to hold prob
          vector[N] prob_d[3]; 
          // to say no
          prob_d[1]=  p_end .* (ar_value * (p_n + random_value * .5 * (1-p_n-p_y)) + 
                              (1-ar_value)*(p_n .* (1-p_y) + .5*p_y .* p_n + random_value*.5*(1-p_n) .* (1-p_y)));
          // to continue
          prob_d[2] = random_value * (1 - p_end) + 
                      (1-random_value) * ar_value*(1-p_y-p_n) + 
                      (1-random_value) * (1-ar_value) * (1-p_n) .* (1-p_y);
          // to say yes
          prob_d[3] = p_end .* (ar_value * (p_y + random_value*.5*(1-p_n-p_y)) + 
                              (1-ar_value)*(p_y .* (1-p_n) + .5*p_y .* p_n + random_value*.5*(1-p_n) .* (1-p_y)));
              // if (terminate == 0 ){//continue
              //     if (random_value == 1){prob_d = 1 - p_end;}
              //     else if (ar_value == 0){prob_d = (1-p_n) * (1-p_y);}
              //     else if (ar_value == 1){prob_d = 1 - p_y - p_n;}
              // }
              // else if (terminate == -1){
              //     if (ar_value == 1){prob_d = p_n + random_value*.5*(1-p_n-p_y);}
              //     else if (ar_value == 0){prob_d = p_n*(1-p_y) + .5*p_y*p_n + random_value*.5*(1-p_n)*(1-p_y);}
              //     prob_d *= p_end;
              // }
              // else if (terminate == 1){
              //     if (ar_value == 1){prob_d = p_y + random_value*.5*(1-p_n-p_y);}
              //     else if (ar_value == 0){prob_d = p_y*(1-p_n) + .5*p_y*p_n + random_value*.5*(1-p_n)*(1-p_y);}
              //     prob_d *= p_end;
              // }
               // print(ar_value," ",terminate," ",p_y," ",p_n," ",prob_d);
            return prob_d;
}
}


data{
  // decision
  real ar_value;
  real random_value;
  int deltaD_value;
  int<lower=1> N;
  int tNo[N];
  int rating[N];
  int terminate[N];  
  int rating_y[N];
  int rating_n[N];
}
parameters{
  // decision
  real<lower=0, upper=1.5> deltaD [deltaD_value == 9];
  real<lower=0> threshold_raw;
  real<lower=0> sigma_mult;
  real<lower=0,upper=1> p_end_raw;
}
transformed parameters {
  // decision
  vector[N] p_n;
  vector[N] p_y; 

  real p_end;
  real sigma;
  real threshold;
  if (random_value == 1){
    threshold = 0;
    p_end = p_end_raw;
  }else{
    threshold = threshold_raw;
    p_end = 1;
  }
  if(random_value == 1){
    sigma = sigma_mult*10;
  }else{
    sigma = sigma_mult * threshold;
  }
  // print(sigma,threshold);
{
  real utility_n;
  real utility_y;
  if (deltaD_value == 9){
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
       p_n[n] = 1 - normal_cdf(utility_n, -threshold, sigma);
       p_y[n] = normal_cdf(utility_y, threshold, sigma);
      }
  }else{
      for (n in 1:N){
       p_n[n] = 1 - normal_cdf(rating_n[n], -threshold, sigma);
       p_y[n] = normal_cdf(rating_y[n], threshold, sigma);
      }    
  }
}
}


model{
  // decision
  vector[N] decision_p[3];
  if(deltaD_value == 9){
    deltaD[1] ~ normal(0.8,0.1)T[0,1.5];
  }
  threshold_raw ~ normal(15,6)T[0,];
  sigma_mult ~ normal(0.4,0.2)T[0,];
  decision_p = decision_l (N, ar_value, random_value, to_vector(terminate), p_n, p_y, to_vector(rep_array(p_end,N)));
  for (n in 1:N){
      // decision
      target += log(decision_p[terminate[n]+2,n]);    
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] decision_p[3];
  // vector[N] decision_0;
  // vector[N] decision_no;
  // vector[N] decision_yes;
  real dev_d;
  dev_d = 0;
  decision_p = decision_l (N, ar_value, random_value, to_vector(terminate), p_n, p_y, to_vector(rep_array(p_end,N)));
  for (n in 1:N){
      // decision
      // decision_0[n] = (decision_l (ar_value, random_value, 0, p_n[n], p_y[n] ,p_end));
      // decision_no[n] = (decision_l (ar_value, random_value, -1, p_n[n], p_y[n] ,p_end));
      // decision_yes[n] = (decision_l (ar_value, random_value, 1, p_n[n], p_y[n] ,p_end));
      log_lik[n] = log(decision_p[terminate[n]+2,n]);
      dev_d += - 2*log_lik[n];
  }
}

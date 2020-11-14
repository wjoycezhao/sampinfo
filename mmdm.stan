// with deltaM; when deltaM is constrained to a fixed number x, lpdf = log(1/(1-1))=0
// https://groups.google.com/forum/#!topic/stan-users/pIEJOTEtbAU
// generate prediction for support distribution



functions{
//
      vector[] theta_l (int M, int lastC, 
                      vector sameC_acc, vector sameA_acc, vector dist_acc, 
                      vector[] sameC_m, vector[] sameA_m, vector[] dist_m,
                      real deltaM,
                      vector b_int, real b_sameC, real b_sameA, real b_dist){
            vector[M] sameC_acc1;
            vector[M] sameA_acc1;
            vector[M] dist_acc1;
            vector[M] theta_raw;
            vector[M] theta;
            vector[M] returnList[4];
            sameC_acc1 =  sameC_acc * deltaM + sameC_m[lastC];
            sameA_acc1 =  sameA_acc * deltaM + sameA_m[lastC];
            dist_acc1 = dist_acc * deltaM + dist_m[lastC];
            theta_raw = b_int + sameC_acc1 * b_sameC  + sameA_acc1 * b_sameA + dist_acc1 * b_dist;
            theta = softmax(theta_raw); 
            returnList[1] = sameC_acc1;
            returnList[2] = sameA_acc1;
            returnList[3] = dist_acc1;
            returnList[4] = theta;
            return returnList;
      }
// 
      real memory_l (int MC,  
                      vector theta_cNo1,vector theta_cNo1_no,vector theta_cNo1_yes, 
                      vector theta, 
                      int cNo, int cID, real prime_value){
            real lprob_m; //variable to hold log prob
            lprob_m = 0;
            if (cNo == 1){
              if (prime_value == 0){
                  lprob_m += (theta_cNo1[cID]); 
              }else if (cID < MC+1){
                  lprob_m += (theta_cNo1_no[cID]); 
              }else {
                  lprob_m += (theta_cNo1_yes[cID]); 
              }
            }else{
                lprob_m += (theta[cID]); 
            }
            return lprob_m;
    }
//
      vector pypn_l (int n, int ar_value, real score0_n, real score0_y,
                    int[] score_n, int[] score_y, int[] cNo,
                    real sigma, real threshold, real deltaD){
          real utility_n;
          real utility_y;
          vector[2] p_no_yes;
          utility_n = score0_n * (deltaD^cNo[n]);
          utility_y = score0_y * (deltaD^cNo[n]);
          for (k in 0:(cNo[n]-1)){
            utility_n += score_n[n-k] * (deltaD^k);
            utility_y += score_y[n-k] * (deltaD^k);
          }
          p_no_yes[1] = normal_cdf(-threshold, utility_n, sigma); //say no
          p_no_yes[2] = 1 - normal_cdf(threshold, utility_y, sigma); //say yes
          // print(p_no_yes);
          return p_no_yes;
          }
// 
      real decision_l (int ar_value, int random_value,
                          int decision, real p_n, real p_y, real p_end){
            real prob_d; //variable to hold log prob
              if (decision == 0 ){//continue
                  if (random_value == 1){prob_d = 1 - p_end;}
                  else if (ar_value == 1){prob_d = 1 - p_y - p_n;}
                  else if (ar_value == 0){prob_d = (1-p_n) * (1-p_y);}
              }
              else{
                  if (decision == 1){
                      if (ar_value == 1){prob_d = p_y;}
                      else if (ar_value == 0){prob_d = p_y*(1-p_n) + .5*p_y*p_n + random_value*.5*(1-p_n)*(1-p_y);}
                    }
                  else{
                      if (ar_value == 1){prob_d = p_n;}
                      else if (ar_value == 0){prob_d = p_n*(1-p_y) + .5*p_y*p_n + random_value*.5*(1-p_n)*(1-p_y);}
                    }
                  prob_d = prob_d * p_end;
               } 
               // print(ar_value," ",decision," ",p_y," ",p_n," ",prob_d);
            return prob_d;
}
int sign(real x) {
  int x_sign;
  if(x<0){x_sign=-1;}
  if(x==0){x_sign = 0;}
  if(x>0){x_sign=1;}
  return x_sign;
  }
}

data{
  // shared
  int<lower=2> M;
  int<lower=1> MC;
  int<lower=1> N;
  int<lower=1> S;
  int cNo[N];
  int partNo[N];
  // memory
  real b_sameC_value;
  real b_sameA_value;
  real b_dist_value;
  real deltaM_value;
  real b_int_value;
  real prime_value;
  vector[M] sameC_m[M];
  vector[M] sameA_m[M];
  vector[M] dist_m[M];
  int cID[N];
  //decision
  int binary_value;
  int ar_value;
  int deltaD_value;
  int random_value;
  int score0_value;
  int score[N];
  vector[7] clusterScoreMatrix[M];
  int<lower=-1, upper=1> decision[N]; 
}
transformed data{
  // decision
    int score_y[N];
    int score_n[N];
  if (ar_value == 0){
    for (n in 1:N){
      score_n[n] = min(0,score[n]);
      score_y[n] = max(0,score[n]);
    }
  } else{
    score_y = score;
    score_n = score;
  }
}
parameters{
  // memory
  real<lower=0, upper=1> deltaM_raw;
  // matt trick
  vector[M-1] b_int_mu_raw0;
  real b_sameC_raw;
  real b_sameA_raw;
  real b_dist_raw;
  // decision
  real<lower=0, upper=1> deltaD_raw;
  // matt trick
  real score0_raw;
  real<lower=0> threshold_raw;
  real<lower=0> sigma_raw;
  real<lower=0,upper=1> p_end_raw;
}
transformed parameters {
  // memory
  real deltaM;
  vector[M] theta[N];
  vector[M] theta_cNo1;
  vector[M] theta_cNo1_yes;
  vector[M] theta_cNo1_no;
  real b_sameC;
  real b_sameA;
  real b_dist;
  // decision
  real p_y[N]; real p_n[N];
  real score0;
  real score0_n;
  real score0_y;
  real deltaD;
  real p_end;
  real sigma;
  real threshold;
  // scale memory
  vector[M-1] b_int_mu_raw;
  vector[M] b_int_mu;
  real<lower=0> b_int_sigma;
  if (b_sameC_value==0){
    b_sameC = b_sameC_value;
  }else{
    b_sameC = 5*b_sameC_raw;
  }
  if (b_sameA_value==0){
    b_sameA = b_sameA_value;
  }else{
    b_sameA = 5*b_sameA_raw;
  }
  if (b_dist_value==0){
    b_dist = b_dist_value;
  }else{
    b_dist = 5*b_dist_raw;
  }
  if (deltaM_value<9){
    deltaM=deltaM_value;
  }else{
    deltaM = deltaM_raw;
  }
  b_int_mu_raw = 2*b_int_mu_raw0;
  b_int_mu[1:(M-1)] = b_int_mu_raw;
  b_int_mu[M] = 0;
  b_int_sigma = 0;
  // scale decision
  if (score0_value==0){
    score0 = 0;
  }else{
    score0 = 15*score0_raw;
  }
  if (ar_value == 0){
    score0_n = fmin(0.0,score0);
    score0_y = fmax(0.0,score0);
  }
  else{
    score0_n = score0;
    score0_y = score0;
  }
  sigma = 35*sigma_raw;
  if (random_value == 1){
    threshold = 0;
    p_end = p_end_raw;
  }else{
    threshold = 35*threshold_raw;
    p_end = 1;
  }
  if (deltaD_value<9){
    deltaD=deltaD_value;
  }else{
    deltaD = deltaD_raw;
  }
{
  // memory
  vector[M] sameC_acc;
  vector[M] sameA_acc;
  vector[M] dist_acc;
  vector[M] returnList[4]; 
  vector[2] p_no_yes[N];
  theta_cNo1 = softmax(b_int_mu);
  theta_cNo1_yes = append_row(to_vector(rep_array(0,MC)),softmax(b_int_mu[(MC+1):M]));
  theta_cNo1_no = append_row(softmax(b_int_mu[1:(MC+1)]),to_vector(rep_array(0,MC)));
  for (n in 1:N){
      if (cNo[n] == 1){
            sameC_acc = to_vector(rep_array(0,M));
            sameA_acc = to_vector(rep_array(0,M));
            dist_acc = to_vector(rep_array(0,M));
            theta[n] = to_vector(rep_array(0,M));
     }else if (cNo[n] > 1){
          returnList = theta_l(M, cID[n-1],
                               sameC_acc, sameA_acc, dist_acc,
                               sameC_m, sameA_m, dist_m,
                               deltaM,
                               b_int_mu, b_sameC, b_sameA, b_dist);
          sameC_acc=returnList[1];
          sameA_acc=returnList[2];
          dist_acc=returnList[3];
          theta[n]=returnList[4];
      }
  }
  for (n in 1:N){
      p_no_yes[n] = pypn_l(n, ar_value, score0_n,score0_y, score_n, score_y, cNo,
                    sigma, threshold, deltaD);
      p_n[n] = p_no_yes[n,1];
      p_y[n] = p_no_yes[n,2];
  }
}
}

model{
  //memory
  deltaM_raw ~ uniform(0,1);
  b_sameC_raw ~ normal(0,1);
  b_sameA_raw ~ normal(0,1);
  b_dist_raw ~ normal(0,1);
  b_int_mu_raw0 ~ normal(0,1);
  //decision
  score0_raw ~ normal(0,1);
  deltaD_raw ~ uniform(0,1);
  threshold_raw ~ normal(0,1)T[0,];
  sigma_raw ~ normal(0,1)T[0,];
  //memory
  for (n in 1:N){
      target += log(memory_l(MC, theta_cNo1, theta_cNo1_no,theta_cNo1_yes,theta[n], cNo[n], cID[n], prime_value));
  }
  //decision
  for (n in 1:N){
      // decision
      target += log(decision_l (ar_value, random_value, decision[n], p_n[n], p_y[n] ,p_end));    
  }
}
generated quantities {
  // DICs
  vector[N] log_lik_m;
  real dev_m;
  vector[N] log_lik_d;
  real dev_d;
  // simulations
  int x_prd;
  int decision_prd;
  int score_prd[100];
  int cID_prd[100];
  x_prd = 100;
  // DICs
  dev_m = 0;
  dev_d = 0;
  for (n in 1:N){
      // memory
      log_lik_m[n] = log(memory_l(MC, theta_cNo1, theta_cNo1_no,theta_cNo1_yes,theta[n], cNo[n], cID[n], prime_value));
      dev_m += - 2*log_lik_m[n];
      // decision
      log_lik_d[n] = log(decision_l (ar_value, random_value, decision[n], p_n[n], p_y[n] ,p_end));
      dev_d += - 2*log_lik_d[n];
  }
 // simulations
  cID_prd = rep_array(100,x_prd);
  {   int score_n_prd[x_prd];
      int score_y_prd[x_prd];
      vector[2] p_no_yes_pred;
      real p_n_prd;
      real p_y_prd;
      int decisionMade;
      int i;
      vector[M] sameC_acc_prd;
      vector[M] sameA_acc_prd;
      vector[M] dist_acc_prd;
      vector[M] theta_prd;
      vector[M] returnList[4];
      vector[3] decision_p_prd;
      i= 1;
      decisionMade= 0;
      decision_prd = 3;
      while (!decisionMade){
           if(i==1){
                  cID_prd[i] = categorical_rng(theta_cNo1);
                  sameC_acc_prd = to_vector(rep_array(0,M));
                  sameA_acc_prd = to_vector(rep_array(0,M));
                  dist_acc_prd  = to_vector(rep_array(0,M));
           }else{
                  returnList = theta_l(M, cID_prd[i-1],
                       sameC_acc_prd, sameA_acc_prd, dist_acc_prd,
                       sameC_m, sameA_m, dist_m,
                       deltaM, 
                       b_int_mu, b_sameC, b_sameA, b_dist);
                  sameC_acc_prd=returnList[1];
                  sameA_acc_prd=returnList[2];
                  dist_acc_prd=returnList[3];
                  theta_prd=returnList[4];
                  cID_prd[i] = categorical_rng(theta_prd);
            }
            score_prd[i] = categorical_rng(clusterScoreMatrix[cID_prd[i]])-4;
            // print(score_prd[i]);
            if (binary_value == 1) {score_prd[i] = sign(score_prd[i]);}
            // print(score_prd[i]);
          if (ar_value == 0){
              score_n_prd[i] = min(score_prd[i],0);
              score_y_prd[i] = max(score_prd[i],0);
          }else{
              score_n_prd[i] = score_prd[i];
              score_y_prd[i] = score_prd[i];
          }
          p_no_yes_pred = pypn_l(i, ar_value, score0_n,score0_y,
                                score_n_prd, score_y_prd, cNo,
                               sigma, threshold, deltaD);
          p_n_prd = p_no_yes_pred[1];
          p_y_prd = p_no_yes_pred[2];

          decision_p_prd[1] = decision_l (ar_value,random_value,-1,
                                            p_n_prd, p_y_prd,p_end); //no
          decision_p_prd[2] = decision_l (ar_value,random_value,1,
                                            p_n_prd, p_y_prd,p_end); //yes
          decision_p_prd[3] = 1 - decision_p_prd[1] - decision_p_prd[2]; //ar_value
          decision_prd = categorical_rng(decision_p_prd); // 1:no 2:yes 3:ar_value
          if (decision_prd == 1){
              decision_prd = -1;
          }
          if (decision_prd == 2){
              decision_prd = 1;
          }
          i = i+1;
          if (decision_prd != 3 || i > x_prd){
              decisionMade = 1;
          }
      }
  }
}

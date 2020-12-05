functions{
//
      vector[] decision_l (int N, real ar_value,
                          vector p_n, vector p_y){
            //variable to hold prob
          vector[N] prob_d[3];
          // to say no
          prob_d[1]=  ar_value * p_n +
                      (1 - ar_value) * (p_n .* (1 - p_y) + .5 * p_y .* p_n);
          // to continue
          prob_d[2] = ar_value * (1 - p_y - p_n) +
                      (1 - ar_value) * (1 - p_n) .* (1 - p_y);
          // to say yes
          prob_d[3] = ar_value * p_y +
                      (1 - ar_value) * (p_y .* (1 - p_n) + .5 * p_y .* p_n);
          return prob_d;
      }
      int sign (int x) {
          int x_sign;
          if (x < 0){x_sign = -1;}
          if (x == 0){x_sign = 0;}
          if (x > 0){x_sign = 1;}
          return x_sign;
      }
}


data{
  int ar_value; //0 if absolute accumulation; 1 if relative accumulation
  int deltaM_value; //  8 for flexible memory decay parameter
  int deltaD_value; //  8 for flexible decision decay parameter
  int<lower=2> C; // number of option categories
  int<lower=2> MC; // ID of the neutral category
  int<lower=0> K; // number of features; if 0 then no betas
  matrix[C,K] X_sim[C]; // feature matrix, for categoreis
  int max_tNo_prd; // maximum length of a simulated trial
  real<lower=0, upper=1> deltaM [2]; // memory decay
  vector[C] alpha; // base rate
  vector[K] beta; // group-level mean for beta
  real<lower=0, upper=1> deltaD [2]; // decision decay
  real<lower=0> threshold; // decision threshold
  real<lower=0> sigma; // decision noise
}

parameters{
}

model{
}

generated quantities {
int rating_prd[max_tNo_prd];
int cID_prd[max_tNo_prd];
int terminate_prd;
vector[3] decision_p_prd1[max_tNo_prd];
int tNo_prd;

  // simulations
cID_prd = rep_array(99,max_tNo_prd);
rating_prd = rep_array(99,max_tNo_prd);
for (i in 1:max_tNo_prd){
  decision_p_prd1[i] = to_vector(rep_array(99,3));
}

{
  matrix[C,K] X_acc_prd;
  vector[C] theta_prd;
  int rating_n_prd;
  int rating_y_prd;
  real utility_n_prd;
  real utility_y_prd;
  real p_n_prd;
  real p_y_prd;
  vector[1] decision_p_prd[3];

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

    rating_prd[tNo_prd] = sign(cID_prd[tNo_prd] - MC);

    if (ar_value == 0){
      rating_n_prd = min(rating_prd[tNo_prd],0);
      rating_y_prd = max(rating_prd[tNo_prd],0);
    }else{
      rating_n_prd = rating_prd[tNo_prd];
      rating_y_prd = rating_prd[tNo_prd];
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
    p_n_prd = 1 - normal_cdf(utility_n_prd, -threshold, sigma);
    p_y_prd = normal_cdf(utility_y_prd, threshold, sigma);
    decision_p_prd = decision_l (1, ar_value, to_vector(rep_array(p_n_prd,1)), to_vector(rep_array(p_y_prd,1)));
    for(i in 1:3){decision_p_prd1[tNo_prd, i] = decision_p_prd[i,1];}
    terminate_prd = categorical_rng(decision_p_prd1[tNo_prd]) - 2;
    tNo_prd += 1;
  }
}
  tNo_prd += -1;
}

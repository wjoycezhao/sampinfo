max_tNo = 30;
cID_prd = rep_array(99,max_tNo);
{   int max_tNo;
  vector[C] theta_prd;
  matrix[C,K * (deltaM_value == 8)] X_acc_prd;
  vector[M] theta_prd;
  int rating_prd;
  int rating_n_prd;
  int rating_y_prd;
  real utility_n_prd;
  real utility_y_prd;
  real p_n;
  real p_y;
  vector[3] decision_p_prd;
  int terminate_prd;
  int tNo_prd;

  tNo_prd= 1;
  terminate_prd = 0;
  while (terminate_prd == 0 & tNo_prd < max_tNo + 1){
    if (K > 0 && deltaM_value == 8){
      if (tNo_prd[n] == 1){
        X_acc_prd = to_matrix(rep_array(0,C,K));
      }else{
        X_acc_prd *= deltaM[1];
        X_acc_prd += X_sim[cID_prd[tNo_prd - 1]];
      }
      theta_prd = X_acc_prd * beta;
      theta_prd += rep_matrix(alpha,N);
    }else if (K > 0 && deltaM_value == 0){
      if (tNo_prd[n] == 1){
        X_acc_prd = to_matrix(rep_array(0,C,K));
      }else{
        X_acc_prd = X_sim[cID_prd[tNo_prd - 1]];
      }
      theta_prd = X_acc_prd * beta;
      theta_prd += rep_matrix(alpha,N);
    }else if (K > 0 && deltaM_value == 1){
      if (tNo_prd[n] == 1){
        X_acc_prd = to_matrix(rep_array(0,C,K));
      }else{
        X_acc_prd += X_sim[cID_prd[tNo_prd - 1]];
      }
      theta_prd = X_acc_prd * beta;
      theta_prd += rep_matrix(alpha,N);
    }else{
      theta_prd = rep_matrix(alpha,N);
    }

    cID_prd[tNo_prd] ~ categorical_logit(theta_prd);

    rating_prd[i] = categorical_rng(cluster_rating_m[cID_prd[tNo_prd]])-4;

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
         p_n_prd = 1 - normal_cdf(utility_n_prd, -threshold, sigma);
         p_y_prd = normal_cdf(utility_y_prd, threshold, sigma);
    }else if (deltaD_value == 0){
            utility_n_prd = rating_n_prd;
            utility_y_prd = rating_y_prd;
         p_n_prd = 1 - normal_cdf(utility_n_prd, -threshold, sigma);
         p_y_prd = normal_cdf(utility_y_prd, threshold, sigma);
    }else if (deltaD_value == 1){
          if (tNo_prd == 1){
            utility_n_prd = rating_n_prd;
            utility_y_prd = rating_y_prd;
            }else{
            utility_n_prd += rating_n_prd;
            utility_y_prd += rating_y_prd;
            }
         p_n_prd = 1 - normal_cdf(utility_n_prd, -threshold, sigma);
         p_y_prd = normal_cdf(utility_y_prd, threshold, sigma);
    }
    decision_p_prd = decision_l (N, ar_value, to_vector(terminate), p_n_prd, p_y_prd);
    terminate_prd = categorical_rng(decision_p_prd) - 2;
    if (terminate_prd == 1){
      terminate_prd = -1;
    }
    if (terminate_prd == 2){
      terminate_prd = 1;
    }
    tNo_prd += 1;
    if (terminate_prd != 3 || i > max_tNo){
      terminate_prdMade = 1;
    }
  }
}


#' Simulate option
#'
#' Feature mattrix is dependent on the option sampled at the last time point
#' @param beta_true Values of beta vector (feature * 1)
#' @param delta_true Values of delta
#' @param alpha_true Values of alpha
#' @param feature_data feature matrix (option * feature);
#' dependent on the last chosen option
#' @param question_num totaal simulation number
#' @inheritParams  getStanFit
#' @export
#'
simulateChoice = function(beta_true,
                          delta_true,
                          alpha_true,
                          feature_data,
                          option_num,
                          question_num,
                          seed_no = 1) {
  set.seed(seed_no)
  sample_num = round(truncnorm::rtruncnorm(question_num, a=1, b=20, mean = 5, sd = 5))
  feature_num = length(beta_true)
  N = sum(sample_num)
  sID = unlist(lapply(1:question_num, function(x)
    rep(x, each = sample_num[x])))
  tNo = unlist(lapply(sample_num, function(x)
    1:x))
  cID = rep(NA, each = N)
  base_p = softmax(alpha_true)
  for (n in 1:N) {
    if (tNo[n] == 1) {
      acc = matrix(0, nrow = option_num, ncol = feature_num)
      cID[n] = sample(1:option_num, size = 1, p = base_p)
    } else{
      acc = acc * delta_true + feature_data[[cID[n - 1]]]
      prob = softmax(t(acc %*% beta_true + alpha_true))
      cID[n] = sample(1:option_num, size = 1, p = prob)
    }
  }
  return(data.frame(qID = 1, sID = sID, tNo = tNo, cID = cID))
}



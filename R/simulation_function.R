#' Model recovery
#'
#' \code{getSimulation} simulates data and fit memory and decision model to the simulated dataset.
#'
#' @inheritParams getStanFit
#' @param deltaM_true 1 minus decay for memory
#' @param alpha_true  baseline activations for memory
#' @param beta_true feature coefficients
#' @param deltaD_true  1 minus decay for decision
#' @param threshold_true  threshold for decision
#' @param sigma_true  noise for decision
#'
#' @export
#'
#'
getSimulation = function(cluster_num = 3,
                         sim_total = NULL,
                         X_sim = NULL,
                         deltaM_value = 8,
                         deltaD_value = 8,
                         ar_value = NULL,
                         max_tNo_prd = 30,
                         seed_no = seed_no,
                         deltaM_true = NULL,
                         alpha_true = NULL,
                         beta_true = NULL,
                         deltaD_true = NULL,
                         threshold_true = NULL,
                         sigma_true = NULL,
                         iter_num = NULL,
                         chain_num = NULL,
                         warmup_num = NULL,
                         core_num = NULL,
                         refresh = 500) {
  if (deltaM_value != 8) {
    deltaM = deltaM_value
  }
  if (deltaD_value != 8) {
    deltaD = deltaD_value
  }
  alpha_true = c(alpha_true, 0)

  ## simulate the dataset

  stan_data_sim = list(
    C = 2 * cluster_num + 1,
    MC = cluster_num + 1,
    K = length(beta_true),
    X_sim = X_sim,
    max_tNo_prd = max_tNo_prd,
    deltaM = c(deltaM_true, deltaM_true),
    alpha = alpha_true,
    beta = beta_true,
    deltaD = c(deltaD_true, deltaD_true),
    threshold = threshold_true,
    sigma = sigma_true,
    ar_value = ar_value,
    deltaM_value = deltaM_value,
    deltaD_value = deltaD_value
  )

  stan_fit = rstan::sampling(
    stanmodels$simulate,
    data = stan_data_sim,
    seed = seed_no,
    iter = sim_total,
    chains = 1,
    warmup = 0,
    cores = 1,
    refresh = refresh,
    algorithm = "Fixed_param"
  )

  ## computate median thought length in the dataset; and exclude datasets that are too long/short

  temp = median(unlist(rstan::extract(stan_fit, pars = c('tNo_prd'))))

  if (temp < 4 | temp > 6) {
    para_mm = NULL
    para_dm = NULL
  } else{

    # format the simulated dataset

    sim_data = as.data.frame(rstan::extract(stan_fit, pars = c('cID_prd'))) %>%
      mutate(sID = 1:sim_total) %>%
      mutate(sID = 1:sim_total) %>%
      tidyr::gather(tNo, cID, 'cID_prd.1':'cID_prd.30') %>%
      mutate(
        qID = 1,
        tNo = readr::parse_number(
          tNo,
          locale = readr::locale(grouping_mark = ". ", decimal_mark = ",")
        ),
        rating = sign(cID - cluster_num - 1)
      ) %>%
      filter(cID != 99) %>% arrange(sID, tNo) %>%
      group_by(sID) %>%
      mutate(terminate = ifelse(tNo == max(tNo), 1, 0)) %>% ungroup()
    sim_data$terminate[sim_data$terminate == 1] = unlist(rstan::extract(stan_fit, pars = c('terminate_prd')))
    sim_data = as.data.frame(sim_data)
    dist = data.frame(qID = 1, do.call(cbind, lapply(X_sim, function(x)
      x[, 3])))
    sim_data = getData(cluster_num = 3,
                       data = as.data.frame(sim_data),
                       dist = dist)[[1]]

    option_num = 2 * cluster_num + 1

    # fit the full memory model

    stan_data_fit = getStanFit(
      beta = c('sameC', 'sameA', 'dist'),
      deltaM_value = deltaM_value,
      option_num = option_num,
      format_data = sim_data,
      iter_num = iter_num,
      chain_num = chain_num,
      warmup_num = warmup_num,
      core_num = core_num,
      adapt_delta = 0.9,
      stepsize = 0.5,
      max_treedepth = 10,
      refresh = refresh,
      hier_value = 0
    )
    para_mm = getParaSummary(stan_data_fit = stan_data_fit)

    # fit the primary decision model

    stan_data_fit = getStanFit(
      binary_value = 1,
      ar_value = ar_value,
      random_value = 0,
      deltaD_value = deltaD_value,
      option_num = option_num,
      format_data = sim_data, #sim_data[,1:6],
      iter_num = iter_num,
      chain_num = chain_num,
      warmup_num = warmup_num,
      core_num = core_num,
      adapt_delta = 0.9,
      stepsize = 0.5,
      max_treedepth = 10,
      refresh = refresh,
      hier_value = 0
    )

    para_dm = getParaSummary(stan_data_fit = stan_data_fit)
  }
  return(list(para_mm = para_mm, para_dm = para_dm))
}






#' Fit the model!
#'
#' \code{getStanFit} uses RStan to fit a memory model or a decision model.
#' A memory model will be fit if deltaM_value is specified.
#' A decision model will be fit if deltaD_value is specified.
#' If both specified then the two models will be fit simultaneously, and posterior predictions will be simulated
#' Note that for memory models the option number (C) should be the same across all questions.
#'
#' @param beta Variable names to be included in the memory model, e.g., c('sameC','sameA','dist').
#' Note that you have to include the associated columns in the data.
#' For example, to fit sameC, you have to include columns with names of 'sameC_1'...'sameC_C'.
#' To fit sameA, you have to include 'sameA_1',...,'sameA_C'.
#' In the memory model, we could include sameC (cluster resampling; sameC_value== 1), sameA (decision congruence; sameA_value== 1),
#' and dist (semantic congruence; dist_value == 1)
#' @param deltaM_value 1 minus memory decay; 0 for full decay (marcov property); 1 for no decay;
#'  8 for estimating it as a free parameter
#' @param binary_value 0 for continuous inputs; 1 for discrete inputs
#' @param ar_value, 0 for absolute accumulation; 1 for relative accumulation
#' @param random_value, 0 for threhold-based;  1 for time limit;
#' @param deltaD_value, 1 minus decision decay; 0 for full decay (marcov property); 1 for no decay;
#'  8 for estimating it as a free parameter
#' @param option_num Number of total clusters in each choice (C)
#' @param format_data This should be a formated data.frame.
#' sID: participant IDs.
#' qID: question ID.
#' cID: option ID being sampled (can take 1 to C).
#' tNo: thought number in that trial (from 1 to T).
#' terminate: 0 if continue sampling; 1 if yes; 0 if no (only needed for decision models)
#' @param save_model_file Specify a name if you want to save the RStan samples and diagnostics to a csv file.
#' This is sample_file for rstan::sampling
#' @param init_values Default is 'random'; can be specified using a function or a list
#' @param iter_num Iteration number
#' @param chain_num Chain number
#' @param warmup_num Warm-up or burn-in number
#' @param core_num Number of cores to be used
#' @param seed_no Seed number
#' @param adapt_delta Default 0.999
#' @param stepsize Default 0.01
#' @param max_treedepth Default 20
#' @param refresh Default 500
#' @param hier_value 0 for fitting a non-hierarchical model; 1 for fitting a hierarchical model. Default is 1.
#' Note that the hierachical version only works for decision or memory models seperately.
#' @param X_sim Feature matrix, for different clusters. Needed for simulation.
#' @param cluster_rating_m Base rates of different ratings (-3 to 3) for each cluster. Needed for simulation.
#' @param max_tNo_prd Maximum length of a simulated trial
#'
#' @return Use $stan_fit to access the returned object from Stan;
#' Use $stan_data to see the input data of the Stan model
#' @export

getStanFit = function (beta = NULL,
                       deltaM_value = NULL,
                       binary_value = NULL,
                       ar_value = NULL,
                       random_value = NULL,
                       deltaD_value = NULL,
                       option_num,
                       format_data,
                       save_model_file =NULL,
                       init_values = "random",
                       iter_num,
                       chain_num,
                       warmup_num,
                       core_num,
                       adapt_delta = 0.999,
                       stepsize = 0.01,
                       max_treedepth = 20,
                       save_warmup = FALSE,
                       refresh = 500,
                       init_r = 0.5,
                       hier_value = 1,
                       X_sim = NULL,
                       cluster_rating_m = NULL,
                       max_tNo_prd = NULL) {
  ## prepare data for RStan
  stan_data = getStanData(
    beta = beta,
    deltaM_value = deltaM_value,
    binary_value = binary_value,
    ar_value = ar_value,
    random_value = random_value,
    deltaD_value = deltaD_value,
    option_num = option_num,
    format_data = format_data,
    X_sim = X_sim,
    cluster_rating_m = cluster_rating_m,
    max_tNo_prd = max_tNo_prd
  )

  ## get the appropriate rstan code
  stan_code = getStanCode(deltaM_value, deltaD_value, hier_value, random_value)

  ## fit model
  stan_fit = rstan::sampling(
    stan_code,
    data = stan_data,
    seed = seed_no,
    sample_file = save_model_file,
    init = init_values,
    iter = iter_num,
    chains = chain_num,
    warmup = warmup_num,
    cores = core_num,
    refresh = refresh,
    save_warmup = save_warmup,
    init_r = init_r,
    control = list(
      adapt_delta = adapt_delta,
      stepsize = stepsize,
      max_treedepth = max_treedepth
    )
  )
  return(stan_data_fit = list(
    beta = beta,
    stan_data = stan_data,
    stan_fit = stan_fit
  ))
}


#' Retrieve stan code
#'
#' @inheritParams getStanFit
#' @export
#'
getStanCode = function(deltaM_value = NULL,
                       deltaD_value = NULL,
                       hier_value = NULL,
                       random_value = NULL) {
  if(!is.null(deltaM_value) & is.null(deltaD_value)){
    if (hier_value == 0){
      stan_code = stanmodels$mm
    }else if (hier_value == 1){
      stan_code = stanmodels$mmh
    }
  }else if(is.null(deltaM_value) & !is.null(deltaD_value)){
    if (random_value == 0) {
      if (hier_value == 0){
        stan_code = stanmodels$dm_threshold
      }else if (hier_value == 1){
        stan_code = stanmodels$dmh_threshold
      }
    } else if (random_value == 1){
      if (hier_value == 0){
        stan_code = stanmodels$dm_timelimit
      }else if (hier_value == 1){
        stan_code = stanmodels$dmh_timelimit
      }
    }
  }else if(!is.null(deltaM_value) & !is.null(deltaD_value)){
    if (random_value == 0) {
      if (hier_value == 0){
        stan_code = stanmodels$mmdm_threshold
      }
    } else if (random_value == 1){
      if (hier_value == 0){
        stan_code = stanmodels$mmdm_timelimit
      }
    }
  } else{
    stop('Cannot fit this model. Please check deltaM_value, deltaD_value, random_value, and hier_value.')
  }
  return(stan_code)
}


#' Format data for Stan
#'
#' @inheritParams getStanFit
#'
#' @export

getStanData = function(beta = character(0),
                       deltaM_value = NULL,
                       binary_value = NULL,
                       ar_value = NULL,
                       random_value = NULL,
                       deltaD_value = NULL,
                       option_num = 3 * 2 + 1,
                       format_data,
                       X_sim = NULL,
                       cluster_rating_m = NULL,
                       max_tNo_prd = NULL) {
  ## memory model
  if (min(format_data$sID) < 1 |
      length(unique(format_data$sID)) != max(format_data$sID)) {
    stop(
      'Participant ID should start from 1, and be consecutive; otherwise it will cause error when sampling.'
    )
  }
  if (min(format_data$qID) < 1 |
      length(unique(format_data$qID)) != max(format_data$qID)) {
    stop(
      'Question ID should start from 1, and be consecutive; otherwise it will cause error when sampling'
    )
  }
  format_data = dplyr::arrange(format_data, qID, sID, tNo)
  for (i in 2:nrow(format_data)) {
    if ((format_data$sID[i] == format_data$sID[i - 1] &
         format_data$tNo[i] != format_data$tNo[i - 1] + 1) |
        (format_data$sID[i] != format_data$sID[i - 1] &
         format_data$tNo[i] != 1)) {
      stop('Thought index should start from 1 and be consecutive in each trial; otherwise it will cause error when sampling')
    }
  }


  # feature matrix for fitting models and for simulations

  if (length(beta) > 0) {
    X = getFeatureMatrices(beta, format_data, deltaM_value,
                           option_num, print_names = F)
    X_sim_sel = lapply(X_sim, function(x) x[,beta])
    X_sim_sel = lapply(X_sim_sel, function(x) as.matrix(x))
  } else{
    X = replicate(nrow(format_data), matrix(0, option_num, 0), simplify = FALSE)
    X_sim_sel = lapply(1:option_num, function(x) matrix(0,nrow=option_num,ncol=0))
  }

  # variables used for fitting decision model

  rating = as.numeric(as.character(format_data$rating))
  rating_y = c()
  rating_n = c()
  if (!is.null(binary_value)){
    # make ratings binary for binary_value = 1
    if(binary_value == 1){
      rating = sign(rating)
    }
    # get ratings for yes and no accumulators
    if (ar_value == 0){
      # absolute accumulation
      rating_n = pmin(0, rating);
      rating_y = pmax(0, rating);
    } else{
      # relative accumulation
      rating_y = rating;
      rating_n = rating;
    }

  # accumulte the ratings if the model has no decay

    if (deltaD_value == 1) {
      for (n in 1:length(rating)) {
        if (format_data$tNo[n] > 1) {
          rating[n] = rating[n-1] + rating[n]
          rating_n[n] = rating_n[n-1] + rating_n[n]
          rating_y[n] = rating_y[n-1] + rating_y[n]
        }
      }
    }
  }

  # indicate whether participants keep sampling (0) or respond with yes (1) or no (-1)

  terminate = as.numeric(as.character(format_data$terminate))

  if(!is.null(cluster_rating_m)){
    cluster_rating_m = as.matrix(cluster_rating_m)
    if(ncol(cluster_rating_m) != 7){
      stop('Wrong cluster_rating_m matrix. Should have 7 columns to indicate ratings from -3 to 3.')
    }
  }

  stan_data <- list(
    ## specifications for memory model
    deltaM_value = deltaM_value,
    condition_value = ifelse(is.null(format_data$condition), 0, 1), # 1 for priming experiments
    C = option_num,
    MC = (option_num + 1)/2, # the neutral cluster
    N = nrow(format_data),
    K = length(beta),
    S = length(unique(format_data$sID)),
    Q = length(unique(format_data$qID)),
    sID =  as.numeric(as.character(format_data$sID)),
    qID =  as.numeric(as.character(format_data$qID)),
    tNo = as.numeric(as.character(format_data$tNo)),
    cID = as.numeric(as.character(format_data$cID)),
    X = X,
    ## specifications for decision model
    binary_value = binary_value,
    ar_value = ar_value,
    random_value = random_value,
    deltaD_value = deltaD_value,
    rating = rating,
    rating_n = rating_n,
    rating_y = rating_y,
    terminate = terminate,
    ## used for simulations
    X_sim = X_sim_sel,
    cluster_rating_m = cluster_rating_m,
    max_tNo_prd = max_tNo_prd
  )

  return(stan_data)
}


#' Extract features from formatted data
#'
#' @param print_names whether to include column and row namaes; default FALSE
#' @inheritParams getStanFit
#'
#' @export

getFeatureMatrices = function(beta,
                              format_data,
                              deltaM_value,
                              option_num,
                              print_names = F) {
  temp = lapply(beta, function(x)
    as.matrix(format_data[, paste0(x, '_', 1:option_num)]))
  X = lapply(1:nrow(format_data), function(x) {
    sapply(temp, function(y)
      as.numeric(y[x,]))
  })

  # accumulte the features if the model has no decay

  if (deltaM_value == 1) {
    for (i in 1:length(X)) {
      if (format_data$tNo[i] > 1) {
        X[[i]] = X[[i - 1]] + X[[i]]
      }
    }
  }
  if (print_names) {
    X = lapply(X, function(x) {
      y = x
      colnames(y) = beta
      rownames(y) = paste0('Option', 1:option_num)
      return(y)
    })
  }
  return(X)
}


#' Change variable names
#'
#' @param names A vector of characters. Variables taken from rstan outputs.
#' @inheritParams getStanFit
#' @param Q Number of questions
#' @param S Number of participants
#' @export
changeBetaNames = function(names, beta, Q, S) {
  for (x in 1:length(beta)) {
    names = plyr::mapvalues(
      names,
      warn_missing = FALSE,
      from = c(
        'deltaM.1.',
        'deltaD.1.',
        paste0('alpha.',1:9,'.'),
        paste0('beta_mu.', x,'.'),
        paste0('beta.', x,'.'),
        'beta_mu',
        paste0('beta_q_sd.', x,'.'),
        paste0('beta_s_sd.', x,'.'),
        paste0('beta_qmean.', x, '.', 1:Q,'.'),
        paste0('beta_smean.', x, '.', 1:S,'.'),
        paste0('beta_qdev.', x, '.', 1:Q,'.'),
        paste0('beta_sdev.', x, '.', 1:S,'.'),
        'beta_q_sd',
        'beta_s_sd',
        unlist(lapply(1:Q, function(y)
          paste0('beta.', y, '.', x, '.', 1:S)))
      ),
      to = c(
        'deltaM',
        'deltaD',
        paste0('alpha.',1:9),
        paste0('beta_', beta[x], '_mu'),
        paste0('beta_', beta[x]),
        paste0('beta_', beta[1], '_mu'),
        paste0('beta_', beta[x], '_q_sd'),
        paste0('beta_', beta[x], '_s_sd'),
        paste0('beta_', beta[x], '_qmean_', 1:Q),
        paste0('beta_', beta[x], '_smean_', 1:S),
        paste0('beta_', beta[x], '_qdev_', 1:Q),
        paste0('beta_', beta[x], '_sdev_', 1:S),
        paste0('beta_', beta[1], '_q_sd'),
        paste0('beta_', beta[1], '_s_sd'),
        unlist(lapply(1:Q, function(y)
          paste0('beta_', beta[x], '_', y, '_', 1:S)))
      )
    )
    names = plyr::mapvalues(
      names,
      warn_missing = FALSE,
      from = c(
        paste0('beta_mu[', x, ']'),
        paste0('beta[', x, ']'),
        paste0('beta_q_sd[', x, ']'),
        paste0('beta_s_sd[', x, ']'),
        paste0('beta_qmean[', x, ',', 1:Q, ']'),
        paste0('beta_smean[', x, ',', 1:S, ']'),
        paste0('beta_qdev[', x, ',', 1:Q, ']'),
        paste0('beta_sdev[', x, ',', 1:S, ']'),
        unlist(lapply(1:Q, function(y)
          paste0('beta.', y, '.', x, '.', 1:S)))
      ),
      to = c(
        paste0('beta_', beta[x], '_mu'),
        paste0('beta_', beta[x]),
        paste0('beta_', beta[x], '_q_sd'),
        paste0('beta_', beta[x], '_s_sd'),
        paste0('beta_', beta[x], '_qmean_', 1:Q),
        paste0('beta_', beta[x], '_smean_', 1:S),
        paste0('beta_', beta[x], '_qdev_', 1:Q),
        paste0('beta_', beta[x], '_sdev_', 1:S),
        unlist(lapply(1:Q, function(y)
          paste0('beta_', beta[x], '_', y, '_', 1:S)))
      )
    )
  }
  return(names)
}

#' Parameter summary.
#'
#' \code{getParaSummary} returns parameter summaries for group and individual level parameters,
#' including means, percentiles, and sd.
#' It also saves the results to a .csv file when \code{csv_name} is specified
#'
#' @param stan_data_fit model outputs obtained from getStanFit
#' @param csv_name if specified then outputs will be saved in a .csv file starting with this name
#'
#' @export
getParaSummary = function(stan_data_fit, csv_name = NULL) {
  stan_data = stan_data_fit$stan_data
  stan_fit = stan_data_fit$stan_fit
  para_name = colnames(as.data.frame(stan_data_fit$stan_fit))
  beta = stan_data_fit$beta
  parameters0 = para_name[stringr::str_detect(para_name,
                          'delta|beta|alpha|p_end|sigma|threshold')]
  parameters0 = parameters0[!stringr::str_detect(parameters0,
                          'raw')]
  parameters0 = c(parameters0, para_name[stringr::str_detect(para_name,
                                        'threshold_mu_raw|sigma_mult_mu_raw')])

  df = t(sapply(as.data.frame(
    rstan::extract(stan_fit, pars = parameters0)
  ), quantile1))
  rownames(df) = changeBetaNames(rownames(df), beta, stan_data$Q, stan_data$S)
  rownames(df) = changeBetaNames(rownames(df), beta, stan_data$Q, stan_data$S)
  if (!is.null(csv_name)) {
    write.table(
      df,
      csv_name,
      sep = ",",
      append = FALSE,
      quote = FALSE,
      col.names = NA,
      row.names = TRUE
    )
  }
  return(df)
}


#' Model diagnostics.
#'
#' \code{getModelDiag} returns model diagnostics,
#' and save the results to a .csv file when \code{csv_name} is specified
#'
#' @param dev_name string. name for dev in rstan codes
#' @inheritParams getParaSummary
#' @inheritParams getWAIC
#'
#' @return A dataframe indicating DIC, WAIC, model running time,
#' maximum rhat, minimum effective sample size,
#' number of divergent transitions and number of max treedepth hit.
#' @export

getModelDiag = function(stan_data_fit,
                        csv_name = NULL,
                        dev_name = 'dev_m',
                        logName = 'log_lik') {
  stan_fit = stan_data_fit$stan_fit
  time = getTime(stan_fit)
  rhat = getRhat(stan_fit)
  ess = getESS(stan_fit)
  ess = ess[!grepl('X_acc', names(ess))]
  df = cbind.data.frame(
    getDIC(stan_fit, dev_name),
    getWAIC(stan_fit, logName = logName),
    elapsed_time_min = min(time),
    elapsed_time_max = max(time),
    rhat_max = max(rhat, na.rm = T),
    ess_min = min(ess, na.rm = T),
    divergent = rstan::get_num_divergent(stan_fit),
    max_tree = rstan::get_num_max_treedepth(stan_fit)
  )
  if (!is.null(csv_name)) {
    write.table(
      df,
      csv_name,
      sep = ",",
      append = FALSE,
      quote = FALSE,
      col.names = NA,
      row.names = TRUE
    )
  }
  return(df)
}

#' Parameter summaries
#'
#' \code{getPara} returns parameter mean, quantiles and sd
#' and save the results to a .csv file when \code{csv_name} is specified
#'
#' @param para_name parameter to extract
#' @param quantiles quantiles to extract
#' @inheritParams getParaSummary
#'
#' @return A dataframe with parameter mean, quantiles and sd, can be saved to a csv file
#' @export

getPara = function(para_name, stan_data_fit, csv_name = NULL,
                   quantiles = c(0.025,0.05,0.1,0.2,0.25,0.5,0.75,0.8,0.9,0.95,0.975)){
  pp_all = as.data.frame(rstan::extract(stan_data_fit$stan_fit,pars=para_name))
  temp = lapply(1:ncol(pp_all),function(x){quantile(pp_all[,x],quantiles)})
  pp_all1 = rbind(colMeans(pp_all),do.call(cbind,temp),sapply(pp_all,sd))
  colnames(pp_all1) = colnames(pp_all)
  rownames(pp_all1)[1] = 'mean'
  rownames(pp_all1)[2 + length(quantiles)] = 'sd'
  if (!is.null(csv_name)) {
    write.table(
      pp_all1,
      csv_name,
      sep = ",",
      append = FALSE,
      quote = FALSE,
      col.names = NA,
      row.names = TRUE
    )
  }
  return(pp_all1)
}


#' Fit the model!
#'
#' \code{getStanFit} use RStan to fit your information sampling model.
#' Note that the option number (C) should be the same across all questions.
#'
#' @param beta Variable names to be included in the model, e.g., c('beta','gamma').
#' Note that in the formatted data, columns 'beta_1'...'beta_C' and 'gamma_1',...,'gamma_C' should be specified
#' @param deltaM_value 0 for full decay (marcov property); 1 for no decay;
#'  9 for estimating it as a free parameter
#' @param option_num Number of total options in each choice (C)
#' @param format_data This should be a formated data.frame. sID: participant IDs. qID: question ID.
#' cID: option ID being sampled (can take 1 to C). tNo: thought number in that question.
#' terminate: 0 if continue sampling (not necessary)
#' @param save_model_file Specify a name if you want to save the RStan samples and diagnostics to a csv file
#' @param init_values Default is 'random'; can be specified using a function or a list
#' @param iter_num Iteration number
#' @param chain_num Chain number
#' @param warmup_num Warm-up or burn-in number
#' @param core_num Number of cores to be used. equal or less to chain_num
#' @param adapt_delta Default 0.999
#' @param stepsize Default 0.01
#' @param max_treedepth Default 20
#' @param refresh Default 500
#'
#' @return Use $stan_fit to access the returned object from Stan;
#' Use $stan_data to see the input data of the Stan model
#' @export

getStanFit = function (beta, deltaM_value, option_num, format_data,
                       save_model_file, init_values="random",
                       iter_num,chain_num,warmup_num,core_num,
                       adapt_delta=0.999, stepsize = 0.01, max_treedepth = 20,
                       save_warmup=FALSE, refresh=500, init_r=0.5){
  ## prepare inputs for RStan; getHMRStanData in 'model_functions.R'
  stan_data = getStanData(beta = beta,
                             deltaM_value = deltaM_value,
                             option_num = option_num,
                             format_data = format_data)
  ##fit model
  stan_code = getStanCode(deltaM_value)

  if (init_values!='random'){
    init_values = function() { list(deltaM=0.5) }
  }
  stan_fit = rstan::sampling(stan_code,
                  data = stan_data, seed=1,
                  sample_file = save_model_file,init=init_values,
                  iter=iter_num, chains=chain_num, warmup=warmup_num, cores=core_num,
                  refresh=refresh,save_warmup=save_warmup,init_r=init_r,
                  control=list(adapt_delta=adapt_delta,stepsize = stepsize, max_treedepth = max_treedepth))
  return(stan_data_fit = list(beta=beta,
                              stan_data = stan_data,
                              stan_fit = stan_fit))
}


#' Retrieve stan code
#'
#' There are two Rstan code files. One for flexible delta (deltaM_value=9);
#' the other for fixed deltaa (deltaM_value = 0 or deltaM_value = 1)
#' @inheritParams getStanFit
#' @export
#'
getStanCode = function(deltaM_value){
  if (deltaM_value == 9){
    stan_code = stanmodels$mm10d9
  }else{
    stan_code = stanmodels$mm10d01
  }
  return(stan_code)
}
#' Format data for Stan
#'
#' @inheritParams getStanFit
#'
#' @export

getStanData = function(beta = character(0), deltaM_value = NULL,
                          option_num = 3*2+1, format_data){
  if(min(format_data$sID)<1 | length(unique(format_data$sID)) != max(format_data$sID)){
    stop('Participant ID should start from 0, and be consecutive; otherwise it will cause error when sampling')}
  if(min(format_data$qID)<1 | length(unique(format_data$qID)) != max(format_data$qID)){
    stop('Question ID should start from 0, and be consecutive; otherwise it will cause error when sampling')}
  format_data = dplyr::arrange(format_data, qID, sID,tNo)
  for (i in 2:nrow(format_data)){
    if((format_data$sID[i]==format_data$sID[i-1] & format_data$tNo[i] !=format_data$tNo[i-1] + 1) |
       (format_data$sID[i]!=format_data$sID[i-1] & format_data$tNo[i] != 1)){stop('wrong thought index')}
  }

  # feature matrix
  if(length(beta>0)){
    X = getFeatureMatrices(beta, format_data, deltaM_value,
                           option_num, print_names = F)
  }else{
    X = replicate(nrow(format_data), matrix(0,option_num,0), simplify = FALSE)
  }

  stan_data <- list(deltaM_value = deltaM_value,
                    C = option_num, # cluster numberss for one of the options
                    N = nrow(format_data),
                    K = length(beta),
                    S = length(unique(format_data$sID)), # unique participant numbers
                    Q = length(unique(format_data$qID)), # unique participant numbers
                    sID =  as.numeric(as.character(format_data$sID)),
                    # qID = rep(1,each=nrow(format_data)),##!! to test for no question differences
                    qID =  as.numeric(as.character(format_data$qID)),
                    tNo = as.numeric(as.character(format_data$tNo)),
                    cID = as.numeric(as.character(format_data$cID)),
                    X = X)
  if(!is.null(format_data$terminate)){
    stan_data$terminate = as.numeric(as.character(format_data$terminate))
    }
  return(stan_data)
}


#' Title
#'
#' @param print_names whether to include column and row namaes; default FALSE
#' @inheritParams getStanFit
#'
#' @export

getFeatureMatrices = function(beta, format_data, deltaM_value,
                              option_num, print_names = F){
  temp = lapply(beta, function(x)
    as.matrix(format_data[,paste0(x,'_',1:option_num)]))
  X = lapply(1:nrow(format_data),function(x){
    sapply(temp, function(y) as.numeric(y[x,]))})
  if(deltaM_value == 1){
    for (i in 1:length(X)){
      if(format_data$tNo[i]>1){
        X[[i]] = X[[i-1]] + X[[i]]
      }
    }
  }
  if(print_names){
    X = lapply(X, function(x) {
      y = x
      colnames(y) = beta
      rownames(y) = paste0('Option',1:option_num)
      return(y)})
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
changeBetaNames = function(names, beta,Q,S){
  for (x in 1:length(beta)){
    names = plyr::mapvalues(names,
           warn_missing = FALSE,
           from = c(paste0('beta_mu.',x),
                    'beta_mu',
                    paste0('beta_q_sd.',x),
                    paste0('beta_s_sd.',x),
                    paste0('beta_qmean.',x,'.',1:Q),
                    paste0('beta_smean.',x,'.',1:S),
                    paste0('beta_qdev.',x,'.',1:Q),
                    paste0('beta_sdev.',x,'.',1:S),
                    'beta_q_sd',
                    'beta_s_sd',
                    unlist(lapply(1:Q,function(y) paste0('beta.',y,'.',x,'.',1:S)))),
           to = c(paste0('beta_',beta[x],'_mu'),
                  paste0('beta_',beta[1],'_mu'),
                  paste0('beta_',beta[x],'_q_sd'),
                  paste0('beta_',beta[x],'_s_sd'),
                  paste0('beta_',beta[x],'_qmean_',1:Q),
                  paste0('beta_',beta[x],'_smean_',1:S),
                  paste0('beta_',beta[x],'_qdev_',1:Q),
                  paste0('beta_',beta[x],'_sdev_',1:S),
                  paste0('beta_',beta[1],'_q_sd'),
                  paste0('beta_',beta[1],'_s_sd'),
                  unlist(lapply(1:Q,function(y) paste0('beta_',beta[x],'_',y,'_',1:S)))))
    names = plyr::mapvalues(names,
                           warn_missing = FALSE,
                           from = c(paste0('beta_mu[',x,']'),
                                    paste0('beta_q_sd[',x,']'),
                                    paste0('beta_s_sd[',x,']'),
                                    paste0('beta_qmean[',x,',',1:Q,']'),
                                    paste0('beta_smean[',x,',',1:S,']'),
                                    paste0('beta_qdev[',x,',',1:Q,']'),
                                    paste0('beta_sdev[',x,',',1:S,']'),
                                    unlist(lapply(1:Q,function(y) paste0('beta.',y,'.',x,'.',1:S)))),
                           to = c(paste0('beta_',beta[x],'_mu'),
                                  paste0('beta_',beta[x],'_q_sd'),
                                  paste0('beta_',beta[x],'_s_sd'),
                                  paste0('beta_',beta[x],'_qmean_',1:Q),
                                  paste0('beta_',beta[x],'_smean_',1:S),
                                  paste0('beta_',beta[x],'_qdev_',1:Q),
                                  paste0('beta_',beta[x],'_sdev_',1:S),
                                  unlist(lapply(1:Q,function(y) paste0('beta_',beta[x],'_',y,'_',1:S)))))
  }
  return(names)
}

#' Parameter summary.
#'
#' \code{getParaSummary} return parameter summaries for group and individual level parameters,
#' including means, percentiles, and sd.
#' It also saves the results to a .csv file when \code{csv_name} is specified
#'
#' @param stan_data_fit model outputs obtained from getStanFit
#' @param csv_name if specified then outputs will be saved in a .csv file starting with this name
#'
#' @export
getParaSummary = function(stan_data_fit, csv_name = NULL){
  stan_data = stan_data_fit$stan_data
  stan_fit = stan_data_fit$stan_fit
  beta = stan_data_fit$beta
  parameters0 = c('deltaM','beta_mu','beta_q_sd',
                  'beta_s_sd','beta_qmean','beta_smean','beta_qdev','beta_sdev',
                  'alpha','beta')
  df = t(sapply(as.data.frame(rstan::extract(stan_fit,pars=parameters0)),quantile1))
  rownames(df) = changeBetaNames(rownames(df), beta, stan_data$Q, stan_data$S)
  if (!is.null(csv_name)){
    write.table(df,csv_name,sep = ",",
                append = FALSE,quote = FALSE, col.names = NA, row.names = TRUE)
  }
  return(df)
}


#' Model diagnostics.
#'
#' \code{getModelDiag} returns model diagnostics,
#' and save the results to a .csv file when \code{csv_name} is specified
#'
#' @param devName string. name for dev in rstan codes
#' @inheritParams getParaSummary
#'
#' @return A dataframe indicating DIC, WAIC, model running time,
#' maximum rhat, minimum effective sample size,
#' number of divergent transitions and number of max treedepth hit.
#' @export

getModelDiag = function(stan_data_fit, csv_name = NULL, devName = 'dev_m'){
  stan_fit = stan_data_fit$stan_fit
  time = getTime(stan_fit)
  df = cbind.data.frame(getDIC(stan_fit,devName),
                        getWAIC(stan_fit),
                        elapsed_time_min = min(time),
                        elapsed_time_max = max(time),
                        rhat_max = max(getRhat(stan_fit),na.rm=T),
                        ess_min = min(getESS(stan_fit),na.rm=T),
                        divergent = rstan::get_num_divergent(stan_fit),
                        max_tree = rstan::get_num_max_treedepth(stan_fit))
  if (!is.null(csv_name)){
    write.table(df,csv_name,sep = ",",
                 append = FALSE,quote = FALSE, col.names = NA, row.names = TRUE)
  }
  return(df)
}


# Little function to calculate posterior variances from simulation
#' Title
#'
#' @param a
#'
#' @return
#' @export
#'
#' @examples
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}

# The calculation of Waic!  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
#' Title
#'
#' @param stan_fit
#' @param logName
#'
#' @return
#' @export
#'
#' @examples
getWAIC <- function (stan_fit, logName = "log_lik"){
  log_lik <- rstan::extract (stan_fit, logName)[[1]]
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}

#' Title
#'
#' @param stan_fit
#' @param devName
#'
#' @return
#' @export
#'
#' @examples
getDIC <- function (stan_fit, devName = "dev"){
  dev <- unlist(rstan::extract (stan_fit, devName))
  return (data.frame (dic_1 = mean(dev) + .5 * var(dev),
                      dev_mean = mean(dev),
                      pd_1 = .5 * var(dev),
                      dic_2 = mean(dev) + mean(dev) - min(dev),
                      min_dev = min(dev),
                      pd_2 = mean(dev) - min(dev)))
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
quantile1 = function(x){
  quantile1 = c(mean=mean(x), median=median(x), quantile(x, c(0.025,0.05,0.1,0.2,0.25,0.5,0.75,0.8,0.9,0.95,0.975)), sd=sqrt(var(x)));
  return(quantile1)
}
#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
quantile2 = function(x){
  quantile2 = c(quantile(x, c(0.1,0.3,0.5,0.7,0.9)));
  return(quantile2)
}


#' Title
#'
#' @param stan_fit
#' @param parameters
#'
#' @return
#' @export
#'
#' @examples
getRhat = function(stan_fit,parameters){
  return(rstan::summary(stan_fit)$summary[parameters,"Rhat"] )
}

#' Title
#'
#' @param stan_fit
#' @param parameters
#'
#' @return
#' @export
#'
#' @examples
getESS = function(stan_fit,parameters){
  return(rstan::summary(stan_fit)$summary[parameters,"n_eff"] )
}

#' Title
#'
#' @param stan_fit
#'
#' @return
#' @export
#'
#' @examples
getTime = function(stan_fit){
  return(rowSums(rstan::get_elapsed_time(stan_fit)))
}


#' Title
#'
#' @param v
#'
#' @return
#' @export
#'
#' @examples
getMode <- function(v) {
  uniqv <- unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

#' Title
#'
#' @param data
#' @param measurevar
#' @param groupvars
#' @param na.rm
#' @param conf.interval
#' @param .drop
#'
#' @return
#' @export
#'
#' @examples
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )

  # Rename the "mean" column
  datac <- plyr::rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}


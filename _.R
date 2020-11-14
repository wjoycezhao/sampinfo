# install.packages("rstantools")
# rstantools::rstan_create_package(path = '/Users/joyce/Dropbox/sampinfo',
# stan_files = c('mm10d01.stan','mm10d9.stan'))

# I finally found the solution. It is important to write the
# following lines of code in the NAMESPACE file, otherwise you will have the error mentioned above.
# import(Rcpp)
# import(methods)
# importFrom(rstan,sampling)
# useDynLib(rstanlm, .registration = TRUE)
#
import(Rcpp)
import(methods)
importFrom(rstan, sampling)
useDynLib(sampinfo)

require(devtools)

pkgbuild::compile_dll(force = TRUE)
roxygen2::roxygenise(clean = TRUE)

devtools::install(quick = FALSE) ##MUST BE FALSE THE FIRST TIME

# devtools::build_vignettes()

# load_all()
use_package('plyr')
use_package('dplyr')
use_package('tidyr')
use_package('rstan')
use_package('utils')
use_package('stats')
use_package('truncnorm')
use_package('bayesplot')


require(dplyr)

# md_data = (subset(read.csv('exp1_data_C3.csv'),qID==1))
# md_data_pri = (subset(read.csv('exp2_data_C3.csv'),qID==1))
# dist_data_pri = (subset(read.csv('exp2_dist_C3.csv'),qID==1))
# head(md_data)
# usethis::use_data( md_data_pri, dist_data_pri, overwrite = FALSE)

# raw_data = subset(read.csv('exp2_data_C3.csv'),sID<11&qID<6)[,-1]
# format_data = subset(getDataHier(3,'exp2')$format_data[,-1],sID<11&qID<6)
# usethis::use_data(raw_data, format_data, overwrite = FALSE)

# install.packages("../sampinfo", repos = NULL, type = "source")
# iter_num = 1000,chain_num = 2,warmup_num = 100, core_num=2

# usethis::use_vignette("tutorial")
# rmarkdown::render("vignettes/tutorial.Rmd", "all")

# feature_data_8 = cbind(read.csv('exp2_dist_C3.csv')[,2],
#                        do.call(rbind, replicate(8, cbind(diag(7), getSameA(3)),
#                                                 simplify=FALSE)),
#                        read.csv('exp2_dist_C3.csv')[,-c(1:2)])
# colnames(feature_data_8) = c('qID',unlist(lapply(c('sameC_','sameA_','dist_'), function(x) paste0(x,1:7))))
# head(feature_data_8)
### usethis::use_data(feature_data_8, overwrite = FALSE)
temp_cp1 = subset(getData(2,'exp1','')$format_data,tNo==1)%>%
  mutate(cID_s = sign(cID -3)) %>%
  group_by(qID,cID_s)%>%
  summarize(cp = mean(final_choice == 1),n=n())
temp_cp2 = subset(getData(2,'exp2','')$format_data,tNo==1)%>%
  mutate(cID_s = sign(cID -3)) %>%
  group_by(qID,condition)%>%
  summarize(cp = mean(final_choice == 1),n=n())


qq=5
temp = getData(4,'exp1','')
format_data = subset(temp$format_data, qID == qq)
format_data$qID = 1
X_sim = temp$X_sim[[qq]]
cluster_rating_m = subset(temp$cluster_rating_m, qID == qq)[,-c(1:2)]
stan_data = getStanData(
  beta = c('sameC','sameA', 'dist'),
  deltaM_value = 8,
  binary_value = 0,
  ar_value = 0,
  random_value = 0,
  deltaD_value = 8,
  option_num = 9,
  format_data = format_data,
  X_sim = X_sim,
  cluster_rating_m = cluster_rating_m,
  max_tNo_prd = 30
)
stan_fit = rstan::stan(
  # file = 'mmdm_threshold.stan',
  file = 'mmdm_timelimit.stan',
  data = stan_data,
  seed = 1,
  # sample_file = save_model_file,
  # init = init_values,
  iter = 20,
  chains = 2,
  warmup = 2,
  cores = 2
  # refresh = 100,
  # save_warmup = F,
  # init_r = 1
)
stan_data_fit = list(
  stan_data = stan_data,
  stan_fit = stan_fit
)
rbind(getModelDiag(stan_data_fit,dev_name = 'dev_d',logName = 'log_lik_d'),
      getModelDiag(stan_data_fit,dev_name = 'dev_m',logName = 'log_lik_m'))

getParaSummary(stan_data_fit = stan_data_fit)
as.data.frame(
  rstan::extract(stan_data_fit$stan_fit, pars = c('terminate_prd','cID_prd'))
)

temp = getPara(para_name = 'theta', stan_data_fit)
temp[1:10,1:100]
getData(3,'exp1','')

qq=1
temp = getData(4,'exp1','')
format_data = subset(temp$format_data, qID == qq)
stan_data_fit = getStanFit(beta = c('sameA', 'dist'), deltaM_value = 8,
                           option_num = 9, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 20,chain_num = 1,warmup_num = 5, core_num=8,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE,
                           hier_value = 0)
getParaSummary(stan_data_fit = stan_data_fit)
stan_fit = stan_data_fit$stan_fit
r = getRhat(stan_fit)
head(r)
stan_data_fit = getStanFit(binary_value = 0,
                           ar_value = 0,
                           random_value = 0,
                           deltaD_value = 1,
                           option_num = 7, format_data = md_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 500,chain_num = 1,warmup_num = 200, core_num=1,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE,
                           hier_value = 0)
getParaSummary(stan_data_fit = stan_data_fit)
getModelDiag(stan_data_fit,dev_name = 'dev_d')


parameters1 = c('decision_0','decision_no','decision_yes')


postPred = cbind.data.frame(#model_name = model_name,
                            # binary_value=binary_value[i],ar_value=ar_value[i],random_value=random_value[i],deltaD_value=deltaD_value[i],
                            pp_all1)

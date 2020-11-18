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


qq=3
temp = getData(3,'exp1','')
format_data = subset(temp$format_data, qID == qq)
format_data$qID = 1
X_sim = format_data$X_sim[[qq]]
cluster_rating_m = subset(temp$cluster_rating_m, qID == qq)[,-c(1:2)]
stan_data = getStanData(
  # beta = character(0),
  beta = c('sameA','dist'),
  deltaM_value = 8,
  binary_value = 0,
  ar_value = 0,
  random_value = 0,
  deltaD_value = 8,
  option_num = 7,
  format_data = format_data
  X_sim = X_sim,
  cluster_rating_m = cluster_rating_m,
  max_tNo_prd = 30
)
stan_fit = rstan::stan(
  # file = 'mmdm_threshold.stan',
  file = 'mmdm_threshold.stan',
  data = stan_data,
  seed = 1,
  # sample_file = save_model_file,
  # init = init_values,
  iter = 3000,
  chains = 4,
  warmup = 500,
  cores = 4
  # refresh = 100,
  # save_warmup = F,
  # init_r = 1
)
stan_data_fit = list(
  stan_data = stan_data,
  stan_fit = stan_fit
)
temp = as.data.frame(rstan::extract(stan_data_fit$stan_fit,
                                    pars = c('terminate_prd','cID_prd')))
sel = (c(1:nrow(temp))%%1==0)
temp = temp[sel, ]
write.csv(temp,  paste0(save_folder, 'exp', exp_no, '_', model_name, '_simulate.csv'))

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
                           hier_value = 1)
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


require(sampinfo)
require(readr)
require(plyr)
dist_all = (subset(read.csv('exp1_dist_C3.csv'),qID==3))
data_all = read.csv('exp3_m3_0118_d_1008_simulate.csv')[1:500,] %>%
  gather(tNo, cID,cID_prd.1:cID_prd.30) %>%
  mutate(tNo = parse_number(tNo, locale = locale(
    grouping_mark = ". ", decimal_mark = ","))) %>%
  filter(cID != 99) %>%
  mutate(rating = sign(cID - 4),sID = X, qID = 3) %>%
  mutate(sID = mapvalues(sID, from = unique(sID), to = 1:length(unique(sID)))) %>%
  arrange(X)
option_num = 2 * cluster_num + 1
## add cluster resampling feature for predicting the next trial
data_all = addFeatureMatrix(diag(option_num),
                            data_all, option_num, 'sameC', 'cID')
data_all = addFeatureMatrix(getSameA(cluster_num),
                            data_all, option_num, 'sameA', 'cID')
## add semantic distance for predicting the next trial
# dist_all = read.csv(paste0(data_folder, exp_name, '_dist_c', cluster_num, '.csv'),
#                     row.names = 1)
qq = 3
data_all = addFeatureMatrix(dist_all[dist_all$qID == qq, -c(1:2)],
                   data_all[data_all$qID == qq, ],
                   option_num, 'dist', 'cID')
data_all$qID = 1
stan_data = getStanData(
  # beta = character(0),
  beta = c('sameA','dist'),
  deltaM_value = 8,
  binary_value = 1,
  ar_value = 0,
  random_value = 0,
  deltaD_value = 8,
  option_num = 7,
  format_data = data_all
  # X_sim = X_sim,
  # cluster_rating_m = cluster_rating_m,
  # max_tNo_prd = 30
)
stan_fit = rstan::stan(
  # file = 'mmdm_threshold.stan',
  # file = 'dm_threshold.stan',
  file = 'mm.stan',
  data = stan_data,
  seed = 1,
  # sample_file = save_model_file,
  # init = init_values,
  iter = 2000,
  chains = 4,
  warmup = 500,
  cores = 4
  # refresh = 100,
  # save_warmup = F,
  # init_r = 1
)
stan_data_fit = list(
  stan_data = stan_data,
  stan_fit = stan_fit
)
getParaSummary(stan_data_fit = stan_data_fit)



dist_all = (subset(read.csv('exp1_dist_C3.csv'),qID==3))
temp = read.csv('happy_3_pp.csv')[1:500,]
data_all =  temp%>%
  subset(b_sameC_value = 0, b_sameA_value = 1,
         b_dist_value = 1, ar_value = 0, random_value = 0, deltaD_value = 9)
data_all = data_all[,-c(1:12)] %>%
  mutate(sID = 1:500) %>%
  gather(tNo, cID, cID_prd.1:cID_prd.100) %>%
  mutate(terminate = decision_prd, tNo = parse_number(tNo, locale = locale(
    grouping_mark = ". ", decimal_mark = ","))) %>%
  filter(cID != 100) %>%
  mutate(rating = sign(cID - 4),qID = 3) %>%
  mutate(sID = mapvalues(sID, from = unique(sID), to = 1:length(unique(sID)))) %>%
  arrange(sID)
option_num = 2 * cluster_num + 1
## add cluster resampling feature for predicting the next trial
data_all = addFeatureMatrix(diag(option_num),
                            data_all, option_num, 'sameC', 'cID')
data_all = addFeatureMatrix(getSameA(cluster_num),
                            data_all, option_num, 'sameA', 'cID')
## add semantic distance for predicting the next trial
# dist_all = read.csv(paste0(data_folder, exp_name, '_dist_c', cluster_num, '.csv'),
#                     row.names = 1)
qq = 3
data_all = addFeatureMatrix(dist_all[dist_all$qID == qq, -c(1:2)],
                            data_all[data_all$qID == qq, ],
                            option_num, 'dist', 'cID')
data_all$qID = 1
stan_data = getStanData(
  # beta = character(0),
  beta = c('sameA','dist'),
  deltaM_value = 8,
  binary_value = 1,
  ar_value = 0,
  random_value = 0,
  deltaD_value = 8,
  option_num = 7,
  format_data = data_all
  # X_sim = X_sim,
  # cluster_rating_m = cluster_rating_m,
  # max_tNo_prd = 30
)
stan_fit = rstan::stan(
  # file = 'mmdm_threshold.stan',
  file = 'dm_threshold.stan',
  # file = 'mm.stan',
  data = stan_data,
  seed = 1,
  # sample_file = save_model_file,
  # init = init_values,
  iter = 2000,
  chains = 4,
  warmup = 500,
  cores = 4
  # refresh = 100,
  # save_warmup = F,
  # init_r = 1
)
stan_data_fit = list(
  stan_data = stan_data,
  stan_fit = stan_fit
)
getParaSummary(stan_data_fit = stan_data_fit)

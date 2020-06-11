# install.packages("rstantools")
# rstantools::rstan_create_package(path = '/Users/joyce/Dropbox/sampinfo',
# stan_files = c('mm10d01.stan','mm10d9.stan'))

# I finally found the solution. It is important to write the
# following lines of code in the NAMESPACE file, otherwise you will have the error mentioned above.
#
import(Rcpp)
import(methods)
importFrom(rstan, sampling)
useDynLib(sampinfo)

require(devtools)

pkgbuild::compile_dll(force = TRUE)
roxygen2::roxygenise(clean = TRUE)
devtools::install(quick = F)##MUST BE FALSE THE FIRST TIME
devtools::build_vignettes()

load_all()
use_package('plyr')
use_package('dplyr')
use_package('tidyr')
use_package('rstan')
use_package('utils')
use_package('stats')
use_package('truncnorm')
use_package('bayesplot')

require(dplyr)


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
stan_data_fit = getStanFit(beta = c('sameA', 'dist'), deltaM_value = 9,
                           option_num = 7, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 200,chain_num = 1,warmup_num = 100, core_num=1,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE)
stan_data_fit = getStanFit(beta = c('sameA', 'dist'), deltaM_value = 8,
                           option_num = 7, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 200,chain_num = 1,warmup_num = 100, core_num=1,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE)

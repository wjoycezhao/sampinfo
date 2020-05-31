# install.packages("rstantools")
# rstantools::rstan_create_package(path = '/Users/joyce/Dropbox/sampinfo',
                                 # stan_files = c('mm10d01.stan','mm10d9.stan'))

# I finally found the solution. It is important to write the
# following lines of code in the NAMESPACE file, otherwise you will have the error mentioned above.
#
import(Rcpp)
import(methods)
importFrom(rstan,sampling)
useDynLib(sampinfo)

pkgbuild::compile_dll(force=TRUE)
roxygen2::roxygenise(clean=TRUE)


require(devtools)
load_all()
use_package('plyr')
use_package('dplyr')
use_package('tidyr')
use_package('rstan')
use_package('utils')
use_package('stats')
use_package('bayesplot')

require(dplyr)

# raw_data = subset(read.csv('exp2_data_C3.csv'),sID<11&qID<6)[,-1]
# format_data = subset(getDataHier(3,'exp2')$format_data[,-1],sID<11&qID<6)
# usethis::use_data(raw_data, format_data, overwrite = FALSE)
devtools::install(quick=F)##MUST BE FALSE THE FIRST TIME
# install.packages("../sampinfo", repos = NULL, type = "source")
# iter_num = 1000,chain_num = 2,warmup_num = 100, core_num=2

# usethis::use_vignette("tutorial")
#。。

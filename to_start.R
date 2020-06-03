## Install devtools, so that we can run devtools::install_github

if(require("devtools")){
} else {
    install.packages("devtools")
    if(require("devtools")){
        print("devtools installed and loaded")
    } else {
        stop("devtools not installed")
    }
}

## Install codes for the sampling model
install_github('wjoycezhao/sampinfo')
requre(sampinfo)

## If the code is installed successfully, you should have a very small toy dataset now.
## See the first 6 rows of it
head(format_data,6)

## Now we can try fitting a model (takes about 5s)
## if you can see the data but cannot fit the model
## it is likely that RStan did not install successfully..
stan_data_fit = getStanFit(beta = c('sameA', 'dist'), deltaM_value = 9, 
                           option_num = 7, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 200,chain_num = 1,warmup_num = 100, core_num=1,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE)


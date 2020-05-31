## install devtools

if(require("devtools")){
} else {
    install.packages("devtools")
    if(require("devtools")){
        print("devtools installed and loaded")
    } else {
        stop("devtools not install lme4")
    }
}

## install codes for the sampling model
install_github('wjoycezhao/sampinfo')
requre(sampinfo)

## We have included a very small toy dataset. See the first 6 rows of it
head(format_data,6)

## fit the model
stan_data_fit = getStanFit(beta = c('sameA', 'dist'), deltaM_value = 9, 
                           option_num = 7, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 200,chain_num = 4,warmup_num = 100, core_num=4,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE)

## see tutorial.html for other functions
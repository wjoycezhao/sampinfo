``` r
library(sampinfo)
library(bayesplot)
#> This is bayesplot version 1.7.2
#> - Online documentation and vignettes at mc-stan.org/bayesplot
#> - bayesplot theme set to bayesplot::theme_default()
#>    * Does _not_ affect other ggplot2 plots
#>    * See ?bayesplot_theme_set for details on theme setting
library(rstan)
#> Loading required package: StanHeaders
#> Loading required package: ggplot2
#> rstan (Version 2.19.3, GitRev: 2e1f913d3ca3)
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores()).
#> To avoid recompilation of unchanged Stan programs, we recommend calling
#> rstan_options(auto_write = TRUE)
library(knitr)
```

Prepare your data:
------------------

To begin with, you have a dataset with many participants, each finishing a same set of questions. Within each questions, there are several time points, where they need to sample from a set of C options.
In the toy data, total question number *Q* = 5. The question ID for each sampling time point is stored in qID.

``` r
subset(sampinfo::raw_data,sID==1)
#>    qID sID final_choice tNo cID rating terminate
#> 1    1   1            1   1   5      3         1
#> 2    2   1            1   1   4      0         0
#> 3    2   1            1   2   3     -1         0
#> 4    2   1            1   3   6      3         1
#> 5    3   1            1   1   6      3         0
#> 6    3   1            1   2   2     -3         0
#> 7    3   1            1   3   6      2         1
#> 8    4   1            1   1   5      3         0
#> 9    4   1            1   2   7      3         0
#> 10   4   1            1   3   5      3         0
#> 11   4   1            1   4   5      2         0
#> 12   4   1            1   5   6      3         1
#> 13   5   1            0   1   2     -3         0
#> 14   5   1            0   2   6      3         0
#> 15   5   1            0   3   4      0         0
#> 16   5   1            0   4   1     -3        -1
```

And participant number *S* = 10. The participant ID is listed in sID.

``` r
subset(sampinfo::raw_data,qID==1)
#>     qID sID final_choice tNo cID rating terminate
#> 1     1   1            1   1   5      3         1
#> 35    1   2            1   1   5      3         0
#> 36    1   2            1   2   5      2         0
#> 37    1   2            1   3   4      0         0
#> 38    1   2            1   4   7      3         0
#> 39    1   2            1   5   5      3         1
#> 85    1   3            1   1   5      3         0
#> 86    1   3            1   2   7      3         0
#> 87    1   3            1   3   5      3         0
#> 88    1   3            1   4   7      3         1
#> 115   1   4            1   1   6      3         0
#> 116   1   4            1   2   5      3         0
#> 117   1   4            1   3   5      3         1
#> 152   1   5            1   1   5      3         0
#> 153   1   5            1   2   3     -1         1
#> 200   1   6            1   1   5      3         0
#> 201   1   6            1   2   1     -3         0
#> 202   1   6            1   3   5      3         0
#> 203   1   6            1   4   7      2         0
#> 204   1   6            1   5   7      3         1
#> 242   1   7            1   1   5      2         0
#> 243   1   7            1   2   7      3         0
#> 244   1   7            1   3   7      3         0
#> 245   1   7            1   4   3     -3         1
#> 270   1   8            0   1   3     -3         0
#> 271   1   8            0   2   3     -2         0
#> 272   1   8            0   3   7      1        -1
#> 293   1   9            1   1   7      2         0
#> 294   1   9            1   2   7      2         0
#> 295   1   9            1   3   5      2         1
#> 317   1  10            1   1   3     -2         0
#> 318   1  10            1   2   5      1         0
#> 319   1  10            1   3   5      3         1
```

For each question, there are many time points, where participant sample from different options. Here time points are saved in tNo, and option sampled in cID. We have in total 7 different options (cID from 1 to 7).

Understand the model:
---------------------

Our goal is to predict the probabilities of sampling different features at each time point, using several known features of the options.

In this example, we have three features: sameC, sameA and dist.

-   sameC: whether the option is in the same as the one sampled in the previous time point
-   sameA: whether they are decision congruent
-   dist: the semantic distance between the option and the previously sampled one

Suppose that we use a single feature, and assume the effect can be captured by a linear multiplier *β*. Then for a trial with C options, we write the feature values into a vector
**X**<sub>**C** **×** **1**</sub> = \[*X*<sub>1</sub>, ..., *X*<sub>*C*</sub>\]<sup>*T*</sup>
 The predicted choice *Y* thus follows a categorical distribution.
*Y* ∼ *C**a**t*(*s**o**f**t**m**a**x*(**X**<sub>*C* × 1</sub>*β*))
.

In a case with *K* features, the effects can be captured by a vector **β**<sub>**K** **×** **1**</sub>, and the feature matrix is **X**<sub>**C** **×** **K**</sub>. 

The predicted smapling probability
*Y* ∼ *C**a**t*(*s**o**f**t**m**a**x*(**X**<sub>**C** **×** **K**</sub>**β**<sub>**K** **×** **1**</sub>)

This is the model for predicting the sampling probability in one time point, of one question, completed by one participant.

Now how about heterogeneity across paricipants and questions? We need a *hiearchical model* to address that.

We first assume there is a grand mean for the feature effects, and denote it as **μ**<sub>**K** **×** **1**</sub><sup>**β**</sup>.

The *S* participants deviate from the grand mean on each of the *K* features. Here the deviation can be written as **δ**<sub>**K** **×** **S**</sub><sup>**β****,** **S**</sup>.

Similarly, the deviations of the *Q* questions on the features can be written as **δ**<sub>**K** **×** **Q**</sub><sup>**β****,** **Q**</sup>.

Each row of these deviation matrices forllows a centered normal distribution. Therefore, we need scale parameters to quantify these deviations. For participant and question level deviations, the scale parameters are **σ**<sub>**K** **×** **1**</sub><sup>**β****,** **S**</sup> and **σ**<sub>**K** **×** **1**</sub><sup>**β****,** **Q**</sup> respectively.

Therefore, the parameters we need to fit include: **μ**<sub>**K** **×** **1**</sub><sup>**β**</sup>, **σ**<sub>**K** **×** **1**</sub><sup>**β****,** **S**</sup>, **σ**<sub>**K** **×** **1**</sub><sup>**β****,** **Q**</sup>, **δ**<sub>**K** **×** **S**</sub><sup>**β****,** **S**</sup>, and **δ**<sub>**K** **×** **Q**</sub><sup>**β****,** **Q**</sup>. For participant *s* and question *q*, the predicted samplng choice
*Y* ∼ *C**a**t*(*s**o**f**t**m**a**x*(**X**<sub>**C** **×** **K**</sub>(**μ**<sub>**K** **×** **1**</sub><sup>**β**</sup> + **δ**<sub>**\*****,** **s**</sub><sup>**β****,** **S**</sup> + **δ**<sub>**\*****,** **q**</sub><sup>**β****,** **Q**</sup>))

Beyond these we can also have decay, and different base rates for different questions.

Now let's take a break from maths...

And take a look at the code.

``` r
sampinfo::getStanCode(9)
#> S4 class stanmodel 'mm10d9' coded as follows:
#> 
#> data{
#>   // memory
#>   int<lower=2> C; // number of option categories
#>   int<lower=0> K; // number of features; if 0 then no betas
#>   int<lower=1> N; // number of time points
#>   int<lower=1> S; // number of participants
#>   int<lower=1> Q; // number of questions
#>   int<lower=0> sID[N]; // participant ID variable
#>   int<lower=0> qID[N]; // question ID variable
#>   int<lower=1> tNo[N]; // thought No.x; used to reset decay
#>   matrix[C,K] X[N]; // feature matrix, for different time points
#>   int<lower=1, upper=C> cID[N]; // response cluster ID, 1-7 for 3-cluster solution
#> }
#> 
#> parameters{
#>   real<lower=0> deltaM; // decay
#>   matrix[C-1,Q] alpha_raw; // base rate
#>   vector[K] beta_mu; // group-level mean for beta
#>   vector<lower=0>[K] beta_s_sd; // sd for beta across individuals
#>   vector<lower=0>[K] beta_q_sd; // sd for beta across questions
#>   matrix[K,Q] beta_q_raw; // question-level beta
#>   matrix[K,S] beta_s_raw; // individual-level beta
#> }
#> 
#> transformed parameters {
#>   matrix[C,Q] alpha;
#>   matrix[C,N] theta;
#>   matrix[K,S] beta[Q];
#>   if (K > 0){
#>     for (q in 1:Q){
#>     beta[q] = rep_matrix(beta_mu,S) +
#>               diag_pre_multiply(beta_s_sd,beta_s_raw) +
#>               rep_matrix(beta_q_sd .* beta_q_raw[,q],S);
#>     }
#>   }
#> 
#>   alpha[1:(C-1),] = alpha_raw/2;
#>   alpha[C,] = to_row_vector(rep_array(0,Q));
#> 
#>   {
#>     if (K > 0){
#>       matrix[C,K] X_acc;
#>       for (n in 1:N){
#>         if (tNo[n] == 1){
#>           X_acc = to_matrix(rep_array(0,C,K));
#>         }else{
#>           X_acc *= deltaM;
#>         }
#>         X_acc += X[n];
#>         theta[,n] = X_acc * beta[qID[n],,sID[n]];
#>       }
#>       theta += alpha[,qID];
#>     }else{
#>       theta = alpha[,qID];
#>     }
#>   }
#> }
#> 
#> model{
#>   deltaM ~ normal(0.8,0.8)T[0,];
#>   beta_mu ~ std_normal();
#>   beta_q_sd ~ std_normal();
#>   beta_s_sd ~ std_normal();
#>   to_vector(beta_q_raw) ~ std_normal();
#>   to_vector(beta_s_raw) ~ std_normal();
#>   to_vector(alpha_raw) ~ std_normal();
#> 
#>   for (n in 1:N){
#>     cID[n] ~ categorical_logit(theta[,n]);
#>   }
#> }
#> 
#> generated quantities {
#>   vector[N] log_lik;
#>   real dev_m;
#>   matrix[K,Q] beta_qdev;
#>   matrix[K,S] beta_sdev;
#>   vector[K] beta_q_temp;
#>   vector[K] beta_s_temp;
#>   matrix[K,Q] beta_qmean;
#>   matrix[K,S] beta_smean;
#>   dev_m = 0;
#>   for (n in 1:N){
#>     log_lik[n] = categorical_logit_lpmf(cID[n]|theta[,n]);
#>     dev_m += - 2*log_lik[n];
#>   }
#>   beta_qdev = diag_pre_multiply(beta_q_sd, beta_q_raw);
#>   beta_sdev = diag_pre_multiply(beta_s_sd,beta_s_raw);
#>   if(K>0){
#>    for (k in 1:K) {
#>      beta_q_temp[k] = mean(to_array_1d(beta_qdev[k,]));
#>      beta_s_temp[k] = mean(to_array_1d(beta_sdev[k,]));
#>      }
#>    beta_qmean = rep_matrix(beta_mu+beta_s_temp,Q) + beta_qdev;
#>    beta_smean = rep_matrix(beta_mu+beta_q_temp,S) + beta_sdev;
#>   }
#> }
```

Try fitting your model to the data:
-----------------------------------

We will use a vector of feature names to specify the which features to be included in the model. For example:

``` r
beta = c('sameC', 'sameA', 'dist')
```

If you need a baseline model without any features, set

``` r
beta = c()
```

Next we need to decide whether we want decay in our model.
If we assume there is no memory at all (i.e., Markov property), then we set

``` r
deltaM_value = 0
```

On the other extreme, if no information is ever forgotten within a question (i.e., no decay). We set

``` r
deltaM_value = 1
```

Finally we can estimate delta as a free parameter. Suppose we need a free parameter for delta, then we can specify

``` r
deltaM_value = 9
```

With the data, and the model specifications, we can now run the following to fit the model.
To make sure that the samples from the posterior distributions are not biased, we need to make sure the chains converge. Only then are we sampling the stationary distribution of the Marco chain. Because many models take minutes/hours/days/... weeks to run, we need to choose a warm-up number carefully. Try some warmup number, test convergence, and repeat until things look good.

``` r
beta = c('sameA', 'dist')
stan_data_fit = sampinfo::getStanFit(beta = beta, deltaM_value = 9, 
                           option_num = 7, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 200,chain_num = 4,warmup_num = 100, core_num=4,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE)
```

Check model convergence using traceplots.

``` r
posterior_w = rstan::extract(stan_data_fit$stan_fit, inc_warmup = TRUE, permuted = FALSE)
dim(posterior_w)
#> [1]  200    4 1940
lapply(dimnames(posterior_w), head)
#> $iterations
#> NULL
#> 
#> $chains
#> [1] "chain:1" "chain:2" "chain:3" "chain:4"
#> 
#> $parameters
#> [1] "deltaM"         "alpha_raw[1,1]" "alpha_raw[2,1]" "alpha_raw[3,1]"
#> [5] "alpha_raw[4,1]" "alpha_raw[5,1]"
bayesplot::color_scheme_set("viridis")
pars0 = c("beta_mu[1]","beta_s_sd[1]","beta_q_sd[1]","beta_s_raw[1,10]","deltaM","alpha[3,1]")
bayesplot::mcmc_trace(posterior_w, pars = pars0, 
                      n_warmup=100, size = 0.5,  
                      facet_args = list(nrow = 3))
```

<img src="/private/var/folders/fr/yxz4164j33v3xmtf243926gr0000gn/T/RtmpgyVIoW/preview-11122ac34df1.dir/my-vignette_files/figure-markdown_github/unnamed-chunk-11-1.png" width="85%" style="display: block; margin: auto;" />

Another useful plot, which shows you if the ranks of the samples are mixed well among chains.

``` r
bayesplot::mcmc_rank_overlay(posterior_w, pars = pars0)
```

<img src="/private/var/folders/fr/yxz4164j33v3xmtf243926gr0000gn/T/RtmpgyVIoW/preview-11122ac34df1.dir/my-vignette_files/figure-markdown_github/unnamed-chunk-12-1.png" width="85%" style="display: block; margin: auto;" />

Also check gelman-rubin rhat.

``` r
rhat = sampinfo::getRhat(stan_data_fit$stan_fit)
print(rhat[rhat>1.05&!is.na(rhat)])
#>           deltaM     beta_s_sd[1] beta_s_raw[1,10]      theta[1,60] 
#>         1.055197         1.087298         1.053247         1.051933 
#>      theta[1,61]      theta[2,97]      theta[2,98]      theta[5,97] 
#>         1.052093         1.053726         1.052967         1.053479 
#>      theta[5,98]      theta[6,97] 
#>         1.063804         1.052355
```

Run the model:
--------------

Now take the (conservative) warm-up number, and choose the total iteration number so that the sample size will be enough for your statistics of interests.

``` r
stan_data_fit = sampinfo::getStanFit(beta= c('sameA', 'dist'), deltaM_value = 9, 
                           option_num = 7, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 1200,chain_num = 4,warmup_num = 300, core_num=4,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE)
```

Check the model diagnostics and goodness of fit:
------------------------------------------------

ESS is effective sample size, 1000 looks good in my opinion.

``` r
model_diag = sampinfo::getModelDiag(stan_data_fit = stan_data_fit)
print(model_diag)
#>      dic_1 dev_mean     pd_1    dic_2  min_dev     pd_2     waic   p_waic
#> 1 749.5508 696.5161 53.03465 729.3835 663.6487 32.86743 726.0886 27.93877
#>        lppd p_waic_1 elapsed_time_min elapsed_time_max rhat_max  ess_min
#> 1 -335.1055 26.30509         21.79493          37.8335 1.006804 1046.775
#>   divergent max_tree
#> 1         0        0
```

DIC and WAIC can be used for model comparisons. Lisheng also mentioned bayes factors, LOO-CV.

Good reference:
Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. Statistics and computing, 24(6), 997-1016. [Link](https://arxiv.org/abs/1307.5928)

Summarize and plot parameters
-----------------------------

``` r
para = sampinfo::getParaSummary(stan_data_fit = stan_data_fit)
print(round(para[1:3,],1))
#>               mean median 2.5%   5%  10%  20%  25%  50%  75%  80% 90% 95% 97.5%
#> deltaM         0.5    0.4  0.0  0.1  0.1  0.2  0.2  0.4  0.7  0.7 0.9 1.1   1.2
#> beta_sameA_mu  0.4    0.4 -0.5 -0.3 -0.1  0.1  0.1  0.4  0.6  0.7 0.8 1.0   1.1
#> beta_dist_mu  -0.5   -0.5 -1.4 -1.2 -1.1 -0.9 -0.8 -0.5 -0.3 -0.2 0.0 0.2   0.3
#>                sd
#> deltaM        0.3
#> beta_sameA_mu 0.4
#> beta_dist_mu  0.4
```

Here are the question level betas (grand mean plus question-level deviance)

``` r
print(round(para[stringr::str_detect(rownames(para),'qmean'),c(1:3,8,13:14)],1))
#>                    mean median 2.5%  50% 97.5%  sd
#> beta_sameA_qmean_1  0.5    0.5 -0.2  0.5   1.3 0.4
#> beta_dist_qmean_1  -0.5   -0.5 -1.5 -0.5   0.6 0.5
#> beta_sameA_qmean_2  0.2    0.2 -0.4  0.2   0.8 0.3
#> beta_dist_qmean_2  -0.7   -0.6 -1.6 -0.6   0.3 0.5
#> beta_sameA_qmean_3 -0.4   -0.4 -1.2 -0.4   0.4 0.4
#> beta_dist_qmean_3  -0.9   -0.9 -2.0 -0.9  -0.1 0.5
#> beta_sameA_qmean_4  1.0    1.0  0.3  1.0   1.8 0.4
#> beta_dist_qmean_4  -0.5   -0.5 -1.4 -0.5   0.5 0.5
#> beta_sameA_qmean_5  0.7    0.7  0.1  0.7   1.4 0.3
#> beta_dist_qmean_5  -0.4   -0.4 -1.4 -0.4   0.9 0.6
```

and the participant level betas (grand mean plus participant-level deviance)

``` r
print(round(para[stringr::str_detect(rownames(para),'smean'),c(1:3,8,13:14)],1))
#>                     mean median 2.5%  50% 97.5%  sd
#> beta_sameA_smean_1   0.2    0.2 -0.8  0.2   1.0 0.4
#> beta_dist_smean_1   -0.4   -0.4 -1.4 -0.4   1.1 0.6
#> beta_sameA_smean_2   0.3    0.3 -0.3  0.3   1.0 0.3
#> beta_dist_smean_2   -0.7   -0.7 -1.9 -0.7   0.1 0.5
#> beta_sameA_smean_3   0.3    0.3 -0.6  0.3   1.1 0.4
#> beta_dist_smean_3   -0.5   -0.5 -1.5 -0.5   0.8 0.5
#> beta_sameA_smean_4   0.9    0.9  0.2  0.9   2.1 0.5
#> beta_dist_smean_4   -0.5   -0.5 -1.4 -0.5   0.6 0.5
#> beta_sameA_smean_5   0.3    0.3 -0.3  0.3   1.0 0.3
#> beta_dist_smean_5   -0.7   -0.7 -1.8 -0.7   0.1 0.5
#> beta_sameA_smean_6   0.4    0.4 -0.3  0.4   1.1 0.3
#> beta_dist_smean_6   -0.5   -0.5 -1.5 -0.5   0.5 0.5
#> beta_sameA_smean_7   0.4    0.4 -0.4  0.4   1.4 0.4
#> beta_dist_smean_7   -0.5   -0.5 -1.6 -0.5   0.7 0.5
#> beta_sameA_smean_8   0.2    0.2 -0.7  0.2   1.1 0.4
#> beta_dist_smean_8   -0.6   -0.6 -1.7 -0.6   0.6 0.6
#> beta_sameA_smean_9   1.1    1.0  0.1  1.0   2.8 0.7
#> beta_dist_smean_9   -0.8   -0.7 -2.2 -0.7   0.2 0.6
#> beta_sameA_smean_10 -0.1    0.0 -1.2  0.0   0.7 0.5
#> beta_dist_smean_10  -0.8   -0.7 -2.2 -0.7   0.1 0.6
```

Also option base rates in differnet questions.

``` r
print(round(para[stringr::str_detect(rownames(para),'alpha'),c(1:3,8,13:14)],1)[1:20,])
#>           mean median 2.5%  50% 97.5%  sd
#> alpha.1.1 -0.4   -0.4 -1.2 -0.4   0.4 0.4
#> alpha.2.1 -0.6   -0.6 -1.4 -0.6   0.2 0.4
#> alpha.3.1  0.1    0.1 -0.7  0.1   0.8 0.4
#> alpha.4.1 -0.4   -0.4 -1.2 -0.4   0.3 0.4
#> alpha.5.1  0.7    0.7  0.1  0.7   1.4 0.3
#> alpha.6.1 -0.5   -0.5 -1.3 -0.5   0.2 0.4
#> alpha.7.1  0.0    0.0  0.0  0.0   0.0 0.0
#> alpha.1.2 -0.6   -0.6 -1.4 -0.6   0.1 0.4
#> alpha.2.2 -0.7   -0.7 -1.5 -0.7   0.0 0.4
#> alpha.3.2  0.5    0.5 -0.1  0.5   1.1 0.3
#> alpha.4.2 -0.3   -0.3 -1.0 -0.3   0.4 0.4
#> alpha.5.2 -0.2   -0.2 -0.9 -0.2   0.5 0.4
#> alpha.6.2  0.3    0.3 -0.3  0.3   0.9 0.3
#> alpha.7.2  0.0    0.0  0.0  0.0   0.0 0.0
#> alpha.1.3 -0.3   -0.3 -1.1 -0.3   0.4 0.4
#> alpha.2.3  0.4    0.4 -0.3  0.4   1.0 0.3
#> alpha.3.3  0.1    0.1 -0.6  0.1   0.8 0.3
#> alpha.4.3 -0.7   -0.6 -1.4 -0.6   0.1 0.4
#> alpha.5.3 -0.4   -0.3 -1.1 -0.3   0.4 0.4
#> alpha.6.3  1.2    1.2  0.6  1.2   1.8 0.3
```

There are many ways to visualize the samples. For example, check the bayesplot package [Link](https://mc-stan.org/bayesplot/articles/plotting-mcmc-draws.html)

``` r
posterior = rstan::extract(stan_data_fit$stan_fit, inc_warmup = TRUE, permuted = FALSE)
dim(posterior)
#> [1] 1200    4 1940
head(dimnames(posterior)[[3]])
#> [1] "deltaM"         "alpha_raw[1,1]" "alpha_raw[2,1]" "alpha_raw[3,1]"
#> [5] "alpha_raw[4,1]" "alpha_raw[5,1]"
## run changeBetaNames to make the parameter names understandable (so that we dont need to deal with matrix index)
dimnames(posterior)[[3]] = sampinfo::changeBetaNames(dimnames(posterior)[[3]],
                            beta = beta,
                            Q = max(format_data$qID), 
                            S = max(format_data$sID))

bayesplot::mcmc_intervals(posterior, regex_pars = "^beta_dist_qmean")
```

<img src="/private/var/folders/fr/yxz4164j33v3xmtf243926gr0000gn/T/RtmpgyVIoW/preview-11122ac34df1.dir/my-vignette_files/figure-markdown_github/unnamed-chunk-20-1.png" width="50%" style="display: block; margin: auto;" />

Additional diagnostics in Rstan
-------------------------------

There are other diagnostics to pay attention to if you run your model with Stan (HMC and NUTS)

``` r
model_diag = getModelDiag(stan_data_fit = stan_data_fit)
print(model_diag[,c('divergent', 'max_tree')])
#>   divergent max_tree
#> 1         0        0
```

For example you need to make sure there is no divergent transitions, which is not uncommon for hierarchical models. When it happens you often need to: 1) debug your rstan code. 2) try recovering your parameters 3) reparameterize your model or 4) increase alpha\_delta, increase tree\_maxdepth and decrease stepsize. I'll explain this later.
If the maximum tree depth is hit in NUTS, we need to increase max\_treedepth.
RStan will send warnings regarding these issues (remember to read it if the models are fit on GPC).

Other topics:
-------------

### Neal’s Funnel

Note that there are two ways to parameterize normal distributions (in hierarchical models):

-   Centered parameterization,
    *γ* ∼ *n**o**r**m**l*(*γ*<sub>*m**u*</sub>, *γ*<sub>*s**d*</sub>)
     Non-centered parameterization,
    *γ*<sub>*r**a**w*</sub> ∼ *n**o**r**m**l*(0, 1)
    *γ* = *γ*<sub>*m**u*</sub> + *γ*<sub>*s**d*</sub>*γ*<sub>*r**a**w*</sub>

If your individual level data are not very informative (and thus the individual parameters, *γ*, are influenced by the group level priors *γ*<sub>*m**u*</sub> and *γ*<sub>*s**d*</sub> a lot), it is much better to use the non-centered approach. Detailed explanations can be found here.

-   RStan manual 22.7 [Link](https://mc-stan.org/docs/2_23/stan-users-guide/reparameterization-section.html)
-   Betancourt, M., & Girolami, M. (2015). Hamiltonian Monte Carlo for hierarchical models. Current trends in Bayesian methodology with applications, 79(30), 2-4. [Link](https://arxiv.org/abs/1312.0906)

Intuitively, when *γ*<sub>*s**d*</sub> is smalll, *γ* tends to be very close to gamma\_mu. When *γ*<sub>*s**d*</sub> is large, *γ* varies a lot. That forms a tunnel shape. The stepsize needed for exploring the 'neck' is much much smaller than that required for the 'body'. Hence there's no single stepsize that stan can use to sample efficiently.
If we use the non-centered parameterization, then conditional on the data, the actively sampled parameters becomes independent from each other.

### HMC and Nuts

To understand how HMC works, consider reading the following:

-   Hamiltonian Monte Carlo explained. [Link](http://arogozhnikov.github.io/2016/12/19/Markov_chain_monte_carlo.html)
-   Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo. arXiv preprint arXiv:1701.02434. [Link](https://arxiv.org/abs/1701.02434)
-   Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo. Journal of Machine Learning Research, 15(1), 1593-1623. [Link](http://www.jmlr.org/papers/volume15/hoffman14a/hoffman14a.pdf)

---
title: "Hierarchical model for sampling"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{my-vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup, warning=FALSE, warning=FALSE}
library(sampinfo)
library(bayesplot)
library(rstan)
library(knitr)
```

Prepare your data:
--------

To begin with, you have a dataset with many participants, each finishing a same set of questions. Within each questions, there are several time points, where they need to sample from a set of C options.  
In the toy data, total question number $Q = 5$. The question ID for each sampling time point is stored in qID.
```{r}
subset(sampinfo::raw_data,sID==1)
```

And participant number $S=10$. The participant ID is listed in sID. 

```{r}
subset(sampinfo::raw_data,qID==1)
```

For each question, there are many time points, where participant sample from different options. Here time points are saved in tNo, and option sampled in cID. We have in total 7 different options (cID from 1 to 7).

Understand the model:
--------

Our goal is to predict the probabilities of sampling different features at each time point, using several known features of the options. 

In this example, we have three features: sameC, sameA and dist. 

* sameC: whether the option is in the same as the one sampled in the previous time point  
* sameA: whether they are decision congruent  
* dist: the semantic distance between the option and the previously sampled one  

Suppose that we use a single feature, and assume the effect can be captured by a linear multiplier $\beta$. Then for a trial with C options, we write the feature values into a vector \[\boldsymbol{X_{C\times1}} = [X_1, ..., X_C]^T\]
The predicted choice $Y$ thus follows a categorical distribution. \[Y \sim Cat(softmax(\boldsymbol{X}_{C\times1}\beta))\]. 

In a case with $K$ features, the effects can be captured by a vector $\boldsymbol{\beta_{K\times1}}$, and the feature matrix is $\boldsymbol{X_{C\times K}}$.\  

The predicted smapling probability \[Y \sim Cat(softmax(\boldsymbol{X_{C\times K}}\boldsymbol{\beta_{K\times1}})\] 

This is the model for predicting the sampling probability in one time point, of one question, completed by one participant. 

Now how about heterogeneity across paricipants and questions? We need a *hiearchical model* to address that. 

We first assume there is a grand mean for the feature effects, and denote it as $\boldsymbol{\mu^{\beta}_{K\times1}}$. 

The $S$ participants deviate from the grand mean on each of the $K$ features. Here the deviation can be written as $\boldsymbol{\delta^{\beta,S}_{K\times S}}$. 

Similarly, the deviations of the $Q$ questions on the features can be written as $\boldsymbol{\delta^{\beta,Q}_{K\times Q}}$. 

Each row of these deviation matrices forllows a centered normal distribution. Therefore, we need scale parameters to quantify these deviations. For participant and question level deviations, the scale parameters are $\boldsymbol{\sigma^{\beta,S}_{K\times 1}}$ and $\boldsymbol{\sigma^{\beta,Q}_{K\times 1}}$ respectively.


Therefore, the parameters we need to fit include: $\boldsymbol{\mu^{\beta}_{K\times1}}$, $\boldsymbol{\sigma^{\beta,S}_{K\times 1}}$, $\boldsymbol{\sigma^{\beta,Q}_{K\times 1}}$, $\boldsymbol{\delta^{\beta,S}_{K\times S}}$, and $\boldsymbol{\delta^{\beta,Q}_{K\times Q}}$. For participant $s$ and question $q$, the predicted samplng choice 
\[Y \sim Cat(softmax(\boldsymbol{X_{C\times K}}(\boldsymbol{\mu^{\beta}_{K\times1}} + \boldsymbol{\delta^{\beta,S}_{*,s}} + \boldsymbol{\delta^{\beta,Q}_{*,q}}))\]
 
Beyond these we can also have decay, and different base rates for different questions.

Now let's take a break from maths...

And take a look at the code.

```{r}
sampinfo::getStanCode(9)
```

Try fitting your model to the data:
--------

We will use a vector of feature names to specify the which features to be included in the model. For example:
```{r}
beta = c('sameC', 'sameA', 'dist')
```

If you need a baseline model without any features, set
```{r}
beta = c()
```


Next we need to decide whether we want decay in our model.  
If we assume there is no memory at all (i.e., Markov property), then we set
```{r}
deltaM_value = 0
```


On the other extreme, if no information is ever forgotten within a question (i.e., no decay). We set
```{r}
deltaM_value = 1
```


Finally we can estimate delta as a free parameter. Suppose we need a free parameter for delta, then we can specify
```{r}
deltaM_value = 9
```
With the data, and the model specifications, we can now run the following to fit the model.   
To make sure that the samples from the posterior distributions are not biased, we need to make sure the chains converge. Only then are we sampling the stationary distribution of the Marco chain. Because many models take minutes/hours/days/... weeks to run, we need to choose a warm-up number carefully. Try some warmup number, test convergence, and repeat until things look good. 
```{r, warning=FALSE}
beta = c('sameA', 'dist')
stan_data_fit = sampinfo::getStanFit(beta = beta, deltaM_value = 9, 
                           option_num = 7, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 200,chain_num = 4,warmup_num = 100, core_num=4,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE)
```


Check model convergence using traceplots.

```{r, dpi=300,   fig.height = 6, fig.width = 8, fig.align = "center", out.width = "85%"}
posterior_w = rstan::extract(stan_data_fit$stan_fit, inc_warmup = TRUE, permuted = FALSE)
dim(posterior_w)
lapply(dimnames(posterior_w), head)
bayesplot::color_scheme_set("viridis")
pars0 = c("beta_mu[1]","beta_s_sd[1]","beta_q_sd[1]","beta_s_raw[1,10]","deltaM","alpha[3,1]")
bayesplot::mcmc_trace(posterior_w, pars = pars0, 
                      n_warmup=100, size = 0.5,  
                      facet_args = list(nrow = 3))
```

Another useful plot, which shows you if the ranks of the samples are mixed well among chains.

```{r, dpi=300,   fig.height = 6, fig.width = 8, fig.align = "center", out.width = "85%"}
bayesplot::mcmc_rank_overlay(posterior_w, pars = pars0)
```

Also check gelman-rubin rhat.
```{r, warning=FALSE}
rhat = sampinfo::getRhat(stan_data_fit$stan_fit)
print(rhat[rhat>1.05&!is.na(rhat)])
```

Run the model:
--------

Now take the (conservative) warm-up number, and choose the total iteration number so that the sample size will be enough for your statistics of interests.
```{r, warning=FALSE}
stan_data_fit = sampinfo::getStanFit(beta= c('sameA', 'dist'), deltaM_value = 9, 
                           option_num = 7, format_data = format_data,
                           save_model_file = NULL, init_values="random",
                           iter_num = 1200,chain_num = 4,warmup_num = 300, core_num=4,
                           adapt_delta=0.9, stepsize = 0.1, max_treedepth = 10,
                           refresh=1000, save_warmup = TRUE)
```


Check the model diagnostics and goodness of fit:
--------

ESS is effective sample size, 1000 looks good in my opinion.

```{r}
model_diag = sampinfo::getModelDiag(stan_data_fit = stan_data_fit)
print(model_diag)
```

DIC and WAIC can be used for model comparisons. Lisheng also mentioned bayes factors, LOO-CV. 

Good reference:    
Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. Statistics and computing, 24(6), 997-1016. [Link](https://arxiv.org/abs/1307.5928)

Summarize and plot parameters
--------
```{r}
para = sampinfo::getParaSummary(stan_data_fit = stan_data_fit)
print(round(para[1:3,],1))
```


Here are the question level betas (grand mean plus question-level deviance)
```{r}
print(round(para[stringr::str_detect(rownames(para),'qmean'),c(1:3,8,13:14)],1))
```


and the participant level betas (grand mean plus participant-level deviance)

```{r}
print(round(para[stringr::str_detect(rownames(para),'smean'),c(1:3,8,13:14)],1))
```

Also option base rates in differnet questions.
```{r}
print(round(para[stringr::str_detect(rownames(para),'alpha'),c(1:3,8,13:14)],1)[1:20,])
```


There are many ways to visualize the samples. For example, check the bayesplot package [Link](https://mc-stan.org/bayesplot/articles/plotting-mcmc-draws.html)

```{r, dpi=300,   fig.height = 3.6, fig.width = 3.6, fig.align = "center", out.width = "50%"}
posterior = rstan::extract(stan_data_fit$stan_fit, inc_warmup = TRUE, permuted = FALSE)
dim(posterior)
head(dimnames(posterior)[[3]])
## run changeBetaNames to make the parameter names understandable (so that we dont need to deal with matrix index)
dimnames(posterior)[[3]] = sampinfo::changeBetaNames(dimnames(posterior)[[3]],
                            beta = beta,
                            Q = max(format_data$qID), 
                            S = max(format_data$sID))

bayesplot::mcmc_intervals(posterior, regex_pars = "^beta_dist_qmean")
```


Additional diagnostics in Rstan
--------

There are other diagnostics to pay attention to if you run your model with Stan (HMC and NUTS)
```{r}
model_diag = getModelDiag(stan_data_fit = stan_data_fit)
print(model_diag[,c('divergent', 'max_tree')])
```


For example you need to make sure there is no divergent transitions, which is not uncommon for hierarchical models. When it happens you often need to: 1) debug your rstan code. 2) try recovering your parameters 3) reparameterize your model or 4) increase alpha_delta, increase tree_maxdepth and decrease stepsize. I'll explain this later.   
If the maximum tree depth is hit in NUTS, we need to increase max_treedepth.  
RStan will send warnings regarding these issues (remember to read it if the models are fit on GPC).

Other topics:
--------
### Neal’s Funnel

Note that there are two ways to parameterize normal distributions (in hierarchical models): 

* Centered parameterization,    
\[\gamma \sim norml(\gamma_{mu}, \gamma_{sd})\]
* Non-centered parameterization,    
\[\gamma_{raw} \sim norml(0, 1) \]   
\[\gamma = \gamma_{mu} + \gamma_{sd} * \gamma_{raw} \]


If your individual level data are not very informative (and thus the individual parameters, $\gamma$, are influenced by the group level priors $\gamma_{mu}$ and $\gamma_{sd}$ a lot), it is much better to use the non-centered approach.
Detailed explanations can be found here.  

* RStan manual 22.7 [Link](https://mc-stan.org/docs/2_23/stan-users-guide/reparameterization-section.html)
* Betancourt, M., & Girolami, M. (2015). Hamiltonian Monte Carlo for hierarchical models.   Current trends in Bayesian methodology with applications, 79(30), 2-4. [Link](https://arxiv.org/abs/1312.0906)
  
Intuitively, when $\gamma_{sd}$ is smalll, $\gamma$ tends to be very close to gamma_mu. When $\gamma_{sd}$ is large, $\gamma$ varies a lot. That forms a tunnel shape. The stepsize needed for exploring the 'neck' is much much smaller than that required for the 'body'. Hence there's no single stepsize that stan can use to sample efficiently.  
If we use the non-centered parameterization, then conditional on the data, the actively sampled parameters becomes independent from each other. 

### HMC and Nuts

To understand how HMC works, consider reading the following:  

- Hamiltonian Monte Carlo explained. [Link](http://arogozhnikov.github.io/2016/12/19/Markov_chain_monte_carlo.html)
- Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo. arXiv preprint arXiv:1701.02434. [Link](https://arxiv.org/abs/1701.02434)
- Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo. Journal of Machine Learning Research, 15(1), 1593-1623. [Link](http://www.jmlr.org/papers/volume15/hoffman14a/hoffman14a.pdf)




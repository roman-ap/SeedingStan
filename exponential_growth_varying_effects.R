##### Set-up
# Packages
library(cmdstanr)
library(posterior)
library(bayesplot)
library(extraDistr)
library(withr)
library(ggplot2)
library(latex2exp)

##### Simulation of observations
# Dimensions
T=22
# Objects for storing the simulated data
lar_raw_obs=ar_obs=array(dim=T)
lbr_raw_obs=br_obs=array(dim=T)
lfr_raw_obs=fr_obs=array(dim=T)
M_obs=array(dim=T)
A_obs=array(dim=T)
Npre_obs=array(dim=T) 
Npost_obs=array(dim=T)
### (Hyper)priors - we should start by eliciting (hyper)priors
# Propose your (hyper)priors 
pm_LN0 = 6
psd_LN0 = 0.5
pm_mu_lar = -1
psd_mu_lar = 0.04
pm_sigma_lar = 0.2
psd_sigma_lar = 0.04
pm_mu_lbr = -1.45
psd_mu_lbr = 0.04
pm_sigma_lbr = 0.2
psd_sigma_lbr = 0.04
pm_mu_lfr = 1
psd_mu_lfr = 0.04
pm_sigma_lfr = 0.2
psd_sigma_lfr = 0.04
### Realisation
# Draw from the (hyper)priors such that A > 0
with_seed(87L,{
  N0_obs <- floor(rlnorm(n = 1, pm_LN0, psd_LN0))
  lar_raw_obs <- rnorm(n = T, 0, 1)
  mu_lar_obs <- rnorm(n = 1, pm_mu_lar, psd_mu_lar)
  sigma_lar_obs <- rtnorm(n = 1, pm_sigma_lar, psd_sigma_lar, a=0)
  lbr_raw_obs <- rnorm(n = T, 0, 1)
  mu_lbr_obs <- rnorm(n = 1, pm_mu_lbr, psd_mu_lbr)
  sigma_lbr_obs <- rtnorm(n = 1, pm_sigma_lbr, psd_sigma_lbr, a=0)
  lfr_raw_obs <- rnorm(n = T, 0, 1)
  mu_lfr_obs <- rnorm(n = 1, pm_mu_lfr, psd_mu_lfr)
  sigma_lfr_obs <- rtnorm(n=1, pm_sigma_lfr, psd_sigma_lfr, a=0)
  for(t in 1:T){
    if(t==1){
      Npre_obs[t]=N0_obs
    }
    else{
      Npre_obs[t]=rpois(n = 1, fr_obs[t-1]*Npost_obs[t-1]) 
    }
    ar_obs[t] = exp(mu_lar_obs + sigma_lar_obs*lar_raw_obs[t])
    br_obs[t] = exp(mu_lbr_obs + sigma_lbr_obs*lbr_raw_obs[t])
    fr_obs[t] = exp(mu_lfr_obs + sigma_lfr_obs*lfr_raw_obs[t])
    M_obs[t] = rbinom(n = 1, size = Npre_obs[t], prob = 1-exp(-(ar_obs[t]+br_obs[t])))
    A_obs[t]=rbinom(n = 1, size = M_obs[t], prob = ar_obs[t]/(ar_obs[t]+br_obs[t]))
    Npost_obs[t] = Npre_obs[t] - M_obs[t] 
  }
})

##### Stan model
exp_model <- cmdstan_model("./exponential_growth_varying_effects.stan")
##### Input data for the Stan model
exp_data <- list(T = T,
                 A = A_obs,
                 pm_LN0 = pm_LN0,
                 psd_LN0 = psd_LN0,
                 pm_mu_lar = pm_mu_lar,
                 psd_mu_lar = psd_mu_lar,
                 pm_sigma_lar = pm_sigma_lar,
                 psd_sigma_lar = psd_sigma_lar,
                 pm_mu_lbr = pm_mu_lbr,
                 psd_mu_lbr = psd_mu_lbr,
                 pm_sigma_lbr = pm_sigma_lbr,
                 psd_sigma_lbr = psd_sigma_lbr,
                 pm_mu_lfr = pm_mu_lfr,
                 psd_mu_lfr = psd_mu_lfr,
                 pm_sigma_lfr = pm_sigma_lfr,
                 psd_sigma_lfr = psd_sigma_lfr)

##### This blocks aims at setting the benchmark
### True hyperparameters' values as initial values
tinit <- function(){
  list(
    N0_raw = Npre_obs[1] - A_obs[1],
    NB_raw = Npre_obs[2:T] - A_obs[2:T],
    U = (M_obs - A_obs)/(Npre_obs - A_obs),
    lar_raw = lar_raw_obs,
    mu_lar = mu_lar_obs,
    sigma_lar = sigma_lar_obs,
    lbr_raw = lbr_raw_obs,
    mu_lbr = mu_lbr_obs,
    sigma_lbr = sigma_lbr_obs,
    lfr_raw = lfr_raw_obs,
    mu_lfr = mu_lfr_obs,
    sigma_lfr = sigma_lfr_obs
  )
}
### Stan benchmark sampling
exp_model_bm <- exp_model$sample(data = exp_data,
                                 chains = 4,
                                 parallel_chains = 4,
                                 refresh = 100,
                                 iter_warmup = 1000,
                                 iter_sampling = 1000,
                                 init = list(tinit(), 
                                             tinit(),
                                             tinit(),
                                             tinit()),
                                 save_warmup = TRUE)
### Summary of benchmark
exp_model_bm$summary(variables = c("N0",
                                   "mu_lar", "sigma_lar",
                                   "mu_lbr", "sigma_lbr",
                                   "mu_lfr", "sigma_lfr"))
exp_model_bm$summary(variables = "lp__")

##### Seeding the initial values
### Functions for generating initial seeds 
ginit <- function(D, cutoff){
  # Objects for storing the prior predictive draws
  N0_prip=array(dim=D)
  mu_lar_prip=sigma_lar_prip=array(dim=D)
  lar_raw_prip=ar_prip=array(dim=c(D,T))
  mu_lbr_prip=sigma_lbr_prip=array(dim=D)
  lbr_raw_prip=br_prip=array(dim=c(D,T))
  mu_lfr_prip=sigma_lfr_prip=array(dim=D)
  lfr_raw_prip=fr_prip=array(dim=c(D,T))
  A_prip=array(dim=c(D,T))
  M_prip=array(dim=c(D,T))
  Npre_prip=array(dim=c(D,T)) 
  Npost_prip=array(dim=c(D,T))
  good_draws=numeric()
  d = 1
  # Realisations
  while (d <= D) {
    N0_prip[d] <- floor(rlnorm(n = 1, pm_LN0, psd_LN0))
    mu_lar_prip[d] <- rnorm(n = 1, pm_mu_lar, psd_mu_lar)
    sigma_lar_prip[d] <- rtnorm(n= 1, pm_sigma_lar, psd_sigma_lar, a=0)
    lar_raw_prip[d,] <- rnorm(n = T, 0, 1)
    mu_lbr_prip[d] <- rnorm(n = 1, pm_mu_lbr, psd_mu_lbr)
    sigma_lbr_prip[d] <- rtnorm(n = 1, pm_sigma_lbr, psd_sigma_lbr, a=0)
    lbr_raw_prip[d,] <- rnorm(n = T, 0, 1)
    mu_lfr_prip[d] <- rnorm(n = 1, pm_mu_lfr, psd_mu_lfr)
    sigma_lfr_prip[d] <- rtnorm(n = 1, pm_sigma_lfr, psd_sigma_lfr, a=0)
    lfr_raw_prip[d,] <- rnorm(n = T, 0, 1)
    for(t in 1:T){
      if(t==1){
        Npre_prip[d,t] = N0_prip[d]
      }
      else{
        Npre_prip[d,t] <- rpois(n = 1, fr_prip[d,t-1]*Npost_prip[d,t-1]) 
      }
      ar_prip[d,t] = exp(mu_lar_prip[d] + sigma_lar_prip[d]*lar_raw_prip[d,t])
      br_prip[d,t] = exp(mu_lbr_prip[d] + sigma_lbr_prip[d]*lbr_raw_prip[d,t])
      fr_prip[d,t] = exp(mu_lfr_prip[d] + sigma_lfr_prip[d]*lfr_raw_prip[d,t])
      M_prip[d,t] <- rbinom(n = 1, size = Npre_prip[d,t], prob = 1-exp(-(ar_prip[d,t]+br_prip[d,t])))
      A_prip[d,t] <- rbinom(n = 1, size = M_prip[d,t], prob = ar_prip[d,t]/(ar_prip[d,t]+br_prip[d,t]))
      Npost_prip[d,t] = Npre_prip[d,t] - M_prip[d,t]
    }
    if(dist(rbind(as.vector(A_obs), as.vector(A_prip[d,])), method = "maximum")<cutoff){
      print(d)
      good_draws <- append(good_draws,d)
      d = d+1
    }
  }
  # Objects for storing the mean of prior predictive draws
  N0_pripm=array(dim=1)
  mu_lar_pripm=sigma_lar_pripm=array(dim=1)
  lar_raw_pripm=ar_pripm=array(dim=T)
  mu_lbr_pripm=sigma_lbr_pripm=array(dim=1)
  lbr_raw_pripm=br_pripm=array(dim=T)
  mu_lfr_pripm=sigma_lfr_pripm=array(dim=1)
  lfr_raw_pripm=fr_pripm=array(dim=T)
  A_pripm=array(dim=T)
  M_pripm=array(dim=T)
  Npre_pripm=array(dim=T) 
  Npost_pripm=array(dim=T)
  # Means
  N0_pripm = mean(N0_prip[good_draws])
  mu_lar_pripm = mean(mu_lar_prip[good_draws])
  sigma_lar_pripm = mean(sigma_lar_prip[good_draws])
  mu_lbr_pripm = mean(mu_lbr_prip[good_draws])
  sigma_lbr_pripm = mean(sigma_lbr_prip[good_draws])
  mu_lfr_pripm = mean(mu_lfr_prip[good_draws])
  sigma_lfr_pripm = mean(sigma_lfr_prip[good_draws])
    for(t in 1:T){
      lar_raw_pripm[t] = mean(lar_raw_prip[good_draws,t])
      ar_pripm[t] = mean(ar_prip[good_draws,t])
      lbr_raw_pripm[t] = mean(lbr_raw_prip[good_draws,t])
      br_pripm[t] = mean(br_prip[good_draws,t])
      lfr_raw_pripm[t] = mean(lfr_raw_prip[good_draws,t])
      fr_pripm[t] = mean(fr_prip[good_draws,t])
      M_pripm[t] = mean(M_prip[good_draws,t])
      A_pripm[t] = mean(A_prip[good_draws,t])
      Npost_pripm[t] = mean(Npost_prip[good_draws,t])
      Npre_pripm[t] = mean(Npre_prip[good_draws,t])
    }
  return(
    list(
      N0_raw = Npre_pripm[1] - A_pripm[1],
      NB_raw = Npre_pripm[2:T] - A_pripm[2:T],
      U = (M_pripm - A_pripm)/(Npre_pripm - A_pripm),
      mu_lar = mu_lar_pripm,
      sigma_lar = sigma_lar_pripm,
      lar_raw = lar_raw_pripm,
      mu_lbr = mu_lbr_pripm,
      sigma_lbr = sigma_lbr_pripm,
      lbr_raw = lbr_raw_pripm,
      mu_lfr = mu_lfr_pripm,
      sigma_lfr = sigma_lfr_pripm,
      lfr_raw = lfr_raw_pripm
    )
  )
}
robust_optimize <- function(stan_model, data, D, cutoff) {
  attempt <- 1
  success <- FALSE
  result <- NULL
  while(!success) {
    message("Attempt ", attempt)
    result <- tryCatch(
      expr = {
        out <- capture.output(
          map <- stan_model$optimize(data = data,
                                     init = list(ginit(D, cutoff)),
                                     jacobian = TRUE,
                                     iter = 10000
          ),
          type = "message"
        )
        if (length(out) > 0) message(out[1])
        map
      },
      error = function(e) {
        message("Error on attempt ", attempt, ": ", e$message)
        return(NULL)
      }, 
      warning = function(w) {
        message("Warning on attempt ", attempt, ": ", w$message)
        return(NULL)
      }
    )
    if (!is.null(result)) {
      success <- TRUE
      message("Optimization succeeded on attempt ", attempt)
    } 
    else {
      attempt <- attempt + 1
    }
  }
  return(result)
}

### Deploying approximations to get better initial values
# Repeat until you get a successful run for the optimization algorithm
exp_model_map <- robust_optimize(exp_model,
                                 exp_data, 
                                 D=1, cutoff=10000)
exp_model_map$summary(variables = c("N0", 
                                    "mu_lar", "sigma_lar", 
                                    "mu_lbr", "sigma_lbr", 
                                    "mu_lfr", "sigma_lfr"))
# You may want to use the MAP as initial values for approximations
ecomod_v1_la <- exp_model$laplace(data = ecomod_v1_data, init = ecomod_v1_map, jacobian = TRUE, draws = 10000)
ecomod_v1_la$summary(variables = c("N0", "mu_lhr", "sigma_lhr", "mu_lmr", "sigma_lmr", "mu_lfr", "sigma_lfr"))
ecomod_v1_pf <- exp_model$pathfinder(data = ecomod_v1_data, init = ecomod_v1_map, num_paths = 5, psis_resample = FALSE)
ecomod_v1_pf$summary(variables = c("N0", "mu_lhr", "sigma_lhr", "mu_lmr", "sigma_lmr", "mu_lfr", "sigma_lfr"))
### HMC (Stan) sampling - use any of the previous approximations
exp_model_hmc <- exp_model$sample(data = exp_data,
                                  chains = 4,
                                  parallel_chains = 4,
                                  refresh = 100,
                                  iter_warmup = 1000,
                                  iter_sampling = 1000,
                                  init = exp_model_map,
                                  save_warmup = TRUE)
### Summaries
exp_model_hmc$summary(variables = c("N0",
                                    "mu_lar", "sigma_lar",
                                    "mu_lbr", "sigma_lbr",
                                    "mu_lfr", "sigma_lfr"))
exp_model_hmc$summary(variables = "lp__")


##### Visualizations
### Data frame with draws 
ecomod_v1_draws <- exp_model_hmc$draws(format = "df")
### Initial population sizes 
ggplot(data.frame(draw=ecomod_v1_draws$`N0`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dlnorm, args = list(meanlog=pm_LN0, sdlog=psd_LN0), colour="black", linewidth=1) +
  xlim(0, exp(pm_LN0+0.5*(psd_LN0**2)) + 3*sqrt((exp(psd_LN0**2)-1)*exp(2*pm_LN0+(psd_LN0**2)))) +
  geom_vline(xintercept = N0_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $N^{pre}_{1}$"),
       y = TeX("Density"), x = TeX("$N^{pre}_{1}$$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
### Accident rates
# Mean of log harvest rates for juveniles
ggplot(data.frame(draw=ecomod_v1_draws$`mu_lhr[1]`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lhr[1], sd=psd_mu_lhr[1]), colour="black", linewidth=1) +
  xlim(pm_mu_lhr[1] - 6*psd_mu_lhr[1], pm_mu_lhr[1] + 6*psd_mu_lhr[1]) +
  geom_vline(xintercept = mu_lhr_obs[1], color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lhr_{j}}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lhr_{j}}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log harvest rate for juveniles
ggplot(data.frame(draw=ecomod_v1_draws$`sigma_lhr[1]`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lhr[1], sd=psd_sigma_lhr[1]), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lhr[1] - 6*psd_sigma_lhr[1]), pm_sigma_lhr[1] + 6*psd_sigma_lhr[1]) +
  geom_vline(xintercept = sigma_lhr_obs[1], color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lhr_{j}}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lhr_{j}}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
### Mortality rates
## Juveniles
# Mean of log harvest rates for juveniles
ggplot(data.frame(draw=ecomod_v1_draws$`mu_lbr[1]`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lbr[1], sd=psd_mu_lbr[1]), colour="black", linewidth=1) +
  xlim(pm_mu_lbr[1] - 6*psd_mu_lbr[1], pm_mu_lbr[1] + 6*psd_mu_lbr[1]) +
  geom_vline(xintercept = mu_lbr_obs[1], color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lbr_{j}}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lbr_{j}}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log harvest rate for juveniles
ggplot(data.frame(draw=ecomod_v1_draws$`sigma_lbr[1]`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lbr[1], sd=psd_sigma_lbr[1]), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lbr[1] - 6*psd_sigma_lbr[1]), pm_sigma_lbr[1] + 6*psd_sigma_lbr[1]) +
  geom_vline(xintercept = sigma_lbr_obs[1], color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lbr_{j}}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lbr_{j}}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
### Fertility rates 
# Mean of log fertility rates
ggplot(data.frame(draw=ecomod_v1_draws$`mu_lfr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lfr, sd=psd_mu_lfr), colour="black", linewidth=1) +
  xlim(pm_mu_lfr - 6*psd_mu_lfr, pm_mu_lfr + 6*psd_mu_lfr) +
  geom_vline(xintercept = mu_lfr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lfr}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lfr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log harvest rate for juveniles
ggplot(data.frame(draw=ecomod_v1_draws$`sigma_lfr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lfr, sd=psd_sigma_lfr), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lfr - 6*psd_sigma_lfr), pm_sigma_lfr + 6*psd_sigma_lfr) +
  geom_vline(xintercept = sigma_lfr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lfr}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lfr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

mcmc_pairs(ecomod_v1_draws,
           pars = c("mu_lbr[1]","mu_lbr[2]","mu_lhr[1]","mu_lhr[2]","mu_lfr"), 
           np = nuts_params(ecomod_v1_hmc))

# Harvest varying effects
ggplot(data.frame(draw = c(ecomod_v1_draws$`hr_prior[1,21]`, ecomod_v1_draws$`hr[1,21]`),
                  distribution = rep(c("prior predictive", "posterior"), each = 4000))) +
  geom_histogram(aes(x = draw, y = after_stat(density), fill= distribution), alpha = 0.75, bins = 100) +
  geom_vline(xintercept = hr_obs[1,21], color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior predictive and posterior distribution of $hr_{j,21}$"),
       y = TeX("Density"), x = TeX("$hr_{j,21}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  scale_fill_manual(values = c("prior predictive" = "black", "posterior" = "pink"))

ggplot(data.frame(draw = c(ecomod_v1_draws$`hr_prior[2,21]`, ecomod_v1_draws$`hr[2,21]`),
                  distribution = rep(c("prior predictive", "posterior"), each = 4000))) +
  geom_histogram(aes(x = draw, y = after_stat(density), fill= distribution), alpha = 0.75, bins = 100) +
  geom_vline(xintercept = hr_obs[2,21], color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior predictive and posterior distribution of $hr_{a,21}$"),
       y = TeX("Density"), x = TeX("$hr_{a,21}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  scale_fill_manual(values = c("prior predictive" = "black", "posterior" = "pink"))


# Abundancy
ggplot(data.frame(draw = c(ecomod_v1_draws$`Npre_prior[1,21]`, ecomod_v1_draws$`Npre[1,21]`),
                  distribution = rep(c("prior predictive", "posterior"), each = 4000))) +
  geom_histogram(aes(x = draw, y = after_stat(density), fill= distribution), alpha = 0.75, bins = 100) +
  geom_vline(xintercept = Npre_obs[1,21], color="blue") +
  xlim(0, 50000) +
  theme_minimal() +
  labs(title = TeX("Prior predictive and posterior distribution of $N^{pre}_{j,21}$"),
       y = TeX("Density"), x = TeX("$N^{pre}_{j,21}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  scale_fill_manual(values = c("prior predictive" = "black", "posterior" = "pink"))

ggplot(data.frame(draw = c(ecomod_v1_draws$`Npre_prior[2,21]`, ecomod_v1_draws$`Npre[2,21]`),
                  distribution = rep(c("prior predictive", "posterior"), each = 4000))) +
  geom_histogram(aes(x = draw, y = after_stat(density), fill= distribution), alpha = 0.75, bins = 100) +
  geom_vline(xintercept = Npre_obs[2,21], color="blue") +
  xlim(0, 50000) +
  theme_minimal() +
  labs(title = TeX("Prior predictive and posterior distribution of $N^{pre}_{a,21}$"),
       y = TeX("Density"), x = TeX("$N^{pre}_{a,21}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  scale_fill_manual(values = c("prior predictive" = "black", "posterior" = "pink"))

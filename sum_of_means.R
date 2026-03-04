##### Set-up
# Working directory
setwd("/home/roman/Nextcloud/Research/Stan/")
# Packages
library(mvtnorm)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(extraDistr)
library(withr)
library(ggplot2)

##### Simulation of data
# Dimensions
N = 1000
# Objects for storing the simulated data
mus_obs = array(dim = 2)
sum_mus_obs = array(dim = 1)
y_obs = array(dim = N)
### Priors - we should start by eliciting priors
var = c(1, 1)
cor = 0.5 
mu = c(0, 0) 
S = matrix(c(var[1],cor*sqrt(var[1]*var[2]),cor*sqrt(var[1]*var[2]),var[2]),
           nrow=2, ncol=2, byrow=T)
# Propose your priors
pm_mus = mu 
pcov_mus = S
psd_mu = 1
### Realisation
# Draw from the (hyper)priors
with_seed(8L, {
  mus_obs <- rmvnorm(1, pm_mus, pcov_mus)
  mu_obs = sum(mus_obs)
  y_obs <- rnorm(n = N, mu_obs, psd_mu)
})
##### Stan model
sum_mod <- cmdstan_model("./sum_of_means.stan")
##### Input data for the Stan model
sum_data <- list(
  N = N,
  y = y_obs,
  pm_mus = pm_mus,
  pcov_mus = pcov_mus,
  psd_mu = psd_mu
)
### Stan benchmark sampling
sum_fit <- sum_mod$sample(
  data = sum_data,
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  iter_warmup = 1000,
  iter_sampling = 1000,
  save_warmup = TRUE
)
### Summary of benchmark
sum_fit$summary(variables = c("mu", "mus"))
sum_fit$summary(variables = "lp__")


sum_draws <- sum_fit$draws(format = "df")

mcmc_pairs(sum_draws, pars = c("mu", "mus[1]", "mus[2]"), np = nuts_params(sum_fit))

mcmc_hist(sum_fit$draws("mu")) +
  geom_vline(xintercept = mu_obs)
mcmc_hist(sum_fit$draws("mus[1]")) +
  geom_vline(xintercept = mus_obs[1])
mcmc_hist(sum_fit$draws("mus[2]")) +
  geom_vline(xintercept = mus_obs[2])

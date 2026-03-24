##### Set-up
# Packages
library(cmdstanr)
library(posterior)
library(bayesplot)
library(extraDistr)
library(withr)
library(ggplot2)

##### Simulation of data
# Dimensions
N = 20
# Objects for storing the simulated data
mu_obs = sigma_obs = array(dim = 1)
x_obs = array(dim = N)
y_obs = array(dim = N)
### (Hyper)priors
pm_mu = 5 
psd_mu = 0.1
pm_sigma = 0
psd_sigma = 0.05
sd_y = 0.01
### Realisation
# Draw from the (hyper)priors
with_seed(8L, {
  mu_obs <- rnorm(n = 1, pm_mu, psd_mu)
  sigma_obs <- rtnorm(n = 1, pm_sigma, psd_sigma, a = 0)
  x_obs <- rnorm(n = N, mu_obs, sigma_obs)
  y_obs <- rnorm(n = N, x_obs, sd_y)
})
##### Stan model
randommeans_model <- cmdstan_model("./random_means.stan")
##### Input data for the Stan model
randommeans_data <- list(
  N = N,
  y = y_obs,
  pm_mu = pm_mu,
  psd_mu = psd_mu,
  pm_sigma = pm_sigma,
  psd_sigma = psd_sigma
)
### Stan benchmark sampling 
randommeans_map <- randommeans_model$optimize(data = randommeans_data,
                                              jacobian = TRUE,
                                              iter = 10000)
randommeans_fit <- randommeans_model$sample(
  data = randommeans_data,
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  iter_warmup = 1000,
  iter_sampling = 1000,
  #init = randommeans_map,
  save_warmup = TRUE
)
### Summary of benchmark
randommeans_fit$summary(variables = c("mu", "sigma", "x"))
randommeans_fit$summary(variables = "lp__")


randommeans_draws <- randommeans_fit$draws(format = "df")

mcmc_pairs(randommeans_draws, pars = c("mu", "sigma"), np = nuts_params(randommeans_fit))

mcmc_hist(randommeans_fit$draws("mu")) +
  geom_vline(xintercept = mu_obs)
mcmc_hist(randommeans_fit$draws("sigma")) +
  geom_vline(xintercept = sigma_obs)
mcmc_hist(randommeans_fit$draws("x[1]")) +
  geom_vline(xintercept = x_obs[1])
mcmc_hist(randommeans_fit$draws("x[2]")) +
  geom_vline(xintercept = x_obs[2])
mcmc_hist(randommeans_fit$draws("x[3]")) +
  geom_vline(xintercept = x_obs[3])
mcmc_hist(randommeans_fit$draws("x[4]")) +
  geom_vline(xintercept = x_obs[4])
mcmc_hist(randommeans_fit$draws("x[5]")) +
  geom_vline(xintercept = x_obs[5])
mcmc_hist(randommeans_fit$draws("x[6]")) +
  geom_vline(xintercept = x_obs[6])
mcmc_hist(randommeans_fit$draws("x[7]")) +
  geom_vline(xintercept = x_obs[7])
mcmc_hist(randommeans_fit$draws("x[8]")) +
  geom_vline(xintercept = x_obs[8])
mcmc_hist(randommeans_fit$draws("x[9]")) +
  geom_vline(xintercept = x_obs[9])
mcmc_hist(randommeans_fit$draws("x[10]")) +
  geom_vline(xintercept = x_obs[10])
mcmc_hist(randommeans_fit$draws("x[11]")) +
  geom_vline(xintercept = x_obs[11])
mcmc_hist(randommeans_fit$draws("x[12]")) +
  geom_vline(xintercept = x_obs[12])
mcmc_hist(randommeans_fit$draws("x[13]")) +
  geom_vline(xintercept = x_obs[13])
mcmc_hist(randommeans_fit$draws("x[14]")) +
  geom_vline(xintercept = x_obs[14])
mcmc_hist(randommeans_fit$draws("x[15]")) +
  geom_vline(xintercept = x_obs[15])
mcmc_hist(randommeans_fit$draws("x[16]")) +
  geom_vline(xintercept = x_obs[16])
mcmc_hist(randommeans_fit$draws("x[17]")) +
  geom_vline(xintercept = x_obs[17])
mcmc_hist(randommeans_fit$draws("x[18]")) +
  geom_vline(xintercept = x_obs[18])
mcmc_hist(randommeans_fit$draws("x[19]")) +
  geom_vline(xintercept = x_obs[19])
mcmc_hist(randommeans_fit$draws("x[20]")) +
  geom_vline(xintercept = x_obs[20])

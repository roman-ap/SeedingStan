##### Set-up
# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(nimble)
library(extraDistr)
library(withr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

##### Data simulation
# Dimension
Ty <- 20
# Objects for storing the simulated observations
N0_obs=array(dim=1)
mu_lmr_obs=sigma_lmr_obs=array(dim=1)
lmr_raw_obs=lmr_obs=mr_obs=array(dim=Ty)
mu_lfr_obs=sigma_lfr_obs=array(dim=1)
lfr_raw_obs=lfr_obs=fr_obs=array(dim=Ty)
M_obs=array(dim=Ty)
Npre_obs=array(dim=Ty) 
Npost_obs=array(dim=Ty)
### (Hyper)priors
# Propose your (hyper)priors 
pm_LN0 = 6
pm_mu_lmr = -1
psd_mu_lmr = 0.1
pm_sigma_lmr = 0
psd_sigma_lmr = 0.05
pm_mu_lfr = 0.5
psd_mu_lfr = 0.1
pm_sigma_lfr = 0
psd_sigma_lfr = 0.05
### Realisation
# Draw from the (hyper)priors such that M > 0
with_seed(87L,{
  N0_obs <- rpois(n = 1, exp(pm_LN0))
  lmr_raw_obs <- rnorm(n = Ty, 0, 1)
  mu_lmr_obs <- rnorm(n = 1, pm_mu_lmr, psd_mu_lmr)
  sigma_lmr_obs <- rtnorm(n = 1, pm_sigma_lmr, psd_sigma_lmr, a=0)
  lmr_obs = mu_lmr_obs + sigma_lmr_obs*lmr_raw_obs
  mr_obs = exp(lmr_obs)
  lfr_raw_obs <- rnorm(n = Ty, 0, 1)
  mu_lfr_obs <- rnorm(n = 1, pm_mu_lfr, psd_mu_lfr)
  sigma_lfr_obs <- rtnorm(n = 1, pm_sigma_lfr, psd_sigma_lfr, a=0)
  lfr_obs = mu_lfr_obs + sigma_lfr_obs*lfr_raw_obs
  fr_obs = exp(lfr_obs)
  for(t in 1:Ty){
    if(t==1){
      Npre_obs[t] = N0_obs
    }
    else{
      Npre_obs[t] <- rpois(n = 1, fr_obs[t-1]*Npost_obs[t-1]) 
    }
    M_obs[t] <- rbinom(n = 1, size = Npre_obs[t], prob = 1-exp(-mr_obs[t]))
    Npost_obs[t] = Npre_obs[t] - M_obs[t] 
  }
})

##### Bayesian Inference 
### NIMBLE code 
ss_code <- nimbleCode({
  # ----- Priors -----
  N0 ~ dpois(lambda = exp(pm_LN0))
  mu_lmr ~ dnorm(pm_mu_lmr, sd = psd_mu_lmr)
  sigma_lmr ~ T(dnorm(pm_sigma_lmr, sd = psd_sigma_lmr), 0, ) 
  mu_lfr ~ dnorm(pm_mu_lfr, sd = psd_mu_lfr)
  sigma_lfr ~ T(dnorm(pm_sigma_lfr, sd = psd_sigma_lfr), 0, )
  # ----- Parameters -----
  for(t in 1:T){
    lmr_raw[t] ~ dnorm(0, sd = 1)
    lmr[t] <- mu_lmr + sigma_lmr * lmr_raw[t]
    mr[t] <- exp(lmr[t])
    lfr_raw[t] ~ dnorm(0, sd = 1)
    lfr[t] <- mu_lfr + sigma_lfr * lfr_raw[t]
    fr[t] <- exp(lfr[t])
  }
  # ----- State process -----
  Npre[1] <- N0
  M[1] ~ dbinom(1-exp(-mr[1]), Npre[1])
  Npost[1] <- Npre[1] - M[1]
  for(t in 2:T){
    Npre[t] ~ dpois(fr[t-1]*Npost[t-1])
    M[t] ~ dbinom(1-exp(-mr[t]), Npre[t])
    Npost[t] <- Npre[t] - M[t]
  }
})
### NIMBLE constants
ss_constants <- list(T = Ty,
                     pm_LN0 = pm_LN0,
                     pm_mu_lmr = pm_mu_lmr,
                     psd_mu_lmr = psd_mu_lmr,
                     pm_sigma_lmr = pm_sigma_lmr,
                     psd_sigma_lmr = psd_sigma_lmr,
                     pm_mu_lfr = pm_mu_lfr,
                     psd_mu_lfr = psd_mu_lfr,
                     pm_sigma_lfr = pm_sigma_lfr,
                     psd_sigma_lfr = psd_sigma_lfr)
### NIMBLE data
ss_data <- list(M = M_obs)
### NIMBLE without seeds (default)
ss_inits_van <- list()
# NIMBLE Model
ss_model_van <- nimbleModel(code = ss_code,
                            constants = ss_constants,
                            data = ss_data,
                            inits = ss_inits_van)
ss_cmodel_van <- compileNimble(ss_model_van)
ss_mcmc_van <- configureMCMC(ss_model_van, print = TRUE)
ss_mcmc_van$addMonitors(c("mu_lmr", "sigma_lmr", "mr",
                          "mu_lfr", "sigma_lfr", "fr",
                          "N0", "Npre", "Npost"))
ss_modelMCMC_van <- buildMCMC(ss_mcmc_van)
ss_cmodelMCMC_van <- compileNimble(ss_modelMCMC_van, project = ss_model_van)
ss_samples_van <- runMCMC(ss_cmodelMCMC_van,
                          niter = 10000,
                          nburnin = 5000,
                          nchains = 4)
# Summaries
### Visualizations
# Data frame with posterior draws when not seeding 
ss_van_draws <- do.call(rbind, lapply(ss_samples_van, as.data.frame))
## Initial population size
gg_van_N0 <- ggplot(data.frame(draw=ss_van_draws$`N0`)) +
  geom_bar(aes(x=draw, y=after_stat(prop)), color="pink", fill="pink") +
  stat_function(fun = dpois, args = list(lambda=exp(pm_LN0)), geom = "step", colour="black", linewidth=1) +
  #xlim(0, 1000) +
  geom_vline(xintercept = N0_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $N^{pre}_{1}$"),
       y = TeX("Density"), x = TeX("$N^{pre}_{1}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
grid.arrange(gg_van_N0)
## Hyperparameters
# Mean of log mortality rates
gg_van_mu_lmr <-ggplot(data.frame(draw=ss_van_draws$`mu_lmr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lmr, sd=psd_mu_lmr), colour="black", linewidth=1) +
  xlim(min(pm_mu_lmr - 6*psd_mu_lmr, min(ss_van_draws$`mu_lmr`)), max(pm_mu_lmr + 6*psd_mu_lmr, max(ss_van_draws$`mu_lmr`))) +
  geom_vline(xintercept = mu_lmr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lmr}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lmr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log mortality rates
gg_van_sigma_lmr <- ggplot(data.frame(draw=ss_van_draws$`sigma_lmr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lmr, sd=psd_sigma_lmr), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lmr - 6*psd_sigma_lmr), max(pm_sigma_lmr + 6*psd_sigma_lmr, max(ss_van_draws$`sigma_lmr`))) +
  geom_vline(xintercept = sigma_lmr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lbr}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lbr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Mean of log fertility rates
gg_van_mu_lfr <- ggplot(data.frame(draw=ss_van_draws$`mu_lfr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lfr, sd=psd_mu_lfr), colour="black", linewidth=1) +
  xlim(min(pm_mu_lfr - 6*psd_mu_lfr, min(ss_van_draws$`mu_lfr`)), max(pm_mu_lfr + 6*psd_mu_lfr, max(ss_van_draws$`mu_lfr`))) +
  geom_vline(xintercept = mu_lfr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lfr}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lfr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log fertility rates
gg_van_sigma_lfr <- ggplot(data.frame(draw=ss_van_draws$`sigma_lfr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lfr, sd=psd_sigma_lfr), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lfr - 6*psd_sigma_lfr), max(pm_sigma_lfr + 6*psd_sigma_lfr, max(ss_van_draws$`sigma_lfr`))) +
  geom_vline(xintercept = sigma_lfr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lfr}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lfr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Plots
grid.arrange(gg_van_mu_lfr, gg_van_sigma_lfr,
             gg_van_mu_lmr, gg_van_sigma_lmr,
             nrow = 2, ncol = 2,
             layout_matrix = rbind(c(1,2), 
                                   c(3,4)))
### Abundance
# Figures
pre_abundancy_van <- ss_van_draws |>
  select(starts_with("Npre")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(index = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(parameter, index) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop") |>
  arrange(index) |>
  select(-index)
post_abundancy_van <- ss_van_draws |>
  select(starts_with("Npost")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(index = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(parameter, index) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop") |>
  arrange(index) |>
  select(-index)
# Plots
gg_van_Npre <- ggplot(pre_abundancy_van, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = Npre_obs)) +
  ylim(0, max(pre_abundancy_van$q97.5)) +
  theme_minimal() +
  labs(title = TeX("Posterior of the pre-mortality population"),
       y = TeX("Total population"), x = "Year") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
gg_van_Npost <- ggplot(post_abundancy_van, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = Npost_obs)) +
  ylim(0, max(pre_abundancy_van$q97.5)) +
  theme_minimal() +
  labs(title = TeX("Posterior of the post-mortality population"),
       y = TeX("Total population"), x = "Year") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_van_Npre, gg_van_Npost, nrow = 1, ncol = 2)
### Fertility rates
# Figures
fr_van <- ss_van_draws |>
  select(starts_with("fr")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(index = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(parameter, index) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop") |>
  arrange(index) |>
  select(-index)
# Plots
gg_van_fr <- ggplot(fr_van, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = fr_obs)) +
  ylim(0, 4) +
  theme_minimal() +
  labs(title = TeX("Fertility rates"),
       y = TeX("Fertility rate"), x = "Year") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_van_fr, nrow = 1, ncol = 1)
## Mortality rates
mr_van <- ss_van_draws |>
  select(starts_with("mr")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(index = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(parameter, index) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop") |>
  arrange(index) |>
  select(-index)
# Plots
gg_van_mr <- ggplot(mr_van, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = mr_obs)) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = TeX("Mortality rates"),
       y = TeX("Mortality rate"), x = "Year") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_van_mr, nrow = 1, ncol = 1)

### Seeding Stan
## Functions for generating initial seeds 
ginit <- function(D, cutoff){
  # Objects for storing the prior predictive draws
  N0_prip=array(dim=D)
  mu_lmr_prip=sigma_lmr_prip=array(dim=D)
  lmr_raw_prip=mr_prip=array(dim=c(D,Ty))
  mu_lfr_prip=sigma_lfr_prip=array(dim=D)
  lfr_raw_prip=fr_prip=array(dim=c(D,Ty))
  M_prip=array(dim=c(D,Ty))
  Npre_prip=array(dim=c(D,Ty)) 
  Npost_prip=array(dim=c(D,Ty))
  good_draws=numeric()
  d = 1
  # Realisations
  while (d <= D) {
    N0_prip[d] <- rpois(n = 1, exp(pm_LN0))
    mu_lmr_prip[d] <- rnorm(n = 1, pm_mu_lmr, psd_mu_lmr)
    sigma_lmr_prip[d] <- rtnorm(n = 1, pm_sigma_lmr, psd_sigma_lmr, a=0)
    lmr_raw_prip[d,] <- rnorm(n = Ty, 0, 1)
    mu_lfr_prip[d] <- rnorm(n = 1, pm_mu_lfr, psd_mu_lfr)
    sigma_lfr_prip[d] <- rtnorm(n = 1, pm_sigma_lfr, psd_sigma_lfr, a=0)
    lfr_raw_prip[d,] <- rnorm(n = Ty, 0, 1)
    for(t in 1:Ty){
      if(t==1){
        Npre_prip[d,t] = N0_prip[d]
      }
      else{
        Npre_prip[d,t] <- rpois(n = 1, fr_prip[d,t-1]*Npost_prip[d,t-1]) 
      }
      mr_prip[d,t] = exp(mu_lmr_prip[d] + sigma_lmr_prip[d]*lmr_raw_prip[d,t])
      fr_prip[d,t] = exp(mu_lfr_prip[d] + sigma_lfr_prip[d]*lfr_raw_prip[d,t])
      M_prip[d,t] <- rbinom(n = 1, size = Npre_prip[d,t], prob = 1-exp(-mr_prip[d,t]))
      Npost_prip[d,t] = Npre_prip[d,t] - M_prip[d,t]
    }
    if(dist(rbind(as.vector(M_obs), as.vector(M_prip[d,])), method = "maximum")<cutoff){
      print(d)
      good_draws <- append(good_draws,d)
      d = d+1
    }
  }
  # Objects for storing the mean of prior predictive draws
  N0_pripm=array(dim=1)
  mu_lmr_pripm=sigma_lmr_pripm=array(dim=1)
  lmr_raw_pripm=mr_pripm=array(dim=Ty)
  mu_lfr_pripm=sigma_lfr_pripm=array(dim=1)
  lfr_raw_pripm=fr_pripm=array(dim=Ty)
  M_pripm=array(dim=Ty)
  Npre_pripm=array(dim=Ty) 
  Npost_pripm=array(dim=Ty)
  # Means
  N0_pripm = mean(N0_prip[good_draws])
  mu_lmr_pripm = mean(mu_lmr_prip[good_draws])
  sigma_lmr_pripm = mean(sigma_lmr_prip[good_draws])
  mu_lfr_pripm = mean(mu_lfr_prip[good_draws])
  sigma_lfr_pripm = mean(sigma_lfr_prip[good_draws])
    for(t in 1:Ty){
      lmr_raw_pripm[t] = mean(lmr_raw_prip[good_draws,t])
      mr_pripm[t] = mean(mr_prip[good_draws,t])
      lfr_raw_pripm[t] = mean(lfr_raw_prip[good_draws,t])
      fr_pripm[t] = mean(fr_prip[good_draws,t])
      M_pripm[t] = mean(M_prip[good_draws,t])
      Npost_pripm[t] = mean(Npost_prip[good_draws,t])
      Npre_pripm[t] = mean(Npre_prip[good_draws,t])
    }
  return(
    list(
      N0 = Npre_pripm[1],
      Npre = Npre_pripm[1:Ty],
      mu_lmr = mu_lmr_pripm,
      sigma_lmr = sigma_lmr_pripm,
      lmr_raw = lmr_raw_pripm,
      mu_lfr = mu_lfr_pripm,
      sigma_lfr = sigma_lfr_pripm,
      lfr_raw = lfr_raw_pripm
    )
  )
}
ss_inits_seeded <- list(ginit(D=1, cutoff=1000))
# NIMBLE Model
ss_model_seeded <- nimbleModel(code = ss_code,
                               constants = ss_constants,
                               data = ss_data,
                               inits = ss_inits_seeded)
ss_cmodel_seeded <- compileNimble(ss_model_seeded)
ss_mcmc_seeded <- configureMCMC(ss_model_seeded, print = TRUE)
ss_mcmc_seeded$addMonitors(c("mu_lmr", "sigma_lmr", "mr",
                             "mu_lfr", "sigma_lfr", "fr",
                             "N0", "Npre", "Npost"))
ss_modelMCMC_seeded <- buildMCMC(ss_mcmc_seeded)
ss_cmodelMCMC_seeded <- compileNimble(ss_modelMCMC_seeded, project = ss_model_seeded)
ss_samples_seeded <- runMCMC(ss_cmodelMCMC_seeded,
                             niter = 10000,
                             nburnin = 5000,
                             nchains = 4)
# Summaries

##### Visualizations
# Data frame with draws 
ss_seeded_draws <- do.call(rbind, lapply(ss_samples_seeded, as.data.frame))
### Hyperparameters
# Initial population size
gg_seeded_N0 <- ggplot(data.frame(draw=ss_seeded_draws$`N0`)) +
  geom_bar(aes(x=draw, y=after_stat(prop)), color="pink", fill="pink") +
  stat_function(fun = dpois, args = list(lambda=exp(pm_LN0)), geom = "step", colour="black", linewidth=1) +
  xlim(0, 1000) +
  geom_vline(xintercept = N0_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $N^{pre}_{1}$"),
       y = TeX("Density"), x = TeX("$N^{pre}_{1}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
grid.arrange(gg_seeded_N0)
# Mean of log mortality rates
gg_seeded_mu_lmr <-ggplot(data.frame(draw=ss_seeded_draws$`mu_lmr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lmr, sd=psd_mu_lmr), colour="black", linewidth=1) +
  xlim(min(pm_mu_lmr - 6*psd_mu_lmr, min(ss_seeded_draws$`mu_lmr`)), max(pm_mu_lmr + 6*psd_mu_lmr, max(ss_seeded_draws$`mu_lmr`))) +
  geom_vline(xintercept = mu_lmr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lmr}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lmr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log mortality rates
gg_seeded_sigma_lmr <- ggplot(data.frame(draw=ss_seeded_draws$`sigma_lmr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lmr, sd=psd_sigma_lmr), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lmr - 6*psd_sigma_lmr), max(pm_sigma_lmr + 6*psd_sigma_lmr, max(ss_seeded_draws$`sigma_lmr`))) +
  geom_vline(xintercept = sigma_lmr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lmr}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lmr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Mean of log fertility rates
gg_seeded_mu_lfr <- ggplot(data.frame(draw=ss_seeded_draws$`mu_lfr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lfr, sd=psd_mu_lfr), colour="black", linewidth=1) +
  xlim(min(pm_mu_lfr - 6*psd_mu_lfr, min(ss_seeded_draws$`mu_lfr`)), max(pm_mu_lfr + 6*psd_mu_lfr, max(ss_seeded_draws$`mu_lfr`))) +
  geom_vline(xintercept = mu_lfr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lfr}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lfr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log fertility rates
gg_seeded_sigma_lfr <- ggplot(data.frame(draw=ss_seeded_draws$`sigma_lfr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lfr, sd=psd_sigma_lfr), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lfr - 6*psd_sigma_lfr), max(pm_sigma_lfr + 6*psd_sigma_lfr, max(ss_seeded_draws$`sigma_lfr`))) +
  geom_vline(xintercept = sigma_lfr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lfr}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lfr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Plots
grid.arrange(gg_seeded_mu_lfr, gg_seeded_sigma_lfr,
             gg_seeded_mu_lmr, gg_seeded_sigma_lmr,
             nrow = 2, ncol = 2,
             layout_matrix = rbind(c(1,2), 
                                   c(3,4)))
### Abundance
# Figures
pre_abundancy_seeded <- ss_seeded_draws |>
  select(starts_with("Npre")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(index = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(parameter, index) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop") |>
  arrange(index) |>
  select(-index)
post_abundancy_seeded <- ss_seeded_draws |>
  select(starts_with("Npost")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(index = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(parameter, index) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop") |>
  arrange(index) |>
  select(-index)
# Plots
gg_seeded_Npre <- ggplot(pre_abundancy_seeded, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = Npre_obs)) +
  ylim(0, max(pre_abundancy_seeded$q97.5)) +
  theme_minimal() +
  labs(title = TeX("Posterior of the pre-mortality population"),
       y = TeX("Total population"), x = "Year") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
gg_seeded_Npost <- ggplot(post_abundancy_seeded, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = Npost_obs)) +
  ylim(0, max(pre_abundancy_seeded$q97.5)) +
  theme_minimal() +
  labs(title = TeX("Posterior of the post-mortality population"),
       y = TeX("Total population"), x = "Year") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_seeded_Npre, gg_seeded_Npost, nrow = 1, ncol = 2)
### Fertility rates
# Figures
fr_seeded <- ss_seeded_draws |>
  select(starts_with("fr")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(index = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(parameter, index) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop") |>
  arrange(index) |>
  select(-index)
# Plots
gg_seeded_fr <- ggplot(fr_seeded, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = fr_obs)) +
  ylim(0, 4) +
  theme_minimal() +
  labs(title = TeX("Fertility rates"),
       y = TeX("Fertility rate"), x = "Year") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_seeded_fr, nrow = 1, ncol = 1)
## Mortality rates
mr_seeded <- ss_seeded_draws |>
  select(starts_with("mr")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(index = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(parameter, index) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975),
    .groups = "drop") |>
  arrange(index) |>
  select(-index)
# Plots
gg_seeded_mr <- ggplot(mr_seeded, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = mr_obs)) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = TeX("Mortality rates"),
       y = TeX("Mortality rate"), x = "Year") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_seeded_mr, nrow = 1, ncol = 1)

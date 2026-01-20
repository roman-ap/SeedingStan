##### Set-up
# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(cmdstanr)
library(extraDistr)
library(withr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

##### Data simulation
# Dimension
Ty <- 20
# Objects for storing the simulated observations
mu_lar_obs=sigma_lar_obs=array(dim=1)
lar_raw_obs=ar_obs=array(dim=Ty)
mu_lbr_obs=sigma_lbr_obs=array(dim=1)
lbr_raw_obs=br_obs=array(dim=Ty)
mu_lfr_obs=sigma_lfr_obs=array(dim=1)
lfr_raw_obs=fr_obs=array(dim=Ty)
M_obs=array(dim=Ty)
A_obs=array(dim=Ty)
Npre_obs=array(dim=Ty) 
Npost_obs=array(dim=Ty)
### (Hyper)priors
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
  lar_raw_obs <- rnorm(n = Ty, 0, 1)
  mu_lar_obs <- rnorm(n = 1, pm_mu_lar, psd_mu_lar)
  sigma_lar_obs <- rtnorm(n = 1, pm_sigma_lar, psd_sigma_lar, a=0)
  lbr_raw_obs <- rnorm(n = Ty, 0, 1)
  mu_lbr_obs <- rnorm(n = 1, pm_mu_lbr, psd_mu_lbr)
  sigma_lbr_obs <- rtnorm(n = 1, pm_sigma_lbr, psd_sigma_lbr, a=0)
  lfr_raw_obs <- rnorm(n = Ty, 0, 1)
  mu_lfr_obs <- rnorm(n = 1, pm_mu_lfr, psd_mu_lfr)
  sigma_lfr_obs <- rtnorm(n = 1, pm_sigma_lfr, psd_sigma_lfr, a=0)
  for(t in 1:Ty){
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

##### Bayesian Inference 
### Stan model
expgro_model <- cmdstan_model("./exponential_growth_varying_effects.stan")
### Input data for the Stan model
expgro_data <- list(T = Ty,
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
### Stan without seeds fails
expgro_model$sample(data = expgro_data,
                    chains = 4,    
                    parallel_chains = 4,
                    refresh = 100,
                    iter_warmup = 1000,
                    iter_sampling = 1000,
                    save_warmup = TRUE)
### This blocks aims at setting the benchmark
# True hyperparameters' values as initial values
tinit <- function(){
  list(
    N0_raw = Npre_obs[1] - A_obs[1],
    NB_raw = Npre_obs[2:Ty] - A_obs[2:Ty],
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
# Stan benchmark sampling
expgro_bm <- expgro_model$sample(data = expgro_data,
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
# Summary of benchmark
expgro_bm$summary(variables = c("N0",
                                "mu_lar", "sigma_lar",
                                "mu_lbr", "sigma_lbr",
                                "mu_lfr", "sigma_lfr"))
expgro_bm$summary(variables = "lp__")

### Seeding Stan
## Functions for generating initial seeds 
ginit <- function(D, cutoff){
  # Objects for storing the prior predictive draws
  N0_prip=array(dim=D)
  mu_lar_prip=sigma_lar_prip=array(dim=D)
  lar_raw_prip=ar_prip=array(dim=c(D,Ty))
  mu_lbr_prip=sigma_lbr_prip=array(dim=D)
  lbr_raw_prip=br_prip=array(dim=c(D,Ty))
  mu_lfr_prip=sigma_lfr_prip=array(dim=D)
  lfr_raw_prip=fr_prip=array(dim=c(D,Ty))
  A_prip=array(dim=c(D,Ty))
  M_prip=array(dim=c(D,Ty))
  Npre_prip=array(dim=c(D,Ty)) 
  Npost_prip=array(dim=c(D,Ty))
  good_draws=numeric()
  d = 1
  # Realisations
  while (d <= D) {
    N0_prip[d] <- floor(rlnorm(n = 1, pm_LN0, psd_LN0))
    mu_lar_prip[d] <- rnorm(n = 1, pm_mu_lar, psd_mu_lar)
    sigma_lar_prip[d] <- rtnorm(n= 1, pm_sigma_lar, psd_sigma_lar, a=0)
    lar_raw_prip[d,] <- rnorm(n = Ty, 0, 1)
    mu_lbr_prip[d] <- rnorm(n = 1, pm_mu_lbr, psd_mu_lbr)
    sigma_lbr_prip[d] <- rtnorm(n = 1, pm_sigma_lbr, psd_sigma_lbr, a=0)
    lbr_raw_prip[d,] <- rnorm(n = Ty, 0, 1)
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
  lar_raw_pripm=ar_pripm=array(dim=Ty)
  mu_lbr_pripm=sigma_lbr_pripm=array(dim=1)
  lbr_raw_pripm=br_pripm=array(dim=Ty)
  mu_lfr_pripm=sigma_lfr_pripm=array(dim=1)
  lfr_raw_pripm=fr_pripm=array(dim=Ty)
  A_pripm=array(dim=Ty)
  M_pripm=array(dim=Ty)
  Npre_pripm=array(dim=Ty) 
  Npost_pripm=array(dim=Ty)
  # Means
  N0_pripm = mean(N0_prip[good_draws])
  mu_lar_pripm = mean(mu_lar_prip[good_draws])
  sigma_lar_pripm = mean(sigma_lar_prip[good_draws])
  mu_lbr_pripm = mean(mu_lbr_prip[good_draws])
  sigma_lbr_pripm = mean(sigma_lbr_prip[good_draws])
  mu_lfr_pripm = mean(mu_lfr_prip[good_draws])
  sigma_lfr_pripm = mean(sigma_lfr_prip[good_draws])
    for(t in 1:Ty){
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
      NB_raw = Npre_pripm[2:Ty] - A_pripm[2:Ty],
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
## Deploying approximations to get better initial values
# Maximum a posteriori
expgro_map <- robust_optimize(expgro_model,
                              expgro_data, 
                              D=1, cutoff=10000)
expgro_map$summary(variables = c("N0", 
                                 "mu_lar", "sigma_lar", 
                                 "mu_lbr", "sigma_lbr", 
                                 "mu_lfr", "sigma_lfr"))
### HMC (Stan) sampling
expgro_hmc <- expgro_model$sample(data = expgro_data,
                                  chains = 4,
                                  parallel_chains = 4,
                                  refresh = 100,
                                  iter_warmup = 1000,
                                  iter_sampling = 1000,
                                  init = expgro_map,
                                  save_warmup = TRUE)
# Summaries
expgro_hmc$summary(variables = c("N0",
                                 "mu_lar", "sigma_lar",
                                 "mu_lbr", "sigma_lbr",
                                 "mu_lfr", "sigma_lfr"))
expgro_hmc$summary(variables = "lp__")


##### Visualizations
# Data frame with draws 
expgro_draws <- expgro_hmc$draws(format = "df")
### Hyperparameters
# Initial population size
ggplot(data.frame(draw=expgro_draws$`N0`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dlnorm, args = list(meanlog=pm_LN0, sdlog=psd_LN0), colour="black", linewidth=1) +
  xlim(0, exp(pm_LN0+0.5*(psd_LN0**2)) + 3*sqrt((exp(psd_LN0**2)-1)*exp(2*pm_LN0+(psd_LN0**2)))) +
  geom_vline(xintercept = N0_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $N^{pre}_{1}$"),
       y = TeX("Density"), x = TeX("$N^{pre}_{1}$$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Mean of log accident rates
gg_mu_lar <- ggplot(data.frame(draw=expgro_draws$`mu_lar`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lar, sd=psd_mu_lar), colour="black", linewidth=1) +
  xlim(pm_mu_lar - 6*psd_mu_lar, pm_mu_lar + 6*psd_mu_lar) +
  geom_vline(xintercept = mu_lar_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lar}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lar}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log accident rates
gg_sigma_lar <- ggplot(data.frame(draw=expgro_draws$`sigma_lar`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lar, sd=psd_sigma_lar), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lar - 6*psd_sigma_lar), pm_sigma_lar + 6*psd_sigma_lar) +
  geom_vline(xintercept = sigma_lar_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lar}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lar}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Mean of log background mortality rates
gg_mu_lbr <-ggplot(data.frame(draw=expgro_draws$`mu_lbr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lbr, sd=psd_mu_lbr), colour="black", linewidth=1) +
  xlim(pm_mu_lbr[1] - 6*psd_mu_lbr, pm_mu_lbr + 6*psd_mu_lbr) +
  geom_vline(xintercept = mu_lbr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lbr}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lbr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log background mortality rates
gg_sigma_lbr <- ggplot(data.frame(draw=expgro_draws$`sigma_lbr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lbr, sd=psd_sigma_lbr), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lbr - 6*psd_sigma_lbr), pm_sigma_lbr + 6*psd_sigma_lbr) +
  geom_vline(xintercept = sigma_lbr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lbr}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lbr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Mean of log fertility rates
gg_mu_lfr <- ggplot(data.frame(draw=expgro_draws$`mu_lfr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu_lfr, sd=psd_mu_lfr), colour="black", linewidth=1) +
  xlim(pm_mu_lfr - 6*psd_mu_lfr, pm_mu_lfr + 6*psd_mu_lfr) +
  geom_vline(xintercept = mu_lfr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu_{lfr}$"),
       y = TeX("Density"), x = TeX("$\\mu_{lfr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation of log fertility rates
gg_sigma_lfr <- ggplot(data.frame(draw=expgro_draws$`sigma_lfr`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dtnorm, args = list(mean=pm_sigma_lfr, sd=psd_sigma_lfr), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma_lfr - 6*psd_sigma_lfr), pm_sigma_lfr + 6*psd_sigma_lfr) +
  geom_vline(xintercept = sigma_lfr_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma_{lfr}$"),
       y = TeX("Density"), x = TeX("$\\sigma_{lfr}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Plots
grid.arrange(gg_mu_lfr, gg_sigma_lfr,
             gg_mu_lbr, gg_sigma_lbr,
             gg_mu_lar, gg_sigma_lar,
             nrow = 3, ncol = 2,
             layout_matrix = rbind(c(1,2), 
                                   c(3,4),
                                   c(5,6)))
### Abundance
# Figures
pre_abundancy <- expgro_draws |>
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
post_abundancy <- expgro_draws |>
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
gg_Npre <- ggplot(pre_abundancy, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = Npre_obs)) +
  ylim(0, max(pre_abundancy$q97.5)) +
  theme_minimal() +
  labs(title = TeX("Posterior of the pre-mortality population"),
       y = TeX("Total population"), x = "Year") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
gg_Npost <- ggplot(post_abundancy, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = Npost_obs)) +
  ylim(0, max(pre_abundancy$q97.5)) +
  theme_minimal() +
  labs(title = TeX("Posterior of the post-mortality population"),
       y = TeX("Total population"), x = "Year") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_Npre, gg_Npost, nrow = 1, ncol = 2)
### Fertility rates
# Figures
fr <- expgro_draws |>
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
gg_fr <- ggplot(fr, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = fr_obs)) +
  ylim(0, 6) +
  theme_minimal() +
  labs(title = TeX("Fertility rates"),
       y = TeX("Fertility rate"), x = "Year") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_fr, nrow = 1, ncol = 1)
## Bckground mortality rates
br <- expgro_draws |>
  select(starts_with("br")) |>
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
gg_br <- ggplot(br, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = br_obs)) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = TeX("Background mortality rates of juveniles for 2003-2024"),
       y = TeX("Background mortality rate"), x = "Year") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_br, nrow = 1, ncol = 1)
## Accident rates
ar <- expgro_draws |>
  select(starts_with("ar")) |>
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
gg_ar <- ggplot(ar, aes(x = 1:Ty, y = mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), 
              fill = "pink", alpha = 0.5) +
  geom_line(color = "magenta", linewidth = 1) +
  geom_point(aes(x = 1:Ty, y = ar_obs)) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = TeX("Accident rates of juveniles for 2003-2024"),
       y = TeX("Accident rate"), x = "Year") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 1)
  )
grid.arrange(gg_ar, nrow = 1, ncol = 1)

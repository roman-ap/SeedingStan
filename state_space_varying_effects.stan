
data {
  int<lower=1> T;
  vector<lower=0>[T] M;
  real pm_LN0;
  real<lower=0> psd_LN0;
  real pm_mu_lmr;
  real<lower=0> psd_mu_lmr;
  real<lower=0> pm_sigma_lmr;
  real<lower=0> psd_sigma_lmr;
  real pm_mu_lfr;
  real<lower=0> psd_mu_lfr;
  real<lower=0> pm_sigma_lfr;
  real<lower=0> psd_sigma_lfr;
}

parameters {
  real<lower=0> N0_raw;
  vector<lower=0>[T-1] NB_raw;
  vector[T] lmr_raw;
  real mu_lmr;
  real<lower=0> sigma_lmr;
  vector[T] lfr_raw;
  real mu_lfr;
  real<lower=0> sigma_lfr;
}

transformed parameters{
  real N0;
  vector[T-1] NB;
  vector[T] Npre;
  vector[T] Npost;
  vector[T] lmr;
  vector[T] mr;
  vector[T] lfr;
  vector[T] fr;
  vector<lower=0>[T-1] meanandvar;
  vector[T-1] mulognormal;
  vector[T-1] sdlognormal;
  lmr = mu_lmr + sigma_lmr*lmr_raw;
  mr = exp(lmr);
  lfr = mu_lfr + sigma_lfr*lfr_raw;
  fr = exp(lfr);
  for (t in 1:T) {
    if(t==1){
      N0 = N0_raw + M[1];
      Npre[t] = N0;
    }
    else{
      meanandvar[t-1] = fr[t-1] * Npost[t-1];
      mulognormal[t-1] = 1.5*log(meanandvar[t-1]) - 0.5*log1p(meanandvar[t-1]);
      sdlognormal[t-1] = sqrt(log1p(meanandvar[t-1])-log(meanandvar[t-1]));
      NB[t-1] = NB_raw[t-1] + M[t];
      Npre[t] = NB[t-1];
    }
    Npost[t] = Npre[t] - M[t];
  }
}

model {
  for (t in 1:T){
      M[t]/Npre[t] ~ beta( (Npre[t]-1)*(1-exp(-mr[t])),
                           (Npre[t]-1)*(exp(-mr[t])) );
  }
  NB ~ lognormal(mulognormal, sdlognormal);
  lmr_raw ~ normal(0,1);
  lfr_raw ~ normal(0,1);
  N0 ~ lognormal(pm_LN0,psd_LN0);
  mu_lmr ~ normal(pm_mu_lmr,psd_mu_lmr);
  sigma_lmr ~ normal(pm_sigma_lmr,psd_sigma_lmr);
  mu_lfr ~ normal(pm_mu_lfr,psd_mu_lfr);
  sigma_lfr ~ normal(pm_sigma_lfr,psd_sigma_lfr);
}

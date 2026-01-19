
functions{
  real cont_binomial_lpdf(real k, real p, real N) {
    return lgamma (N+1) - lgamma (k+1) -lgamma (N-k+1) + k*log(p) + (N-k)*log(1-p);
  }
  real scaled_beta_lpdf(real k, real alpha, real beta, real N) {
    return (alpha-1)*log(k) + (beta-1)*log(N-k) - (alpha+beta-1)*log(N) - lbeta(alpha,beta);
  }
}

data {
  int<lower=1> T;
  vector<lower=0>[T] A;
  real pm_LN0;
  real<lower=0> psd_LN0;
  real pm_mu_lar;
  real<lower=0> psd_mu_lar;
  real<lower=0> pm_sigma_lar;
  real<lower=0> psd_sigma_lar;
  real pm_mu_lbr;
  real<lower=0> psd_mu_lbr;
  real<lower=0> pm_sigma_lbr;
  real<lower=0> psd_sigma_lbr;
  real pm_mu_lfr;
  real<lower=0> psd_mu_lfr;
  real<lower=0> pm_sigma_lfr;
  real<lower=0> psd_sigma_lfr;
}

parameters {
  real<lower=0> N0_raw;
  vector<lower=0>[T-1] NB_raw;
  vector<lower=0, upper=1>[T] U;
  vector[T] lar_raw;
  real mu_lar;
  real<lower=0> sigma_lar;
  vector[T] lbr_raw;
  real mu_lbr;
  real<lower=0> sigma_lbr;
  vector[T] lfr_raw;
  real mu_lfr;
  real<lower=0> sigma_lfr;
}

transformed parameters{
  real N0;
  vector[T-1] NB;
  vector[T] M;
  vector[T] Npre;
  vector[T] Npost;
  vector[T] lar;
  vector[T] ar;
  vector[T] lbr;
  vector[T] br;
  vector[T] lfr;
  vector[T] fr;
  vector<lower=0>[T-1] meanandvar;
  vector[T-1] mulognormal;
  vector[T-1] sdlognormal;
  lar = mu_lar + sigma_lar*lar_raw;
  ar = exp(lar);
  lbr = mu_lbr + sigma_lbr*lbr_raw;
  br = exp(lbr);
  lfr = mu_lfr + sigma_lfr*lfr_raw;
  fr = exp(lfr);
  for (t in 1:T) {
    if(t==1){
      N0 = N0_raw + A[1];
      Npre[t] = N0;
    }
    else{
      meanandvar[t-1] = fr[t-1] * Npost[t-1];
      mulognormal[t-1] = 1.5*log(meanandvar[t-1]) - 0.5*log1p(meanandvar[t-1]);
      sdlognormal[t-1] = sqrt(log1p(meanandvar[t-1])-log(meanandvar[t-1]));
      NB[t-1] = NB_raw[t-1] + A[t];
      Npre[t] = NB[t-1];
    }
    M[t] = A[t] + U[t]*(Npre[t] - A[t]);
    Npost[t] = Npre[t] - M[t];
  }
}

model {
  for (t in 1:T){
      M[t]/Npre[t] ~ beta( (Npre[t]-1)*(1-exp(-(ar[t]+br[t]))),
                           (Npre[t]-1)*(exp(-(ar[t]+br[t]))) );
      A[t]/M[t] ~ beta( (M[t]-1)*(ar[t]/(ar[t]+br[t])),
                        (M[t]-1)*(br[t]/(ar[t]+br[t])) );
  }
  NB ~ lognormal(mulognormal, sdlognormal);
  lar_raw ~ normal(0,1);
  lbr_raw ~ normal(0,1);
  lfr_raw ~ normal(0,1);
  N0 ~ lognormal(pm_LN0,psd_LN0);
  mu_lar ~ normal(pm_mu_lar,psd_mu_lar);
  sigma_lar ~ normal(pm_sigma_lar,psd_sigma_lar);
  mu_lbr ~ normal(pm_mu_lbr,psd_mu_lbr);
  sigma_lbr ~ normal(pm_sigma_lbr,psd_sigma_lbr);
  mu_lfr ~ normal(pm_mu_lfr,psd_mu_lfr);
  sigma_lfr ~ normal(pm_sigma_lfr,psd_sigma_lfr);
}

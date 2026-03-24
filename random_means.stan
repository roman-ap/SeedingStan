
data {
  int<lower=1> N;
  vector[N] y;
  real pm_mu;
  real<lower=0> psd_mu;
  real<lower=0> pm_sigma;
  real<lower=0> psd_sigma;
  real<lower=0> sd_y;
}

parameters{
  real mu;
  real<lower=0> sigma;
  vector[N] x;
}

model {
    target += normal_lpdf(y | x, sd_y);
    target += normal_lpdf(x | mu, sigma);
    target += normal_lpdf(mu | pm_mu, psd_mu);
    target += normal_lpdf(sigma | pm_sigma, psd_sigma);
}

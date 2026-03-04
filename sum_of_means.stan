
data {
  int<lower=1> N;
  vector[N] y;
  vector[2] pm_mus;
  cov_matrix[2] pcov_mus;
  real<lower=0> psd_mu;
}

parameters {
  vector[2] mus;
}

transformed parameters{
  real mu;
  mu = sum(mus);
}

model {
  target += normal_lpdf(y | mu, psd_mu);
  //target += multi_normal_lpdf(mus | pm_mus, pcov_mus);
}

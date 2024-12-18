data {
  int<lower=1> obs_n;
  int<lower=1> syn_n;
  int<lower=1> col_n;
  matrix[obs_n, col_n] obs_x;
  matrix[syn_n, col_n] syn_x;
  int<lower=0,upper=1> obs_y[obs_n];
  int<lower=0,upper=1> syn_y[syn_n];
  real kappa;
  real tau_alpha;
  real tau_gamma;
}

parameters {
  vector[col_n+1] coefs;
  real<lower=0> theta;
}

model {
  vector[obs_n] obs_mu = coefs[1] + obs_x * coefs[2:(col_n+1)];
  vector[syn_n] syn_mu = coefs[1] + syn_x * coefs[2:(col_n+1)];
  real syn_weight = 1/(theta*syn_n);

  // Likelihood
  for (i in 1:obs_n)
    target += bernoulli_logit_lpmf(obs_y[i] | obs_mu[i]);

  // Catalytic prior
  for (i in 1:syn_n)
      target += syn_weight * bernoulli_logit_lpmf(syn_y[i] | syn_mu[i]);

  // Inverse-tau_gamma function
  target += (col_n+tau_alpha+1)*log(1/theta)-(1/theta)*(kappa+1/tau_gamma);
}


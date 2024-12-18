data {
  int<lower=1> obs_n;
  int<lower=1> syn_n;
  int<lower=1> col_n;
  matrix[obs_n, col_n] obs_x;
  matrix[syn_n, col_n] syn_x;
  vector[obs_n] obs_y;
  vector[syn_n] syn_y;
  real<lower=0> tau;
  real<lower=0> variance_alpha;
  real<lower=0> variance_beta;
}

parameters {
  vector[col_n+1] coefs;
  real<lower=0> sigma;
}

model {
  vector[obs_n] obs_mu = coefs[1] + obs_x * coefs[2:(col_n+1)];
  vector[syn_n] syn_mu = coefs[1] + syn_x * coefs[2:(col_n+1)];
  real syn_weight = tau/syn_n;

  // Variance Prior
  //pow(sigma, 2) ~ inv_gamma(variance_alpha, variance_beta);

  // Likelihood
  for (i in 1:obs_n)
    target += normal_lpdf(obs_y[i] | obs_mu[i], sigma);

  // Catalytic prior
  for (i in 1:syn_n)
      target += syn_weight * normal_lpdf(syn_y[i] | syn_mu[i], sigma);

  // Variance Prior
  target += inv_gamma_lpdf(pow(sigma, 2) | variance_alpha, variance_beta);
}

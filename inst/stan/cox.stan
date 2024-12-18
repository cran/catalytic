functions {
  /**
   * Calculates the log prior for the coefficients (`coefs`) in a Cox model with a customized hazard constant.
   *
   * This function computes the log prior using the synthetic data and a specified hazard constant.
   * The contribution from each synthetic observation is calculated, and then scaled by a downweighting factor.
   *
   * @param syn_size Number of synthetic samples.
   * @param coefs Vector of model coefficients.
   * @param syn_x Matrix of synthetic covariates.
   * @param syn_time Vector of synthetic survival times.
   * @param syn_status Vector of synthetic status indicators.
   * @param hazard_constant Customized hazard constant values evaluated at each synthetic survival time.
   * @param tau Positive number that downweighting the synthetic dataset in the estimation process.
   *
   * @return A scalar representing the log prior for the given `coefs` using the synthetic data.
   */
  real log_of_cat_prior(int syn_size,
                        matrix syn_x,
                        vector syn_time,
                        vector syn_status,
                        vector coefs,
                        real tau,
                        real hazard_constant
                        ) {

    vector[syn_size] syn_lp = syn_x * coefs;
    vector[syn_size] exp_syn_lp = exp(syn_lp);

    real syn_prior = 0;
    for (i in 1:syn_size) {
      syn_prior = syn_prior + (syn_status[i] * syn_lp[i] - exp_syn_lp[i] * syn_time[i] * hazard_constant);
    }
    return syn_prior * tau / syn_size;
  }
}


data {
  int<lower=0> obs_size;
  int<lower=0> syn_size;
  int<lower=0> col_n;
  matrix[obs_size, col_n] obs_x;
  matrix[syn_size, col_n] syn_x;
  vector[syn_size] syn_time;
  vector[syn_size] syn_status;
  real<lower=0> tau;

  int<lower=0> time_interval_num;  // Number of time intervals.
  matrix[obs_size, time_interval_num] adj_risk_set;  // Matrix indicating risk set minus failure set for each observation across intervals.
  matrix[obs_size, time_interval_num] failure_set;  // Matrix indicating which intervals an observation has an event.

  real<lower=0> hazard_constant;
  vector[time_interval_num] hazard_alpha_list;  // Shape parameters for the gamma prior distribution of the baseline hazard.
  real hazard_beta;  // Rate parameter for the gamma prior distribution of the baseline hazard.
}


parameters {
  vector[col_n] coefs;
  vector<lower=0>[time_interval_num] h_j;  // parameter with a gamma prior
}

model {
  /**
   * This model block represents a Bayesian Cox Proportional Hazards model with the inclusion of a customized
   * hazard constant and a log catalytic prior for the coefficients based on synthetic data.
   *
   * The likelihood component of the model is constructed using observed data, while the prior is constructed using
   * both observed and synthetic data.
   *
   * The total target log posterior (unnormalized) is updated with contributions from both the likelihood and the prior.
   */
  vector[time_interval_num] adj_risk_set_sum;  // Vector to store the summation terms for the risk set minus event set.
  vector[time_interval_num] event_set_sum;  // Vector to store the summation terms for the event set.
  matrix[obs_size, time_interval_num] exp_obs_lp = rep_matrix(exp(obs_x * coefs), time_interval_num);  // Matrix where each column is the exponential of obs_x multiplied by coefs.
  matrix[obs_size, time_interval_num] h_seq = rep_matrix(h_j', obs_size);  // Replicating the hazard sequence across `obs_size` rows.
  matrix[obs_size, time_interval_num] h_seq_exp_obs_lp = -h_seq .* exp_obs_lp;  // Matrix storing product of hazard sequence and the exponential transformation of obs_x and coefs.

  for (j in 1:time_interval_num) {
    // Summing over the adjusted risk set for the `j-th` interval.
    adj_risk_set_sum[j] = sum(exp_obs_lp[, j] .* adj_risk_set[, j]);
    // Summing over the event set for the `j-th` interval using the log1m_exp transformation.
    event_set_sum[j] = sum(log1m_exp(h_seq_exp_obs_lp[, j]) .* failure_set[, j]);
  }

  // Update the target log posterior with the likelihood component.
  target += sum(-h_j .* adj_risk_set_sum + event_set_sum);

  // Update the target log posterior with the log of catalytic prior component.
  target += log_of_cat_prior(syn_size,
                             syn_x,
                             syn_time,
                             syn_status,
                             coefs,
                             tau,
                             hazard_constant);

  // Prior distribution for the hazard sequence based on a Gamma distribution.
  h_j ~ gamma(hazard_alpha_list, hazard_beta);

}


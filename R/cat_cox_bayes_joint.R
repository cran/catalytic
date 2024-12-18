#' Bayesian Estimation for Catalytic Cox Proportional-Hazards Model (COX) with adaptive tau
#'
#' This function performs Bayesian estimation for a catalytic Cox proportional-hazards model (COX) using RStan
#' by using adaptive tau. It allows users to estimate the coefficients and cumulative baseline hazard increments
#' over specified time intervals with Bayesian sampling.
#'
#' @param formula A formula specifying the Cox model. Should at least include response variables (e.g. \code{~ .}).
#' @param cat_init A list generated from `cat_cox_initialization`.
#' @param hazard_beta Numeric, default `2`. Shape parameter for the Gamma distribution in the hazard model.
#' @param tau_alpha Numeric, defaults `2`. Scalar for the shape parameter of the Gamma-like function for tau.
#' @param tau_gamma Numeric, defaults `1`. Scalar for the scale parameter of the Gamma-like function for tau.
#' @param chains Integer, default `4`. Number of Markov chains to be run in the RStan sampling.
#' @param iter Integer, default `2000`. Number of iterations per chain in the RStan sampling.
#' @param warmup Integer, default `1000`. Number of warm-up (or burn-in) iterations for each chain.
#' @param ... Additional arguments passed to RStan’s `rstan::sampling` function.
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{stan_data}{A data list used for fitting RStan sampling model.}
#' \item{stan_model}{Compiled RStan model object for Cox regression.}
#' \item{stan_sample_model}{Fitted RStan sampling model containing posterior samples.}
#' \item{tau}{Mean posterior estimates of tau value from `stan_sample_model`.}
#' \item{coefficients}{Mean posterior estimates of model coefficients from `stan_sample_model`.}
#' \item{increment_cumulative_baseline_hazard}{Mean posterior estimates of Estimated
#' cumulative hazard increments over time intervals from `stan_sample_model`.}
#'
#' @examples
#' \donttest{
#' library(survival)
#' data("cancer")
#' cancer$status[cancer$status == 1] <- 0
#' cancer$status[cancer$status == 2] <- 1
#'
#' cat_init <- cat_cox_initialization(
#'   formula = Surv(time, status) ~ 1, # formula for simple model
#'   data = cancer,
#'   syn_size = 100, # Synthetic data size
#'   hazard_constant = 0.1, # Hazard rate value
#'   entry_points = rep(0, nrow(cancer)), # Entry points of each observation
#'   x_degree = rep(1, ncol(cancer) - 2), # Degrees for polynomial expansion of predictors
#'   resample_only = FALSE, # Whether to perform resampling only
#'   na_replace = stats::na.omit # How to handle NA values in data
#' )
#'
#' cat_model <- cat_cox_bayes_joint(
#'   formula = ~., # Should at least include response variables
#'   cat_init = cat_init, # Only accept object generated from `cat_cox_initialization`
#'   hazard_beta = 2, # Shape parameter for the Gamma distribution in the hazard model
#'   tau_alpha = 2, # Shape parameter of the Gamma-like function for tau
#'   tau_gamma = 1, # Scale parameter of the Gamma-like function for tau
#'   chains = 1, # Number of Markov chains to be run in the RStan sampling
#'   iter = 10, # Number of iterations per chain in the RStan sampling
#'   warmup = 5 # Number of warm-up (or burn-in) iterations for each chain
#' )
#' cat_model
#' }
#' @export
cat_cox_bayes_joint <- function(
    formula,
    cat_init,
    hazard_beta = 2,
    tau_alpha = 2,
    tau_gamma = 1,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    ...) {
  # Validate Input Parameters
  validate_cox_input(
    formula = formula,
    cat_init = cat_init,
    hazard_beta = hazard_beta,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  # Get unique times from sorted observed status times and a value larger than the maximum observed time
  unique_observed_status_times_and_beyond <- unique(c(
    # Sort the times where event occurred, status = 1
    sort(cat_init$obs_time[cat_init$obs_status == 1]),
    # Create a value larger than the maximum observed time
    2 * max(cat_init$obs_time) -
      max(cat_init$obs_time[-which(cat_init$obs_time == max(cat_init$obs_time))])
  ))

  # get risk set and failure set
  risk_and_failure_sets <- get_cox_risk_and_failure_sets(
    time_vector = cat_init$obs_time,
    status_vector = cat_init$obs_status,
    indicator_vector = unique_observed_status_times_and_beyond
  )

  # Get the initial guess of H_0 from the intercept-only weibull regression model based on observation data
  weibull_model <- do.call(
    survival::survreg,
    list(
      formula = survival::Surv(cat_init$obs_time, cat_init$obs_status) ~ 1,
      dist = "weibull"
    )
  )
  # Extracts the shape and scale parameters of the Weibull distribution
  weibull_shape <- as.numeric(1 / weibull_model$scale)
  weibull_scale <- as.numeric(exp(stats::coef(weibull_model)))

  # Calculates the cumulative hazard values for the risk_and_failure_sets
  H_star_list <- (unique_observed_status_times_and_beyond / weibull_scale)^(weibull_shape)

  # List of tau_alpha value from each h_j ~ Gamma(alpha_j, hazard_beta), where alpha_j = alpha_0[j] - alpha_0[j-1]
  hazard_alpha_list <- diff(c(0, hazard_beta * H_star_list))

  stan_data <- list(
    obs_size = cat_init$obs_size,
    syn_size = cat_init$syn_size,
    col_n = ncol(cat_init$adj_x),
    obs_x = cat_init$adj_obs_x,
    syn_x = cat_init$adj_syn_x,
    syn_time = cat_init$syn_time,
    syn_status = cat_init$syn_status,
    time_interval_num = length(unique_observed_status_times_and_beyond),
    adj_risk_set = risk_and_failure_sets$risk_set - risk_and_failure_sets$failure_set,
    failure_set = risk_and_failure_sets$failure_set,
    hazard_constant = cat_init$hazard_constant,
    hazard_alpha_list = hazard_alpha_list,
    hazard_beta = hazard_beta,
    kappa = get_cox_kappa(
      X = cat_init$syn_x,
      time = cat_init$syn_time,
      status = cat_init$syn_status,
      hazard_constant = cat_init$hazard_constant
    ),
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma
  )

  stan_model <- get_stan_model(
    type = "cox",
    joint_tau = TRUE
  )

  # Perform Bayesian Sampling via RStan
  stan_sample_model <- rstan::sampling(
    object = stan_model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    ...
  )

  ## Extract the mean of coefs from the stan sampling model
  coefs <- rstan::summary(stan_sample_model)$summary[
    grep("coefs", rownames(rstan::summary(stan_sample_model)$summary)),
    "mean"
  ]

  if (!is.null(coefs)) {
    names(coefs) <- colnames(cat_init$adj_x)
  }

  cat_model <- list(
    function_name = "cat_cox_bayes_joint",
    ## Input/Processed parameters
    formula = stats::as.formula(paste0(
      "survival::Surv(",
      cat_init$time_col_name,
      ",",
      cat_init$status_col_name,
      ")",
      f_rhs
    )),
    cat_init = cat_init,
    hazard_beta = hazard_beta,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma,
    chains = chains,
    iter = iter,
    warmup = warmup,
    ...,
    ## Result
    stan_data = stan_data,
    stan_model = stan_model,
    stan_sample_model = stan_sample_model,
    tau = rstan::summary(stan_sample_model)$summary[
      grep("tau", rownames(rstan::summary(stan_sample_model)$summary)),
      "mean"
    ],
    coefficients = coefs,
    increment_cumulative_baseline_hazard = rstan::summary(stan_sample_model)$summary[
      grep("h_j", rownames(rstan::summary(stan_sample_model)$summary)),
      "mean"
    ]
  )

  class(cat_model) <- c(cat_model$class, "cat_bayes", "cat_cox")

  return(cat_model)
}

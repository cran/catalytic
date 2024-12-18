# ----- GLMs -----
#' Generate Sample Data for GLM
#'
#' This function generates sample data for a specified GLM family. It can generate
#' binomial or Gaussian distributed data based on the provided parameters.
#'
#' @param family_string Character. The family of the GLM. Options are `"binomial"` or `"gaussian"`.
#' @param n Integer. The number of samples to generate.
#' @param mean Numeric. The mean of the distribution (used for both binomial and Gaussian).
#' @param sd Numeric. The standard deviation of the distribution (used only for Gaussian).
#'
#' @return Numeric vector of generated sample data.
get_glm_sample_data <- function(family_string,
                                n = 10,
                                mean = 0,
                                sd = 1) {
  switch(family_string,
    "binomial" = stats::rbinom(n, 1, mean),
    "gaussian" = stats::rnorm(n, mean, sd)
  )
}

#' Run Hamiltonian Monte Carlo to Get MCMC Sample Result
#'
#' This function uses Hamiltonian Monte Carlo (HMC) to generate samples for Markov Chain Monte Carlo (MCMC)
#' sampling from a target distribution specified by `neg_log_den_func`. Each iteration performs a full HMC update
#' to generate a new sample position.
#'
#' @param neg_log_den_func A function that computes the negative log-density of the target distribution at a given position.
#' @param neg_log_den_grad_func A function that computes the gradient of `neg_log_den_func` at a given position.
#' @param coefs_0 A numeric vector specifying the initial position of the chain in the parameter space.
#' @param iter An integer specifying the number of HMC sampling iterations. Defaults to 5.
#' @param hmc_scale A numeric value representing the scale factor for the leapfrog step size in the HMC update. Defaults to 0.01.
#' @param hmc_steps An integer specifying the number of leapfrog steps in each HMC update. Defaults to 5.
#'
#' @return A numeric vector representing the final position in the parameter space after the specified number of iterations.
get_hmc_mcmc_result <- function(
    neg_log_den_func,
    neg_log_den_grad_func,
    coefs_0,
    iter = 5,
    hmc_scale = 0.01,
    hmc_steps = 5) {
  coefs <- coefs_0

  for (i in 1:iter) {
    coefs <- hmc_neal_2010(
      neg_log_den_func = neg_log_den_func,
      neg_log_den_grad_func = neg_log_den_grad_func,
      leapfrog_stepsize = stats::rbeta(1, 0.5, 0.5) * hmc_scale / hmc_steps, # Add random value to leapfrog_stepsize to overcome poor performance
      leapfrog_step = hmc_steps,
      current_pos = coefs
    )$position
  }

  return(coefs)
}

#' Hamiltonian Monte Carlo (HMC) Implementation
#'
#' This function implements the Hamiltonian Monte Carlo algorithm as
#' described by Radford M. Neal (2010) in "MCMC using Hamiltonian dynamics",
#' which is a part of the Handbook of Markov Chain Monte Carlo. The method uses
#' Hamiltonian dynamics to propose new positions and then applies the Metropolis
#' criterion to decide whether to accept or reject the new position.
#'
#' @param neg_log_den_func A function that evaluates the negative log of the density
#'                          (potential energy) of the distribution to be sampled,
#'                          including any constants.
#' @param neg_log_den_grad_func A function that computes the gradient of
#'                              \code{neg_log_den_func}.
#' @param leapfrog_stepsize A numeric value specifying the step size for the
#'                          leapfrog integration method.
#' @param leapfrog_step A numeric value specifying the number of leapfrog steps
#'                      to take to propose a new state.
#' @param current_pos A numeric vector representing the current position (state)
#'                    of the system.
#'
#' @return A list containing the following elements:
#'  - `position`: The position of the system after the leapfrog steps, which is
#'                 the proposed new position if accepted, or the current position
#'                 if rejected.
#'  - `potential_energy`: The potential energy of the proposed position.
#'  - `accepted`: A logical value indicating whether the proposal was accepted
#'                 (TRUE) or rejected (FALSE).
#'
#' @details This function was written for illustrative purposes. More elaborate
#  implementations of this basic HMC method and various variants of HMC are available
#' on Radford M. Neal's personal webpage (http://www.cs.utoronto.ca/~radford/).
#'
#' @references Neal, R. M. (2012). MCMC using Hamiltonian dynamics.
#' arXiv:1206.1901. Available at: \url{https://arxiv.org/pdf/1206.1901}
hmc_neal_2010 <- function(
    neg_log_den_func,
    neg_log_den_grad_func,
    leapfrog_stepsize,
    leapfrog_step,
    current_pos) {
  pos <- current_pos

  # Sample initial momentum from standard normal distribution
  momentum <- stats::rnorm(length(pos), 0, 1)
  current_momentum <- momentum

  # Initial half-step for momentum
  momentum <- momentum - leapfrog_stepsize * neg_log_den_grad_func(pos) / 2

  # Full steps for position and momentum over the leapfrog trajectory (except at end of trajectory)
  for (i in 1:(leapfrog_step - 1)) {
    # Full step for position
    pos <- pos + leapfrog_stepsize * momentum

    # Full step for momentum
    momentum <- momentum - leapfrog_stepsize * neg_log_den_grad_func(pos)
  }

  # Final position update
  pos <- pos + leapfrog_stepsize * momentum

  # Final half-step for momentum and reverse momentum to maintain proposal symmetry
  momentum <- -(momentum - leapfrog_stepsize * neg_log_den_grad_func(pos) / 2)

  # Compute the potential and kinetic energy at the start and end of the trajectory
  current_potential_energy <- neg_log_den_func(current_pos)
  current_kinetic_energy <- sum(current_momentum^2) / 2
  proposed_potential_energy <- neg_log_den_func(pos)
  proposed_kinetic_energy <- sum(momentum^2) / 2

  # Metropolis acceptance criterion
  if (stats::runif(1) < exp(current_potential_energy - current_potential_energy +
    current_kinetic_energy - proposed_kinetic_energy)) {
    # Accept the state at end of trajectory, returning the position at the end of the trajectory
    return(list(
      position = pos,
      potential_energy = proposed_potential_energy,
      accepted = TRUE
    ))
  } else {
    # Reject the state at end of trajectory, returning the current position
    return(list(
      position = current_pos,
      potential_energy = current_potential_energy,
      accepted = FALSE
    ))
  }
}

#' Compute Gradient of Log Density for GLM Families
#'
#' This function calculates the gradient of the log density with respect to the coefficients
#' for a given GLM family based on the provided predictors, response variable, and weights.
#'
#' @param family_string Character. The GLM family to use. Options are `"binomial"` or `"gaussian"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param Y Vector. The response variable.
#' @param coefs Numeric vector. The coefficients for the GLM.
#' @param weights Numeric vector. The weights for the GLM. Default is 1.
#'
#' @return Numeric vector. The gradient of the log density with respect to the coefficients
get_glm_log_density_grad <- function(family_string,
                                     X,
                                     Y,
                                     coefs,
                                     weights = 1) {
  # Calculate the gradient of log density
  gradient <- weights * (Y - get_glm_mean(family_string, X, coefs))

  return(c(t(gradient) %*% as.matrix(cbind(1, X))))
}

#' Compute Diagonal Approximate Covariance Matrix
#'
#' This function computes the diagonal elements of the approximate covariance matrix
#' for the coefficients in a generalized linear model (GLM). The covariance is derived
#' from the second derivative (Hessian) of the log-likelihood function.
#'
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param model A fitted GLM model object. The object should contain the fitted values
#'        and prior weights necessary for computing the Hessian.
#'
#' @return Numeric vector. The diagonal elements of the approximate covariance matrix.
get_glm_diag_approx_cov <- function(X,
                                    model) {
  # Nested function to compute the second derivative (Hessian) matrix
  get_second_derivative <- function(X,
                                    fitted_values,
                                    prior_weights) {
    return((t(cbind(1, X))
    ) %*% diag(fitted_values * (1 - fitted_values) * prior_weights) %*% cbind(1, X))
  }

  # Calculate the diagonal elements of the covariance approximation
  hessian <- get_second_derivative(
    X = as.matrix(X),
    fitted_values = model$fitted.values,
    prior_weights = model$prior.weights
  )

  tryCatch(
    return(diag(solve(hessian))),
    error = function(.) {
      # Regularization for dealing with computationally singular issue
      return(diag(solve(hessian + diag(ncol(hessian)) * 1e-5)))
    }
  )
}



#' Generate Suggestions for Bayesian Joint Binomial GLM Parameter Estimation
#'
#' This function provides suggestions for improving the parameter estimation process
#' in Bayesian joint Binomial GLM modeling based on the diagnostic output from a Stan model.
#' It evaluates the results and suggests adjustments to improve model fit and convergence.
#'
#' @param alpha Numeric. The alpha parameter used in the prior distribution.
#' @param stan_iter Integer. The number of iterations used in Stan sampling.
#' @param stan_sample_model Stan model object containing sampling results.
#' @param binomial_joint_theta Logical. Whether to use theta in the binomial model.
#' @param binomial_joint_alpha Logical. Whether to use joint alpha in the binomial model.
#' @param binomial_tau_lower Numeric. The lower bound for tau in the binomial model.
#'
#' @return NULL. The function prints suggestions to the console based on the model diagnostics.
print_glm_bayes_joint_binomial_suggestion <- function(
    alpha,
    stan_iter,
    stan_sample_model,
    binomial_joint_theta,
    binomial_joint_alpha,
    binomial_tau_lower) {
  # Helper function to calculate the error between estimated and theoretical percentiles
  get_est_params <- function(params, percentiles, mean_, distribution) {
    alpha <- exp(params[1]) # Convert log-alpha to alpha
    beta <- exp(params[2]) # Convert log-beta to beta

    # Calculate percentiles based on the distribution
    calc_percentiles <- switch(distribution,
      "invgamma" = invgamma::qinvgamma(c(0.25, 0.5, 0.75), shape = alpha, scale = beta),
      "gamma" = stats::qgamma(c(0.25, 0.5, 0.75), shape = alpha, scale = beta)
    )

    # Compute error between calculated percentiles and those from the Stan model
    error <- sum((percentiles - calc_percentiles)^2) + (mean_ - (beta / (alpha - 1)))^2

    return(error)
  }

  if (binomial_joint_theta) {
    # Retrieve theta summary from Stan model
    theta_summary <- rstan::summary(stan_sample_model)$summary["theta", ]

    # Optimize parameters for the inverse gamma distribution
    est_params <- stats::optim(
      par = c(1, 1),
      fn = get_est_params,
      percentiles = theta_summary[c("25%", "50%", "75%")],
      mean_ = theta_summary["mean"],
      distribution = "invgamma"
    )

    # Check if estimated parameters fall outside the 97.5% quantile range
    if (invgamma::qinvgamma(c(0.975),
      shape = exp(est_params$par[1]),
      scale = exp(est_params$par[2])
    ) < theta_summary["97.5%"]) {
      # Provide suggestions based on model diagnostics
      if (any(rstan::summary(stan_sample_model)$summary[, "Rhat"] > 1.2)) {
        cat("\nSUGGESTION: Please consider using more iterations for Stan (e.g.,", stan_iter + 1000, ")\n")
      } else if (alpha < 5 & !binomial_joint_alpha) {
        cat("\nSUGGESTION: Please consider using a larger `alpha` (e.g.,", min(alpha + 0.1, 5), ")\n")
      } else if (alpha >= 5 & !binomial_joint_alpha) {
        cat("\nSUGGESTION: Please consider setting `binomial_joint_alpha = TRUE`.\n")
      }
    }
  } else {
    # Retrieve tau summary from Stan model
    tau_summary <- rstan::summary(stan_sample_model)$summary["tau", ]

    # Optimize parameters for the gamma distribution
    est_params <- stats::optim(
      par = c(1, 1),
      fn = get_est_params,
      percentiles = tau_summary[c("25%", "50%", "75%")],
      mean_ = tau_summary["mean"],
      distribution = "gamma"
    )

    # Check if estimated parameters fall outside the 2.5% quantile range
    if (stats::qgamma(c(0.025),
      shape = exp(est_params$par[1]),
      scale = exp(est_params$par[2])
    ) > tau_summary["2.5%"]) {
      # Provide suggestions based on model diagnostics
      if (any(rstan::summary(stan_sample_model)$summary[, "Rhat"] > 1.2)) {
        cat("\nSUGGESTION: Please consider using more iterations for Stan (e.g.,", stan_iter + 1000, ")\n")
      } else if (binomial_tau_lower < 1) {
        cat("\nSUGGESTION: Please consider using a larger `binomial_tau_lower` (e.g.,", min(binomial_tau_lower + 0.05, 1), ")\n")
      } else {
        cat("\nSUGGESTION: Please consider setting `binomial_joint_theta = TRUE`.\n")
      }
    }
  }
}

#' Get Custom Variance for Generalized Linear Model (GLM)
#'
#' This function calculates a custom variance for a Generalized Linear Model (GLM) based on the specified formula,
#' the model initialization object, and a scaling factor `tau`. The custom variance is computed by adjusting
#' the residuals of the fitted model and returning a weighted sum of squared residuals.
#'
#' @param formula A formula object specifying the GLM model to be fitted, such as `response ~ predictors`.
#' @param cat_init A list object containing the initialization data for the model. Generated from `cat_initialization`
#' @param tau A numeric value representing a scaling factor for down-weighting synthetic data
#'
#' @return A numeric value representing the custom variance for the GLM model.
get_glm_custom_var <- function(
    formula,
    cat_init,
    tau) {
  ncol_model <- do.call(
    stats::glm,
    list(
      formula = formula,
      family = cat_init$family,
      data = cat_init$data,
      weights = c(
        rep(1, cat_init$obs_size),
        rep(tau / cat_init$syn_size, cat_init$syn_size)
      )
    )
  )


  return(
    sum(stats::residuals(ncol_model)^2) / stats::df.residual(ncol_model) +
      sum(stats::residuals(ncol_model)[1:cat_init$obs_size]^2) / cat_init$obs_size
  )
}

#' Compute Lambda Based on Discrepancy Method
#'
#' This function calculates a lambda value based on the selected discrepancy method
#' for a generalized linear model (GLM). The discrepancy method determines the
#' type of error or deviance used in the calculation.
#'
#' @param discrepancy_method Character. A string specifying the type of discrepancy method
#'        to use. Options are `"mean_square_error"`, `"mean_classification_error"`, or `"logistic_deviance"`.
#'        Default is `"mean_square_error"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param coefs Numeric vector. The coefficients for the GLM.
#'
#' @return Numeric. The computed lambda value based on the selected discrepancy method.
get_glm_lambda <- function(
    discrepancy_method = c(
      "mean_square_error",
      "mean_classification_error",
      "logistic_deviance"
    ),
    X,
    coefs) {
  # Validate and match the discrepancy method
  discrepancy_method <- match.arg(discrepancy_method)

  # Compute lambda based on the selected discrepancy method
  return(switch(discrepancy_method,
    "mean_square_error" = 2 * get_glm_mean(family_string = "gaussian", X = X, coefs = coefs),
    "mean_classification_error" = 2 * get_glm_mean(family_string = "binomial", X = X, coefs = coefs),
    "logistic_deviance" = get_linear_predictor(X = X, coefs = coefs)
  ))
}

#' Compute Log Density Based on GLM Family
#'
#' This function calculates the log density of the response variable given a generalized linear model (GLM)
#' based on the specified family. The log density is computed differently for binomial and gaussian families.
#'
#' @param family_string Character. The GLM family to use. Options are `"binomial"` or `"gaussian"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param Y Vector or data frame. The response variable for the GLM. If a data frame, it is converted to a numeric vector.
#' @param coefs Numeric vector. The coefficients for the GLM.
#' @param weights Numeric vector. Weights for the observations. Default is `1` (no weighting).
#'
#' @return Numeric. The computed log density of the response variable based on the specified family.
get_glm_log_density <- function(family_string,
                                X,
                                Y,
                                coefs,
                                weights = 1) {
  # Compute the linear predictor
  lp <- get_linear_predictor(X = X, coefs = coefs)

  return((switch(family_string,
    "binomial" = sum(weights * (Y * lp - log(1 + exp(lp)))),
    "gaussian" = -0.5 * sum(weights * (log(2 * pi * max(mean((Y - lp)^2), 1e-5)) + 1))
  )))
}


#' Compute Mean Based on GLM Family
#'
#' This function calculates the mean of the response variable for a generalized linear model (GLM)
#' based on the specified family. The calculation depends on whether the family is binomial or gaussian.
#'
#' @param family_string Character. The GLM family to use. Options are `"binomial"` or `"gaussian"`.
#' @param X Matrix. The design matrix (predictors) for the GLM.
#' @param coefs Numeric vector. The coefficients for the GLM.
#'
#' @return Numeric vector. The computed mean of the response variable based on the specified family.
get_glm_mean <- function(family_string, X, coefs) {
  # Compute the linear predictor
  lp <- get_linear_predictor(X = X, coefs = coefs)

  # Return the mean based on the specified family
  return(switch(family_string,
    "binomial" = 1 / (1 + exp(-lp)),
    "gaussian" = lp
  ))
}

#' Retrieve GLM Family Name or Name with Link Function
#'
#' This function retrieves the name of a GLM family or, optionally, the family name with the associated link function.
#'
#' @param family Character or function. The name of the GLM family (as a string) or a function that returns a GLM family object.
#' @param with_link Logical. If TRUE, returns the family name along with the link function in the format "family \\[link\\]".
#' If FALSE, only the family name is returned. Default is FALSE.
#'
#' @return A character string. The name of the GLM family, or the name with the link function if `with_link` is TRUE.
get_glm_family_string <- function(family, with_link = FALSE) {
  # Validate and retrieve the family object
  validate_glm_family <- function(f) {
    if (methods::is(f, "character")) {
      tryCatch(
        {
          f <- get(f, mode = "function", envir = parent.frame(2))
        },
        error = function(e) {
          stop("'family' must be a valid GLM family.",
            call. = FALSE
          )
        }
      )
    }
    if (methods::is(f, "function")) {
      f <- f()
    }
    if (!methods::is(f, "family")) {
      stop("'family' must be a valid GLM family.",
        call. = FALSE
      )
    }
    return(f)
  }

  # Retrieve and validate the family object
  family <- validate_glm_family(family)

  # Return the family name or family name with link function
  if (with_link) {
    return(paste0(family$family, " [", family$link, "]"))
  } else {
    return(family$family)
  }
}

#' Validate Inputs for Catalytic Generalized Linear Models (GLMs)
#'
#' This function validates the input parameters for initializing a catalytic Generalized
#' Linear Models (GLMs). It ensures that the provided model formula, family, and
#' additional parameters are suitable for further analysis. The function performs
#' various checks on the input values to confirm they meet expected criteria.
#'
#' @param formula A formula object specifying the GLM to be fitted. The left-hand
#'        side of the formula should at least contains the response variable.
#' @param cat_init An object of class `cat_initialization` generated by
#'        `cat_glm_initialization`. It contains model initialization details,
#'        such as the response variable name and the GLM family.
#' @param tau A positive numeric value for the tau parameter in the model.
#'        It represents a regularization or scaling factor and must be greater
#'        than zero.
#' @param tau_seq A numeric vector specifying a sequence of tau values.
#'        This is used for parameter tuning and must contain positive values.
#' @param tau_0 A positive numeric value for the initial tau parameter,
#'        which must be greater than zero.
#' @param parametric_bootstrap_iteration_times An integer specifying the
#'        number of iterations for the parametric bootstrap method. It must be
#'        greater than zero.
#' @param cross_validation_fold_num An integer for the number of folds in
#'        cross-validation. It must be greater than 1 and less than or equal to
#'        the number of observations.
#' @param risk_estimate_method A character string specifying the method for
#'        estimating risk, such as "parametric_bootstrap" or other options,
#'        depending on the family of the GLM.
#' @param discrepancy_method A character string specifying the method for
#'        calculating discrepancy. The valid options depend on the GLM family
#'        and risk estimation method.
#' @param binomial_joint_theta Logical; if TRUE, uses joint theta (theta = 1/tau) in Binomial models.
#' @param binomial_joint_alpha Logical; if TRUE, uses joint alpha (adaptive tau_alpha) in Binomial models.
#' @param binomial_tau_lower A positive numeric value specifying the lower
#'        bound for tau in binomial GLMs. It must be greater than zero.
#' @param tau_alpha A positive numeric value for the tau alpha parameter.
#' @param tau_gamma A positive numeric value for the tau gamma parameter.
#' @param gibbs_iter An integer for the number of Gibbs iterations in the
#'        sampling process. It must be greater than zero.
#' @param gibbs_warmup An integer for the number of warm-up iterations in the
#'        Gibbs sampling. It must be positive and less than the total number
#'        of iterations.
#' @param coefs_iter An integer specifying the number of iterations for the
#'        coefficient update in the Gibbs sampling. It must be positive.
#' @param gaussian_variance_alpha The shape parameter for the inverse-gamma prior on
#'        variance if the variance is unknown in Gaussian models. It must be positive.
#' @param gaussian_variance_beta  The scale parameter for the inverse-gamma prior on
#'        variance if the variance is unknown in Gaussian models. It must be positive.
#' @details
#' This function performs several checks to ensure the validity of the input
#' parameters:
#'   - Ensures that `tau`, `tau_0`, `parametric_bootstrap_iteration_times`,
#'         `binomial_tau_lower`, `tau_alpha`, `tau_gamma`, `gibbs_iter`,
#'         `gibbs_warmup`, and `coefs_iter` are positive values.
#'   - Verifies that `cat_init` is an object generated by
#'         `cat_glm_initialization`.
#'   - Checks that the `formula` response name matches the response name
#'         used in the `cat_init` object.
#'   - Verifies that `risk_estimate_method` and `discrepancy_method` are
#'         compatible with the GLM family and that no invalid combinations are
#'         used.
#'   - Warns if the dataset size is too large for the specified risk
#'         estimation method.
#' If any of these conditions are not met, the function raises an error or warning to guide the user.
#'
#' @return Returns nothing if all checks pass; otherwise, raises an error or warning.
validate_glm_input <- function(formula,
                               cat_init,
                               tau = NULL,
                               tau_seq = NULL,
                               tau_0 = NULL,
                               parametric_bootstrap_iteration_times = NULL,
                               cross_validation_fold_num = NULL,
                               risk_estimate_method = NULL,
                               discrepancy_method = NULL,
                               binomial_joint_theta = FALSE,
                               binomial_joint_alpha = FALSE,
                               binomial_tau_lower = NULL,
                               tau_alpha = NULL,
                               tau_gamma = NULL,
                               gibbs_iter = NULL,
                               gibbs_warmup = NULL,
                               coefs_iter = NULL,
                               gaussian_variance_alpha = NULL,
                               gaussian_variance_beta = NULL) {
  # Check if below parameters is positive
  validate_positive("tau", tau, incl_0 = TRUE)
  validate_positive("tau_0", tau_0)
  validate_positive("parametric_bootstrap_iteration_times", parametric_bootstrap_iteration_times)
  validate_positive("cross_validation_fold_num", cross_validation_fold_num)
  validate_positive("binomial_tau_lower", binomial_tau_lower, incl_0 = TRUE)
  validate_positive("tau_alpha", tau_alpha)
  validate_positive("tau_gamma", tau_gamma)
  validate_positive("gibbs_iter", gibbs_iter)
  validate_positive("gibbs_warmup", gibbs_warmup)
  validate_positive("coefs_iter", coefs_iter)
  validate_positive("gaussian_variance_alpha", gaussian_variance_alpha)
  validate_positive("gaussian_variance_beta", gaussian_variance_beta)
  validate_positive("tau_seq", tau_seq, incl_0 = TRUE, is_vector = TRUE)

  # Check if cat_init is generated from cat_glm_initialization
  if (!(inherits(cat_init, "cat_initialization") &
    cat_init$function_name == "cat_glm_initialization")) {
    stop(
      paste0(
        "The provided `cat_init` must be generated from `cat_glm_initialization`. <",
        cat_init$function_name,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if formula has the same response name as that from cat_init
  f_lhs <- get_formula_lhs(formula)
  if (f_lhs != "NULL" & f_lhs != cat_init$y_col_name) {
    stop(
      paste0(
        "The response name from provided `formula` <",
        get_formula_string(formula),
        "> should be same as the response name used in the provided `cat_init` object."
      ),
      call. = FALSE
    )
  }


  # Specific check for cat_glm_bayes_joint_gibbs
  if (!is.null(gibbs_iter) & !is.null(gibbs_warmup)) {
    # Check if user apply gibbs to linear regression
    if (cat_init$family != "binomial") {
      stop(
        "Please note that only logistic regression can use `cat_glm_bayes_joint_gibbs`.",
        call. = FALSE
      )
    }

    # Check if warmup is larger than iter
    if (gibbs_warmup >= gibbs_iter) {
      stop(
        paste0(
          "The provided `warmup` should be less than `iter`. <",
          "iter = ",
          gibbs_iter,
          "; warmup =",
          gibbs_warmup,
          "> is invalid."
        ),
        call. = FALSE
      )
    }

    # Check if coefs_iter too large in gibbs
    if (coefs_iter > 10) {
      response <- utils::menu(
        c("YES", "NO"),
        title = paste0(
          "\nPlease note that large `coefs_iter` can slow down gibb's sampling process.\n",
          "Do you want to continue?\n"
        )
      )

      if (response == 2) {
        stop("\nProcess terminated by user.", call. = FALSE)
      }
    }
  }

  # Specific check for cat_glm_tune
  if (!is.null(risk_estimate_method) & !is.null(discrepancy_method)) {
    # Initialize an empty vector for validation errors
    validation_errors <- c()

    # Append messages based on conditions
    if (cat_init$family == "gaussian" & risk_estimate_method == "steinian_estimate") {
      validation_errors <- c(
        validation_errors,
        paste0(
          "Only logistic regression can use `risk_estimate_method = steinian_estimate`.\n",
          "SUGGESTION: Use `mallowian_estimate` for linear regression."
        )
      )
    }

    if (cat_init$family == "binomial" & risk_estimate_method == "mallowian_estimate") {
      validation_errors <- c(
        validation_errors,
        paste0(
          "Only linear regression can use `risk_estimate_method = mallowian_estimate`.\n",
          "SUGGESTION: Use `steinian_estimate` for logistic regression."
        )
      )
    }

    if (cat_init$family == "gaussian" & discrepancy_method != "mean_square_error") {
      validation_errors <- c(
        validation_errors,
        paste0(
          "Only logistic regression can use `discrepancy_method = mean_classification_error` or `logistic_deviance`.\n",
          "SUGGESTION: Use `mean_square_error` for linear regression."
        )
      )
    }

    if (cat_init$family == "binomial" & discrepancy_method == "mean_square_error") {
      validation_errors <- c(
        validation_errors,
        paste0(
          "Only linear regression can use `discrepancy_method = mean_square_error`.\n",
          "SUGGESTION: Use `logistic_deviance` for logistic regression."
        )
      )
    }

    # Stop the process if there are any validation errors
    if (length(validation_errors) > 0) {
      stop(paste(validation_errors, collapse = "\n"), call. = FALSE)
    }

    # Set data size thresholds based on the risk estimation method
    size_threshold <- switch(risk_estimate_method,
      "parametric_bootstrap" = 1000,
      "steinian_estimate" = 800,
      Inf
    )

    # Check if the observational dataset size exceeds the threshold
    if (cat_init$obs_size > size_threshold) {
      response <- utils::menu(
        choices = c("YES", "NO"),
        title = paste0(
          "The large size of the dataset may lead to longer estimation times.\n",
          "Do you want to continue?"
        )
      )
      # Terminate process if the user chooses "NO"
      if (response == 2) {
        stop("\nProcess terminated by user.", call. = FALSE)
      }
    }

    # Check if cross_validation_fold_num is less than 2 or larger than size of observation data
    if (risk_estimate_method == "cross_validation" && (
      cross_validation_fold_num < 2 | cross_validation_fold_num > cat_init$obs_size)) {
      stop(
        paste0(
          "The provided `cross_validation_fold_num` must be larger than 1 and smaller than the observation data size. <",
          cross_validation_fold_num,
          "> is invalid."
        ),
        call. = FALSE
      )
    }

    # Warn User that will apply "logistic_deviance" for "steinian_estimate"
    if (cat_init$family == "binomial" &
      risk_estimate_method == "steinian_estimate" &
      discrepancy_method != "logistic_deviance") {
      warning(
        paste0(
          "Only `discrepancy_method = logistic_deviance` can be applied to `risk_estimate_method = steinian_estimate`.\n",
          "Will continue use `discrepancy_method = logistic_deviance` instead of <",
          discrepancy_method,
          ">."
        ),
        call. = FALSE
      )
    }

    # Warn user if dataset size is small for cross-validation
    if (risk_estimate_method == "cross_validation" && cat_init$obs_size < (5 * cross_validation_fold_num)) {
      warning(
        paste0(
          "The size of the observation dataset is relatively small for cross-validation to yield optimal results.\n",
          "SUGGESTION: Please consider utilizing alternative methods (e.g., risk_estimate_method = parametric_bootstrap) ",
          "or reduce `cross_validation_fold_num` for better estimation."
        ),
        call. = FALSE
      )
    }
  }

  if ((!is.null(gaussian_variance_alpha) | !is.null(gaussian_variance_beta)) & (cat_init$gaussian_known_variance | cat_init$family != "gaussian")) {
    warning(
      paste0(
        "The provided `gaussian_variance_alpha` or\\and `gaussian_variance_beta` only be applicable when ",
        "`gaussian_known_variance = FALSE` and `family = gaussian`."
      ),
      call. = FALSE
    )
  }

  if (cat_init$family != "binomial" & (binomial_joint_theta | binomial_joint_alpha)) {
    warning(
      paste0(
        "The provided `binomial_joint_theta` or\\and `binomial_joint_alpha` only be applicable when ",
        "`family = binomial`."
      ),
      call. = FALSE
    )
  }

  if (cat_init$family == "binomial" & binomial_joint_alpha & !binomial_joint_theta) {
    warning(
      paste0(
        "To activate `binomial_joint_alpha = TRUE` feature, ",
        "both `binomial_joint_theta` and `binomial_joint_alpha` should be assigned as `TRUE`."
      ),
      call. = FALSE
    )
  }
}

#' Validate Inputs for Catalytic Generalized Linear Models (GLMs) Initialization
#'
#' This function validates the input parameters required for initializing a catalytic Generalized Linear Model (GLM).
#' It ensures the appropriate structure and compatibility of the formula, family, data, and additional parameters
#' before proceeding with further modeling.
#'
#' @param formula A formula object specifying the \code{stats::glm} model to be fitted. It must not contain random effects or survival terms.
#' @param family A character or family object specifying the error distribution and link function. Valid values are "binomial" and "gaussian".
#' @param data A \code{data.frame} containing the data to be used in the GLM.
#' @param syn_size A positive integer specifying the sample size used for the synthetic data.
#' @param custom_variance A positive numeric value for the custom variance used in the model (only applicable for Gaussian family).
#' @param gaussian_known_variance A logical indicating whether the variance is known for the Gaussian family.
#' @param x_degree A numeric vector specifying the degree of the predictors. Its length should match the number of predictors (excluding the response variable).
#'
#' @details
#' This function performs the following checks:
#'   - Ensures that \code{syn_size}, \code{custom_variance}, and \code{x_degree} are positive values.
#'   - Verifies that the provided \code{formula} is suitable for GLMs, ensuring no random effects or survival terms.
#'   - Checks that the provided \code{data} is a \code{data.frame}.
#'   - Confirms that the \code{formula} does not contain too many terms relative to the number of columns in \code{data}.
#'   - Ensures that the \code{family} is either "binomial" or "gaussian".
#'   - Validates that \code{x_degree} has the correct length relative to the number of predictors in \code{data}.
#'   - Warns if \code{syn_size} is too small relative to the number of columns in \code{data}.
#'   - Issues warnings if \code{custom_variance} or \code{gaussian_known_variance} are used with incompatible families.
#' If any of these conditions are not met, the function raises an error or warning to guide the user.
#'
#' @return Returns nothing if all checks pass; otherwise, raises an error or warning.
validate_glm_initialization_input <- function(formula,
                                              family,
                                              data,
                                              syn_size,
                                              custom_variance,
                                              gaussian_known_variance,
                                              x_degree) {
  # Check if below parameters is positive
  validate_positive("syn_size", syn_size)
  validate_positive("custom_variance", custom_variance)
  validate_positive("x_degree", x_degree, is_vector = TRUE)

  # Check if formula is for GLMs
  if (grepl("\\(.*\\|.*\\)", deparse(formula)) || grepl("Surv", deparse(formula))) {
    stop(
      paste0(
        "The provided formula <",
        get_formula_string(formula),
        "> is not appropriate for Generalized Linear Models (GLMs)."
      ),
      call. = FALSE
    )
  }

  # Check if data is data.frame
  if (!is.data.frame(data)) {
    stop(
      paste0(
        "The provided `data` must be a `data.frame`. <`",
        class(data)[1],
        "`> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if the formula is overly complex
  if (length(attr(stats::terms(formula), "term.labels")) >= (ncol(data) - 1)) {
    stop(
      paste0(
        "The provided `formula` must contain fewer terms than the number of columns in `data`. <",
        get_formula_string(formula),
        "> is too complex."
      ),
      call. = FALSE
    )
  }

  # Check if family is either binomial or gaussian
  family_string <- get_glm_family_string(family)
  if (!(family_string %in% c("binomial", "gaussian"))) {
    stop(
      paste0(
        "The provided `family` must be either `binomial` or `gaussian`. <",
        family,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if x_degree has same size to the data covariance
  if (!is.null(x_degree) & (length(x_degree) != (ncol(data) - 1))) {
    stop(
      paste0(
        "The provided `x_degree` should be a vector of positive values ",
        "with a length less than the number of covariates. <",
        paste(x_degree, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  if (!is.null(syn_size) && syn_size < (4 * ncol(data))) {
    warning(
      "The provided `syn_size` is too small, recommand to choose a larger value.",
      call. = FALSE
    )
  }


  if (family_string != "gaussian" & (!is.null(custom_variance) | gaussian_known_variance)) {
    warning(
      "The `custom_variance` or `gaussian_known_variance` is not applicable when the family is not `gaussian`.",
      call. = FALSE
    )
  }

  if (family_string == "gaussian" & !is.null(custom_variance) & !gaussian_known_variance) {
    warning(
      "The `custom_variance` is not applicable when `gaussian_known_variance = FALSE`.",
      call. = FALSE
    )
  }
}

# ----- COX -----

#' Calculate Risk and Failure Sets for Cox Proportional Hazards Model
#'
#' This function calculates the risk and failure sets for subjects in a Cox proportional hazards model based on their time-to-event data, status, and an indicator vector.
#'
#' @param time_vector A numeric vector of time-to-event data for each subject.
#' @param status_vector A numeric vector indicating event occurrence (1 = event, 0 = censored).
#' @param indicator_vector A numeric vector representing the indicator times used to define risk and failure sets.
#'
#' @return A list containing two elements:
#' - `risk_set`: A matrix indicating which subjects are at risk at each time point.
#' - `failure_set`: A matrix indicating which subjects experienced an event at each time point.
get_cox_risk_and_failure_sets <- function(time_vector,
                                          status_vector,
                                          indicator_vector) {
  time_vector <- as.vector(time_vector)
  status_vector <- as.vector(status_vector)
  indicator_vector <- as.vector(indicator_vector)

  # Function to find the index where an event occurs or reaches censoring time
  find_occur_idx <- function(time) {
    min(which(indicator_vector >= time))
  }


  # Initialize risk and failure sets as zero matrices
  risk_set <- failure_set <- matrix(
    0,
    length(time_vector),
    length(indicator_vector)
  )

  # Process subjects who experienced an event (status = 1)
  for (i in which(status_vector == 1)) {
    # Mark the exact time of the event in the failure set
    failure_set[i, find_occur_idx(time_vector[i])] <- 1
    # Mark all times up to the event time in the risk set
    risk_set[i, 1:find_occur_idx(time_vector[i])] <- 1
  }

  # Mark all times up to the censoring time in the risk set
  risk_set[which(status_vector == 0), 1:find_occur_idx(time_vector[i])] <- 1

  # Return the risk and failure sets as a list
  return(list(
    risk_set = risk_set,
    failure_set = failure_set
  ))
}

#' Estimate the kappa value for the synthetic Cox proportional hazards model
#'
#' This function iterative estimates the kappa value for the synthetic Cox proportional hazards model using a vectorized approach for efficiency.
#'
#' @param X A matrix of covariates with rows representing observations and columns representing features.
#' @param time A vector of time-to-event data.
#' @param status A vector indicating event occurrence (1 = event, 0 = censored).
#' @param hazard_constant A scalar representing the hazard constant. Defaults to NULL, in which case it's calculated internally.
#'
#' @return A numeric value representing the estimated kappa for the synthetic Cox model.
get_cox_kappa <- function(X,
                          time,
                          status,
                          hazard_constant) {
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Initialize coefficients
  coefs <- rep(0, ncol(X))

  diff_new_old <- TRUE
  iter <- 0

  while ((iter < 50) & diff_new_old) {
    coefs_old <- coefs

    # Calculate linear predictor
    exp_lp <- exp(c(get_linear_predictor(X, coefs_old)))

    # Calculate gradient and Hessian
    gd_old <- get_cox_syn_gradient(X, time, coefs_old, hazard_constant)
    hess_old <- get_cox_syn_hessian(X, time, coefs_old, hazard_constant)

    # Update coefficients using QR solve
    coefs <- coefs_old - get_cox_qr_solve(hess_old, gd_old)

    iter <- iter + 1

    # Check for convergence
    diff_new_old <- norm(as.numeric(coefs_old - coefs), "2") > -5 * norm(as.numeric(coefs_old), "2")
  }

  # Calculate the final linear predictor
  lp_new <- c(get_linear_predictor(X, coefs))

  return(sum(lp_new - exp(lp_new) * hazard_constant * time) / nrow(X))
}

#' Compute the Partial Likelihood for the Cox Proportional Hazards Model
#'
#' This function calculates the partial likelihood for the Cox proportional hazards model. The partial likelihood is computed for the censored observations in the dataset.
#'
#' @param X A matrix of covariates with rows representing observations and columns representing features.
#' @param time A vector of time-to-event data.
#' @param status A vector indicating the event status (1 for event occurred, 0 for censored).
#' @param coefs A vector of regression coefficients.
#' @param entry_points A vector of entry points (optional). Defaults to NULL, in which case a vector of zeros is used.
#'
#' @return A numeric scalar representing the partial likelihood of the Cox model.
get_cox_partial_likelihood <- function(X,
                                       time,
                                       status,
                                       coefs,
                                       entry_points) {
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  pl <- 0

  # Calculate partial likelihood for each censored observation
  for (i in which(status == 1)) {
    risk_set_idx <- get_cox_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      status_vector = status
    )

    # Compute the linear predictors for the risk set
    exp_risk_lp <- exp(
      c(get_linear_predictor(
        X[risk_set_idx, , drop = FALSE],
        coefs
      ))
    )

    # Update partial likelihood
    pl <- pl + sum(X[i, , drop = FALSE] * coefs) - log(sum(exp_risk_lp))
  }

  return(pl)
}


#' Solve Linear System using QR Decomposition
#'
#' This function solves the linear system defined by \code{hessian_matrix} and \code{gradient_vector}
#' using QR decomposition. Any NA values in the resulting solution vector are replaced with 0.0001. If there
#' is an error during the solution process, a vector of default values (0.0001) is returned instead.
#'
#' @param hessian_matrix A matrix of coefficients representing the system of linear equations.
#' @param gradient_vector A numeric vector representing the constants in the linear system.
#'
#' @return A numeric vector representing the solution to the linear system. NA values in the solution
#' are replaced with a small value (0.0001). If an error occurs during solving, a vector of default values
#' (0.0001) is returned.
get_cox_qr_solve <- function(hessian_matrix,
                             gradient_vector) {
  # Calculate solution using QR decomposition with error handling
  qr_solution <- tryCatch(
    {
      # Attempt to solve the linear system
      solve(hessian_matrix, gradient_vector)
    },
    error = function(e) {
      # In case of error, return a vector of default values
      rep(0.0001, length(gradient_vector))
    }
  )

  # Replace any NA values in the solution with 0.0001
  qr_solution[is.na(qr_solution)] <- 0.0001

  return(qr_solution)
}


#' Compute the Gradient for Cox Proportional Hazards Model
#'
#' This function computes the gradient for the Cox proportional hazards model. The gradient
#' is calculated by considering the contributions of each observation to the gradient based on the
#' risk set at each event time.
#'
#' @param X A matrix of covariates (design matrix) for the Cox model.
#' @param time A numeric vector of event times.
#' @param status A numeric vector of event indicators (1 for event, 0 for censored).
#' @param coefs A numeric vector of coefficients for the Cox model.
#' @param entry_points A numeric vector of entry times for the subjects. Defaults to 0.
#'
#' @return A numeric vector representing the gradient of the Cox proportional hazards model.
get_cox_gradient <- function(X,
                             time,
                             status,
                             coefs,
                             entry_points) {
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Initialize gradient matrix
  gradient <- matrix(NA, nrow = ncol(X), ncol = nrow(X))

  # Loop over each observation to compute the gradient
  for (i in 1:(nrow(X))) {
    # Get the indices of the risk set for the current observation
    risk_set_idx <- get_cox_risk_set_idx(
      time[i],
      entry_points,
      time,
      status
    )

    # Skip if the risk set is empty
    if (length(risk_set_idx) == 0) next

    # Extract risk set covariates and compute exponential risk
    risk_set <- as.matrix(X[risk_set_idx, , drop = FALSE])
    exp_risk_set <- exp(risk_set %*% coefs)

    # Compute the gradient for the current observation
    gradient[, i] <- c(X[i, , drop = FALSE]) -
      c(t(risk_set) %*% exp_risk_set / sum(exp_risk_set))
  }

  # Return the gradient vector by summing over all observations weighted by the status vector
  return(c(gradient %*% status))
}

#' Compute the gradient of the synthetic Cox proportional hazards model
#'
#' This function calculates the gradient of the synthetic Cox proportional hazards model using a vectorized approach.
#'
#' @param X A matrix of covariates with rows representing observations and columns representing features.
#' @param time A vector of time-to-event data.
#' @param coefs A vector of regression coefficients.
#' @param hazard_constant A scalar representing the hazard constant.
#'
#' @return A numeric vector representing the gradient of the synthetic Cox model.
get_cox_syn_gradient <- function(X,
                                 time,
                                 coefs,
                                 hazard_constant) {
  X <- as.matrix(X)
  time <- as.vector(time)

  # Calculate linear predictors
  syn_lp <- c(get_linear_predictor(X, coefs))

  # Calculate the gradient in a vectorized manner
  gradient <- c(colSums((1 - hazard_constant * time * exp(syn_lp)) * X))

  return(gradient)
}

#' Identify the risk set indices for Cox proportional hazards model
#'
#' This function returns the indices of the risk set for a given time of interest in the Cox proportional hazards model.
#'
#' @param time_of_interest A numeric value representing the time at which the risk set is calculated.
#' @param entry_vector A numeric vector representing the entry times of subjects.
#' @param time_vector A numeric vector representing the time-to-event or censoring times of subjects.
#' @param status_vector A numeric vector indicating event occurrence (1) or censoring (0) for each subject.
#'
#' @return A vector of indices representing the subjects at risk at the specified time of interest.
get_cox_risk_set_idx <- function(time_of_interest,
                                 entry_vector,
                                 time_vector,
                                 status_vector) {
  time_vector <- as.vector(time_vector)
  entry_vector <- as.vector(entry_vector)
  status_vector <- as.vector(status_vector)

  # Find indices where subjects are at risk at the given time of interest
  return(which((time_of_interest >= entry_vector) & (
    (time_vector == time_of_interest & status_vector == 1) | (
      time_vector + 1e-08 > time_of_interest))))
}

#' Compute the Hessian Matrix for Cox Proportional Hazards Model
#'
#' This function computes the Hessian matrix of the Cox proportional hazards model,
#' which is used for estimating the covariance matrix of the coefficients. The Hessian is calculated
#' by summing contributions from each event time in the risk set.
#'
#' @param X A matrix of covariates (design matrix) for the Cox model.
#' @param time A numeric vector of event times.
#' @param status A numeric vector of event indicators (1 for event, 0 for censored).
#' @param coefs A numeric vector of coefficients for the Cox model.
#' @param entry_points A numeric vector of entry times for the subjects. Defaults to 0.
#'
#' @return A matrix representing the negative Hessian of the Cox model.
get_cox_hessian <- function(X,
                            time,
                            status,
                            coefs,
                            entry_points) {
  X <- as.matrix(X)
  time <- as.vector(time)
  status <- as.vector(status)

  # Initialize the Hessian matrix
  hessian <- matrix(0, nrow = ncol(X), ncol = ncol(X))

  # Calculate the exponential of the linear predictor
  exp_risk_lp <- exp(c(get_linear_predictor(X, coefs)))

  # Compute the Hessian matrix contributions for each censored event
  for (i in which(status == 1)) {
    # Determine the risk set for the current time of interest
    risk_set_idx <- get_cox_risk_set_idx(
      time_of_interest = time[i],
      entry_vector = entry_points,
      time_vector = time,
      status_vector = status
    )

    # Skip if the risk set is empty
    if (length(risk_set_idx) == 0) next

    # Extract the risk set and compute relevant quantities
    risk_set <- as.matrix(X[risk_set_idx, , drop = FALSE])
    exp_risk_lp_risk_set <- exp_risk_lp[risk_set_idx]
    sum_exp_risk_lp <- sum(exp_risk_lp_risk_set)

    # Compute the Hessian matrix contribution
    risk_exp_risk_lp <- t(risk_set) %*% exp_risk_lp_risk_set
    risk_sqrt_exp_risk_lp <- t(sweep(risk_set, 1, sqrt(exp_risk_lp_risk_set), "*"))

    hessian <- hessian +
      tcrossprod(risk_sqrt_exp_risk_lp) / sum_exp_risk_lp -
      tcrossprod(risk_exp_risk_lp) / sum_exp_risk_lp^2
  }

  return(-hessian)
}

#' Compute the Synthetic Hessian Matrix for Cox Proportional Hazards Model
#'
#' This function computes the synthetic Hessian matrix for the Cox proportional hazards model.
#' The Hessian is calculated by summing the contributions from each individual observation, scaled by
#' the hazard constant and the time of the event.
#' @param X A matrix of covariates (design matrix) for the Cox model.
#' @param time A numeric vector of event times.
#' @param coefs A numeric vector of coefficients for the Cox model.
#' @param hazard_constant A numeric value representing the hazard constant.
#'
#' @return A matrix representing the synthetic Hessian of the Cox model.
get_cox_syn_hessian <- function(X,
                                time,
                                coefs,
                                hazard_constant) {
  X <- as.matrix(X)
  time <- as.vector(time)

  # Calculate the linear predictors
  syn_lp <- c(get_linear_predictor(X, coefs))

  # Compute the Hessian using a vectorized approach
  hessian <- hazard_constant *
    t(X) %*% (sweep(X, 1, time * exp(syn_lp), "*"))

  return(-hessian)
}

#' Validate Inputs for Catalytic Cox Model
#'
#' This function validates the parameters provided for setting up a catalytic Cox proportional hazards model
#' with an initialization object created by \code{cat_cox_initialization}.
#'
#' @param formula An object of class \code{formula}. Specifies the model structure for the Cox model, including a \code{Surv} object for survival analysis. Should at least include response variance.
#' @param cat_init An initialization object generated by \code{cat_cox_initialization}. This object should contain necessary information about the dataset, including the time and status column names.
#' @param tau Optional. A numeric scalar, the regularization parameter for the Cox model. Must be positive.
#' @param tau_seq Optional. A numeric vector for specifying a sequence of regularization parameters. Must be positive.
#' @param init_coefficients Optional. A numeric vector of initial coefficients for the Cox model. Should match the number of predictors in the dataset.
#' @param tol Optional. A positive numeric value indicating the tolerance level for convergence in iterative fitting.
#' @param max_iter Optional. A positive integer indicating the maximum number of iterations allowed in the model fitting.
#' @param cross_validation_fold_num Optional. A positive integer specifying the number of folds for cross-validation. Should be greater than 1 and less than or equal to the size of the dataset.
#' @param hazard_beta Optional. A positive numeric value representing a constant for adjusting the hazard rate in the Cox model.
#' @param tau_alpha Optional. A positive numeric value controlling the influence of \code{tau}.
#' @param tau_gamma Optional. A positive numeric value controlling the influence of \code{tau_seq}.
#'
#' @details
#' This function checks:
#'  - That \code{tau}, \code{tol}, \code{max_iter}, \code{cross_validation_fold_num}, \code{hazard_beta}, \code{tau_alpha}, and \code{tau_gamma} are positive.
#'  - That \code{tau_seq} is a non-negative vector.
#'  - That \code{cat_init} is generated from \code{cat_cox_initialization}.
#'  - That \code{formula} uses the same time and status column names as those in \code{cat_init}.
#'  - That \code{init_coefficients} has the correct length for the number of predictors.
#'  - That \code{cross_validation_fold_num} is between 2 and the dataset size.
#'  - That the dataset is sufficiently large for cross-validation, recommending fewer folds if it is not.
#' If any conditions are not met, the function will raise an error or warning.
#'
#' @return Returns nothing if all checks pass; otherwise, raises an error or warning.
validate_cox_input <- function(formula,
                               cat_init,
                               tau = NULL,
                               tau_seq = NULL,
                               init_coefficients = NULL,
                               tol = NULL,
                               max_iter = NULL,
                               cross_validation_fold_num = NULL,
                               hazard_beta = NULL,
                               tau_alpha = NULL,
                               tau_gamma = NULL) {
  # Check if below parameters is positive
  validate_positive("tau", tau)
  validate_positive("tol", tol)
  validate_positive("max_iter", max_iter)
  validate_positive("cross_validation_fold_num", cross_validation_fold_num)
  validate_positive("hazard_beta", hazard_beta)
  validate_positive("tau_alpha", tau_alpha)
  validate_positive("tau_gamma", tau_gamma)
  validate_positive("tau_seq", tau_seq, is_vector = TRUE)

  # Check if cat_init is generated from cat_cox_initialization
  if (!(inherits(cat_init, "cat_initialization") &
    cat_init$function_name == "cat_cox_initialization")) {
    stop(
      paste0(
        "The provided `cat_init` must be generated from `cat_cox_initialization`. <",
        cat_init$function_name,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if formula has the same time and status column name as that from cat_init
  f_lhs <- get_formula_lhs(formula)
  time_col_name <- gsub(
    "`",
    "",
    regmatches(
      f_lhs,
      regexec(
        "\\(([^,]+),",
        f_lhs
      )
    )[[1]][2]
  )

  status_col_name <- gsub(
    "`",
    "",
    regmatches(
      f_lhs,
      regexec(
        ",\\s*([^\\)]+)\\)",
        f_lhs
      )
    )[[1]][2]
  )

  if (f_lhs != "NULL" & (time_col_name != cat_init$time_col_name | status_col_name != cat_init$status_col_name)) {
    stop(
      paste0(
        "The time column or/and status column name from provided `formula` <",
        get_formula_string(formula),
        "> should be same as the response name used in the provided `cat_init` object."
      ),
      call. = FALSE
    )
  }

  # Check if the length of input coefs is match to the data.
  if (!is.null(init_coefficients) && length(init_coefficients) != ncol(cat_init$adj_x)) {
    stop(
      paste0(
        "The provided `init_coefficients` should be same length of number of columns from data."
      ),
      call. = FALSE
    )
  }

  # Check if cross_validation_fold_num is less than 2 or larger than size of observation data
  if (!is.null(cross_validation_fold_num) && (
    cross_validation_fold_num < 2 | cross_validation_fold_num > cat_init$obs_size)
  ) {
    stop(
      paste0(
        "The provided `cross_validation_fold_num` must be larger than 1 and smaller than the observation data size. <",
        cross_validation_fold_num,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if dataset size is small for cross-validation
  if (!is.null(cross_validation_fold_num) && (cat_init$obs_size < (5 * cross_validation_fold_num))) {
    warning(
      paste0(
        "The size of the observation data is relatively small for cross-validation to yield results.\n",
        "SUGGESTION: Please consider reducing `cross_validation_fold_num` for better estimation."
      ),
      call. = FALSE
    )
  }
}

#' Validate Inputs for Catalytic Cox proportional hazards model (COX) Initialization
#'
#' This function performs validation checks on input parameters for initializing a catalytic Cox proportional hazards model.
#' It ensures that essential parameters meet requirements, such as being of the correct type, appropriate length, and having valid values.
#'
#' @param formula An object of class \code{formula}. The model formula specifying the Cox model structure. It must contain a \code{Surv} object to indicate survival analysis.
#' @param data A \code{data.frame} containing the dataset to be used for model fitting. It should include all variables referenced in \code{formula}.
#' @param syn_size A positive integer indicating the size of the synthetic dataset. It is recommended that this value is at least four times the number of columns in \code{data}.
#' @param hazard_constant A positive numeric value representing the hazard constant for the Cox model.
#' @param entry_points A numeric vector representing entry times for observations. This vector should be non-negative and have a length equal to the number of rows in \code{data}.
#' @param x_degree A numeric vector indicating degrees for each covariate. It should be non-negative and match the number of covariates (i.e., \code{ncol(data) - 2}).
#'
#' @details
#' This function checks:
#'   - That \code{syn_size}, \code{hazard_constant}, \code{entry_points}, and \code{x_degree} are positive values.
#'   - That \code{formula} includes a \code{Surv} object to be suitable for Cox models.
#'   - That \code{data} is a \code{data.frame}.
#'   - The complexity of \code{formula} to ensure it has fewer terms than the number of columns in \code{data}.
#'   - The length of \code{x_degree} and \code{entry_points} to match the dimensions of \code{data}.
#' If the conditions are not met, descriptive error messages are returned.
#'
#' @return Returns nothing if all checks pass; otherwise, raises an error.
validate_cox_initialization_input <- function(formula,
                                              data,
                                              syn_size,
                                              hazard_constant,
                                              entry_points,
                                              x_degree) {
  # Check if below parameters is positive
  validate_positive("syn_size", syn_size)
  validate_positive("hazard_constant", hazard_constant)
  validate_positive("entry_points", entry_points, incl_0 = TRUE, is_vector = TRUE)
  validate_positive("x_degree", x_degree, is_vector = TRUE)

  # Check if formula is for COX
  if (!grepl("Surv", deparse(formula))) {
    stop(
      paste0(
        "The provided formula <",
        get_formula_string(formula),
        "> is not appropriate for a Cox Regression Models (COX)."
      ),
      call. = FALSE
    )
  }

  # Check if data is data.frame
  if (!is.data.frame(data)) {
    stop(
      paste0(
        "The provided `data` must be a `data.frame`. <`",
        class(data)[1],
        "`> is invalid."
      ),
      call. = FALSE
    )
  }
  # Check if the formula is overly complex
  if (length(attr(stats::terms(formula), "term.labels")) >= (ncol(data) - 2)) {
    stop(
      paste0(
        "The provided `formula` must contain fewer terms than the number of columns in `data`. <",
        get_formula_string(formula),
        "> is too complex."
      ),
      call. = FALSE
    )
  }

  # Check if x_degree is positive and same size to the data covariance
  if (!is.null(x_degree) && (length(x_degree) != (ncol(data) - 2))) {
    stop(
      paste0(
        "The provided `x_degree` should be a vector of non-negtive values ",
        "with a length same as the number of covariates. <",
        paste(x_degree, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if entry_points is positive and same size to the data covariance
  if (!is.null(entry_points) && (length(entry_points) != nrow(data))) {
    stop(
      paste0(
        "The provided `entry_points` should be a vector of values ",
        "with a length same as the number of covariates. <",
        paste(entry_points, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  if (!is.null(syn_size) && syn_size < (4 * ncol(data))) {
    warning(
      "The provided `syn_size` is too small, recommand to choose a larger value.",
      call. = FALSE
    )
  }
}

# ----- LMM -----
#' Calculates the log-likelihood for linear mixed models (LMMs) by combining
#' observed and synthetic log-likelihoods based on the variance parameters.
#'
#' This function evaluates the log-likelihood of observed and synthetic data,
#' using residual and random-effect variance terms to determine the fit
#' of variance parameters in the mixed model context.
#'
#' @param residual_variance Numeric, the variance associated with the residual errors.
#' @param random_effect_variance Numeric, the variance associated with random effects.
#' @param obs_z_eigenvalues Vector, eigenvalues of the observed Z matrix of data.
#' @param syn_z_eigenvalues Vector, eigenvalues of the synthetic Z matrix of data.
#' @param obs_adjusted_residuals Vector, adjusted residuals of observed data.
#' @param syn_adjusted_residuals Vector, adjusted residuals of synthetic data.
#' @param tau Numeric, weight factor for the synthetic data.
#'
#' @return The sum of observed and synthetic log-likelihoods.
update_lmm_variance <- function(residual_variance,
                                random_effect_variance,
                                obs_z_eigenvalues,
                                syn_z_eigenvalues,
                                obs_adjusted_residuals,
                                syn_adjusted_residuals,
                                tau) {
  obs_log_likelihood <- sum(
    log(random_effect_variance * obs_z_eigenvalues + residual_variance) +
      obs_adjusted_residuals^2 / (random_effect_variance * obs_z_eigenvalues + residual_variance)
  )

  syn_log_likelihood <- tau / length(syn_z_eigenvalues) * sum(
    log(random_effect_variance * syn_z_eigenvalues + residual_variance) +
      syn_adjusted_residuals^2 / (random_effect_variance * syn_z_eigenvalues + residual_variance)
  )

  return(obs_log_likelihood + syn_log_likelihood)
}

#' Validate Inputs for Catalytic Linear Mixed Model (LMM)
#'
#' This function validates the parameters needed for fitting a catalytic Linear Mixed Model (LMM) or Generalized Linear Model (GLM),
#' specifically for the use with the categorical initialization from `cat_lmm_initialization`.
#'
#' @param cat_init An object of class \code{cat_initialization}, typically generated from the \code{cat_lmm_initialization} function.
#' @param tau A positive numeric value specifying the penalty parameter for the model.
#' @param residual_variance_0 A positive numeric value for the initial residual variance estimate.
#' @param random_effect_variance_0 A positive numeric value for the initial random effect variance estimate.
#' @param coefs_0 A numeric vector of length equal to the number of columns in the observation matrix. This represents the initial values for the model coefficients.
#' @param optimize_domain A numeric vector of length 2 specifying the domain for the optimization procedure.
#' @param max_iter A positive integer specifying the maximum number of iterations for the optimization.
#' @param tol A positive numeric value indicating the tolerance level for convergence.
#' @param tau_seq A numeric vector representing a sequence of values for the penalty parameter.
#' @param cross_validation_fold_num A positive integer specifying the number of folds for cross-validation.
#'
#' @details
#' This function performs the following checks:
#'   - Ensures that \code{tau}, \code{tau_seq}, \code{residual_variance_0}, \code{random_effect_variance_0}, \code{optimize_domain}, \code{max_iter}, and \code{tol} are positive values.
#'   - Verifies that \code{cat_init} is an object generated by \code{cat_lmm_initialization}.
#'   - Checks if \code{coefs_0} has the same length as the number of columns in the observation matrix of \code{cat_init}.
#'   - Ensures \code{optimize_domain} is a numeric vector of length 2.
#'   - Confirms that \code{cross_validation_fold_num} is greater than 1 and less than the number of observations in \code{cat_init}.
#' If any of these conditions are not met, the function raises an error to guide the user.
#'
#' @return Returns nothing if all checks pass; otherwise, raises an error.
validate_lmm_input <- function(
    cat_init,
    tau = NULL,
    residual_variance_0 = NULL,
    random_effect_variance_0 = NULL,
    coefs_0 = NULL,
    optimize_domain = NULL,
    max_iter = NULL,
    tol = NULL,
    tau_seq = NULL,
    cross_validation_fold_num = NULL) {
  # Check if below parameters is positive
  validate_positive("tau", tau, incl_0 = TRUE)
  validate_positive("tau_seq", tau_seq, incl_0 = TRUE, is_vector = TRUE)
  validate_positive("residual_variance_0", residual_variance_0)
  validate_positive("random_effect_variance_0", random_effect_variance_0)
  validate_positive("optimize_domain", optimize_domain, incl_0 = TRUE, is_vector = TRUE)
  validate_positive("max_iter", max_iter)
  validate_positive("tol", tol)

  # Check if cat_init is generated from cat_lmm_initialization
  if (!(inherits(cat_init, "cat_initialization") &
    cat_init$function_name == "cat_lmm_initialization")) {
    stop(
      paste0(
        "The provided `cat_init` must be generated from `cat_lmm_initialization`. <",
        cat_init$function_name,
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if length of coefs_0 equal to number of columns of observation x
  if (!is.null(coefs_0) && length(coefs_0) != ncol(cat_init$obs_x)) {
    stop(
      paste0(
        "The provided `coefs_0` must have same length as the number of columns of `cat_init$obs_x`. <",
        paste(coefs_0, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if optimize_domain is length of 2
  if (!is.null(optimize_domain) && length(optimize_domain) != 2) {
    stop(
      paste0(
        "The provided `optimize_domain` must have only length 2. <",
        paste(optimize_domain, collapse = ","),
        "> is invalid."
      ),
      call. = FALSE
    )
  }


  # Check if cross_validation_fold_num is less than 2 or larger than size of observation data
  if (!is.null(cross_validation_fold_num) && (
    cross_validation_fold_num < 2 | cross_validation_fold_num > cat_init$obs_size)
  ) {
    stop(
      paste0(
        "The provided `cross_validation_fold_num` must be larger than 1 and smaller than the observation data size. <",
        cross_validation_fold_num,
        "> is invalid."
      ),
      call. = FALSE
    )
  }
}

#' Validate Inputs for Catalytic Linear Mixed Model (LMM) Initialization
#'
#' This function validates the parameters needed for initializing a catalytic Linear Mixed Model (LMM) or Generalized Linear Model (GLM)
#' based on the input formula, data, and column specifications.
#'
#' @param formula An object of class \code{formula} representing the model formula, typically including fixed and random effects for LMMs or for GLMs.
#' @param data A \code{data.frame} containing the data for model fitting. This should include all columns specified in \code{x_cols}, \code{y_col}, \code{z_cols}, and \code{group_col}.
#' @param x_cols A character vector of column names to be used as predictor variables in the model.
#' @param y_col A single character string specifying the name of the response variable column.
#' @param z_cols A character vector of column names to be used as additional predictors or grouping factors, depending on the model structure.
#' @param group_col A single character string specifying the name of the grouping variable for random effects.
#' @param syn_size Optional. A positive integer indicating the synthetic data size, typically for use in data augmentation or model diagnostics.
#'
#' @details
#' This function performs the following checks:
#'   - Ensures \code{syn_size} is a positive integer.
#'   - Verifies that \code{formula} is not for survival analysis (e.g., does not contain \code{Surv} terms).
#'   - Checks that the formula is not overly complex by confirming it has fewer terms than the total columns in \code{data}.
#'   - Ensures \code{y_col} and \code{group_col} each contain only one column name.
#'   - Confirms \code{data} is a \code{data.frame}.
#'   - Validates that all specified columns in \code{x_cols}, \code{y_col}, \code{z_cols}, and \code{group_col} exist in \code{data} without overlap or missing values.
#'   - Warns if \code{syn_size} is set too small relative to the data dimensions, recommending a larger value.
#' If any of these conditions are not met, the function raises an error or warning to guide the user.
#'
#' @return Returns nothing if all checks pass; otherwise, raises an error or warning.
validate_lmm_initialization_input <- function(
    formula,
    data,
    x_cols,
    y_col,
    z_cols,
    group_col,
    syn_size) {
  # Check if below parameters is positive
  validate_positive("syn_size", syn_size)

  # Check if formula is for GLMs or LMM
  if (grepl("Surv", deparse(formula))) {
    stop(
      paste0(
        "The provided formula <",
        get_formula_string(formula),
        "> is not appropriate for Linear Mixed Model (LMM) or Generalized Linear Model (GLMs)."
      ),
      call. = FALSE
    )
  }

  # Check if the formula is overly complex
  if (length(attr(stats::terms(formula), "term.labels")) >= (ncol(data) - 1)) {
    stop(
      paste0(
        "The provided `formula` must contain fewer terms than the number of columns in `data`. <",
        get_formula_string(formula),
        "> is too complex."
      ),
      call. = FALSE
    )
  }

  # Check if y_col is only one
  if (length(y_col) != 1) {
    stop(
      paste0(
        "The provided `y_col` must be a single column. <",
        paste(y_col, collapse = ", "),
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if group_col is only one
  if (length(group_col) != 1) {
    stop(
      paste0(
        "The provided `group_col` must be a single column. <",
        paste(group_col, collapse = ", "),
        "> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if data is data.frame
  if (!is.data.frame(data)) {
    stop(
      paste0(
        "The provided `data` must be a `data.frame`. <`",
        class(data)[1],
        "`> is invalid."
      ),
      call. = FALSE
    )
  }

  # Check if all columns has been identified and inlcuded in `data`
  colnames_subset <- c(x_cols, y_col, z_cols, group_col)
  if ((!all(colnames(data) %in% colnames_subset)
  ) | (!all(colnames_subset %in% colnames(data))
  ) | (length(unique(colnames_subset)) != length(colnames_subset))) {
    stop(
      paste0(
        "Please enter `x_cols`, `y_col`, `z_cols` and `group_col` with no overlap, ",
        "and include all column names showed in `data`. ",
        "Or, please enter right group column name in the formula"
      ),
      call. = FALSE
    )
  }

  if (!is.null(syn_size) && syn_size < (4 * ncol(data))) {
    warning(
      "The provided `syn_size` is too small, recommand to choose a larger value.",
      call. = FALSE
    )
  }
}

# ----- General -----
#' Validate Positive or Non-negative Parameter
#'
#' This function checks whether a parameter value is positive (or non-negative if `incl_0` is set to `TRUE`).
#' It can handle both single numeric values and vectors, and it raises an error with an informative
#' message if the validation fails.
#'
#' @param param_name A string representing the name of the parameter. Used in the error message.
#' @param param_value The parameter value to validate, either a single numeric or a numeric vector.
#' @param incl_0 Logical, if `TRUE`, allows non-negative values (larger or equal to 0); if `FALSE`, requires positive values (larger than 0).
#' @param is_vector Logical, if `TRUE`, treats `param_value` as a vector; otherwise, expects a single numeric value.
#'
#' @return `NULL` if validation passes; otherwise, an error is raised.
validate_positive <- function(param_name,
                              param_value,
                              incl_0 = FALSE,
                              is_vector = FALSE) {
  if (!is.null(param_value)) {
    condition <- ifelse(incl_0, any(param_value < 0), any(param_value <= 0))
    stop_add_on <- ifelse(incl_0, "non-negative", "positive")
    if (is_vector) {
      if (condition | length(param_value) == 1) {
        stop(
          paste0(
            "The provided `",
            param_name,
            "` should be a vector of ",
            stop_add_on,
            " numbers. <",
            paste(param_value, collapse = ", "),
            "> is invalid."
          ),
          call. = FALSE
        )
      }
    } else if (condition | length(param_value) != 1) {
      stop(
        paste0(
          "The provided `",
          param_name,
          "` should be a ",
          stop_add_on,
          " number. <",
          paste(param_value, collapse = ", "),
          "> is invalid."
        ),
        call. = FALSE
      )
    }
  }
}

#' Generate Stan Model Based on Specified Parameters
#'
#' This function retrieves a Stan model file based on a combination of input parameters,
#' constructs the file path, and loads the Stan model.
#'
#' @param type Character, either `"glm"` or `"cox"`, specifying the type of model to load.
#' @param glm_family_string Character, specifying the family for GLM models, either `"binomial"` or `"gaussian"`.
#'        Required for `"glm"` models, ignored for `"cox"` models.
#' @param joint_tau Logical, if `TRUE`, includes "joint" in the file name to indicate a joint model with tau parameter.
#' @param glm_binomial_joint_theta Logical, if `TRUE` and `glm_family_string` is `"binomial"`, includes "theta" in the
#'        file name for joint theta parameter.
#' @param glm_binomial_joint_alpha Logical, if `TRUE` and `glm_family_string` is `"binomial"`, includes "alpha" in the
#'        file name for joint alpha parameter.
#' @param glm_gaussian_known_variance Logical, if `TRUE` and `glm_family_string` is `"gaussian"`, includes
#'        "known_variance" in the file name to specify known variance.
#'
#' @return A compiled Stan model loaded by `rstan::stan_model`.
get_stan_model <- function(type = c("glm", "cox"),
                           glm_family_string = c("gaussian", "binomial"),
                           joint_tau = FALSE,
                           glm_binomial_joint_theta = FALSE,
                           glm_binomial_joint_alpha = FALSE,
                           glm_gaussian_known_variance = FALSE) {
  # Ensure type and glm family is valid
  type <- match.arg(type)
  glm_family_string <- match.arg(glm_family_string)

  # Construct the filename based on parameters
  stan_file <- system.file(
    "stan",
    paste0(
      paste0(
        c(
          type,
          if (type == "glm") glm_family_string,
          if (any(joint_tau, glm_binomial_joint_theta, glm_binomial_joint_alpha)) "joint",
          if (type == "glm" & glm_family_string == "binomial" & glm_binomial_joint_theta) "theta",
          if (type == "glm" & glm_family_string == "binomial" & glm_binomial_joint_theta & glm_binomial_joint_alpha) "alpha",
          if (type == "glm" & glm_family_string == "gaussian" & glm_gaussian_known_variance) "known_variance"
        ),
        collapse = "_"
      ),
      ".stan"
    ),
    package = "catalytic"
  )

  # Load and return the Stan model
  return(rstan::stan_model(file = stan_file))
}

#' Compute Linear Predictor
#'
#' This function computes the linear predictor from a matrix of predictor variables and a vector of coefficients.
#' It handles cases with and without an intercept term.
#'
#' @param X A matrix of predictor variables.
#' @param coefs A vector of coefficients. It should be either the same length as the number of columns in X
#'               (for models without an intercept) or one more than the number of columns in X (for models with an intercept).
#'
#' @return A vector of linear predictor values.
get_linear_predictor <- function(X, coefs) {
  X <- as.matrix(X)
  coefs <- as.vector(coefs)

  # Case 1: Coefficients vector length matches number of columns in X (no intercept)
  if (length(coefs) == ncol(X)) {
    return(X %*% coefs)

    # Case 2: Coefficients vector length is one more than number of columns in X (with intercept)
  } else if (length(coefs) == (ncol(X) + 1)) {
    return(coefs[1] + X %*% coefs[-1])
  }

  # Error case: dimensions do not match
  stop(
    paste0(
      "Coefficient and covariance dimensions do not match.\n",
      "Coefficients Length: ", length(coefs),
      "\nCovariance Dimention: nrow = ", nrow(X), "; ncol = ", ncol(X)
    ),
    call. = FALSE
  )
}

#' Compute Discrepancy Measures
#'
#' This function computes various discrepancy measures between observed and estimated values.
#' It supports different methods including logarithmic error, square error, classification error, and logistic deviance.
#'
#' @param discrepancy_method A character string specifying the discrepancy method to use. Options are:
#'   \describe{
#'     \item{"logarithmic_error"}{Logarithmic error, suitable for probabilities.}
#'     \item{"mean_square_error"}{Mean squared error.}
#'     \item{"mean_classification_error"}{Mean of classification error, suitable for binary outcomes.}
#'     \item{"logistic_deviance"}{Logistic deviance, computed using a GLM model.}
#'   }
#' @param family_string A GLM family in string (e.g., "binomial") used to compute logistic deviance.
#' @param X A matrix of predictor variables.
#' @param Y A vector or data frame of observed values.
#' @param coefs A vector of coefficients for the GLM model.
#' @param est_Y A vector of estimated values. If not provided, it will be computed using `get_glm_mean` with the specified `family`.
#'
#' @return A numeric value representing the discrepancy between observed and estimated values.
get_discrepancy <- function(
    discrepancy_method = c(
      "mean_logarithmic_error",
      "mean_square_error",
      "mean_classification_error",
      "logistic_deviance"
    ),
    family_string = NULL,
    X = NULL,
    Y = NULL,
    coefs = NULL,
    est_Y = NULL) {
  # Match the provided discrepancy method to one of the allowed methods
  discrepancy_method <- match.arg(discrepancy_method)

  # Compute estimated values (est_Y) if not provided, using get_glm_mean if family is specified
  if (is.null(est_Y) & !is.null(family_string)) {
    est_Y <- get_glm_mean(
      family_string = family_string,
      X = X,
      coefs = coefs
    )
  }

  # Apply specific preprocessing for logarithmic error
  if (discrepancy_method == "mean_logarithmic_error") {
    Y <- pmax(0.0001, pmin(0.9999, Y))
    est_Y <- pmax(0.0001, pmin(0.9999, est_Y))
  }

  # Compute and return the discrepancy measure based on the chosen method
  return(
    switch(discrepancy_method,
      "mean_logarithmic_error" = mean(Y * log(Y / est_Y) + (1 - Y) * log((1 - Y) / (1 - est_Y))),
      "mean_square_error" = mean((Y - est_Y)^2),
      "mean_classification_error" = mean(Y - 2 * Y * est_Y + est_Y),
      "logistic_deviance" = -get_glm_log_density(family_string = "binomial", X = X, Y = Y, coefs = coefs)
    )
  )
}

#' Adjusted Cat Initialization
#'
#' This function adjusts the categorical initialization by creating a model frame
#' for the predictors specified in the right-hand side of the formula and splits
#' the adjusted data into observed and synthetic parts.
#'
#' @param cat_init The object generated from `cat_glm_initialization`, `cat_cox_initialization` or `cat_lmm_initialization`
#' @param formula_rhs A formula specifying the right-hand side of the model for predictors.
#'
#' @return A list containing the original `cat_init` with added components:
#'     - `adj_x`: The adjusted model frame for the predictors.
#'     - `adj_obs_x`: The observed part of the adjusted predictors.
#'     - `adj_syn_x`: The synthetic part of the adjusted predictors.
get_adjusted_cat_init <- function(cat_init,
                                  formula_rhs) {
  # Generate model frame for predictors using the right-hand side of the formula
  cat_init$adj_x <- stats::model.frame(
    formula_rhs,
    as.data.frame(cat_init$x)
  )

  # Split adjusted data into observed and synthetic parts
  cat_init$adj_obs_x <- cat_init$adj_x[1:cat_init$obs_size, , drop = FALSE]
  cat_init$adj_syn_x <- cat_init$adj_x[(cat_init$obs_size + 1):cat_init$size, , drop = FALSE]

  return(cat_init)
}

#' Resampling Methods for Data Processing
#'
#' This function includes various resampling methods applied to input data for each column to prepare it
#' for analysis. These methods help to transform the data distribution and improve model fitting.
#'
#' @param data A data frame to be resampled.
#' @param resample_size An integer specifying the size of the resample.
#' @param data_degree A numeric vector indicating the degree of each column in the data (optional).
#' @param resample_only A logical value indicating whether to return only the resampled data (default is FALSE).
#'
#' @return A list containing:
#'   \item{resampled_df}{A data frame of resampled data.}
#'   \item{resampled_df_log}{A data frame recording the resampling process for each column.}
#'
#' @details
#' - **Coordinate**: This method refers to the preservation of the original data values as reference coordinates during processing.
#' It ensures that the transformations applied are based on the initial structure of the data.
#'
#' - **Deskewing**:Deskewing is the process of adjusting the data distribution to reduce skewness, making it more symmetric.
#' If the absolute value of skewness is greater than or equal to 1, deskewing techniques will be applied
#' to normalize the distribution, which can enhance model performance.
#'
#' - **Smoothing**: Smoothing techniques reduce noise in the data by averaging or modifying data points.
#' This is especially useful when there are many unique values in the original data column, as it helps
#' to stabilize the dataset and prevent overfitting during model training.
#'
#' - **Flattening**: Flattening modifies the data to create a more uniform distribution across its range.
#' This method is employed when the frequency of certain categories in categorical variables is low,
#' replacing some original values with randomly selected unique values from the dataset to reduce sparsity.
#'
#' - **Symmetrizing**: Symmetrizing adjusts the data so that it becomes more balanced around its mean.
#' This is crucial for achieving better statistical properties and improving the robustness of the model
#' fitting process.
get_resampled_df <- function(data,
                             resample_size,
                             data_degree = NULL,
                             resample_only = FALSE) {
  # Resample a column
  resample_column <- function(col, size) {
    sample(col, size = size, replace = TRUE)
  }

  # Calculate skewness
  get_skewness <- function(lst) {
    n <- length(lst)
    (sum((lst - mean(lst))^3) / n) / (sum((lst - mean(lst))^2) / n)^(3 / 2)
  }

  # Resampling methods
  get_resampling_methods <- function(x,
                                     original_x,
                                     method = c(
                                       "symmetrizing",
                                       "smoothing",
                                       "flattening",
                                       "deskewing"
                                     )) {
    # Interquartile range (IQR) calculation
    get_iqr <- function(x, invl = 0.5) {
      right_quantile <- (2 - invl) / 2
      iqr_value <- 0
      while (iqr_value == 0 & right_quantile < 1) {
        iqr_value <- diff(
          stats::quantile(
            as.numeric(x),
            c(1 - right_quantile, right_quantile),
            na.rm = FALSE,
            names = FALSE,
            type = 7
          )
        )
        right_quantile <- right_quantile + 0.05
      }

      return(iqr_value)
    }

    # Apply the selected resampling method
    switch(match.arg(method),
      "symmetrizing" = mean(original_x) + sample(
        c(-1, 1),
        length(x),
        replace = TRUE
      ) * (x - mean(original_x)),
      "smoothing" = x + stats::rnorm(
        length(x),
        0,
        min(diff(sort(unique(original_x)))) / 10
      ),
      "flattening" = replace(
        x,
        sample(length(x), length(x) / 2),
        sample(unique(original_x),
          length(x) / 2,
          replace = TRUE
        )
      ),
      "deskewing" = replace(
        x,
        sample(length(x), length(x) / 2),
        truncnorm::rtruncnorm(length(x) / 2,
          a = min(original_x),
          b = max(original_x),
          mean(original_x, trim = 0.1),
          -get_iqr(original_x) / (2 * stats::qnorm(1 / 4))
        )
      )
    )
  }

  # Process each column based on resampling methods
  process_column <- function(column, original_column, degree) {
    process <- "Coordinate"
    data_type <- ifelse(is.continuous(original_column), "Continuous", "Category")

    if (data_type == "Continuous") {
      if (degree >= (length(unique(original_column)) - 1)) {
        process <- c(process, "Smoothing")
        column <- get_resampling_methods(
          column,
          original_column,
          "smoothing"
        )
      }

      if (abs(get_skewness(original_column)) >= 1) {
        process <- c(process, "Deskewing")
        column <- get_resampling_methods(
          column,
          original_column,
          "deskewing"
        )
      }
    } else {
      if (any(table(original_column) < length(original_column) / (2 * length(table(original_column))))) {
        process <- c(process, "Flattening")
        column <- get_resampling_methods(
          column,
          original_column,
          "flattening"
        )
      }
    }

    return(list(
      column = column,
      data_type = data_type,
      process = paste(process, collapse = "->")
    ))
  }

  # Resample the data by columns
  resampled_df <- as.data.frame(lapply(data, resample_column, size = resample_size))

  # Return resampled data only if resample_only is TRUE
  if (resample_only) {
    return(list(
      resampled_df = resampled_df,
      resampled_df_log = NULL
    ))
  }

  if (is.null(data_degree)) data_degree <- rep(1, (ncol(data)))
  # Record resampling processes for each column
  resampled_df_log <- do.call(
    rbind,
    lapply(
      1:ncol(resampled_df),
      function(n) {
        processed <- process_column(resampled_df[, n], data[, n], data_degree[n])
        resampled_df[, n] <<- processed$column
        data.frame(
          Covariate = colnames(resampled_df)[n],
          Type = processed$data_type,
          Process = processed$process
        )
      }
    )
  )

  return(list(
    resampled_df = resampled_df,
    resampled_df_log = resampled_df_log
  ))
}

#' Convert Formula to String
#'
#' This function converts a formula object to a character string. It removes extra whitespace and formats the formula as a single line.
#'
#' @param formula A formula object to be converted to a string.
#'
#' @return A character string representing the formula.
get_formula_string <- function(formula) {
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = " "))
  return(char)
}

#' Extract Left-Hand Side of Formula as String
#'
#' This function extracts the left-hand side (LHS) of a formula object and converts it to a character string. It uses `get_formula_string` to ensure consistent formatting.
#'
#' @param formula A formula object from which the LHS will be extracted.
#'
#' @return A character string representing the left-hand side of the formula.
get_formula_lhs <- function(formula) {
  return(get_formula_string(rlang::f_lhs(formula)))
}

#' Extract the Right-Hand Side of a Formula
#'
#' This function extracts the right-hand side (RHS) of a formula and returns it as a character string.
#' Optionally, it can include a tilde (`~`) at the beginning of the RHS.
#'
#' @param formula A formula object from which to extract the RHS.
#' @param with_tilde Logical, indicating whether to include a tilde (`~`) at the beginning of the RHS.
#'   Defaults to `FALSE`.
#'
#' @return A character string representing the right-hand side of the formula. If `with_tilde` is `TRUE`,
#'   the string includes a leading tilde.
get_formula_rhs <- function(formula, with_tilde = FALSE) {
  if (with_tilde) {
    return(paste("~", deparse(rlang::f_rhs(formula)), collapse = " "))
  }
  return(get_formula_string(rlang::f_rhs(formula)))
}

#' Check if a Variable is Continuous
#'
#' This function checks whether a given vector represents a continuous variable. A continuous variable is numeric and has more than two unique values.
#'
#' @param lst A vector to be checked.
#'
#' @return A logical value indicating whether the input vector is considered continuous. Returns `TRUE` if the vector is numeric and has more than two unique values; otherwise, returns `FALSE`.
is.continuous <- function(lst) {
  # Check if the input 'lst' is numeric and has more than 2 unique values
  return(is.numeric(lst) & length(unique(lst)) > 2)
}

#' Standardize Data
#'
#' This function standardizes a dataset by converting columns to numeric or factor types and replacing NA values.
#' For continuous variables, NA values are replaced with either a specific numeric value or a computed statistic.
#' For categorical variables, NA values are replaced with the mode of the column.
#'
#' @param data A data frame to be standardized.
#' @param na_replace A function or numeric value used to replace NA values. If a function, it should take a vector and return a replacement value.
#'   If a numeric value, it is used directly to replace NA values in continuous columns.
#'   The default is `stats::na.omit`, which omits rows with NA values (used as an indicator here, not the actual replacement value).
#'
#' @return A data frame where columns have been converted to numeric or factor types, and NA values have been replaced according to the method specified.
get_standardized_data <- function(
    data,
    na_replace = stats::na.omit) {
  # Function to replace NA values based on data type
  replace_na <- function(x) {
    if (anyNA(x)) {
      replace_value <- if (is.continuous(stats::na.omit(x))) {
        if (is.numeric(na_replace)) na_replace else na_replace(x, na.rm = TRUE)
      } else {
        unique_vals <- unique(x)
        unique_vals[which.max(tabulate(match(x, unique_vals)))]
      }

      x[is.na(x)] <- replace_value
    }
    x
  }

  # Convert columns and replace NA values
  data <- as.data.frame(if (identical(na_replace, stats::na.omit)) {
    stats::na.omit(data)
  } else {
    lapply(data, replace_na)
  })

  return(data)
}

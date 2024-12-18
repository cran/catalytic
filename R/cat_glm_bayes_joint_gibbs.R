#' Bayesian Estimation with Gibbs Sampling for Catalytic Generalized Linear Models (GLMs) Binomial Family for Coefficients and tau
#'
#' This function uses Gibbs sampling to estimate a Bayesian GLMs Binomial Family, where both the coefficients
#' and tau parameter are jointly sampled. tau is updated via a gamma distribution, while coefficients are
#' updated using Hamiltonian Monte Carlo (HMC) sampling. The model allows for progress updates,
#' warm-up iterations, and initial coefficient estimation based on initial tau value.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables (e.g. \code{~.}).
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param iter Integer; the number of Gibbs sampling iterations (default = 1000).
#' @param warmup Integer; the number of initial iterations for warm-up (default = 500).
#' @param coefs_iter Integer; the number of iterations for the HMC step to update coefficients.
#' @param tau_0 Initial value for tau; defaults to the number of predictors / 4 if NULL.
#' @param tau_alpha Shape parameter for the gamma distribution when updating tau. Default is 2.
#' @param tau_gamma Scale parameter for the gamma distribution when updating tau. Default is 1.
#' @param refresh Logical; if TRUE, displays sampling progress. Default is TRUE.
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{gibbs_iteration_log}{Matrix containing the coefficients and tau values from each Gibbs iteration.}
#' \item{inform_df}{Summary statistics of each parameter, including mean, standard error, quantiles, and effective sample size.}
#' \item{tau}{Mean of sampled tau values.}
#' \item{coefficients}{Mean of sampled coefficient values.}
#'
#' @examples
#' \donttest{
#' binomial_data <- data.frame(
#'   X1 = stats::rnorm(10),
#'   X2 = stats::rnorm(10),
#'   Y = stats::rbinom(10, 1, 0.5)
#' )
#'
#' cat_init <- cat_glm_initialization(
#'   formula = Y ~ 1, # formula for simple model
#'   data = binomial_data,
#'   family = binomial,
#'   syn_size = 100, # Synthetic data size
#'   custom_variance = NULL, # User customized variance value
#'   gaussian_known_variance = FALSE, # Indicating whether the data variance is unknown
#'   x_degree = c(1, 1), # Degrees for polynomial expansion of predictors
#'   resample_only = FALSE, # Whether to perform resampling only
#'   na_replace = stats::na.omit # How to handle NA values in data
#' )
#'
#' cat_model <- cat_glm_bayes_joint_gibbs(
#'   formula = ~.,
#'   cat_init = cat_init, # Only accept object generated from `cat_glm_initialization`
#'   iter = 10, # Number of Gibbs sampling iterations
#'   warmup = 5, # Number of warm-up (or burn-in) iterations for initial iterations
#'   coefs_iter = 2, # Number of iterations for the HMC step to update coefficients
#'   tau_alpha = 1, # Shape parameter for the gamma distribution when updating tau
#'   tau_gamma = 2, # Scale parameter for the gamma distribution when updating tau
#'   refresh = TRUE # Indicator for displaying sampling progress
#' )
#' cat_model
#' }
#' @export
cat_glm_bayes_joint_gibbs <- function(formula,
                                      cat_init,
                                      iter = 1000,
                                      warmup = 500,
                                      coefs_iter = 5,
                                      tau_0 = NULL,
                                      tau_alpha = 2,
                                      tau_gamma = 1,
                                      refresh = TRUE) {
  # Validate Input Parameters
  validate_glm_input(
    formula = formula,
    cat_init = cat_init,
    gibbs_iter = iter,
    gibbs_warmup = warmup,
    coefs_iter = coefs_iter,
    tau_0 = tau_0,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma,
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  full_formula <- stats::as.formula(paste0(cat_init$y_col_name, f_rhs))

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  tau <- ifelse(is.null(tau_0), ncol(cat_init$adj_x) / 4, tau_0)

  ## Suppress warning for `binomial` family when weights contains value < 1
  suppressWarnings(
    kappa_ <- get_glm_log_density(
      family_string = cat_init$family,
      X = cat_init$adj_syn_x,
      Y = cat_init$syn_y,
      coefs = stats::coef(do.call(
        stats::glm,
        list(
          formula = full_formula,
          family = cat_init$family,
          data = cat_init$syn_data
        )
      )),
      weights = 1 / cat_init$syn_size
    )
  )

  # Fit GLM to obtain initial coefficients
  ## Suppress warning for `binomial` family when weights contains value < 1
  suppressWarnings(
    coefs <- stats::coef(do.call(
      stats::glm,
      list(
        formula = full_formula,
        family = cat_init$family,
        data = cat_init$data,
        weights = c(
          rep(1, cat_init$obs_size),
          rep(tau / cat_init$syn_size, cat_init$syn_size)
        )
      )
    ))
  )

  # Initialize a matrix to collect samples
  gibbs_iteration_log <- matrix(
    nrow = iter,
    ncol = (1 + length(coefs))
  )
  colnames(gibbs_iteration_log) <- c("(Intercept)", colnames(cat_init$adj_x), "tau")
  gibbs_iteration_log[1, ] <- c(coefs, tau)

  if (refresh) {
    cat(
      "\nSAMPLING FROM GIBB'S NOW",
      "\nGibb's sampling:",
      "\nGibb's sampling: \n"
    )
  }

  # Define print boundaries for progress reporting
  print_boundary <- if (iter >= 10) c(1, round(seq(0.1, 1, 0.1) * iter)) else seq(1, iter)

  # Initialize flags and timers for tracking progress
  sampling_start_flag <- FALSE
  warmup_start_time <- Sys.time()

  for (gibbs_i in 2:iter) {
    # Update tau by using the gamma distribution
    tau <- stats::rgamma(
      n = 1,
      shape = ncol(cat_init$adj_x) + tau_alpha,
      scale = kappa_ + 1 / tau_gamma - get_glm_log_density(
        family_string = cat_init$family,
        X = cat_init$adj_syn_x,
        Y = cat_init$syn_y,
        coefs = coefs,
        weights = 1 / cat_init$syn_size
      )
    )

    ## Suppress warning for `binomial` family when weights contains value < 1
    tau_model <- suppressWarnings(
      do.call(
        stats::glm,
        list(
          formula = full_formula,
          family = cat_init$family,
          data = cat_init$data,
          weights = c(
            rep(1, cat_init$obs_size),
            rep(tau / cat_init$syn_size, cat_init$syn_size)
          )
        )
      )
    )

    # Compute diagonal approximation for HMC
    hmc_scale <- sqrt(get_glm_diag_approx_cov(X = cat_init$adj_x, model = tau_model))

    tau_weights <- c(
      rep(1, cat_init$obs_size),
      rep(tau / cat_init$syn_size, cat_init$syn_size)
    )

    # Update coefs by performing HMC sampling
    coefs <- get_hmc_mcmc_result(
      neg_log_den_func = function(coefs) {
        return(
          -get_glm_log_density(
            family_string = cat_init$family,
            X = cat_init$adj_x,
            Y = cat_init$y,
            coefs = coefs,
            weights = tau_weights
          )
        )
      },
      neg_log_den_grad_func = function(coefs) {
        return(
          -get_glm_log_density_grad(
            family_string = cat_init$family,
            X = cat_init$adj_x,
            Y = cat_init$y,
            coefs = coefs,
            weights = tau_weights
          )
        )
      },
      coefs_0 = coefs,
      iter = coefs_iter,
      hmc_scale = hmc_scale
    )

    gibbs_iteration_log[gibbs_i, ] <- c(coefs, tau)

    # Print progress updates
    if (((gibbs_i %in% print_boundary) | gibbs_i == warmup) & refresh) {
      print_str <- paste("Gibb's sampling:", gibbs_i, " / ", iter, "[", round(gibbs_i / iter * 100), "%]")

      if (gibbs_i < warmup) {
        cat(print_str, " (Warmup)\n")
      } else {
        if (!sampling_start_flag) {
          sampling_start_time <- Sys.time()
          sampling_start_flag <- TRUE
        }
        cat(print_str, " (Sampling)\n")
      }
    }
  }

  sampling_end_time <- Sys.time()

  # Print elapsed times for warmup and sampling
  if (refresh) {
    cat(paste0(
      "Gibb's sampling: \n",
      "Gibb's sampling:",
      "  Elapsed Time: ",
      round(difftime(
        sampling_start_time,
        warmup_start_time,
        units = "secs"
      ), 3),
      " seconds ",
      "(Warm-up)\n",
      "Gibb's sampling:                ",
      round(difftime(
        sampling_end_time,
        sampling_start_time,
        units = "secs"
      ), 3),
      " seconds ",
      "(Sampling)\n",
      "Gibb's sampling:                ",
      round(difftime(
        sampling_end_time,
        warmup_start_time,
        units = "secs"
      ), 3),
      " seconds ",
      "(Total)\n",
      "Gibb's sampling: \n"
    ))
  }

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_glm_bayes_joint_gibbs",
    ## Input/Processed parameters
    formula = formula,
    cat_init = cat_init,
    iter = iter,
    warmup = warmup,
    coefs_iter = coefs_iter,
    tau_0 = tau_0,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma,
    refresh = refresh,
    sys_time = sampling_end_time,
    ## Result
    gibbs_iteration_log = gibbs_iteration_log,
    inform_df = tryCatch(
      {
        t(apply(
          gibbs_iteration_log,
          2,
          function(x) {
            c(
              mean = mean(x),
              se_mean = stats::sd(x) / sqrt(length(x)),
              sd = stats::sd(x),
              stats::quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
              n_eff = as.double(coda::effectiveSize(coda::as.mcmc(x)))
            )
          }
        ))
      },
      error = function(.) {
        t(sapply(
          gibbs_iteration_log,
          function(x) {
            c(
              mean = mean(x), se_mean = NA, sd = NA,
              stats::quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), n_eff = 0
            )
          }
        ))
      }
    )
  )

  cat_model$tau <- cat_model$inform_df[(length(coefs) + 1), "mean"]
  cat_model$coefficients <- cat_model$inform_df[1:length(coefs), "mean"]

  class(cat_model) <- c(cat_model$class, "cat_gibbs", "cat_glm")

  return(cat_model)
}

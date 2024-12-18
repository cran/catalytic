#' Bayesian Estimation for Catalytic Generalized Linear Models (GLMs)  with adaptive tau
#'
#' This function performs Bayesian estimation for a catalytic Generalized Linear Models (GLMs) using RStan
#' by using adaptive tau. It supports both Gaussian and Binomial family models, enabling
#' flexibility in prior specifications and algorithm configurations.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables (e.g. \code{~.}).
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param chains Number of Markov chains to run. Default is 4.
#' @param iter Total number of iterations per chain. Default is 2000.
#' @param warmup Number of warm-up iterations per chain (discarded from final analysis). Default is 1000.
#' @param algorithm The sampling algorithm to use in \code{rstan::sampling}. Default is "NUTS"
#' (No-U-Turn Sampler).
#' @param tau_alpha Shape parameter of the prior for tau. Default is 2.
#' @param tau_gamma Scale parameter of the prior for tau. Default is 1.
#' @param binomial_tau_lower A numeric lower bound for tau in Binomial models. Default is 0.05.
#' @param binomial_joint_theta Logical; if TRUE, uses joint theta (theta = 1/tau) in Binomial models. Default is FALSE.
#' @param binomial_joint_alpha Logical; if TRUE, uses joint alpha (adaptive tau_alpha)
#'   in Binomial models. Default is FALSE. To activate this feature, both
#'   `binomial_joint_theta = TRUE` and `binomial_joint_alpha = TRUE` must be set.
#' @param gaussian_variance_alpha The shape parameter for the inverse-gamma prior on
#' variance if the variance is unknown in Gaussian models. Defaults to the number of predictors.
#' @param gaussian_variance_beta The scale parameter for the inverse-gamma prior on
#' variance if the variance is unknown in Gaussian models. Defaults to the number of predictors
#' times variance of observation response.
#' @param ... Additional parameters to pass to \code{rstan::sampling}.
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{stan_data}{A data list used for fitting RStan sampling model.}
#' \item{stan_model}{Compiled RStan model object for GLMs.}
#' \item{stan_sample_model}{Fitted RStan sampling model containing posterior samples.}
#' \item{coefficients}{Mean posterior estimates of model coefficients from `stan_sample_model`.}
#' \item{tau}{Mean posterior of tau (or transformed theta if applicable).}
#'
#' @examples
#' \donttest{
#' gaussian_data <- data.frame(
#'   X1 = stats::rnorm(10),
#'   X2 = stats::rnorm(10),
#'   Y = stats::rnorm(10)
#' )
#'
#' cat_init <- cat_glm_initialization(
#'   formula = Y ~ 1, # formula for simple model
#'   data = gaussian_data,
#'   syn_size = 100, # Synthetic data size
#'   custom_variance = NULL, # User customized variance value
#'   gaussian_known_variance = FALSE, # Indicating whether the data variance is unknown
#'   x_degree = c(1, 1), # Degrees for polynomial expansion of predictors
#'   resample_only = FALSE, # Whether to perform resampling only
#'   na_replace = stats::na.omit # How to handle NA values in data
#' )
#'
#' cat_model <- cat_glm_bayes_joint(
#'   formula = ~.,
#'   cat_init = cat_init, # Only accept object generated from `cat_glm_initialization`
#'   chains = 1, # Number of Markov chains to be run in the RStan sampling
#'   iter = 10, # Number of iterations per chain in the RStan sampling
#'   warmup = 5, # Number of warm-up (or burn-in) iterations for each chain
#'   algorithm = "NUTS", # Sampling algorithm to use in \code{rstan::sampling}
#'   tau_alpha = 1, # Shape parameter of the prior for tau
#'   tau_gamma = 2, # Scale parameter of the prior for tau
#'   binomial_tau_lower = 0.05, # Lower bound for tau in Binomial models.
#'   binomial_joint_theta = FALSE, # Indicator for using joint theta for Binomial models
#'   binomial_joint_alpha = FALSE, # Indicator for using oint alpha for Binomial models
#'   gaussian_variance_alpha = 1, # The shape parameter for the inverse-gamma prior for variance
#'   gaussian_variance_beta = 2 # The scale parameter for the inverse-gamma prior for variance
#' )
#' cat_model
#' }
#' @export
cat_glm_bayes_joint <- function(
    formula,
    cat_init,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    algorithm = "NUTS",
    tau_alpha = 2,
    tau_gamma = 1,
    binomial_tau_lower = 0.05,
    binomial_joint_theta = FALSE,
    binomial_joint_alpha = FALSE,
    gaussian_variance_alpha = NULL,
    gaussian_variance_beta = NULL,
    ...) {
  # Validate Input Parameters
  validate_glm_input(
    formula = formula,
    cat_init = cat_init,
    binomial_tau_lower = binomial_tau_lower,
    binomial_joint_theta = binomial_joint_theta,
    binomial_joint_alpha = binomial_joint_alpha,
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma,
    gaussian_variance_alpha = gaussian_variance_alpha,
    gaussian_variance_beta = gaussian_variance_beta
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  full_formula <- stats::as.formula(paste0(cat_init$y_col_name, f_rhs))

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  stan_data <- list(
    obs_n = cat_init$obs_size,
    syn_n = cat_init$syn_size,
    col_n = ncol(cat_init$adj_x),
    obs_x = cat_init$adj_obs_x,
    syn_x = cat_init$adj_syn_x,
    obs_y = cat_init$obs_y,
    syn_y = cat_init$syn_y,
    kappa = get_glm_log_density(
      family_string = cat_init$family,
      X = cat_init$adj_syn_x,
      Y = cat_init$syn_y,
      coefs = stats::coef(
        do.call(
          stats::glm,
          list(
            formula = full_formula,
            family = cat_init$family,
            data = cat_init$syn_data
          )
        )
      ),
      weights = 1 / cat_init$syn_size
    ),
    tau_alpha = tau_alpha,
    tau_gamma = tau_gamma
  )

  if (cat_init$family == "gaussian") {
    binomial_joint_theta <- FALSE
    binomial_joint_alpha <- FALSE
    # Include Gaussian-specific parameters if the family is Gaussian
    if (cat_init$gaussian_known_variance) {
      # Include sigma value if user set the variance is known
      add_on_list <- list(
        sigma = sqrt(
          ifelse(
            is.null(cat_init$custom_variance),
            get_glm_custom_var(
              formula = full_formula,
              cat_init = cat_init,
              tau = ncol(cat_init$adj_x) # Specifying tau as the number of predictors for variance adjustment
            ),
            cat_init$custom_variance
          )
        )
      )
      gaussian_variance_alpha <- NULL
      gaussian_variance_beta <- NULL
    } else {
      if (is.null(gaussian_variance_alpha)) gaussian_variance_alpha <- ncol(cat_init$adj_x)
      if (is.null(gaussian_variance_beta)) gaussian_variance_beta <- ncol(cat_init$adj_x) * stats::var(cat_init$obs_y)

      # Include the alpha and beta for the inverse-gamma prior for unknown variance
      add_on_list <- list(
        variance_alpha = gaussian_variance_alpha,
        variance_beta = gaussian_variance_beta
      )
    }
  } else {
    gaussian_variance_alpha <- NULL
    gaussian_variance_beta <- NULL

    # Include Binomial-specific parameters
    add_on_list <- list(tau_lower = binomial_tau_lower)
  }

  stan_data <- c(stan_data, add_on_list)

  stan_model <- get_stan_model(
    type = "glm",
    glm_family_string = cat_init$family,
    joint_tau = TRUE,
    glm_binomial_joint_theta = binomial_joint_theta,
    glm_binomial_joint_alpha = binomial_joint_alpha,
    glm_gaussian_known_variance = cat_init$gaussian_known_variance
  )

  # Perform Bayesian Sampling via RStan
  stan_sample_model <- rstan::sampling(
    object = stan_model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    algorithm = algorithm,
    ...
  )

  ## Extract the mean of coefs from the stan sampling model
  coefs <- rstan::summary(stan_sample_model)$summary[
    grep("coefs", rownames(rstan::summary(stan_sample_model)$summary)),
    "mean"
  ]

  if (!is.null(coefs)) {
    names(coefs) <- c("(Intercept)", colnames(cat_init$adj_x))
  }

  if ((cat_init$family == "binomial") && (nrow(stan_sample_model) != 0)) {
    print_glm_bayes_joint_binomial_suggestion(
      alpha = tau_alpha,
      stan_iter = iter,
      stan_sample_model = stan_sample_model,
      binomial_joint_theta = binomial_joint_theta,
      binomial_joint_alpha = binomial_joint_alpha,
      binomial_tau_lower = binomial_tau_lower
    )
  }

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_glm_bayes_joint",
    ## Input/Processed parameters
    formula = full_formula,
    cat_init = cat_init,
    tau_alpha = ifelse(
      (binomial_joint_alpha & binomial_joint_theta),
      rstan::summary(stan_sample_model)$summary[
        grep("alpha", rownames(rstan::summary(stan_sample_model)$summary)),
        "mean"
      ],
      tau_alpha
    ),
    tau_gamma = tau_gamma,
    binomial_tau_lower = binomial_tau_lower,
    binomial_joint_theta = binomial_joint_theta,
    binomial_joint_alpha = binomial_joint_alpha,
    chains = chains,
    iter = iter,
    warmup = warmup,
    algorithm = algorithm,
    gaussian_variance_alpha = gaussian_variance_alpha,
    gaussian_variance_beta = gaussian_variance_beta,
    ...,
    ## Result
    stan_data = stan_data,
    stan_model = stan_model,
    stan_sample_model = stan_sample_model,
    coefficients = coefs,
    tau = if (binomial_joint_theta) {
      1 / rstan::summary(stan_sample_model)$summary[
        "theta", "mean"
      ]
    } else {
      rstan::summary(stan_sample_model)$summary[
        "tau", "mean"
      ]
    }
  )

  class(cat_model) <- c(cat_model$class, "cat_bayes", "cat_glm")

  return(cat_model)
}

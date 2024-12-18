#' Bayesian Estimation for Catalytic Generalized Linear Models (GLMs) with Fixed tau
#'
#' Fits a Bayesian generalized linear model using synthetic and observed data based on an initial
#' data structure, formula, and other model specifications. Supports only Gaussian and Binomial
#' distributions in the GLM family.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables (e.g. \code{~.}).
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param tau Optional numeric scalar controlling the weight of the synthetic data in the coefficient estimation.
#'  Defaults to the number of predictors / 4 for Gaussian models or the number of predictors otherwise.
#' @param chains Number of Markov chains to run. Default is 4.
#' @param iter Total number of iterations per chain. Default is 2000.
#' @param warmup Number of warm-up iterations per chain (discarded from final analysis). Default is 1000.
#' @param algorithm The sampling algorithm to use in \code{rstan::sampling}. Default is "NUTS"
#' (No-U-Turn Sampler).
#' @param gaussian_variance_alpha The shape parameter for the inverse-gamma prior on
#' variance if the variance is unknown in Gaussian models. Defaults to the number of predictors.
#' @param gaussian_variance_beta The scale parameter for the inverse-gamma prior on
#' variance if the variance is unknown in Gaussian models. Defaults to the number of predictors
#' times variance of observation response.
#' @param ... Additional parameters to pass to \code{rstan::sampling}.
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{stan_data}{The data list used for fitting RStan sampling model.}
#' \item{stan_model}{Compiled RStan model object for GLMs.}
#' \item{stan_sample_model}{Fitted RStan sampling model containing posterior samples.}
#' \item{coefficients}{Mean posterior estimates of model coefficients from `stan_sample_model`.}
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
#' cat_model <- cat_glm_bayes(
#'   formula = ~.,
#'   cat_init = cat_init, # Only accept object generated from `cat_glm_initialization`
#'   tau = 1, # Weight for synthetic data
#'   chains = 1, # Number of Markov chains to be run in the RStan sampling
#'   iter = 10, # Number of iterations per chain in the RStan sampling
#'   warmup = 5, # Number of warm-up (or burn-in) iterations for each chain
#'   algorithm = "NUTS", # Sampling algorithm to use in \code{rstan::sampling}
#'   gaussian_variance_alpha = 1, # The shape parameter for the inverse-gamma prior for variance
#'   gaussian_variance_beta = 2 # The scale parameter for the inverse-gamma prior for variance
#' )
#' cat_model
#' }
#' @export
cat_glm_bayes <- function(
    formula,
    cat_init,
    tau = NULL,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    algorithm = "NUTS",
    gaussian_variance_alpha = NULL,
    gaussian_variance_beta = NULL,
    ...) {
  # Validate Input Parameters
  validate_glm_input(
    formula = formula,
    cat_init = cat_init,
    tau = tau,
    gaussian_variance_alpha = gaussian_variance_alpha,
    gaussian_variance_beta = gaussian_variance_beta,
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  full_formula <- stats::as.formula(paste0(cat_init$y_col_name, f_rhs))

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  if (is.null(tau)) {
    tau <- ifelse(
      cat_init$family == "gaussian",
      ncol(cat_init$adj_x) / 4,
      ncol(cat_init$adj_x)
    )
  }

  stan_data <- list(
    obs_n = cat_init$obs_size,
    syn_n = cat_init$syn_size,
    col_n = ncol(cat_init$adj_x),
    obs_x = cat_init$adj_obs_x,
    syn_x = cat_init$adj_syn_x,
    obs_y = cat_init$obs_y,
    syn_y = cat_init$syn_y,
    tau = tau
  )

  if (cat_init$family == "gaussian") {
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
    add_on_list <- NULL
  }


  stan_data <- c(stan_data, add_on_list)

  stan_model <- get_stan_model(
    type = "glm",
    glm_family_string = cat_init$family,
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

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_glm_bayes",
    ## Input/Processed parameters
    formula = full_formula,
    cat_init = cat_init,
    tau = tau,
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
    coefficients = coefs
  )

  class(cat_model) <- c(cat_model$class, "cat_bayes", "cat_glm")

  return(cat_model)
}

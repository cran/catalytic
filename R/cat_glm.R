#' Catalytic Generalized Linear Models (GLMs) Fitting Function with Fixed Tau
#'
#' Fits a Catalytic Generalized Linear Models (GLMs) by using observed and synthetic data.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables (e.g. \code{~ .}).
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param tau Optional numeric scalar controlling the weight of the synthetic data in the coefficient estimation.
#'  Defaults to the number of predictors / 4 for Gaussian models or the number of predictors otherwise.
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{coefficients}{Estimated coefficient vector.}
#' \item{model}{Fitted GLMs object (`stats::glm`).}
#'
#' @examples
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
#'   gaussian_known_variance = TRUE, # Indicating whether the data variance is known
#'   x_degree = c(1, 1), # Degrees for polynomial expansion of predictors
#'   resample_only = FALSE, # Whether to perform resampling only
#'   na_replace = stats::na.omit # How to handle NA values in data
#' )
#'
#' cat_model <- cat_glm(
#'   formula = ~.,
#'   cat_init = cat_init, # Only accept object generated from `cat_glm_initialization`
#'   tau = 1 # Weight for synthetic data
#' )
#' cat_model
#' @export
cat_glm <- function(formula,
                    cat_init,
                    tau = NULL) {
  # Validate Input Parameters
  validate_glm_input(
    formula = formula,
    cat_init = cat_init,
    tau = tau
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

  if (tau == 0 & (dim(cat_init$adj_obs_x)[1] < dim(cat_init$adj_obs_x)[2])) {
    warning(
      paste0(
        "The number of columns in the observation data exceeds its data size. \n",
        "To avoid issues, the value of `tau` will be assigned with `0.01` instead of `0`."
      ),
      call. = FALSE
    )
    tau <- 0.01
  }

  # Perform Catalytic Fitting
  ## Suppress warning for `binomial` family when weights contains value < 1
  glm_model <- suppressWarnings(do.call(
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


  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_glm",
    ## Input/Processed parameters
    formula = full_formula,
    cat_init = cat_init,
    tau = tau,
    ## Result
    model = glm_model,
    coefficients = stats::coef(glm_model)
  )

  class(cat_model) <- c(cat_model$class, "cat", "cat_glm")

  return(cat_model)
}

#' Initialization for Catalytic Generalized Linear Models (GLMs)
#'
#' This function prepares and initializes a catalytic Generalized Linear Models (GLMs) by processing input data,
#' extracting necessary variables, generating synthetic datasets, and fitting a model.
#'
#' @param formula A formula specifying the GLMs. Should include response and predictor variables.
#' @param family The type of GLM family. Defaults to Gaussian.
#' @param data A data frame containing the data for modeling.
#' @param syn_size An integer specifying the size of the synthetic dataset to be generated. Default is four times the number of predictor columns.
#' @param custom_variance A custom variance value to be applied if using a Gaussian model. Defaults to `NULL`.
#' @param gaussian_known_variance A logical value indicating whether the data variance is known. Defaults to `FALSE`. Only applicable to Gaussian family.
#' @param x_degree A numeric vector indicating the degree for polynomial expansion of predictors. Default is 1 for each predictor.
#' @param resample_only A logical indicating whether to perform resampling only. Default is FALSE.
#' @param na_replace A function to handle NA values in the data. Default is `stats::na.omit`.
#'
#' @return A list containing the values of all the input arguments and the following components:
#'
#' - **Function Information**
#'   - `function_name`: The name of the function, "cat_glm_initialization".
#'   - `y_col_name`: The name of the response variable in the dataset.
#'   - `simple_model`: An object of class `stats::glm`, representing the fitted model for generating synthetic response from the original data.
#'
#' - **Observation Data Information**
#'   - `obs_size`: Number of observations in the original dataset.
#'   - `obs_data`: Data frame of standardized observation data.
#'   - `obs_x`: Predictor variables for observed data.
#'   - `obs_y`: Response variable for observed data.
#'
#' - **Synthetic Data Information**
#'   - `syn_size`: Number of synthetic observations generated.
#'   - `syn_data`: Data frame of synthetic predictor and response variables.
#'   - `syn_x`: Synthetic predictor variables.
#'   - `syn_y`: Synthetic response variable.
#'   - `syn_x_resample_inform`: Information about resampling methods for synthetic predictors:
#'     - Coordinate: Preserves the original data values as reference coordinates during processing.
#'     - Deskewing: Adjusts the data distribution to reduce skewness and enhance symmetry.
#'     - Smoothing: Reduces noise in the data to stabilize the dataset and prevent overfitting.
#'     - Flattening: Creates a more uniform distribution by modifying low-frequency categories in categorical variables.
#'     - Symmetrizing: Balances the data around its mean to improve statistical properties for model fitting.
#'
#' - **Whole Data Information**
#'   - `size`: Total number of combined original and synthetic observations.
#'   - `data`: Data frame combining original and synthetic datasets.
#'   - `x`: Combined predictor variables from original and synthetic data.
#'   - `y`: Combined response variable from original and synthetic data.
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
#' cat_init
#' @export
cat_glm_initialization <- function(formula,
                                   family = "gaussian",
                                   data,
                                   syn_size = NULL,
                                   custom_variance = NULL,
                                   gaussian_known_variance = FALSE,
                                   x_degree = NULL,
                                   resample_only = FALSE,
                                   na_replace = stats::na.omit) {
  # Validate Input Parameters
  validate_glm_initialization_input(
    formula = formula,
    family = family,
    data = data,
    syn_size = syn_size,
    custom_variance = custom_variance,
    gaussian_known_variance = gaussian_known_variance,
    x_degree = x_degree
  )

  # Extract Information from Input Parameters
  ## Standardize observation covariate and response
  obs_data <- get_standardized_data(
    data = data,
    na_replace = na_replace
  )

  ## Extract observation covariate and response
  y_col_name <- get_formula_lhs(formula)
  obs_x <- as.data.frame(obs_data[, c(setdiff(names(obs_data), c(y_col_name))), drop = FALSE])
  obs_y <- obs_data[, y_col_name]

  # Generate Synthetic Dataset
  if (is.null(x_degree)) x_degree <- rep(1, (ncol(obs_x)))
  if (is.null(syn_size)) syn_size <- ncol(obs_x) * 4

  ## Generate synthetic covariate by resampling observation covariate
  syn_x_inform <- get_resampled_df(
    data = obs_x,
    resample_size = syn_size,
    data_degree = x_degree,
    resample_only = resample_only
  )
  syn_x <- as.data.frame(syn_x_inform$resampled_df)

  ## Generate synthetic response by using simple model
  simple_model <- stats::glm(
    formula = formula,
    family = family,
    data = obs_data
  )

  family <- get_glm_family_string(family)

  syn_y <- as.numeric(
    if (family == "binomial") {
      stats::rbinom(syn_size, 1, stats::predict(simple_model, syn_x, type = "response"))
    } else {
      stats::predict(simple_model, syn_x)
    }
  )

  # Finalize Setup and Output
  cat_init <- list(
    function_name = "cat_glm_initialization",
    ## Input/Processed parameters
    formula = formula,
    family = family,
    syn_size = syn_size,
    custom_variance = if ((family == "gaussian") & gaussian_known_variance) custom_variance else NULL,
    gaussian_known_variance = if (family == "gaussian") gaussian_known_variance else FALSE,
    x_degree = x_degree,
    resample_only = resample_only,
    na_replace = na_replace,
    y_col_name = y_col_name,
    simple_model = simple_model,
    ## Observation data information
    obs_size = nrow(obs_data),
    obs_data = obs_data,
    obs_x = obs_x,
    obs_y = obs_y,
    ## Synthetic data information
    syn_size = syn_size,
    syn_data = cbind(syn_x, stats::setNames(as.data.frame(syn_y), y_col_name)),
    syn_x = syn_x,
    syn_y = syn_y,
    syn_x_resample_inform = syn_x_inform$resampled_df_log
  )

  cat_init <- c(
    cat_init,
    list(
      ## Whole data information
      size = cat_init$obs_size + cat_init$syn_size,
      data = rbind(cat_init$obs_data, cat_init$syn_data),
      x = rbind(cat_init$obs_x, cat_init$syn_x),
      y = c(cat_init$obs_y, cat_init$syn_y)
    )
  )

  class(cat_init) <- c(cat_init$class, "cat_initialization")

  return(cat_init)
}

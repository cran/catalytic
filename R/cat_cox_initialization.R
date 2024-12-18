#' Initialization for Catalytic Cox proportional hazards model (COX)
#'
#' This function prepares and initializes a catalytic Cox proportional hazards model by processing input data,
#' extracting necessary variables, generating synthetic datasets, and fitting a model.
#'
#' @param formula A formula specifying the Cox model. Should include response and predictor variables.
#' @param data A data frame containing the data for modeling.
#' @param syn_size An integer specifying the size of the synthetic dataset to be generated. Default is four times the number of predictor columns.
#' @param hazard_constant A constant hazard rate for generating synthetic time data if not using a fitted Cox model. Default is NULL and will calculate in function.
#' @param entry_points A numeric vector for entry points of each observation. Default is NULL.
#' @param x_degree A numeric vector indicating the degree for polynomial expansion of predictors. Default is 1 for each predictor.
#' @param resample_only A logical indicating whether to perform resampling only. Default is FALSE.
#' @param na_replace A function to handle NA values in the data. Default is `stats::na.omit`.
#'
#' @return A list containing the values of all the input arguments and the following components:
#'
#' - **Function Information**:
#'   - `function_name`: The name of the function, "cat_cox_initialization".
#'   - `time_col_name`: The name of the time variable in the dataset.
#'   - `status_col_name`: The name of the status variable (event indicator) in the dataset.
#'   - `simple_model`: If the formula has no predictors, a constant hazard rate model is used; otherwise, a fitted Cox model object.
#'
#' - **Observation Data Information**:
#'   - `obs_size`: Number of observations in the original dataset.
#'   - `obs_data`: Data frame of standardized observation data.
#'   - `obs_x`: Predictor variables for observed data.
#'   - `obs_time`: Observed survival times.
#'   - `obs_status`: Event indicator for observed data.
#'
#' - **Synthetic Data Information**:
#'   - `syn_size`: Number of synthetic observations generated.
#'   - `syn_data`: Data frame of synthetic predictor and response variables.
#'   - `syn_x`: Synthetic predictor variables.
#'   - `syn_time`: Synthetic survival times.
#'   - `syn_status`: Event indicator for synthetic data (defaults to 1).
#'   - `syn_x_resample_inform`: Information about resampling methods for synthetic predictors:
#'     - Coordinate: Preserves the original data values as reference coordinates during processing.
#'     - Deskewing: Adjusts the data distribution to reduce skewness and enhance symmetry.
#'     - Smoothing: Reduces noise in the data to stabilize the dataset and prevent overfitting.
#'     - Flattening: Creates a more uniform distribution by modifying low-frequency categories in categorical variables.
#'     - Symmetrizing: Balances the data around its mean to improve statistical properties for model fitting.
#'
#' - **Whole Data Information**:
#'   - `size`: Total number of combined original and synthetic observations.
#'   - `data`: Data frame combining original and synthetic datasets.
#'   - `x`: Combined predictor variables from original and synthetic data.
#'   - `time`: Combined survival times from original and synthetic data.
#'   - `status`: Combined event indicators from original and synthetic data.
#'
#' @examples
#' library(survival)
#' data("cancer")
#' cancer$status[cancer$status == 1] <- 0
#' cancer$status[cancer$status == 2] <- 1
#'
#' cat_init <- cat_cox_initialization(
#'   formula = Surv(time, status) ~ 1, # formula for simple model
#'   data = cancer,
#'   syn_size = 100, # Synthetic data size
#'   hazard_constant = NULL, # Hazard rate value
#'   entry_points = rep(0, nrow(cancer)), # Entry points of each observation
#'   x_degree = rep(1, ncol(cancer) - 2), # Degrees for polynomial expansion of predictors
#'   resample_only = FALSE, # Whether to perform resampling only
#'   na_replace = stats::na.omit # How to handle NA values in data
#' )
#' cat_init
#' @export
cat_cox_initialization <- function(formula,
                                   data,
                                   syn_size = NULL,
                                   hazard_constant = NULL,
                                   entry_points = NULL,
                                   x_degree = NULL,
                                   resample_only = FALSE,
                                   na_replace = stats::na.omit) {
  # Validate Input Parameters
  validate_cox_initialization_input(
    formula = formula,
    data = data,
    syn_size = syn_size,
    hazard_constant = hazard_constant,
    entry_points = entry_points,
    x_degree = x_degree
  )

  # Extract Information from Input Parameters
  ## Standardize observation covariate and response
  obs_data <- get_standardized_data(
    data = data,
    na_replace = na_replace
  )
  entry_points <- entry_points[as.numeric(rownames(obs_data))]

  ## Extract observation covariate and response
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

  obs_time <- obs_data[, time_col_name]
  obs_status <- obs_data[, status_col_name]
  obs_x <- obs_data[, setdiff(colnames(data), c(time_col_name, status_col_name)), drop = FALSE]
  obs_size <- nrow(obs_x)

  # Generate Synthetic Dataset
  if (is.null(x_degree)) x_degree <- rep(1, ncol(obs_x))
  if (is.null(syn_size)) syn_size <- 4 * ncol(obs_x)
  if (is.null(entry_points)) entry_points <- rep(0, obs_size)

  ## Generate synthetic covariate by resampling observation covariate
  syn_x_inform <- get_resampled_df(
    data = obs_x,
    resample_size = syn_size,
    data_degree = x_degree,
    resample_only = resample_only
  )
  syn_x <- as.data.frame(syn_x_inform$resampled_df)

  # Generate synthetic time by using simple model
  if (get_formula_rhs(formula) == 1) {
    simple_model <- sum(obs_status) / mean(obs_time) / obs_size
    syn_time <- stats::rexp(
      syn_size,
      rate = simple_model
    )
  } else {
    simple_model <- survival::coxph(
      formula = formula,
      data = data
    )

    syn_time <- summary(
      survival::survfit(
        simple_model,
        newdata = as.data.frame(syn_x)
      )
    )$table[, "median"]
  }

  # Delete synthetic data that syn_time is NA
  syn_time_na_idx <- which(is.na(syn_time))

  if (length(syn_time_na_idx) > 0) {
    warning(paste(
      "Found and removed",
      length(syn_time_na_idx),
      "NA when calculate syn_y."
    ))

    syn_time <- syn_time[-syn_time_na_idx]
    syn_x <- syn_x[-syn_time_na_idx, ]
  }
  syn_status <- rep(1, syn_size - length(syn_time_na_idx))

  cat_init <- list(
    function_name = "cat_cox_initialization",
    ## Input/Processed parameters
    formula = formula,
    hazard_constant = if (is.null(hazard_constant)) sum(obs_status) / mean(obs_time) / obs_size else hazard_constant,
    entry_points = entry_points,
    x_degree = x_degree,
    resample_only = resample_only,
    na_replace = na_replace,
    time_col_name = time_col_name,
    status_col_name = status_col_name,
    simple_model = simple_model,
    # Observation data Information
    obs_size = obs_size,
    obs_data = obs_data,
    obs_x = obs_x,
    obs_time = obs_time,
    obs_status = obs_status,
    # Synthetic data Information
    syn_size = syn_size - length(syn_time_na_idx),
    syn_data = cbind(
      syn_x,
      stats::setNames(as.data.frame(syn_time), time_col_name),
      stats::setNames(as.data.frame(syn_status), status_col_name)
    ),
    syn_x = syn_x,
    syn_time = syn_time,
    syn_status = syn_status,
    syn_x_resample_inform = syn_x_inform$resampled_df_log
  )

  cat_init <- c(
    cat_init,
    list(
      ## Whole data information
      size = cat_init$obs_size + cat_init$syn_size,
      data = rbind(cat_init$obs_data, cat_init$syn_data),
      x = rbind(cat_init$obs_x, cat_init$syn_x),
      time = c(cat_init$obs_time, cat_init$syn_time),
      status = c(cat_init$obs_status, cat_init$syn_status)
    )
  )

  class(cat_init) <- c(cat_init$class, "cat_initialization")

  return(cat_init)
}

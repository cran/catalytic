#' Initialization for Catalytic Linear Mixed Model (LMM)
#'
#' This function prepares and initializes a catalytic linear mixed model by processing input data,
#' extracting necessary variables, generating synthetic datasets, and fitting a model.
#' (Only consider one random effect variance)
#'
#' @param formula A formula specifying the model. Should include response and predictor variables.
#' @param data A data frame containing the data for modeling.
#' @param x_cols A character vector of column names for fixed effects (predictors).
#' @param y_col A character string for the name of the response variable.
#' @param z_cols A character vector of column names for random effects.
#' @param group_col A character string for the grouping variable (optional). If not given (NULL), it is extracted from the formula.
#' @param syn_size An integer specifying the size of the synthetic dataset to be generated, default is length(x_cols) * 4.
#' @param resample_by_group A logical indicating whether to resample by group, default is FALSE.
#' @param resample_only A logical indicating whether to perform resampling only, default is FALSE.
#' @param na_replace A function to replace NA values in the data, default is mean.
#'
#' @return A list containing the values of all the input arguments and the following components:
#'
#' - **Function Information**:
#'   - `function_name`: A character string representing the name of the function, "cat_lmm_initialization".
#'   - `simple_model`: An object of class `lme4::lmer` or `stats::lm`, representing the fitted model for generating synthetic response from the original data.
#'
#' - **Observation Data Information**:
#'   - `obs_size`: An integer representing the number of observations in the original dataset.
#'   - `obs_data`: The original data used for fitting the model, returned as a data frame.
#'   - `obs_x`: A data frame containing the standardized predictor variables from the original dataset.
#'   - `obs_y`: A numeric vector of the standardized response variable from the original dataset.
#'   - `obs_z`: A data frame containing the standardized random effect variables from the original dataset.
#'   - `obs_group`: A numeric vector representing the grouping variable for the original observations.
#'
#' - **Synthetic Data Information**:
#'   - `syn_size`: An integer representing the number of synthetic observations generated.
#'   - `syn_data`: A data frame containing the synthetic dataset, combining synthetic predictor and response variables.
#'   - `syn_x`: A data frame containing the synthetic predictor variables.
#'   - `syn_y`: A numeric vector of the synthetic response variable values.
#'   - `syn_z`: A data frame containing the synthetic random effect variables.
#'   - `syn_group`: A numeric vector representing the grouping variable for the synthetic observations.
#'   - `syn_x_resample_inform`: A data frame containing information about the resampling process for synthetic predictors:
#'     - Coordinate: Preserves the original data values as reference coordinates during processing.
#'     - Deskewing: Adjusts the data distribution to reduce skewness and enhance symmetry.
#'     - Smoothing: Reduces noise in the data to stabilize the dataset and prevent overfitting.
#'     - Flattening: Creates a more uniform distribution by modifying low-frequency categories in categorical variables.
#'     - Symmetrizing: Balances the data around its mean to improve statistical properties for model fitting.
#'   - `syn_z_resample_inform`: A data frame containing information about the resampling process for synthetic random effects. The resampling methods are the same as those from `syn_x_resample_inform`.
#'
#' - **Whole Data Information**:
#'   - `size`: An integer representing the total size of the combined original and synthetic datasets.
#'   - `data`: A combined data frame of the original and synthetic datasets.
#'   - `x`: A combined data frame of the original and synthetic predictor variables.
#'   - `y`: A combined numeric vector of the original and synthetic response variables.
#'   - `z`: A combined data frame of the original and synthetic random effect variables.
#'   - `group`: A combined numeric vector representing the grouping variable for both original and synthetic datasets.
#'
#' @examples
#' data(mtcars)
#' cat_init <- cat_lmm_initialization(
#'   formula = mpg ~ wt + (1 | cyl), # formula for simple model
#'   data = mtcars,
#'   x_cols = c("wt"), # Fixed effects
#'   y_col = "mpg", # Response variable
#'   z_cols = c("disp", "hp", "drat", "qsec", "vs", "am", "gear", "carb"), # Random effects
#'   group_col = "cyl", # Grouping column
#'   syn_size = 100, # Synthetic data size
#'   resample_by_group = FALSE, # Resampling option
#'   resample_only = FALSE, # Resampling method
#'   na_replace = mean # NA replacement method
#' )
#' cat_init
#' @export
cat_lmm_initialization <- function(formula,
                                   data,
                                   x_cols,
                                   y_col,
                                   z_cols,
                                   group_col = NULL,
                                   syn_size = NULL,
                                   resample_by_group = FALSE,
                                   resample_only = FALSE,
                                   na_replace = mean) {
  if (is.null(group_col)) {
    group_col <- gsub(
      "`",
      "",
      gsub(
        ".*\\((.*?)\\|\\s*(\\S+)\\)",
        "\\2",
        deparse(formula)
      )
    )
  }

  # Validate Input Parameters
  validate_lmm_initialization_input(
    formula = formula,
    data = data,
    x_cols = x_cols,
    y_col = y_col,
    z_cols = z_cols,
    group_col = group_col,
    syn_size = syn_size
  )

  # Extract Information from Input Parameters
  ## Standardize observation covariate and response
  obs_data <- get_standardized_data(
    data = data,
    na_replace = na_replace
  )

  ## Extract observation fixed effect, group, response and random effect
  obs_x <- as.data.frame(obs_data[, x_cols, drop = FALSE])
  obs_group <- obs_data[, group_col]
  obs_y <- obs_data[, y_col]
  obs_z <- as.data.frame(obs_data[, z_cols, drop = FALSE])

  # Generate Synthetic Dataset
  if (is.null(syn_size)) syn_size <- length(x_cols) * 4

  ## Generating synthetic group by resampling observation group
  syn_group_inform <- get_resampled_df(
    data = as.data.frame(obs_group),
    resample_size = syn_size,
    resample_only = TRUE
  )

  syn_group <- as.numeric(syn_group_inform$resampled_df[[1]])

  ## Generate synthetic fixed effect and random effect
  if (resample_by_group) {
    syn_x <- as.data.frame(matrix(
      nrow = syn_size,
      ncol = dim(obs_x)[2],
      dimnames = list(NULL, c(x_cols))
    ))

    syn_z <- as.data.frame(matrix(
      nrow = syn_size,
      ncol = dim(obs_z)[2],
      dimnames = list(NULL, c(z_cols))
    ))

    for (g in unique(c(obs_group)[[1]])) {
      group_idx <- which(syn_group == g)

      syn_x_inform <- get_resampled_df(
        data = obs_x[which(obs_group == g), ],
        resample_size = length(group_idx)
      )
      syn_x[group_idx, ] <- syn_x_inform$resampled_df

      syn_z_inform <- get_resampled_df(
        data = obs_z[which(obs_group == g), ],
        resample_size = length(group_idx)
      )
      syn_z[group_idx, ] <- syn_z_inform$resampled_df
    }

    syn_x_inform <- NULL
    syn_z_inform <- NULL
  } else {
    syn_x_inform <- get_resampled_df(
      data = obs_x,
      resample_size = syn_size
    )

    syn_z_inform <- get_resampled_df(
      data = obs_z,
      resample_size = syn_size
    )

    syn_x <- syn_x_inform$resampled_df
    syn_z <- syn_z_inform$resampled_df
  }

  ## Generate synthetic response by using simple model
  simple_model <- if (grepl("\\(.*\\|.*\\)", deparse(formula))) {
    do.call(
      lme4::lmer,
      list(
        formula = formula,
        data = obs_data
      )
    )
  } else {
    do.call(
      stats::lm, list(
        formula = formula,
        data = data
      )
    )
  }

  syn_data <- cbind(
    syn_x,
    syn_z,
    stats::setNames(as.data.frame(syn_group), group_col)
  )
  syn_y <- as.numeric(stats::predict(simple_model, newdata = syn_data))

  # Finalize Setup and Output
  cat_init <- list(
    function_name = "cat_lmm_initialization",
    ## Input/Processed parameters
    formula = formula,
    x_cols = x_cols,
    y_col = y_col,
    z_cols = z_cols,
    group_col = group_col,
    syn_size = syn_size,
    resample_by_group = resample_by_group,
    resample_only = resample_only,
    na_replace = na_replace,
    simple_model = simple_model,
    ## Observation data information
    obs_size = nrow(obs_x),
    obs_data = data,
    obs_x = obs_x,
    obs_y = obs_y,
    obs_z = obs_z,
    obs_group = obs_group,
    ## Synthetic data information
    syn_size = syn_size,
    syn_data = cbind(syn_data, stats::setNames(as.data.frame(syn_y), y_col)),
    syn_x = syn_x,
    syn_y = syn_y,
    syn_z = syn_z,
    syn_group = syn_group,
    syn_x_resample_inform = syn_x_inform$resampled_df_log,
    syn_z_resample_inform = syn_z_inform$resampled_df_log
  )

  cat_init <- c(
    cat_init,
    list(
      ## Whole data information
      size = cat_init$obs_size + cat_init$syn_size,
      data = rbind(cat_init$obs_data, cat_init$syn_data),
      x = rbind(cat_init$obs_x, cat_init$syn_x),
      y = c(cat_init$obs_y, cat_init$syn_y),
      z = rbind(cat_init$obs_z, cat_init$syn_z),
      group = c(cat_init$obs_group, cat_init$syn_group)
    )
  )

  class(cat_init) <- c(cat_init$class, "cat_initialization")

  return(cat_init)
}

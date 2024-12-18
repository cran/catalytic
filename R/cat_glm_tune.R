#' Catalytic Generalized Linear Models (GLMs) Fitting Function by Tuning tau from a Sequence of tau Values
#'
#' This function tunes a catalytic catalytic Generalized Linear Models (GLMs) by performing specified risk estimate method
#' to estimate the optimal value of the tuning parameter `tau`.  The resulting `cat_glm_tune` object
#' encapsulates the fitted model, including estimated coefficients and family information, facilitating further analysis.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables (e.g. \code{~ .}).
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param risk_estimate_method Method for risk estimation, chosen from "parametric_bootstrap",
#' "cross_validation", "mallows_estimate", "steinian_estimate". Depends on the size of the data if not provided.
#' @param discrepancy_method Method for discrepancy calculation, chosen from "mean_square_error",
#' "mean_classification_error", "logistic_deviance". Depends on the family if not provided.
#' @param tau_seq Vector of numeric values for down-weighting synthetic data.
#' Defaults to a sequence around one fourth of the number of predictors for gaussian and
#' the number of predictors for binomial.
#' @param tau_0 Initial `tau` value used for discrepancy calculation in risk estimation.
#' Defaults to one fourth of the number of predictors for binomial and 1 for gaussian.
#' @param parametric_bootstrap_iteration_times Number of bootstrap iterations for "parametric_bootstrap" risk estimation. Defaults to 100.
#' @param cross_validation_fold_num Number of folds for "cross_validation" risk estimation.. Defaults to 5.
#'
#' @return A list containing the values of all the arguments and the following components:
#'   \item{tau}{Optimal `tau` value determined through tuning.}
#'   \item{model}{Fitted GLM model object with the optimal `tau` value.}
#'   \item{coefficients}{Estimated coefficients from the `model` fitted by the optimal `tau` value.}
#'   \item{risk_estimate_list}{Collected risk estimates for each `tau`.}
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
#' cat_model <- cat_glm_tune(
#'   formula = ~.,
#'   cat_init = cat_init, # Only accept object generated from `cat_glm_initialization`
#'   risk_estimate_method = "parametric_bootstrap",
#'   discrepancy_method = "mean_square_error",
#'   tau_seq = c(1, 2), # Weight for synthetic data
#'   tau_0 = 2,
#'   parametric_bootstrap_iteration_times = 20, # Number of bootstrap iterations
#'   cross_validation_fold_num = 5 # Number of folds
#' )
#' cat_model
#' @export
cat_glm_tune <- function(
    formula,
    cat_init,
    risk_estimate_method = c(
      "parametric_bootstrap",
      "cross_validation",
      "mallowian_estimate",
      "steinian_estimate"
    ),
    discrepancy_method = c(
      "mean_square_error",
      "mean_classification_error",
      "logistic_deviance"
    ),
    tau_seq = NULL,
    tau_0 = NULL,
    parametric_bootstrap_iteration_times = 100,
    cross_validation_fold_num = 5) {
  if (missing(risk_estimate_method)) {
    risk_estimate_method <- ifelse(
      cat_init$obs_size > 200,
      "cross_validation",
      "parametric_bootstrap"
    )
  }
  if (missing(discrepancy_method)) {
    discrepancy_method <- ifelse(
      cat_init$family == "gaussian",
      "mean_square_error",
      "logistic_deviance"
    )
  }

  risk_estimate_method <- match.arg(risk_estimate_method)
  discrepancy_method <- match.arg(discrepancy_method)

  # Validate Input Parameters
  validate_glm_input(
    formula = formula,
    cat_init = cat_init,
    tau_seq = tau_seq,
    tau_0 = tau_0,
    parametric_bootstrap_iteration_times = parametric_bootstrap_iteration_times,
    cross_validation_fold_num = cross_validation_fold_num,
    risk_estimate_method = risk_estimate_method,
    discrepancy_method = discrepancy_method
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  full_formula <- stats::as.formula(paste0(cat_init$y_col_name, f_rhs))

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  if (is.null(tau_seq)) {
    tau_seq <- if (cat_init$family == "gaussian") {
      seq(max(0.01, ncol(cat_init$adj_x) / 4 - 2), ncol(cat_init$adj_x) / 4 + 2, 0.5)
    } else {
      seq(max(0.01, ncol(cat_init$adj_x) - 2), ncol(cat_init$adj_x) + 2, 0.5)
    }
  }

  tau_seq <- sort(unique(tau_seq))

  if ((0 %in% tau_seq) & (dim(cat_init$adj_obs_x)[1] < dim(cat_init$adj_obs_x)[2])) {
    warning(
      paste0(
        "Found 0 in input `tau_seq`, and the number of columns in the observation data exceeds its data size. \n",
        "To avoid issues, 0 in `tau_seq` will be replaced with a small positive value."
      ),
      call. = FALSE
    )

    tau_seq[1] <- tau_seq[2] / 2
  }

  if (is.null(tau_0)) {
    tau_0 <- ifelse(
      cat_init$family == "gaussian",
      1,
      ncol(cat_init$adj_x) / 4
    )
  }

  # Estimate Risk Using the Specified Method and Processed Parameters
  risk_estimate_list <- do.call(
    risk_estimate_method,
    list(
      formula = full_formula,
      cat_init = cat_init,
      discrepancy_method = discrepancy_method,
      tau_seq = tau_seq,
      tau_0 = tau_0,
      parametric_bootstrap_iteration_times = parametric_bootstrap_iteration_times,
      cross_validation_fold_num = cross_validation_fold_num
    )
  )

  # Fit the Optimal Model Using the Optimal tau
  optimal_tau <- tau_seq[which.min(risk_estimate_list)]
  ## Suppress warning for `binomial` family when weights contains value < 1
  optimal_model <- suppressWarnings(
    do.call(
      stats::glm,
      list(
        formula = full_formula,
        family = cat_init$family,
        data = cat_init$data,
        weights = c(
          rep(1, cat_init$obs_size),
          rep(optimal_tau / cat_init$syn_size, cat_init$syn_size)
        )
      )
    )
  )

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_glm_tune",
    ## Input/Processed parameters
    formula = full_formula,
    cat_init = cat_init,
    risk_estimate_method = risk_estimate_method,
    discrepancy_method = if (risk_estimate_method == "mallowian_estimate") {
      "mean_square_error"
    } else if (risk_estimate_method == "steinian_estimate") {
      "logistic_deviance"
    } else {
      discrepancy_method
    },
    tau_seq = tau_seq,
    tau_0 = tau_0,
    parametric_bootstrap_iteration_times = parametric_bootstrap_iteration_times,
    cross_validation_fold_num = cross_validation_fold_num,
    ## Result
    tau = optimal_tau,
    model = optimal_model,
    coefficients = stats::coef(optimal_model),
    risk_estimate_list = risk_estimate_list
  )

  class(cat_model) <- c(cat_model$class, "cat_tune", "cat_glm")

  return(cat_model)
}

#' Perform Parametric Bootstrap for Model Risk Estimation
#'
#' This function performs parametric bootstrapping to estimate model risk. It fits a sequence
#' of Generalized Linear Models (GLMs) with different values of `tau`, calculates the in-sample
#' prediction error, and incorporates deviations from the bootstrap response samples. The final
#' risk estimate is obtained by combining the in-sample error and the covariance penalty derived
#' from the bootstrap samples.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables.
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param tau_seq A sequence of tuning parameter values (`tau`) over which
#' the model risk will be estimated. Each `tau` value is used to weight the synthetic data during model fitting.
#' @param tau_0 A reference value for `tau` used in the preliminary estimate model and variance calculation.
#' @param discrepancy_method The method used to calculate the discrepancy (e.g., logistic deviance).
#' @param parametric_bootstrap_iteration_times The number of bootstrap iterations to perform.
#' @param ... Other arguments passed to other internal functions.
#'
#' @details
#' 1. **Preliminary Estimate Model**: The function first fits a GLM model using the observed
#' and synthetic data with an initial value of `tau_0` for the synthetic data weights.
#'
#' 2. **Bootstrap Samples**: The function generates bootstrap response samples based on the
#' mean and standard deviation of the preliminary estimate model, using parametric bootstrapping.
#'
#' 3. **In-sample Prediction Error**: For each value of `tau` in `tau_seq`, the function computes
#' the in-sample prediction error (e.g., using logistic deviance).
#'
#' 4. **Bootstrap Models**: For each bootstrap iteration, the function fits a GLM using the
#' bootstrap response samples and calculates the corresponding lambda values.
#'
#' 5. **Covariance Penalty**: The function approximates the covariance penalty using the weighted
#' deviations across all bootstrap iterations.
#'
#' 6. **Final Risk Estimate**: The final model risk estimate is calculated by summing the in-sample
#' prediction error and the average weighted deviations from the bootstrap response samples.
#'
#' @return
#' A numeric vector containing the risk estimates for each `tau` in `tau_seq`.
parametric_bootstrap <- function(
    formula,
    cat_init,
    tau_seq,
    tau_0,
    discrepancy_method,
    parametric_bootstrap_iteration_times,
    ...) {
  # The Preliminary Estimate Model from Using Initial Value tau_0
  ## Suppress warning for `binomial` family when weights contains value < 1
  tau_0_model <- suppressWarnings(
    do.call(
      stats::glm,
      list(
        formula = formula,
        family = cat_init$family,
        data = cat_init$data,
        weights = c(
          rep(1, cat_init$obs_size),
          rep(tau_0 / cat_init$syn_size, cat_init$syn_size)
        )
      )
    )
  )

  # The Mean from the Preliminary Estimate Model for Generating Bootstrap Response Samples
  boots_mean_Y <- get_glm_mean(
    family_string = cat_init$family,
    X = cat_init$adj_obs_x,
    coefs = stats::coef(tau_0_model)
  )

  # The Standard Deviance for Generating Bootstrap Response Samples
  boots_sd_Y <- if (cat_init$family == "gaussian") {
    sqrt(
      ifelse(
        (is.null(cat_init$custom_variance)),
        get_glm_custom_var(
          formula = formula,
          cat_init = cat_init,
          tau = ncol(cat_init$adj_x)
        ),
        cat_init$custom_variance
      )
    )
  } else {
    0
  }

  # The data.frame that Contains All Vectors of Bootstrap Response Samples for Each Bootstrap Iterations
  boots_Y_df <- as.data.frame(
    lapply(
      1:parametric_bootstrap_iteration_times,
      function(.) {
        get_glm_sample_data(
          family_string = cat_init$family,
          n = cat_init$obs_size,
          mean = boots_mean_Y,
          sd = boots_sd_Y
        )
      }
    )
  )

  # The Deviations of Bootstrap Response Samples from Their Means
  boots_dev_df <- sweep(boots_Y_df, 1, rowMeans(boots_Y_df))

  risk_estimate_list <- c()

  for (tau in tau_seq) {
    # The tau Model
    ## Suppress warning for `binomial` family when weights contains value < 1
    tau_model <- suppressWarnings(
      do.call(
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
    )


    # The In-sample Prediction Error from the tau Model
    D <- get_discrepancy(
      discrepancy_method = discrepancy_method,
      family_string = cat_init$family,
      X = cat_init$adj_obs_x,
      Y = cat_init$obs_y,
      coefs = stats::coef(tau_model)
    )

    boots_lambda_mu_df <- do.call(
      cbind,
      lapply(
        1:parametric_bootstrap_iteration_times,
        function(boots_i) {
          # The Bootstrap Model from Replacing the Observation Response by Bootstrap Response Samples
          ## Suppress warning for `binomial` family when weights contains value < 1
          boots_model <- suppressWarnings(
            do.call(
              stats::glm, list(
                formula = formula,
                family = cat_init$family,
                data = cbind(
                  cat_init$x,
                  stats::setNames(
                    data.frame(c(
                      boots_Y_df[[boots_i]],
                      cat_init$syn_y
                    )),
                    cat_init$y_col_name
                  )
                ),
                weights = c(
                  rep(1, cat_init$obs_size),
                  rep(tau / cat_init$syn_size, cat_init$syn_size)
                )
              )
            )
          )

          # Return Lambda(mu) values
          return(get_glm_lambda(
            discrepancy_method = discrepancy_method,
            X = cat_init$adj_obs_x,
            coefs = stats::coef(boots_model)
          ))
        }
      )
    )

    risk_estimate_list <- c(
      risk_estimate_list,
      # Covariance Penalty Approximation Using Weighted Deviations Across Bootstrap Iterations
      D + sum(boots_dev_df * boots_lambda_mu_df) / parametric_bootstrap_iteration_times
    )
  }

  return(risk_estimate_list)
}

#' Perform Steinian Estimate for Model Risk (Only Applicable for Binomial Family)
#'
#' This function computes the Steinian estimate for model risk by fitting a sequence of
#' Generalized Linear Models (GLMs) with varying values of `tau`. It combines the preliminary
#' estimate from a model fitted with an initial `tau_0` value with a penalty term that incorporates
#' the in-sample prediction error and a covariance penalty, which is based on models fitted by inverting
#' the response of individual observations.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables.
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param tau_seq A sequence of tuning parameter values (`tau`) over which
#' the Steinian estimate will be computed. Each value of `tau` is used to weight the
#' synthetic data during model fitting.
#' @param tau_0 A reference value for `tau` that is used in the calculation of the
#' preliminary estimate model and the variance term.
#' @param ... Other arguments passed to other internal functions.
#'
#' @details
#' 1. **Preliminary Estimate Model**: The function first fits a GLM model using the observed and
#' synthetic data with an initial value of `tau_0` for the synthetic data weights.
#'
#' 2. **In-sample Prediction Error**: For each value of `tau` in `tau_seq`, the function computes
#' the in-sample prediction error (logistic deviance).
#'
#' 3. **Steinian Penalty**: The function calculates the Steinian covariance penalty for each observation
#' by fitting a modified model that inverts one observation at a time. The penalty is added to the
#' in-sample prediction error to obtain the final risk estimate.
#'
#' 4. **Steinian Risk Estimate**: The final Steinian risk estimate is calculated by summing the
#' in-sample prediction error and the Steinian penalty term for each value of `tau` in `tau_seq`.
#'
#' @return
#' A numeric vector of Steinian risk estimates, one for each value of `tau` in `tau_seq`.
steinian_estimate <- function(
    formula,
    cat_init,
    tau_seq,
    tau_0,
    ...) {
  # The Preliminary Estimate Model from Using Initial Value tau_0
  ## Suppress warning for `binomial` family when weights contains value < 1
  tau_0_model <- suppressWarnings(
    do.call(
      stats::glm,
      list(
        formula = formula,
        family = "binomial",
        data = cat_init$data,
        weights = c(
          rep(1, cat_init$obs_size),
          rep(tau_0 / cat_init$syn_size, cat_init$syn_size)
        )
      )
    )
  )

  # The Mean from the Preliminary Estimate Model for Generating Bootstrap Response Samples
  tau_0_mean <- get_glm_mean(
    family_string = "binomial",
    X = cat_init$adj_obs_x,
    coefs = stats::coef(tau_0_model)
  )

  # The Variance Estimates for the Preliminary Estimate Model
  tau_0_varience <- tau_0_mean * (1 - tau_0_mean)

  risk_estimate_list <- c()

  for (tau in tau_seq) {
    # The tau Model
    ## Suppress warning for `binomial` family when weights contains value < 1
    tau_model <- suppressWarnings(
      do.call(
        stats::glm,
        list(
          formula = formula,
          family = "binomial",
          data = cat_init$data,
          weights = c(
            rep(1, cat_init$obs_size),
            rep(tau / cat_init$syn_size, cat_init$syn_size)
          )
        )
      )
    )

    # The In-sample Prediction Error from the tau Model
    D <- get_discrepancy(
      discrepancy_method = "logistic_deviance",
      family_string = "binomial",
      X = cat_init$adj_obs_x,
      Y = cat_init$obs_y,
      coefs = stats::coef(tau_model)
    )

    lambda_mu_tau_df <- get_glm_lambda(
      discrepancy_method = "logistic_deviance",
      X = cat_init$adj_x,
      coefs = stats::coef(tau_model)
    )

    steinian_penalty_sum <- 0

    for (obs_i in 1:cat_init$obs_size) {
      # The Triangle Response by Inverting One Observation
      triangle_obs_y <- cat_init$obs_y
      triangle_obs_y[obs_i] <- 1 - cat_init$obs_y[obs_i]

      # The Triangle Model Fitted from Triangle Response
      ## Suppress warning for `binomial` family when weights contains value < 1
      triangle_tau_model <- suppressWarnings(
        do.call(
          stats::glm,
          list(
            formula = formula,
            family = "binomial",
            data = cbind(
              cat_init$adj_x,
              stats::setNames(
                data.frame(c(
                  triangle_obs_y,
                  cat_init$syn_y
                )),
                cat_init$y_col_name
              )
            ),
            weights = c(
              rep(1, cat_init$obs_size),
              rep(tau / cat_init$syn_size, cat_init$syn_size)
            )
          )
        )
      )

      #  Steinian Covariance Penalty
      steinian_penalty_sum <- steinian_penalty_sum +
        tau_0_varience[obs_i] * (
          2 * cat_init$obs_y[obs_i] - 1) * (
          lambda_mu_tau_df[obs_i] - get_glm_lambda(
            discrepancy_method = "logistic_deviance",
            X = cat_init$adj_x[obs_i, ],
            coefs = stats::coef(triangle_tau_model)
          )
        )
    }

    risk_estimate_list <- c(risk_estimate_list, D + steinian_penalty_sum)
  }

  return(risk_estimate_list)
}

#' Perform Mallowian Estimate for Model Risk (Only Applicable for Gaussian Family)
#'
#' This function calculates the Mallowian estimate for model risk by fitting a sequence of
#' Generalized Linear Models (GLMs) with varying values of `tau`. It uses the in-sample prediction
#' error along with a regularized projection matrix to estimate the model risk. The `tau` parameter
#' influences the weighting of synthetic data during model fitting.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables.
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param tau_seq A sequence of tuning parameter values (`tau`) over which
#' the Mallowian estimate will be computed. Each value of `tau` is used to weight the
#' synthetic data during model fitting.
#' @param ... Other arguments passed to other internal functions.
#'
#' @details
#' 1. **Model Fitting**: For each value of `tau` in `tau_seq`, the function fits a GLM model
#' using the observed and synthetic data. The synthetic data is weighted by the corresponding
#' `tau` value during the fitting process.
#'
#' 2. **In-sample Prediction Error**: After fitting the model, the function computes the
#' in-sample prediction error (Mean Squared Error) to assess the model's performance.
#'
#' 3. **Regularized Projection Matrix**: The function calculates a regularized projection matrix
#' using the observed and synthetic data, which influences the covariance matrix used in risk estimation.
#'
#' 4. **Mallowian Risk Estimate**: The final Mallowian risk estimate is computed by combining the
#' in-sample prediction error with a penalty term involving the projection matrix and a variance term.
#' This estimate is calculated for each value of `tau` in `tau_seq`.
#'
#' @return
#' A numeric vector of Mallowian risk estimates, one for each value of `tau` in `tau_seq`.
mallowian_estimate <- function(
    formula,
    cat_init,
    tau_seq,
    ...) {
  var_Y <- ifelse(
    (is.null(cat_init$custom_variance)),
    get_glm_custom_var(
      formula = formula,
      cat_init = cat_init,
      tau = ncol(cat_init$adj_x)
    ),
    cat_init$custom_variance
  )

  risk_estimate_list <- c()

  for (tau in tau_seq) {
    # The tau Model
    tau_model <- do.call(
      stats::glm,
      list(
        formula = formula,
        family = "gaussian",
        data = cat_init$data,
        weights = c(
          rep(1, cat_init$obs_size),
          rep(tau / cat_init$syn_size, cat_init$syn_size)
        )
      )
    )

    # The In-sample Prediction Error from the tau Model
    D <- get_discrepancy(
      discrepancy_method = "mean_square_error",
      family_string = "gaussian",
      X = cat_init$adj_obs_x,
      Y = cat_init$obs_y,
      coefs = stats::coef(tau_model)
    )

    obs_x_matrix <- as.matrix(cat_init$adj_obs_x)
    syn_x_matrix <- as.matrix(cat_init$adj_syn_x)

    solved_obs_syn_matrix <- solve(
      crossprod(obs_x_matrix) + (tau / cat_init$syn_size) * crossprod(syn_x_matrix)
    )

    weighted_covariance_matrix <- t(obs_x_matrix) + (
      tau / ((cat_init$syn_size * cat_init$obs_size))
    ) * t(syn_x_matrix) %*% matrix(1, nrow = cat_init$syn_size, ncol = cat_init$obs_size)

    # The Regularized Projection Matrix
    H_tau_matrix <- obs_x_matrix %*% solved_obs_syn_matrix %*% weighted_covariance_matrix

    risk_estimate_list <- c(
      risk_estimate_list,
      D + (2 * var_Y * sum(diag(H_tau_matrix))) / cat_init$obs_size
    )
  }

  return(risk_estimate_list)
}

#' Perform Cross-Validation for Model Estimation
#'
#' This function performs cross-validation for estimating risk over a sequence
#' of tuning parameters (`tau_seq`) by fitting a Generalized Linear Model (GLM) to the data.
#' It evaluates model performance by splitting the dataset into multiple folds, training
#' the model on a subset of the data, and testing it on the remaining portion.
#'
#' @param formula A formula specifying the GLMs. Should at least include response variables.
#' @param cat_init A list generated from `cat_glm_initialization`.
#' @param tau_seq A sequence of tuning parameter values (`tau`) over which
#' cross-validation will be performed. Each value of `tau` is used to weight the
#' synthetic data during model fitting.
#' @param discrepancy_method A function used to calculate the discrepancy (error) between
#' model predictions and actual values.
#' @param cross_validation_fold_num The number of folds to use in cross-validation.
#' The dataset will be randomly split into this number of subsets, and the model will be trained and tested on different combinations of these subsets.
#' @param ... Other arguments passed to other internal functions.
#'
#' @details
#' 1. **Randomization of the Data**: The data is randomly shuffled into `cross_validation_fold_num`
#' subsets to ensure that the model is evaluated across different splits of the dataset.
#'
#' 2. **Model Training and Prediction**: For each fold, a training set is used to fit
#' a GLM with varying values of `tau` (from `tau_seq`), and the model is evaluated on a test set.
#' The training data consists of both the observed and synthetic data, with synthetic data weighted by `tau`.
#'
#' 3. **Risk Estimation**: After fitting the model, the `discrepancy_method` is used to calculate the
#' prediction error for each combination of fold and `tau`. These errors are accumulated for each `tau`.
#'
#' 4. **Average Risk Estimate**: After completing all folds, the accumulated prediction errors
#' are averaged over the number of folds to provide a final risk estimate for each value of `tau`.
#'
#' @return
#' A numeric vector of averaged risk estimates, one for each value of `tau` in `tau_seq`.
cross_validation <- function(
    formula,
    cat_init,
    tau_seq,
    discrepancy_method,
    cross_validation_fold_num,
    ...) {
  random_idxs <- sample(1:cat_init$obs_size)
  risk_estimate_list <- rep(0, length(tau_seq))

  # The perform of cross-validation
  for (fold in seq_len(cross_validation_fold_num)) {
    test_idx <- random_idxs[
      (cat_init$obs_size * (fold - 1) / cross_validation_fold_num + 1):
      (cat_init$obs_size * fold / cross_validation_fold_num)
    ]

    test_x <- cat_init$adj_obs_x[test_idx, , drop = FALSE]
    test_y <- cat_init$obs_y[test_idx]

    train_x <- cat_init$adj_obs_x[-test_idx, , drop = FALSE]
    train_y <- cat_init$obs_y[-test_idx]

    risk_estimate_list_fold <- c()

    # The iteration over the tau sequence for estimation
    for (tau in tau_seq) {
      ## Suppress warning for `binomial` family when weights contains value < 1
      cv_tau_model <- suppressWarnings(
        do.call(
          stats::glm,
          list(
            formula = formula,
            family = cat_init$family,
            data = cbind(
              rbind(train_x, cat_init$adj_syn_x),
              stats::setNames(
                data.frame(c(train_y, cat_init$syn_y)),
                cat_init$y_col_name
              )
            ),
            weights = c(
              rep(1, nrow(train_x)),
              rep(tau / cat_init$syn_size, cat_init$syn_size)
            )
          )
        )
      )

      prediction_error <- get_discrepancy(
        discrepancy_method = discrepancy_method,
        family_string = cat_init$family,
        X = test_x,
        Y = test_y,
        coefs = stats::coef(cv_tau_model)
      )

      risk_estimate_list_fold <- c(
        risk_estimate_list_fold,
        prediction_error
      )
    }

    risk_estimate_list <- risk_estimate_list + risk_estimate_list_fold
  }

  return(risk_estimate_list / cross_validation_fold_num)
}

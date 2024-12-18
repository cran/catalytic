#' Catalytic Cox Proportional-Hazards Model (COX) Fitting Function by Tuning tau from a Sequence of tau Values
#'
#' This function tunes a catalytic Cox proportional-hazards model (COX) by performing cross-validation
#' to estimate the optimal value of the tuning parameter `tau`. It finally uses the optimal tau value
#' in the `cat_cox` function for model fitting.
#'
#' @param formula A formula specifying the Cox model. Should at least include response variables (e.g. \code{~.}).
#' @param cat_init A list generated from `cat_cox_initialization`.
#' @param method The estimation method, either `"CRE"` (Catalytic-regularized Estimator) or `"WME"` (Weighted Mixture Estimator).
#' @param tau_seq A numeric vector specifying the sequence of `tau` values to be tested.
#'                If NULL, a default sequence is generated based on the number of predictors.
#' @param cross_validation_fold_num An integer representing the number of folds for cross-validation.
#'                                   Defaults to 5.
#' @param ... Additional arguments passed to the `cat_cox` function for model fitting.
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{tau}{he optimal `tau` value determined from cross-validation.}
#' \item{model}{The fitted lmer model object by using the optimal `tau` value.}
#' \item{coefficients}{Coefficients of the fitted model by using the optimal `tau` value.}
#' \item{likelihood_list}{Average likelihood value for each `tau `value.}
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
#'   hazard_constant = 0.1, # Hazard rate value
#'   entry_points = rep(0, nrow(cancer)), # Entry points of each observation
#'   x_degree = rep(1, ncol(cancer) - 2), # Degrees for polynomial expansion of predictors
#'   resample_only = FALSE, # Whether to perform resampling only
#'   na_replace = stats::na.omit # How to handle NA values in data
#' )
#'
#' cat_model <- cat_cox_tune(
#'   formula = ~., # Should at least include response variables
#'   cat_init = cat_init, # Only accept object generated from `cat_cox_initialization`
#'   tau_seq = c(1, 2), # Vector of weights for synthetic data
#'   cross_validation_fold_num = 5 # number of folds for cross-validation
#' )
#' cat_model
#' @export
cat_cox_tune <- function(
    formula,
    cat_init,
    method = c("CRE", "WME"),
    tau_seq = NULL,
    cross_validation_fold_num = 5,
    ...) {
  method <- match.arg(method)

  # Validate Input Parameters
  validate_cox_input(
    formula = formula,
    cat_init = cat_init,
    tau_seq = tau_seq,
    cross_validation_fold_num = cross_validation_fold_num,
    ...
  )

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  if (is.null(tau_seq)) {
    tau_seq <- seq(
      max(0.01, ncol(cat_init$adj_obs_x) - 2),
      ncol(cat_init$adj_obs_x) + 2, 0.5
    )
  }

  tau_seq <- sort(unique(tau_seq))

  # Calcualte Likelihood Using the Specified Method (Cross-validation) and Processed Parameters
  likelihood_list <- do.call(
    "cross_validation_cox",
    list(
      formula = formula,
      cat_init = cat_init,
      method = method,
      tau_seq = tau_seq,
      cross_validation_fold_num = cross_validation_fold_num,
      ...
    )
  )

  # Fit the Optimal Model Using the Optimal tau
  optimal_tau <- tau_seq[which.max(likelihood_list)]
  optimal_model <- do.call(
    cat_cox,
    list(
      formula = formula,
      cat_init = cat_init,
      method = method,
      tau = optimal_tau,
      ...
    )
  )

  cat_model <- list(
    function_name = "cat_cox_tune",
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
    tau_seq = tau_seq,
    method = method,
    cross_validation_fold_num = cross_validation_fold_num,
    ...,
    ## Result
    tau = optimal_tau,
    model = optimal_model,
    coefficients = stats::coef(optimal_model),
    likelihood_list = likelihood_list
  )

  class(cat_model) <- c(cat_model$class, "cat_tune", "cat_cox")

  return(cat_model)
}

#' Perform Cross-Validation for Catalytic Cox Proportional-Hazards Model (COX) to Select Optimal tau
#'
#' This function performs cross-validation for the catalytic Cox proportional-hazards model (COX) to
#' estimate the likelihood associated with different values of tau. It splits the data into
#' training and testing sets and computes prediction errors for model evaluation.
#'
#' @param formula A formula specifying the Cox model. Should at least include response variables.
#' @param cat_init A list containing initialized parameters for the catalytic COX.
#' @param method Character string specifying the optimization method used in the Cat-Cox model fitting.
#' @param tau_seq A numeric vector of tau values for which to estimate likelihood.
#' @param cross_validation_fold_num An integer indicating the number of folds for cross-validation.
#' @param ... Additional arguments passed to the `cat_cox` function for model fitting.
#'
#' @return A numeric vector containing the average likelihood estimates for each tau value.
cross_validation_cox <- function(formula,
                                 cat_init,
                                 method,
                                 tau_seq,
                                 cross_validation_fold_num,
                                 ...) {
  random_idxs <- sample(1:cat_init$obs_size)
  likelihood_list <- rep(0, length(tau_seq))

  # The perform of cross-validation
  for (fold in seq_len(cross_validation_fold_num)) {
    not_incl_idx <- random_idxs[
      (cat_init$obs_size * (fold - 1) / cross_validation_fold_num + 1):
      (cat_init$obs_size * fold / cross_validation_fold_num)
    ]

    train_x <- cat_init$obs_x[-not_incl_idx, , drop = FALSE]
    train_time <- cat_init$obs_time[-not_incl_idx]
    train_status <- cat_init$obs_status[-not_incl_idx]

    # Create a specific cat_init from train data
    cv_cat_init <- cat_init
    cv_cat_init$obs_size <- nrow(train_x)
    cv_cat_init$obs_x <- train_x
    cv_cat_init$obs_time <- train_time
    cv_cat_init$obs_status <- train_status
    cv_cat_init$obs_data <- cbind(
      cv_cat_init$obs_x,
      stats::setNames(
        data.frame(cv_cat_init$obs_time),
        cat_init$time_col_name
      ),
      stats::setNames(
        data.frame(cv_cat_init$obs_status),
        cat_init$status_col_name
      )
    )

    cv_cat_init$data <- rbind(cv_cat_init$obs_data, cv_cat_init$syn_data)
    cv_cat_init$x <- rbind(cv_cat_init$obs_x, cv_cat_init$syn_x)
    cv_cat_init$time <- c(cv_cat_init$obs_time, cv_cat_init$syn_time)
    cv_cat_init$status <- c(cv_cat_init$obs_status, cv_cat_init$syn_status)
    cv_cat_init$size <- cv_cat_init$obs_size + cv_cat_init$syn_size
    cv_cat_init$entry_points <- cv_cat_init$entry_points[-not_incl_idx]

    likelihood_list_fold <- c()

    # The iteration over the tau sequence for estimation
    for (tau in tau_seq) {
      cv_model <- do.call(
        cat_cox,
        list(
          formula = formula,
          cat_init = cv_cat_init,
          method = method,
          tau = tau,
          ...
        )
      )

      cox_partial_likelihood <- get_cox_partial_likelihood(
        X = cat_init$obs_x,
        time = cat_init$obs_time,
        status = cat_init$obs_status,
        coefs = cv_model$coefficients,
        entry_points = cat_init$entry_points
      ) - get_cox_partial_likelihood(
        X = train_x,
        time = train_time,
        status = train_status,
        coefs = cv_model$coefficients,
        entry_points = cv_cat_init$entry_points
      )

      likelihood_list_fold <- c(
        likelihood_list_fold,
        cox_partial_likelihood
      )
    }

    likelihood_list <- likelihood_list + likelihood_list_fold
  }

  return(likelihood_list / cross_validation_fold_num)
}

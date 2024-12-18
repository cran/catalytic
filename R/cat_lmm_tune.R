#' Catalytic Linear Mixed Model (LMM) Fitting Function by Tuning tau from a Sequence of tau Values
#'
#' This function tunes a catalytic linear mixed model by performing cross-validation
#' to estimate the optimal value of the tuning parameter tau. It finally uses the optimal tau value
#' in the `lmer` function from the `lme4` package for model fitting. (Only consider one random effect variance)
#'
#' @param cat_init A list generated from `cat_lmm_initialization`.
#' @param tau_seq A numeric vector specifying the sequence of tau values to be tested.
#'                If NULL, a default sequence is generated based on the number of predictors.
#' @param cross_validation_fold_num An integer representing the number of folds for cross-validation.
#'                                   Defaults to 5.
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{tau}{The optimal tau value determined from cross-validation.}
#' \item{model}{The fitted lmer model object by using the optimal tau value.}
#' \item{coefficients}{Coefficients of the fitted model by using the optimal tau value.}
#' \item{risk_estimate_list}{Average prediction errors for each tau value.}
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
#'
#' cat_model <- cat_lmm_tune(
#'   cat_init = cat_init, # Only accept object generated from cat_lmm_initialization
#'   tau_seq = c(1, 2), # Vector of weights for synthetic data
#'   cross_validation_fold_num = 3 # number of folds for cross-validation
#' )
#' cat_model
#'
#' @export
cat_lmm_tune <- function(
    cat_init,
    tau_seq = NULL,
    cross_validation_fold_num = 5) {
  # Validate Input Parameters
  validate_lmm_input(
    cat_init = cat_init,
    tau_seq = tau_seq,
    cross_validation_fold_num = cross_validation_fold_num,
  )

  if (is.null(tau_seq)) {
    tau_seq <- seq(
      max(0, ncol(cat_init$x) / 4 - 2),
      ncol(cat_init$x) / 4 + 2, 0.5
    )
  }

  tau_seq <- sort(unique(tau_seq))

  risk_estimate_list <- do.call(
    "cross_validation_lmm",
    list(
      cat_init = cat_init,
      tau_seq = tau_seq,
      cross_validation_fold_num = cross_validation_fold_num
    )
  )

  optimal_tau <- tau_seq[which.min(risk_estimate_list)]

  Y <- cat_init$y
  X <- as.matrix(cat_init$x)
  group <- cat_init$group

  optimal_model <- do.call(
    lme4::lmer,
    list(
      formula = Y ~ X + (1 | group),
      weight = c(
        rep(1, cat_init$obs_size),
        rep(optimal_tau / cat_init$syn_size, cat_init$syn_size)
      ),
      control = lme4::lmerControl(check.conv.singular = "ignore")
    )
  )

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_lmm_tune",
    ## Input/Processed parameters
    cat_init = cat_init,
    tau_seq = tau_seq,
    cross_validation_fold_num = cross_validation_fold_num,
    discrepancy_method = "mean_square_error",
    ## Result
    tau = optimal_tau,
    model = optimal_model,
    coefficients = as.numeric(lme4::fixef(optimal_model)),
    risk_estimate_list = risk_estimate_list
  )

  class(cat_model) <- c(cat_model$class, "cat_tune", "cat_lmm")

  return(cat_model)
}

#' Perform Cross-Validation for Catalytic Linear Mixed Model (LMM) to Select Optimal tau
#'
#' This function performs cross-validation for the catalytic linear mixed model (LMM) to
#' estimate the risk associated with different values of tau. It splits the data into
#' training and testing sets and computes prediction errors for model evaluation.
#'
#' @param cat_init A list containing initialized parameters for the catalytic LMM.
#' @param tau_seq A numeric vector of tau values for which to estimate risk.
#' @param cross_validation_fold_num An integer indicating the number of folds for cross-validation.
#'
#' @return A numeric vector containing the average risk estimates for each tau value.
cross_validation_lmm <- function(
    cat_init,
    tau_seq,
    cross_validation_fold_num = 5) {
  random_sample <- sample(1:cat_init$obs_size)
  risk_estimate_list <- rep(0, length(tau_seq))

  # The perform of cross-validation
  for (fold in seq_len(cross_validation_fold_num)) {
    test_idx <- random_sample[
      (cat_init$obs_size * (fold - 1) / cross_validation_fold_num + 1):
      (cat_init$obs_size * fold / cross_validation_fold_num)
    ]

    test_x <- cat_init$obs_x[test_idx, , drop = FALSE]
    test_y <- cat_init$obs_y[test_idx]
    test_z <- cat_init$obs_z[test_idx, , drop = FALSE]
    test_group <- cat_init$obs_group[test_idx]

    train_x <- cat_init$obs_x[-test_idx, , drop = FALSE]
    train_y <- cat_init$obs_y[-test_idx]
    train_z <- cat_init$obs_z[-test_idx, , drop = FALSE]
    train_group <- cat_init$obs_group[-test_idx]

    risk_estimate_list_fold <- c()

    # The iteration over the tau sequence for estimation
    for (tau in tau_seq) {
      cv_y <- c(train_y, cat_init$syn_y)
      cv_x <- as.matrix(rbind(train_x, cat_init$syn_x))
      cv_group <- c(train_group, cat_init$syn_group)

      cv_model <- do.call(
        lme4::lmer,
        list(
          formula = cv_y ~ cv_x + (1 | cv_group),
          weight = c(
            rep(1, nrow(train_x)),
            rep(tau / cat_init$syn_size, cat_init$syn_size)
          ),
          control = lme4::lmerControl(check.conv.singular = "ignore")
        )
      )

      prediction_error <- get_discrepancy(
        discrepancy_method = "mean_square_error",
        family_string = "gaussian",
        X = test_x,
        Y = test_y,
        coefs = lme4::fixef(cv_model)
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

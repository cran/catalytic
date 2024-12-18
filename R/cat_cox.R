#' Catalytic Cox Proportional Hazards Model (COX) Fitting Function with Fixed Tau
#'
#' Fits a Catalytic Cox proportional hazards model for survival data with specified variance parameters
#' and iterative coefficient estimation, with either `CRE` (Catalytic-regularized Estimator) or `WME` (Weighted Mixture Estimator) methods.
#'
#' @param formula A formula specifying the Cox model. Should at least include response variables (e.g. \code{~ .}).
#' @param cat_init A list generated from `cat_cox_initialization`.
#' @param tau Optional numeric scalar controlling the weight of the synthetic data in the coefficient estimation, defaults to the number of predictors.
#' @param method The estimation method, either `"CRE"` (Catalytic-regularized Estimator) or `"WME"` (Weighted Mixture Estimator).
#' @param init_coefficients Initial coefficient values before iteration. Defaults to zero if not provided (if using `CRE`).
#' @param tol Convergence tolerance for iterative methods. Default is `1e-5` (if using `CRE`).
#' @param max_iter Maximum number of iterations allowed for convergence. Default is `25` (if using `CRE`).
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{coefficients}{Estimated coefficient vector.}
#' \item{model}{Fitted Cox model object (if using `WME`).}
#' \item{iteration_log}{Matrix logging variance and coefficient values for each iteration(if using `CRE`).}
#' \item{iter}{Number of iterations (if using `CRE`).}
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
#' cat_model_cre <- cat_cox(
#'   formula = ~.,
#'   cat_init = cat_init, # Only accept object generated from `cat_cox_initialization`
#'   tau = 1, # Weight for synthetic data
#'   method = "CRE", # Choose from `"CRE"` or `"WME"`
#'   init_coefficients = rep(0, ncol(cat_init$x)), # Initial coefficient values (Only for `CRE`)
#'   tol = 1e-1, # Tolerance for convergence criterion  (Only for `CRE`)
#'   max_iter = 3 # Maximum number of iterations for convergence  (Only for `CRE`)
#' )
#' cat_model_cre
#'
#' cat_model_wme <- cat_cox(
#'   formula = ~.,
#'   cat_init = cat_init, # Only accept object generated from `cat_cox_initialization`
#'   tau = 1, # Weight for synthetic data
#'   method = "WME"
#' )
#' cat_model_wme
#' @export
cat_cox <- function(formula,
                    cat_init,
                    tau = NULL,
                    method = c("CRE", "WME"),
                    init_coefficients = NULL,
                    tol = 1e-5,
                    max_iter = 25) {
  method <- match.arg(method)

  # Update cat_init with adjusted data based on the formula's right-hand side
  f_rhs <- get_formula_rhs(formula, with_tilde = TRUE)
  full_formula <- stats::as.formula(paste0(
    "survival::Surv(",
    cat_init$time_col_name,
    ",",
    cat_init$status_col_name,
    ")",
    f_rhs
  ))

  cat_init <- get_adjusted_cat_init(
    cat_init = cat_init,
    formula_rhs = f_rhs
  )

  # Validate Input Parameters
  validate_cox_input(
    formula = formula,
    cat_init = cat_init,
    tau = tau,
    init_coefficients = init_coefficients,
    tol = tol,
    max_iter = max_iter
  )

  if (is.null(tau)) tau <- ncol(cat_init$adj_x)

  coefs_model <- NULL
  iteration_log <- NULL
  iter <- NULL

  if (method == "CRE") {
    coefs <- if (is.null(init_coefficients)) rep(0, ncol(cat_init$adj_obs_x)) else init_coefficients

    diff_new_old <- TRUE

    iter <- 0
    iteration_log <- matrix(
      data = c(coefs),
      ncol = (length(coefs)),
      dimnames = list(NULL, colnames(cat_init$adj_x))
    )

    while ((iter < max_iter) & diff_new_old) {
      coefs_old <- coefs

      hessian_matrix <- get_cox_hessian(
        X = cat_init$adj_obs_x,
        time = cat_init$obs_time,
        status = cat_init$obs_status,
        coefs = coefs_old,
        entry_points = cat_init$entry_points
      ) + tau / cat_init$syn_size * get_cox_syn_hessian(
        X = cat_init$adj_syn_x,
        time = cat_init$syn_time,
        coefs = coefs_old,
        hazard_constant = cat_init$hazard_constant
      )

      gradient_vector <- get_cox_gradient(
        X = cat_init$adj_obs_x,
        time = cat_init$obs_time,
        status = cat_init$obs_status,
        coefs = coefs_old,
        entry_points = cat_init$entry_points
      ) + tau / cat_init$syn_size * get_cox_syn_gradient(
        X = cat_init$adj_syn_x,
        time = cat_init$syn_time,
        coefs = coefs_old,
        hazard_constant = cat_init$hazard_constant
      )

      coefs <- coefs_old - get_cox_qr_solve(hessian_matrix, gradient_vector)

      iter <- iter + 1
      iteration_log <- rbind(iteration_log, c(coefs))

      diff_new_old <- norm(as.numeric(coefs_old - coefs), "2") > (tol * norm(as.numeric(coefs_old), "2"))
    }
  } else {
    coefs_model <- do.call(
      survival::coxph,
      list(
        formula = full_formula,
        data = cat_init$data,
        weights = c(
          rep(1, cat_init$obs_size),
          rep(tau / cat_init$syn_size, cat_init$syn_size)
        ),
        ties = "breslow"
      )
    )

    coefs <- stats::coef(coefs_model)
  }

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_cox",
    ## Input/Processed parameters
    formula = full_formula,
    cat_init = cat_init,
    tau = tau,
    method = method,
    init_coefficients = init_coefficients,
    tol = tol,
    max_iter = max_iter,
    ## Result
    iter = iter,
    iteration_log = iteration_log,
    model = coefs_model,
    coefficients = coefs
  )

  class(cat_model) <- c(cat_model$class, "cat", "cat_cox")

  return(cat_model)
}

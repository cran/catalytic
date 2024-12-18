#' Catalytic Linear Mixed Model (LMM) Fitting Function with fixed tau
#'
#' Fits a Catalytic linear mixed model (LMM) for observation and synthetic data with specified variance parameters
#' and iterative coefficient estimation. This function initializes model parameters,
#' sorts synthetic data, calculates Eigen-decomposition, and iterative optimizes
#' variance and coefficient values to convergence, by a single given tau value. (Only consider one random effect variance)
#'
#' @param cat_init A list generated from `cat_lmm_initialization`.
#' @param tau Optional numeric scalar controlling the weight of the synthetic data in the coefficient estimation, defaults to \code{ncol(cat_init$obs_x) / 4}.
#' @param residual_variance_0 Initial value for residual variance, default is 1.
#' @param random_effect_variance_0 Initial value for random effect variance, default is 1.
#' @param coefs_0 Optional initial coefficient vector, default is \code{NULL} which initializes randomly.
#' @param optimize_domain Numeric vector of length 2 defining optimization range for variance parameters, default is \code{c(0, 30)}.
#' @param max_iter Integer specifying maximum number of iterations for convergence, default is 500.
#' @param tol Tolerance for convergence criterion, default is 1e-08.
#'
#' @return A list containing the values of all the arguments and the following components:
#' \item{coefficients}{Estimated coefficient vector.}
#' \item{iteration_log}{Matrix logging variance and coefficient values for each iteration.}
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
#' cat_model <- cat_lmm(
#'   cat_init = cat_init, # Only accept object generated from cat_lmm_initialization
#'   tau = 1, # Weight for synthetic data
#'   residual_variance_0 = 1, # Initial value for residual variance
#'   random_effect_variance_0 = 1, # Initial value for random effect variance
#'   coefs_0 = c(1), # Initial coefficient vector
#'   optimize_domain = c(0, 10), # Optimization range for residual and random effect variance
#'   max_iter = 2, # Maximum number of iterations for convergence
#'   tol = 1e-01 # Tolerance for convergence criterion
#' )
#' cat_model
#' @export
cat_lmm <- function(
    cat_init,
    tau = NULL,
    residual_variance_0 = 1,
    random_effect_variance_0 = 1,
    coefs_0 = NULL,
    optimize_domain = c(0, 30),
    max_iter = 500,
    tol = 1e-08) {
  # Validate Input Parameters
  validate_lmm_input(
    cat_init = cat_init,
    tau = tau,
    residual_variance_0 = residual_variance_0,
    random_effect_variance_0 = random_effect_variance_0,
    coefs_0 = coefs_0,
    optimize_domain = optimize_domain,
    max_iter = max_iter,
    tol = tol
  )

  if (is.null(tau)) tau <- ncol(cat_init$obs_x) / 4

  # Sort synthetic data by group
  sorted_syn_data <- cbind(
    cat_init$syn_x,
    stats::setNames(as.data.frame(cat_init$syn_y), cat_init$y_col),
    cat_init$syn_z,
    stats::setNames(as.data.frame(cat_init$syn_group), cat_init$group_col)
  )[order(cat_init$syn_group), ]

  cat_init$syn_data <- sorted_syn_data
  cat_init$syn_x <- sorted_syn_data[, cat_init$x_cols, drop = FALSE]
  cat_init$syn_y <- sorted_syn_data[, cat_init$y_col]
  cat_init$syn_z <- sorted_syn_data[, cat_init$z_cols, drop = FALSE]
  cat_init$syn_group <- sorted_syn_data[[cat_init$group_col]]

  # Calculates the eigenvalues and eigenvectors of obs_z
  obs_z_eigen <- eigen(tcrossprod(as.matrix(cat_init$obs_z)))
  obs_z_eigenvalues <- obs_z_eigen$values
  obs_z_eigenvectors <- obs_z_eigen$vectors

  # Calculates the eigenvalues and eigenvectors of syb_z
  syn_z_eigenvalues <- c(rep(0, cat_init$syn_size))
  syn_z_eigenvectors <- matrix(0, ncol = cat_init$syn_size, nrow = cat_init$syn_size)
  row_offset <- 0

  for (g in unique(cat_init$group)) {
    group_idx <- which(cat_init$syn_group == g)
    syn_z_eigen <- eigen(tcrossprod(as.matrix(cat_init$syn_z[group_idx, ])))
    syn_z_eigenvalues[group_idx] <- syn_z_eigen$values
    syn_z_eigenvectors[((1:nrow(syn_z_eigen$vectors)) + row_offset), group_idx] <- syn_z_eigen$vectors

    row_offset <- row_offset + nrow(syn_z_eigen$vectors)
  }

  coefs <- if (is.null(coefs_0)) stats::rnorm(ncol(cat_init$obs_x)) else coefs_0
  residual_variance <- residual_variance_0
  random_effect_variance <- random_effect_variance_0

  iteration_log <- matrix(
    data = c(residual_variance, random_effect_variance, coefs),
    ncol = (2 + length(coefs))
  )

  colnames(iteration_log) <- c(
    "residual_variance",
    "random_effect_variance",
    cat_init$x_cols
  )

  iter <- 0
  parameters_converged <- FALSE

  while ((iter < max_iter) & (!parameters_converged)) {
    prev_coefs <- coefs
    prev_residual_variance <- residual_variance
    prev_random_effect_variance <- random_effect_variance

    obs_z_inv <- quadform::quad.tform(
      diag(1 / (residual_variance + random_effect_variance * obs_z_eigenvalues)),
      obs_z_eigenvectors
    )

    syn_z_inv <- quadform::quad.tform(
      diag(1 / (residual_variance + random_effect_variance * syn_z_eigenvalues)),
      syn_z_eigenvectors
    )

    # Update coefficients using generalized inverse and current parameter estimates
    coefs <- MASS::ginv(
      (tau / cat_init$syn_size) * quadform::quad.form(syn_z_inv, as.matrix(cat_init$syn_x)) +
        quadform::quad.form(obs_z_inv, as.matrix(cat_init$obs_x))
    ) %*% (
      quadform::quad3.form(obs_z_inv, as.matrix(cat_init$obs_x), cat_init$obs_y) +
        tau / cat_init$syn_size * quadform::quad3.form(syn_z_inv, as.matrix(cat_init$syn_x), cat_init$syn_y))

    obs_adjusted_residuals <- t(obs_z_eigenvectors) %*% (cat_init$obs_y - as.matrix(cat_init$obs_x) %*% coefs)
    syn_adjusted_residuals <- t(syn_z_eigenvectors) %*% (cat_init$syn_y - as.matrix(cat_init$syn_x) %*% coefs)

    random_effect_variance <- stats::optimize(update_lmm_variance,
      optimize_domain,
      residual_variance = prev_residual_variance,
      obs_z_eigenvalues = obs_z_eigenvalues,
      syn_z_eigenvalues = syn_z_eigenvalues,
      obs_adjusted_residuals = obs_adjusted_residuals,
      syn_adjusted_residuals = syn_adjusted_residuals,
      tau = tau
    )$minimum

    residual_variance <- stats::optimize(update_lmm_variance,
      optimize_domain,
      random_effect_variance = random_effect_variance,
      obs_z_eigenvalues = obs_z_eigenvalues,
      syn_z_eigenvalues = syn_z_eigenvalues,
      obs_adjusted_residuals = obs_adjusted_residuals,
      syn_adjusted_residuals = syn_adjusted_residuals,
      tau = tau
    )$minimum

    if ((norm(
      c(
        residual_variance,
        random_effect_variance
      ) - c(
        prev_residual_variance,
        prev_random_effect_variance
      ),
      "2"
    ) < tol) & (norm(coefs - prev_coefs, "2") < tol)) {
      parameters_converged <- TRUE
    }

    iter <- iter + 1
    iteration_log <- rbind(
      iteration_log,
      c(residual_variance, random_effect_variance, coefs)
    )
  }

  # Finalize Setup and Output
  cat_model <- list(
    function_name = "cat_lmm",
    ## Input/Processed parameters
    cat_init = cat_init,
    tau = tau,
    residual_variance_0 = residual_variance_0,
    random_effect_variance_0 = random_effect_variance_0,
    coefs_0 = coefs_0,
    optimize_domain = optimize_domain,
    max_iter = max_iter,
    tol = tol,
    ## Result
    coefficients = as.numeric(coefs),
    iteration_log = iteration_log
  )

  class(cat_model) <- c(cat_model$class, "cat", "cat_lmm")

  return(cat_model)
}

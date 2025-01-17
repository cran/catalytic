% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cat_glm_tune.R
\name{cat_glm_tune}
\alias{cat_glm_tune}
\title{Catalytic Generalized Linear Models (GLMs) Fitting Function by Tuning tau from a Sequence of tau Values}
\usage{
cat_glm_tune(
  formula,
  cat_init,
  risk_estimate_method = c("parametric_bootstrap", "cross_validation",
    "mallowian_estimate", "steinian_estimate"),
  discrepancy_method = c("mean_square_error", "mean_classification_error",
    "logistic_deviance"),
  tau_seq = NULL,
  tau_0 = NULL,
  parametric_bootstrap_iteration_times = 100,
  cross_validation_fold_num = 5
)
}
\arguments{
\item{formula}{A formula specifying the GLMs. Should at least include response variables (e.g. \code{~ .}).}

\item{cat_init}{A list generated from \code{cat_glm_initialization}.}

\item{risk_estimate_method}{Method for risk estimation, chosen from "parametric_bootstrap",
"cross_validation", "mallows_estimate", "steinian_estimate". Depends on the size of the data if not provided.}

\item{discrepancy_method}{Method for discrepancy calculation, chosen from "mean_square_error",
"mean_classification_error", "logistic_deviance". Depends on the family if not provided.}

\item{tau_seq}{Vector of numeric values for down-weighting synthetic data.
Defaults to a sequence around one fourth of the number of predictors for gaussian and
the number of predictors for binomial.}

\item{tau_0}{Initial \code{tau} value used for discrepancy calculation in risk estimation.
Defaults to one fourth of the number of predictors for binomial and 1 for gaussian.}

\item{parametric_bootstrap_iteration_times}{Number of bootstrap iterations for "parametric_bootstrap" risk estimation. Defaults to 100.}

\item{cross_validation_fold_num}{Number of folds for "cross_validation" risk estimation.. Defaults to 5.}
}
\value{
A list containing the values of all the arguments and the following components:
\item{tau}{Optimal \code{tau} value determined through tuning.}
\item{model}{Fitted GLM model object with the optimal \code{tau} value.}
\item{coefficients}{Estimated coefficients from the \code{model} fitted by the optimal \code{tau} value.}
\item{risk_estimate_list}{Collected risk estimates for each \code{tau}.}
}
\description{
This function tunes a catalytic catalytic Generalized Linear Models (GLMs) by performing specified risk estimate method
to estimate the optimal value of the tuning parameter \code{tau}.  The resulting \code{cat_glm_tune} object
encapsulates the fitted model, including estimated coefficients and family information, facilitating further analysis.
}
\examples{
gaussian_data <- data.frame(
  X1 = stats::rnorm(10),
  X2 = stats::rnorm(10),
  Y = stats::rnorm(10)
)

cat_init <- cat_glm_initialization(
  formula = Y ~ 1, # formula for simple model
  data = gaussian_data,
  syn_size = 100, # Synthetic data size
  custom_variance = NULL, # User customized variance value
  gaussian_known_variance = TRUE, # Indicating whether the data variance is known
  x_degree = c(1, 1), # Degrees for polynomial expansion of predictors
  resample_only = FALSE, # Whether to perform resampling only
  na_replace = stats::na.omit # How to handle NA values in data
)

cat_model <- cat_glm_tune(
  formula = ~.,
  cat_init = cat_init, # Only accept object generated from `cat_glm_initialization`
  risk_estimate_method = "parametric_bootstrap",
  discrepancy_method = "mean_square_error",
  tau_seq = c(1, 2), # Weight for synthetic data
  tau_0 = 2,
  parametric_bootstrap_iteration_times = 20, # Number of bootstrap iterations
  cross_validation_fold_num = 5 # Number of folds
)
cat_model
}

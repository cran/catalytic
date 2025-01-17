% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{validate_glm_initialization_input}
\alias{validate_glm_initialization_input}
\title{Validate Inputs for Catalytic Generalized Linear Models (GLMs) Initialization}
\usage{
validate_glm_initialization_input(
  formula,
  family,
  data,
  syn_size,
  custom_variance,
  gaussian_known_variance,
  x_degree
)
}
\arguments{
\item{formula}{A formula object specifying the \code{stats::glm} model to be fitted. It must not contain random effects or survival terms.}

\item{family}{A character or family object specifying the error distribution and link function. Valid values are "binomial" and "gaussian".}

\item{data}{A \code{data.frame} containing the data to be used in the GLM.}

\item{syn_size}{A positive integer specifying the sample size used for the synthetic data.}

\item{custom_variance}{A positive numeric value for the custom variance used in the model (only applicable for Gaussian family).}

\item{gaussian_known_variance}{A logical indicating whether the variance is known for the Gaussian family.}

\item{x_degree}{A numeric vector specifying the degree of the predictors. Its length should match the number of predictors (excluding the response variable).}
}
\value{
Returns nothing if all checks pass; otherwise, raises an error or warning.
}
\description{
This function validates the input parameters required for initializing a catalytic Generalized Linear Model (GLM).
It ensures the appropriate structure and compatibility of the formula, family, data, and additional parameters
before proceeding with further modeling.
}
\details{
This function performs the following checks:
\itemize{
\item Ensures that \code{syn_size}, \code{custom_variance}, and \code{x_degree} are positive values.
\item Verifies that the provided \code{formula} is suitable for GLMs, ensuring no random effects or survival terms.
\item Checks that the provided \code{data} is a \code{data.frame}.
\item Confirms that the \code{formula} does not contain too many terms relative to the number of columns in \code{data}.
\item Ensures that the \code{family} is either "binomial" or "gaussian".
\item Validates that \code{x_degree} has the correct length relative to the number of predictors in \code{data}.
\item Warns if \code{syn_size} is too small relative to the number of columns in \code{data}.
\item Issues warnings if \code{custom_variance} or \code{gaussian_known_variance} are used with incompatible families.
If any of these conditions are not met, the function raises an error or warning to guide the user.
}
}

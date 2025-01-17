% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cat_glm_tune.R
\name{parametric_bootstrap}
\alias{parametric_bootstrap}
\title{Perform Parametric Bootstrap for Model Risk Estimation}
\usage{
parametric_bootstrap(
  formula,
  cat_init,
  tau_seq,
  tau_0,
  discrepancy_method,
  parametric_bootstrap_iteration_times,
  ...
)
}
\arguments{
\item{formula}{A formula specifying the GLMs. Should at least include response variables.}

\item{cat_init}{A list generated from \code{cat_glm_initialization}.}

\item{tau_seq}{A sequence of tuning parameter values (\code{tau}) over which
the model risk will be estimated. Each \code{tau} value is used to weight the synthetic data during model fitting.}

\item{tau_0}{A reference value for \code{tau} used in the preliminary estimate model and variance calculation.}

\item{discrepancy_method}{The method used to calculate the discrepancy (e.g., logistic deviance).}

\item{parametric_bootstrap_iteration_times}{The number of bootstrap iterations to perform.}

\item{...}{Other arguments passed to other internal functions.}
}
\value{
A numeric vector containing the risk estimates for each \code{tau} in \code{tau_seq}.
}
\description{
This function performs parametric bootstrapping to estimate model risk. It fits a sequence
of Generalized Linear Models (GLMs) with different values of \code{tau}, calculates the in-sample
prediction error, and incorporates deviations from the bootstrap response samples. The final
risk estimate is obtained by combining the in-sample error and the covariance penalty derived
from the bootstrap samples.
}
\details{
\enumerate{
\item \strong{Preliminary Estimate Model}: The function first fits a GLM model using the observed
and synthetic data with an initial value of \code{tau_0} for the synthetic data weights.
\item \strong{Bootstrap Samples}: The function generates bootstrap response samples based on the
mean and standard deviation of the preliminary estimate model, using parametric bootstrapping.
\item \strong{In-sample Prediction Error}: For each value of \code{tau} in \code{tau_seq}, the function computes
the in-sample prediction error (e.g., using logistic deviance).
\item \strong{Bootstrap Models}: For each bootstrap iteration, the function fits a GLM using the
bootstrap response samples and calculates the corresponding lambda values.
\item \strong{Covariance Penalty}: The function approximates the covariance penalty using the weighted
deviations across all bootstrap iterations.
\item \strong{Final Risk Estimate}: The final model risk estimate is calculated by summing the in-sample
prediction error and the average weighted deviations from the bootstrap response samples.
}
}

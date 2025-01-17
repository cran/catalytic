% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{get_cox_partial_likelihood}
\alias{get_cox_partial_likelihood}
\title{Compute the Partial Likelihood for the Cox Proportional Hazards Model}
\usage{
get_cox_partial_likelihood(X, time, status, coefs, entry_points)
}
\arguments{
\item{X}{A matrix of covariates with rows representing observations and columns representing features.}

\item{time}{A vector of time-to-event data.}

\item{status}{A vector indicating the event status (1 for event occurred, 0 for censored).}

\item{coefs}{A vector of regression coefficients.}

\item{entry_points}{A vector of entry points (optional). Defaults to NULL, in which case a vector of zeros is used.}
}
\value{
A numeric scalar representing the partial likelihood of the Cox model.
}
\description{
This function calculates the partial likelihood for the Cox proportional hazards model. The partial likelihood is computed for the censored observations in the dataset.
}

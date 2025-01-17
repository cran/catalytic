% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{get_cox_hessian}
\alias{get_cox_hessian}
\title{Compute the Hessian Matrix for Cox Proportional Hazards Model}
\usage{
get_cox_hessian(X, time, status, coefs, entry_points)
}
\arguments{
\item{X}{A matrix of covariates (design matrix) for the Cox model.}

\item{time}{A numeric vector of event times.}

\item{status}{A numeric vector of event indicators (1 for event, 0 for censored).}

\item{coefs}{A numeric vector of coefficients for the Cox model.}

\item{entry_points}{A numeric vector of entry times for the subjects. Defaults to 0.}
}
\value{
A matrix representing the negative Hessian of the Cox model.
}
\description{
This function computes the Hessian matrix of the Cox proportional hazards model,
which is used for estimating the covariance matrix of the coefficients. The Hessian is calculated
by summing contributions from each event time in the risk set.
}

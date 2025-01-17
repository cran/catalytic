% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{get_cox_syn_gradient}
\alias{get_cox_syn_gradient}
\title{Compute the gradient of the synthetic Cox proportional hazards model}
\usage{
get_cox_syn_gradient(X, time, coefs, hazard_constant)
}
\arguments{
\item{X}{A matrix of covariates with rows representing observations and columns representing features.}

\item{time}{A vector of time-to-event data.}

\item{coefs}{A vector of regression coefficients.}

\item{hazard_constant}{A scalar representing the hazard constant.}
}
\value{
A numeric vector representing the gradient of the synthetic Cox model.
}
\description{
This function calculates the gradient of the synthetic Cox proportional hazards model using a vectorized approach.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{get_cox_risk_set_idx}
\alias{get_cox_risk_set_idx}
\title{Identify the risk set indices for Cox proportional hazards model}
\usage{
get_cox_risk_set_idx(
  time_of_interest,
  entry_vector,
  time_vector,
  status_vector
)
}
\arguments{
\item{time_of_interest}{A numeric value representing the time at which the risk set is calculated.}

\item{entry_vector}{A numeric vector representing the entry times of subjects.}

\item{time_vector}{A numeric vector representing the time-to-event or censoring times of subjects.}

\item{status_vector}{A numeric vector indicating event occurrence (1) or censoring (0) for each subject.}
}
\value{
A vector of indices representing the subjects at risk at the specified time of interest.
}
\description{
This function returns the indices of the risk set for a given time of interest in the Cox proportional hazards model.
}

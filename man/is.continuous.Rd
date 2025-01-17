% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{is.continuous}
\alias{is.continuous}
\title{Check if a Variable is Continuous}
\usage{
is.continuous(lst)
}
\arguments{
\item{lst}{A vector to be checked.}
}
\value{
A logical value indicating whether the input vector is considered continuous. Returns \code{TRUE} if the vector is numeric and has more than two unique values; otherwise, returns \code{FALSE}.
}
\description{
This function checks whether a given vector represents a continuous variable. A continuous variable is numeric and has more than two unique values.
}

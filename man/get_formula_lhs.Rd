% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{get_formula_lhs}
\alias{get_formula_lhs}
\title{Extract Left-Hand Side of Formula as String}
\usage{
get_formula_lhs(formula)
}
\arguments{
\item{formula}{A formula object from which the LHS will be extracted.}
}
\value{
A character string representing the left-hand side of the formula.
}
\description{
This function extracts the left-hand side (LHS) of a formula object and converts it to a character string. It uses \code{get_formula_string} to ensure consistent formatting.
}

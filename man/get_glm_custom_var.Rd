% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{get_glm_custom_var}
\alias{get_glm_custom_var}
\title{Get Custom Variance for Generalized Linear Model (GLM)}
\usage{
get_glm_custom_var(formula, cat_init, tau)
}
\arguments{
\item{formula}{A formula object specifying the GLM model to be fitted, such as \code{response ~ predictors}.}

\item{cat_init}{A list object containing the initialization data for the model. Generated from \code{cat_initialization}}

\item{tau}{A numeric value representing a scaling factor for down-weighting synthetic data}
}
\value{
A numeric value representing the custom variance for the GLM model.
}
\description{
This function calculates a custom variance for a Generalized Linear Model (GLM) based on the specified formula,
the model initialization object, and a scaling factor \code{tau}. The custom variance is computed by adjusting
the residuals of the fitted model and returning a weighted sum of squared residuals.
}

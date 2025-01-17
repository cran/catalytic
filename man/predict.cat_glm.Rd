% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.cat_glm}
\alias{predict.cat_glm}
\title{Predict Outcome for New Data Using a Fitted GLM Model}
\usage{
\method{predict}{cat_glm}(object, newdata = NULL, ...)
}
\arguments{
\item{object}{A fitted model object of class \code{cat_glm}, containing the GLM fit and model details.}

\item{newdata}{An optional data frame containing new predictor values. If \code{NULL},
the function uses the observation data from the model's initialization object.}

\item{...}{Additional arguments passed to \code{stats::predict.glm}, if applicable.
User should input \code{type = c("link", "response", "terms")} for different regression models.}
}
\value{
A vector of predicted values for the specified new data.
}
\description{
This function generates predictions for new data points based on a fitted categorical Generalized Linear Model (GLM) object.
Depending on the type of model, it either uses \code{stats::predict.glm} or calculates predictions based on the model coefficients.
}

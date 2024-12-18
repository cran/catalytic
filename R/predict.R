#' Predict Outcome for New Data Using a Fitted GLM Model
#'
#' This function generates predictions for new data points based on a fitted categorical Generalized Linear Model (GLM) object.
#' Depending on the type of model, it either uses `stats::predict.glm` or calculates predictions based on the model coefficients.
#'
#' @param object A fitted model object of class `cat_glm`, containing the GLM fit and model details.
#' @param newdata An optional data frame containing new predictor values. If `NULL`,
#' the function uses the observation data from the model's initialization object.
#' @param ... Additional arguments passed to `stats::predict.glm`, if applicable.
#' User should input `type = c("link", "response", "terms")` for different regression models.
#'
#' @return A vector of predicted values for the specified new data.
#' @export
predict.cat_glm <- function(
    object,
    newdata = NULL,
    ...) {
  newdata <- if (is.null(newdata)) object$cat_init$adj_obs_x else newdata[, colnames(object$cat_init$adj_x), drop = FALSE]

  if (object$function_name %in% c("cat_glm", "cat_glm_tune")) {
    return(
      stats::predict.glm(
        object$model,
        newdata = as.data.frame(newdata),
        ...
      )
    )
  }

  return(
    c(get_glm_mean(
      family_string = object$cat_init$family,
      X = newdata,
      coefs = stats::coef(object)
    ))
  )
}

#' Predict Linear Predictor for New Data Using a Fitted Cox Model
#'
#' This function calculates the linear predictor (LP) for new data points based on a fitted Cox proportional hazards model.
#'
#' @param object A fitted model object of class `cat_cox`, containing the COX fit and model details.
#' @param newdata An optional data frame with new predictor values. If `NULL`,
#' the function uses the observation data from the model's initialization object.
#' @param ... Additional arguments passed to other functions.
#'
#' @return A vector of linear predictor values for the specified new data.
#' @export
predict.cat_cox <- function(
    object,
    newdata = NULL,
    ...) {
  # Prepare the covariates
  newdata <- if (is.null(newdata)) object$cat_init$adj_obs_x else newdata[, colnames(object$cat_init$adj_x), drop = FALSE]

  return(c(get_linear_predictor(
    X = newdata,
    coefs = stats::coef(object)
  )))
}

#' Predict Linear Predictor for New Data Using a Fitted Linear Mixed Model
#'
#' This function calculates the linear predictor (LP) for new data points based on a fitted linear mixed model (LMM) stored in `object`.
#'
#' @param object A fitted model object of class `cat_lmm`, containing the LMM fit and model details.
#' @param newdata An optional data frame with new predictor values. If `NULL`,
#' the function uses the observation data from the model's initialization object.
#' @param ... Additional arguments passed to other functions.
#'
#' @return A vector of linear predictor values for the specified new data.
#' @export
predict.cat_lmm <- function(
    object,
    newdata = NULL,
    ...) {
  # Prepare the covariates
  newdata <- if (is.null(newdata)) object$cat_init$obs_x else newdata[, colnames(object$cat_init$x), drop = FALSE]

  return(c(get_linear_predictor(
    X = newdata,
    coefs = stats::coef(object)
  )))
}

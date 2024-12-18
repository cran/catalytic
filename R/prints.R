#' Print Summary of `cat_gibbs` Model
#'
#' This function prints a summary of the `cat_gibbs` model, displaying details about the formula, covariate dimensions,
#' family, coefficients, and Gibbs sampling settings.
#'
#' @param x A `cat_gibbs` model object containing the results of a Bayesian GLM fitted using Gibbs sampling.
#' @param digit An integer indicating the number of decimal places for printing
#'   coefficient estimates. Default is 3.
#' @param detail A logical value indicating whether to include additional detailed output at the end of the summary.
#' If `TRUE`, it will print additional interpretation help.
#' @param ... Additional parameters to pass to other functions.
#'
#' @details
#' The summary includes:
#' - The function name and formula used in the model.
#' - Dimensions of the covariate matrix.
#' - Family and link function details.
#' - Sampling information, including the total iterations, warm-up iterations, and effective Gibbs sampling post-warmup.
#' - Coefficients with summary statistics and effective sample size.
#'
#' If `detail` is set to TRUE, additional guidance for interpreting the printed output is provided.
#'
#' @return The `x` object is returned invisibly.
#' @export
print.cat_gibbs <- function(
    x,
    digit = 3,
    detail = TRUE,
    ...) {
  cat(x$function_name)

  cat("\n formula:                ", get_formula_string(x$formula))
  cat("\n covariates dimention:   ", extract_dim(x$cat_init))
  cat("\n family:                 ", get_glm_family_string(x$cat_init$family, with_link = TRUE))

  cat("\n------\n")
  cat("coefficients' information:\n\n")
  cat(paste0(
    "Inference for Gibb's sampling.\n",
    "gibbs_iter=", x$iter,
    "; gibbs_warmup=", x$warmup, "; ",
    "gibb's sampling post=", x$iter - x$warmup,
    ".\n\n"
  ))

  print(round(x$inform_df, digit))

  cat(paste0(
    "\nSamples were drawn using Gibb's Sampling at ",
    formatted_time <- format(x$sys_time, "%a %b %d %H:%M:%S %Y"),
    ".\nFor each parameter, n_eff is a crude measure of effective sample size."
  ))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_gibbs\n")
  }

  invisible(x)
}

#' Print Summary of `cat_bayes` Model
#'
#' This function prints a formatted summary of a `cat_bayes` model object, displaying
#' key parameters and settings of the fitted model, including the formula,
#' covariate dimensions, tau (if applicable), family, and algorithm settings,
#' as well as the coefficients' summary.
#'
#' @param x An object of class `cat_tune`, typically resulting from a tuning process, including
#' `cat_glm_bayes`, `cat_glm_bayes_joint`, `cat_cox_bayes` and `cat_cox_bayes_joint`.
#' @param digit An integer indicating the number of decimal places for printing
#'   coefficient estimates. Default is 3.
#' @param detail A logical value indicating whether to include additional detailed output at the end of the summary.
#' If `TRUE`, it will print additional interpretation help.
#' @param ... Additional parameters to pass to other functions.
#'
#' @details
#' This function provides an organized printout of essential details from a Bayesian
#' model fit. It includes the model formula, dimensionality of covariates, model
#' family, Stan algorithm settings, and summary of the model coefficients. If
#' \code{detail} is set to \code{TRUE}, additional information on interpreting the
#' output is included.
#'
#' @return The `x` object is returned invisibly.
#' @export
print.cat_bayes <- function(
    x,
    digit = 3,
    detail = TRUE,
    ...) {
  cat(x$function_name)

  cat("\n formula:               ", get_formula_string(x$formula))
  cat("\n covariates dimention:  ", extract_dim(x$cat_init))
  if (!grepl("joint", x$function_name)) {
    cat("\n tau:                   ", x$tau)
  }
  if (grepl("glm", x$function_name)) {
    cat("\n family:                ", get_glm_family_string(x$cat_init$family, with_link = TRUE))
    cat("\n stan algorithm:        ", x$algorithm)
  }
  cat("\n stan chain:            ", x$chain)
  cat("\n stan iter:             ", x$iter)
  cat("\n stan warmup:           ", x$warmup)

  cat("\n------\n")
  cat("coefficients' information:\n\n")
  print(extract_stan_summary(x,
    digit,
    with_intercept = (!grepl("cox", x$function_name))
  ))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_bayes\n")
  }

  invisible(x)
}

#' Print Method for `cat_tune` Object
#'
#' This function prints a summary of the `cat_tune` object, displaying key details such as
#' the function name, dimensions of covariates, tau sequence, optimal tau, likelihood or risk
#' estimate, and the model's coefficients.
#'
#' @param x An object of class `cat_tune`, typically resulting from a tuning process, including
#' `cat_glm_tune`, `cat_cox_tune` and `cat_lmm_tune`.
#' @param digit An integer indicating the number of decimal places for printing
#'   coefficient estimates. Default is 3.
#' @param detail A logical value indicating whether to include additional detailed output at the end of the summary.
#' If `TRUE`, it will print additional interpretation help.
#' @param ... Additional parameters to pass to other functions.
#'
#' @details
#' This method provides a comprehensive overview of the tuning process for the model,
#' including the tau sequence and optimal tau, along with either the maximum likelihood (for Cox models)
#' or minimum risk estimate (for other models). It also displays the coefficients of the model.
#'
#' The function also checks if the `x` is a Cox model (`cat_cox_tune`) to adjust the interpretation
#' of the output.
#'
#' @return The `x` object is returned invisibly.
#' @export
print.cat_tune <- function(
    x,
    digit = 3,
    detail = TRUE,
    ...) {
  cat(x$function_name)

  cat("\n formula:                ", get_formula_string(x$formula))
  cat("\n covariates dimention:   ", extract_dim(x$cat_init))
  cat("\n tau sequnce:            ", extract_tau_seq(x, digit))
  if (grepl("glm", x$function_name)) {
    cat("\n family:                 ", x$cat_init$family)
    cat("\n risk estimate method:   ", x$risk_estimate_method)
    cat("\n discrepancy method:     ", x$discrepancy_method)
  }
  cat("\n")
  cat("\n optimal tau:            ", x$tau)
  if (grepl("cox", x$function_name)) {
    cat("\n method:                 ", x$method)
    cat("\n maximun likelihood:     ", round(max(x$likelihood_list), digit))
  } else {
    cat("\n minimun risk estimate:  ", round(min(x$risk_estimate_list), digit))
  }

  cat("\n------\n")
  cat("coefficients' information:\n\n")

  print(extract_coefs(
    x = x,
    digit = digit
  ))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_tune\n")
  }

  invisible(x)
}

#' Print Method for `cat` Object
#'
#' The `print.cat` function provides a detailed summary of the `cat` object, displaying key information about the model and its settings,
#' including model type, covariates, formula, tau values, and relevant coefficients.
#'
#' This function customizes the output based on the model type stored within the `x` object, such as GLM, Cox, or other types of models.
#'
#' @param x An object of class `cat`, representing a fitted model.
#' @param digit An integer indicating the number of decimal places for printing
#'   coefficient estimates. Default is 3.
#' @param detail A logical indicating whether to print additional details for interpreting the output. Default is TRUE.
#' @param ... Additional parameters to pass to other functions.
#'
#' @details
#' The `print.cat` function prints a summary of the model stored in the `x` object. It will display different information
#' depending on the model's type (GLM, Cox, etc.). It will show:
#' - The model's function name.
#' - The dimensions of the covariates used in the model.
#' - The tau values.
#' - Model-specific details such as family for GLMs or method and iteration info for Cox models.
#' - Coefficients related to the model.
#'
#' @return The `x` object is returned invisibly.
#' @export
print.cat <- function(
    x,
    digit = 3,
    detail = TRUE,
    ...) {
  cat(x$function_name)

  cat("\n formula:               ", get_formula_string(x$formula))
  cat("\n covariates dimention:  ", extract_dim(x$cat_init))
  cat("\n tau:                   ", x$tau)

  if (grepl("glm", x$function_name)) {
    cat("\n family:                ", get_glm_family_string(x$cat_init$family, with_link = TRUE))
  } else if (grepl("cox", x$function_name)) {
    cat("\n method:                ", x$method)
    if (x$method == "CRE") {
      cat(
        "\n iter:                  ",
        paste(
          x$iter,
          "/",
          x$max_iter
        )
      )
    }
  } else {
    cat("\n iteration:             ", nrow(x$iteration_log) - 1, "/", x$max_iter)
    cat("\n tolerance:             ", x$tol)
    cat("\n optimize domain:       ", paste0("[", paste(x$optimize_domain, collapse = ", "), "]"))
  }
  cat("\n------\n")
  cat("coefficients' information:\n")

  print(extract_coefs(
    x = x,
    digit = digit
  ))

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat\n")
  }

  invisible(x)
}

#' Print Summary for Catalytic Initialization Model
#'
#' This function provides a comprehensive summary of a catalytic initialization model object (`cat_init`),
#' including formula details, data dimensions, and sample data previews.
#'
#' @param x A catalytic initialization model object containing formula, family, data dimensions, and sampling details.
#' @param show_data Logical, default `TRUE`. If `TRUE`, previews the head of both observation and synthetic data (up to the first 10 columns).
#' @param detail Logical, default `TRUE`. If `TRUE`, adds guidance for interpreting the output.
#' @param ... Additional parameters to pass to other functions.
#'
#' @details
#' The function provides a detailed overview of the initialization process for the `cat_initialization` model,
#' including:
#' - The formula used in the model.
#' - The type of model (if Gaussian, the known or unknown variance is specified).
#' - The family of the Generalized Linear Model (GLM), along with the associated link function.
#' - The dimensions of the observation and synthetic data sets, with an option to display the first few rows.
#' - Information on the data generation process if available.
#'
#' The `show_data` parameter controls whether the first few rows of the data are printed, while the `detail`
#' parameter controls whether additional help for interpreting the printed output is displayed.
#'
#' @return Invisibly returns the `x` object.
#' @export
print.cat_initialization <- function(
    x,
    show_data = TRUE,
    detail = TRUE,
    ...) {
  cat(x$function_name)

  cat("\n formula:              ", get_formula_string(x$formula))
  if (grepl("glm", x$function_name)) {
    if (x$family == "gaussian") {
      if (x$gaussian_known_variance) {
        model_type <- if (is.null(x$custom_variance)) {
          "Known Variance: Generates variance based on data and function."
        } else {
          "Known Variance: Uses the specified `custom_variance`."
        }
      } else {
        model_type <- if (!is.null(x$custom_variance)) {
          "Unknown Variance: The input `custom_variance` will be ignored."
        } else {
          "Unknown Variance"
        }
      }
      cat("\n model:                ", model_type)
      cat("\n custom variance:      ", ifelse((is.null(x$custom_variance)),
        "NULL",
        x$custom_variance
      ))
    }

    cat("\n family:               ", get_glm_family_string(x$family, with_link = TRUE))
  }
  cat("\n covariates dimention: ", extract_dim(x))

  cat("\n------\n")
  cat("Observation Data Information:")
  cat("\n covariates dimention: ", x$obs_size, "rows with ", ncol(x$obs_x), "column(s)")
  if (show_data) {
    cat("\n head(data) :              \n")
    cat(" [only show the head of first 10 columns] \n\n")
    if (ncol(x$obs_data) > 10) {
      print(utils::head(x$obs_data)[, 1:10])
    } else {
      print(utils::head(x$obs_data))
    }
  }

  cat("\n------\n")

  cat("Synthetic Data Information:")
  cat("\n covariates dimention: ", x$syn_size, "rows with ", ncol(x$syn_x), " column(s)")
  if (show_data) {
    cat("\n head(data):              \n")
    cat(" [only show the head of first 10 columns] \n\n")
    if (ncol(x$syn_data) > 10) {
      print(utils::head(x$syn_data)[, 1:10])
    } else {
      print(utils::head(x$syn_data))
    }
  }

  if (!is.null(x$syn_x_resample_inform) || !is.null(x$syn_z_resample_inform)) {
    print_data_gen_pro_df <- rbind(
      x$syn_x_resample_inform,
      x$syn_z_resample_inform
    )

    cat("\n data generation process:")
    cat("\n [only show the first 10 columns] \n\n")

    print(
      format(
        if (nrow(print_data_gen_pro_df) > 10) {
          print_data_gen_pro_df[1:10, ]
        } else {
          print_data_gen_pro_df
        },
        justify = "left",
        width = max(nchar(as.character(unlist(print_data_gen_pro_df)))) + 1
      ),
      quote = FALSE
    )
  }

  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.cat_initialization\n")
  }

  invisible(x)
}



# ----- Helper -----
#' Extract Dimension Information from Model Initialization
#'
#' This function retrieves and formats the dimensions of the dataset used in the model,
#' including the number of observed and synthetic data points and the total number of rows
#' and columns.
#'
#' @param cat_init A list containing model initialization data, expected to include
#' `obs_size` (observed data size), `syn_size` (synthetic data size), `size` (total data size),
#' and `x` (the covariate matrix).
#'
#' @return A character string summarizing the dimensions of the dataset used in the model.
extract_dim <- function(cat_init) {
  return(paste0(
    cat_init$obs_size,
    " (Observation) + ",
    cat_init$syn_size,
    " (Synthetic) = ",
    cat_init$size,
    " rows with ",
    ncol(cat_init$x),
    " column(s)"
  ))
}

#' Extract and Format Model Coefficients
#'
#' This function retrieves the coefficients from a `x` object, formats them
#' with appropriate names, and rounds each coefficient to the specified number of decimal places.
#' Optionally, the intercept can be included or excluded from the output.
#'
#' @param x A model object generated from `catalytic` that containing model coefficients.
#' @param digit An integer specifying the number of decimal places for rounding coefficients. Default is 3.
#'
#' @return A named numeric vector of model coefficients, rounded to the specified number of decimal places.
extract_coefs <- function(
    x,
    digit = 3) {
  print_coef <- stats::coef(x)
  col_names <- if (is.null(x$cat_init$adj_x)) colnames(x$cat_init$x) else colnames(x$cat_init$adj_x)

  if (length(print_coef) == length(col_names)) {
    names(print_coef) <- col_names
  } else {
    names(print_coef) <- c("(Intercept)", col_names)
  }

  return(round(print_coef, digit))
}

#' Extract and Format Sequence of Tau Values
#'
#' This function retrieves the sequence of tau values from a `x` object, rounds each
#' value to the specified number of decimal places, and formats the output as a concise string.
#' If the sequence contains more than 10 values, only the first 3 and last 3 values are shown,
#' with ellipsis ("...") in between.
#'
#' @param x A model object generated from `catalytic` that containing a sequence of tau values.
#' @param digit An integer specifying the number of decimal places to which tau values should
#' be rounded. Default is 3.
#'
#' @return A character string representing the rounded tau values, formatted for readability.
extract_tau_seq <- function(
    x,
    digit = 3) {
  tau_seq <- round(x$tau_seq, digit)
  return(if (length(tau_seq) > 10) {
    paste(
      paste(tau_seq[1:3], collapse = ", "),
      "...",
      paste(
        tau_seq[(length(tau_seq) - 2):length(tau_seq)],
        collapse = ", "
      )
    )
  } else {
    paste(tau_seq, collapse = ", ")
  })
}

#' Extract and Format Summary of Stan Model Results
#'
#' This function extracts the summary statistics from a fitted Stan model stored within
#' a `x` object, formats the parameter names, and rounds values to a specified number
#' of decimal places. By default, the function includes an intercept term in the summary if
#' present.
#'
#' @param x A model object generated from `catalytic` that containing a fitted Stan model.
#' @param digit An integer specifying the number of decimal places to which the summary
#' statistics should be rounded. Default is 3.
#' @param with_intercept A logical value indicating whether the intercept should be included
#' in the summary. If `TRUE`, the intercept is labeled and included in the formatted output.
#' Default is `TRUE`.
#'
#' @return A matrix of rounded summary statistics from the Stan model, with row names
#' representing parameter labels and columns containing summary values.
extract_stan_summary <- function(x,
                                 digit = 3,
                                 with_intercept = TRUE) {
  stan_summary <- rstan::summary(x$stan_sample_model)$summary

  # Remove h_j from cat_cox_bayes or cat_cox_bayes_joint while printing
  stan_summary <- stan_summary[!grepl("h_j", rownames(stan_summary)), ]

  if (with_intercept) {
    rownames(stan_summary)[1:(ncol(x$cat_init$adj_x) + 1)] <- c(
      "(Intercept)",
      colnames(x$cat_init$adj_x)
    )
  } else {
    rownames(stan_summary)[1:(ncol(x$cat_init$adj_x))] <- colnames(x$cat_init$adj_x)
  }

  return(round(stan_summary, digit))
}

#' Print Data Frame with Head and Tail Rows
#'
#' This function displays the first 5 and last 5 rows of a data frame.
#' Column names are displayed only for the first 5 rows, with ellipses (`...`)
#' in the middle to indicate additional rows.
#'
#' @param df A data frame to display.
#' @param digit An integer specifying the number of decimal places to which the summary
#' statistics should be rounded. Default is 3.
#'
#' @return Invisibly returns the original data frame.
print_df_head_tail <- function(df, digit = 3) {
  total_rows <- nrow(df)

  if (total_rows <= 10) {
    # If 10 or fewer rows, print the whole data frame
    print(df)
  } else {
    # Print first 5 rows
    print(utils::head(df, 5), digit = digit)

    # Print ellipses
    cat("...\n")

    # Print last 5 rows
    print(utils::tail(df, 5), digit = digit)
  }

  # Return the data frame invisibly
  invisible(df)
}

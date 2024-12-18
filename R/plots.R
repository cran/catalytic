#' Plot Likelihood or Risk Estimate vs. Tau for Tuning Model
#'
#' This function generates a plot showing the relationship between
#' the tuning parameter \code{tau} and either the likelihood score (for a \code{cat_cox_tune} model)
#' or the risk estimate (for other models) during cross-validation or other model evaluation methods.
#' The plot highlights the optimal \code{tau} value and provides visual cues for
#' the best tuning parameter based on the specified method.
#'
#' @param x A fitted model object of class `cat_tune` that contains the
#' results of the tuning process. This object includes the likelihood or risk estimate lists,
#' the tuning sequence (\code{tau_seq}), and the selected optimal \code{tau}.
#' @param digit An integer specifying the number of decimal places to round the
#' displayed values (default is 2).
#' @param legend_pos A character string specifying the position of the legend on
#' the plot (default is \code{"topright"}).
#' @param text_pos An integer specifying the position of the text label on the
#' plot (default is 3, which places the text above the point).
#' @param ... Additional parameters to pass to other functions.
#'
#' @details
#' The function generates a line plot with \code{tau_seq} on the x-axis and either
#' the likelihood score or risk estimate on the y-axis. If the model is of class
#' \code{cat_cox_tune}, the plot shows the likelihood score, while for other models,
#' it shows the risk estimates. The optimal \code{tau} is marked with a red cross,
#' and red dashed lines are drawn to highlight the optimal point on the plot.
#'
#' @return A plot with the specified y-values plotted against \code{tau_seq},
#' including a highlighted optimal \code{tau} point.
#' @export
plot.cat_tune <- function(x,
                          digit = 2,
                          legend_pos = "topright",
                          text_pos = 3,
                          ...) {
  if (x$function_name == "cat_cox_tune") {
    y_lab <- "likelihood score"
    y_plot <- x$likelihood_list
    main <- paste0("Plot of Likelihood Scores from cat_cox_tune \n(cross_validation)")
    mark_point <- max(x$likelihood_list)
  } else {
    y_lab <- paste0("discrepancy error (", x$discrepancy_method, ")")
    y_plot <- x$risk_estimate_list
    main <- paste0(
      "Plot of Risk Estimates from ",
      x$function_name,
      "\n(",
      ifelse(is.null(x$risk_estimate_method),
        "cross_validation",
        x$risk_estimate_method
      ),
      ")"
    )
    mark_point <- min(x$risk_estimate_list)
  }

  plot(x$tau_seq,
    y_plot,
    type = "l",
    col = "black",
    lwd = 2,
    xlab = "tau",
    ylab = y_lab,
    main = main
  )

  # Add text labels to the highlighted point
  graphics::text(x$tau,
    mark_point,
    labels = paste0("(", round(x$tau, 1), ", ", round(mark_point, 1), ")"),
    pos = text_pos
  )

  # Draw red dotted lines from the axes to the highlighted point
  graphics::segments(x0 = graphics::par("usr")[1], y0 = mark_point, x1 = x$tau, y1 = mark_point, col = "red", lty = 2) # Horizontal line
  graphics::segments(x0 = x$tau, y0 = graphics::par("usr")[3], x1 = x$tau, y1 = mark_point, col = "red", lty = 2) # Vertical line

  graphics::points(x$tau,
    mark_point,
    pch = 4,
    col = "red"
  )

  graphics::legend(legend_pos,
    legend = "Optimal tau",
    col = "red",
    pch = 4,
    text.col = "black"
  )
}

#' Traceplot for Bayesian Sampling Model
#'
#' This function generates a traceplot for the Bayesian sampling model fitted
#' using `rstan`. It utilizes the `traceplot` function from the `rstan` package
#' to visualize the sampling progress and convergence of the Markov Chain Monte Carlo (MCMC) chains
#' for the given model.
#'
#' @param object A fitted model object of class `cat_bayes` that contains
#' the Stan sampling model. The object should include the `stan_sample_model`,
#' which is the result of fitting the model using the `rstan` package.
#' @param ... Additional arguments passed to the `rstan::traceplot` function.
#' These can include customization options for the traceplot, such as `pars`
#' for selecting specific parameters or `inc_warmup` for including or excluding warmup iterations.
#'
#' @details
#' The function calls `rstan::traceplot` on the `stan_sample_model` contained
#' within the `x` object. The resulting plot displays the trace of each parameter
#' across MCMC iterations, helping to assess the convergence and mixing of the chains.
#'
#' @return A traceplot displaying the MCMC chains' trace for each parameter,
#' helping to assess convergence.
#' @export
traceplot.cat_bayes <- function(object, ...) {
  rstan::traceplot(object$stan_sample_model, ...)
}

#' Traceplot for Gibbs Sampling Model
#'
#' This function generates a traceplot for the Gibbs sampling model, which is
#' typically used for posterior sampling in a Bayesian context. The traceplot
#' visualizes the evolution of parameter values across Gibbs sampling iterations.
#' It helps to diagnose the convergence and mixing of the chains.
#'
#' @param object A fitted model object of class `cat_gibbs` that contains
#' Gibbs sampling results. The object must include `gibbs_iteration_log`,
#' which holds the iteration logs for all sampled parameters, and `warmup` and `iter`
#' which indicate the warmup and total iteration counts, respectively.
#' @param pars A character vector specifying the parameter names to plot.
#' If `NULL`, the function will select the first 9 parameters automatically.
#' @param inc_warmup A logical value indicating whether to include warmup iterations
#' in the traceplot. If `TRUE`, warmup iterations are included, otherwise they are excluded.
#' Defaults to `FALSE`.
#' @param ... Additional parameters to pass to other functions.
#'
#' @details
#' The function generates a series of line plots for the selected parameters,
#' displaying their values over the iterations of the Gibbs sampling process.
#' If `inc_warmup` is set to `TRUE`, the traceplot includes the warmup period,
#' otherwise, it starts after the warmup. The traceplots are arranged in a 3x3 grid,
#' and no more than 9 parameters can be selected for plotting at once.
#'
#' @return A series of traceplots for the selected parameters,
#' showing their evolution over the Gibbs sampling iterations.
#' @export
traceplot.cat_gibbs <- function(
    object,
    pars = NULL,
    inc_warmup = FALSE,
    ...) {
  oldpar <- graphics::par(no.readonly = TRUE) # Save original graphical parameters
  on.exit(graphics::par(oldpar)) # Ensure they are restored on exit

  print_df <- if (inc_warmup) object$gibbs_iteration_log else object$gibbs_iteration_log[object$warmup:object$iter, ]

  if (is.null(pars)) {
    pars <- colnames(print_df)[1:min(9, ncol(print_df))]
  }

  if (length(pars) > 9) {
    stop("Please select less than 10 features to plot.")
  }

  graphics::par(mar = c(2, 2, 2, 2))
  graphics::par(mfrow = c(3, 3))

  for (col_name in pars) {
    plot(
      x = if (inc_warmup) seq(1, object$iter) else seq(object$warmup, object$iter),
      y = print_df[, col_name],
      main = col_name,
      xlab = "",
      ylab = "",
      type = "l"
    )
  }
}

#' Traceplot for Bayesian Model Sampling
#'
#' The `traceplot` function is a generic function used to generate traceplots for Bayesian model sampling,
#' primarily for assessing the convergence and mixing of Markov Chain Monte Carlo (MCMC) chains. This function
#' dispatches specific traceplot methods depending on the class of the `object` object.
#'
#' @param object An object representing a Bayesian model, typically generated by the
#' \code{cat_glm_bayes} or \code{cat_cox_bayes} functions, or similar models with Bayesian sampling results.
#' The function uses S3 method dispatch to apply the appropriate \code{traceplot} method based on
#' the class of \code{object}.
#' @param ... Additional arguments passed to specific \code{traceplot} methods for customization, such as
#' selecting parameters to plot or setting display options.
#'
#' @details
#' This generic \code{traceplot} function allows for flexible visualization of MCMC chains across
#' different types of Bayesian models. Specific \code{traceplot} methods, such as \code{traceplot.cat_bayes},
#' are dispatched based on the \code{object} class to produce tailored traceplots, providing insights
#' into the sampling progress and convergence diagnostics of each chain.
#'
#' @return A traceplot displaying the MCMC sampling chains for each parameter, assisting in convergence analysis.
#' The exact output format depends on the specific \code{traceplot} method applied.
#' @export
traceplot <- function(object, ...) {
  UseMethod("traceplot")
}

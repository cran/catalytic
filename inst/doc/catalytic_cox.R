## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 7,
  fig.height = 6,
  cache = TRUE
)

## ----eval=FALSE---------------------------------------------------------------
#  library(survival)
#  library(catalytic)
#  
#  set.seed(1)
#  
#  # Compute the Partial Likelihood for the Cox Proportional Hazards Model
#  get_cox_partial_likelihood <- function(X,
#                                         time,
#                                         status,
#                                         coefs,
#                                         entry_points) {
#    # Identify the risk set indices for Cox proportional hazards model
#    get_cox_risk_set_idx <- function(time_of_interest,
#                                     entry_vector,
#                                     time_vector,
#                                     status_vector) {
#      time_vector <- as.vector(time_vector)
#      entry_vector <- as.vector(entry_vector)
#      status_vector <- as.vector(status_vector)
#  
#      # Find indices where subjects are at risk at the given time of interest
#      return(which((time_of_interest >= entry_vector) & (
#        (time_vector == time_of_interest & status_vector == 1) | (
#          time_vector + 1e-08 > time_of_interest))))
#    }
#  
#    X <- as.matrix(X)
#    time <- as.vector(time)
#    status <- as.vector(status)
#  
#    pl <- 0
#  
#    # Calculate partial likelihood for each censored observation
#    for (i in which(status == 1)) {
#      risk_set_idx <- get_cox_risk_set_idx(
#        time_of_interest = time[i],
#        entry_vector = entry_points,
#        time_vector = time,
#        status_vector = status
#      )
#  
#      # Compute the linear predictors for the risk set
#      exp_risk_lp <- exp(X[risk_set_idx, , drop = FALSE] %*% coefs)
#  
#      # Update partial likelihood
#      pl <- pl + sum(X[i, , drop = FALSE] * coefs) - log(sum(exp_risk_lp))
#    }
#  
#    return(pl)
#  }
#  
#  # Load pbc dataset for Cox proportional hazards model
#  data("pbc")
#  
#  pbc <- pbc[stats::complete.cases(pbc), 2:20] # Remove NA values and `id` column
#  
#  pbc$trt <- ifelse(pbc$trt == 1, 1, 0) # Convert trt to binary
#  pbc$sex <- ifelse(pbc$sex == "f", 1, 0) # Convert sex to binary
#  pbc$edema2 <- ifelse(pbc$edema == 0.5, 1, 0) # Convert edema2 to binary
#  pbc$edema <- ifelse(pbc$edema == 1, 1, 0) # Convert edema to binary
#  pbc$time <- pbc$time / 365 # Convert time to years
#  pbc$status <- (pbc$status == 2) * 1 # Convert status to binary
#  
#  # Identify columns to be standardized
#  not_standarlized_index <- c(1, 2, 3, 5, 6, 7, 8, 9, 20)
#  
#  # Standardize the columns
#  pbc[, -not_standarlized_index] <- scale(pbc[, -not_standarlized_index])
#  
#  # Seperate observation data into train and test data
#  train_n <- (ncol(pbc) - 2) * 5
#  train_idx <- sample(1:nrow(pbc), train_n)
#  train_data <- pbc[train_idx, ]
#  test_data <- pbc[-train_idx, ]
#  
#  test_x <- test_data[, -which(names(test_data) %in% c("time", "status"))]
#  test_time <- test_data$time
#  test_status <- test_data$status
#  test_entry_points <- rep(0, nrow(test_x))
#  
#  dim(train_data)

## ----eval=FALSE---------------------------------------------------------------
#  # Calculate MPLE (Marginal Partial Likelihood Estimate) prediction score with 0 coefficients
#  MPLE_0 <- get_cox_partial_likelihood(
#    X = test_x,
#    time = test_time,
#    status = test_status,
#    coefs = rep(0, (ncol(test_x))),
#    entry_points = test_entry_points
#  )
#  
#  # Fit a COX regression model
#  cox_model <- survival::coxph(
#    formula = survival::Surv(time, status) ~ .,
#    data = train_data
#  )
#  
#  # Calculate MPLE (Marginal Partial Likelihood Estimate) prediction score of estimated coefficients minus 0 coefficients
#  cat(
#    "MLE COX Model - Predition Score:",
#    get_cox_partial_likelihood(
#      X = test_x,
#      time = test_time,
#      status = test_status,
#      coefs = coef(cox_model),
#      entry_points = test_entry_points
#    ) - MPLE_0
#  )

## ----eval=FALSE---------------------------------------------------------------
#  cat_init <- cat_cox_initialization(
#    formula = survival::Surv(time, status) ~ 1,
#    data = train_data,
#    syn_size = 100, # Default: 4 times the number of predictors
#    hazard_constant = 1, # Default: calculated from observed data
#    entry_points = rep(0, nrow(train_data)), # Default: Null, which is vector of zeros
#    x_degree = NULL, # Default: NULL, which means degrees of all predictors are 1
#    resample_only = FALSE, # Default: FALSE
#    na_replace = mean # Default: stats::na.omit
#  )
#  
#  cat_init

## ----eval=FALSE---------------------------------------------------------------
#  names(cat_init)

## ----eval=FALSE---------------------------------------------------------------
#  # The number of observations (rows) in the original dataset (`obs_data`)
#  cat_init$obs_size
#  
#  # The information detailing the process of resampling synthetic data
#  cat_init$syn_x_resample_inform

## ----eval=FALSE---------------------------------------------------------------
#  cat_cox_model <- cat_cox(
#    formula = survival::Surv(time, status) ~ ., # Same as `~ .`
#    cat_init = cat_init,
#    tau = 10, # Default: number of predictors
#    method = "WME", # Default: CRE
#    init_coefficients = NULL, # Default: all zeros
#    tol = 0.01, # Default: 1e-5
#    max_iter = 5 # Default: 25
#  )
#  
#  cat_cox_model

## ----eval=FALSE---------------------------------------------------------------
#  cat(
#    "Catalytic `cat_cox` - Predition Score:",
#    get_cox_partial_likelihood(
#      X = test_x,
#      time = test_time,
#      status = test_status,
#      coefs = coef(cat_cox_model),
#      entry_points = test_entry_points
#    ) - MPLE_0
#  )

## ----eval=FALSE---------------------------------------------------------------
#  names(cat_cox_model)

## ----eval=FALSE---------------------------------------------------------------
#  # The formula used for modeling
#  cat_cox_model$formula
#  
#  # The fitted model object obtained from `survival::coxph` with `tau`
#  cat_cox_model$coefficients

## ----eval=FALSE---------------------------------------------------------------
#  cat_cox_tune_model <- cat_cox_tune(
#    formula = survival::Surv(time, status) ~ ., # Same as `~ .`
#    cat_init = cat_init,
#    method = "CRE", # Default: "CRE"
#    tau_seq = seq(1, 5, 0.5), # Default is a numeric sequence around the number of predictors
#    cross_validation_fold_num = 3, # Default: 5
#    tol = 1 # Additional arguments passed to the `cat_cox` function
#  )
#  
#  cat_cox_tune_model

## ----eval=FALSE---------------------------------------------------------------
#  plot(cat_cox_tune_model, text_pos = 2, legend_pos = "bottomright")

## ----eval=FALSE---------------------------------------------------------------
#  cat(
#    "Catalytic `cat_cox_tune` - Predition Score:",
#    get_cox_partial_likelihood(
#      X = test_x,
#      time = test_time,
#      status = test_status,
#      coefs = coef(cat_cox_tune_model),
#      entry_points = test_entry_points
#    ) - MPLE_0
#  )

## ----eval=FALSE---------------------------------------------------------------
#  names(cat_cox_tune_model)

## ----eval=FALSE---------------------------------------------------------------
#  # The estimated coefficients from the fitted model
#  cat_cox_tune_model$coefficients
#  
#  # Selected optimal tau value from `tau_seq` that maximize the likelihood
#  cat_cox_tune_model$tau

## ----eval = FALSE-------------------------------------------------------------
#  cat_cox_bayes_model <- cat_cox_bayes(
#    formula = survival::Surv(time, status) ~ ., # Same as `~ .
#    cat_init = cat_init,
#    tau = 1, # Default: NULL
#    chains = 1, # Default: 4
#    iter = 100, # Default: 2000
#    warmup = 50, # Default: 1000
#    hazard_beta = 2 # Default: 2
#  )
#  
#  cat_cox_bayes_model

## ----eval=FALSE---------------------------------------------------------------
#  traceplot(cat_cox_bayes_model)

## ----eval=FALSE---------------------------------------------------------------
#  traceplot(cat_cox_bayes_model, inc_warmup = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  cat(
#    "Catalytic `cat_cox_bayes` - Predition Score:",
#    get_cox_partial_likelihood(
#      X = test_x,
#      time = test_time,
#      status = test_status,
#      coefs = coef(cat_cox_bayes_model),
#      entry_points = test_entry_points
#    ) - MPLE_0
#  )

## ----eval = FALSE-------------------------------------------------------------
#  names(cat_cox_bayes_model)

## ----eval = FALSE-------------------------------------------------------------
#  # The number of Markov chains used during MCMC sampling in `rstan`.
#  cat_cox_bayes_model$chain
#  
#  # The mean estimated coefficients from the Bayesian model
#  cat_cox_bayes_model$coefficients

## ----eval = FALSE-------------------------------------------------------------
#  cat_cox_bayes_joint_model <- cat_cox_bayes_joint(
#    formula = survival::Surv(time, status) ~ ., # Same as `~ .
#    cat_init = cat_init,
#    chains = 1, # Default: 4
#    iter = 100, # Default: 2000
#    warmup = 50, # Default: 1000
#    tau_alpha = 2, # Default: 2
#    tau_gamma = 1, # Default: 1
#    hazard_beta = 2 # Default: 2
#  )
#  
#  cat_cox_bayes_joint_model

## ----eval = FALSE-------------------------------------------------------------
#  traceplot(cat_cox_bayes_joint_model)

## ----eval = FALSE-------------------------------------------------------------
#  traceplot(cat_cox_bayes_joint_model, inc_warmup = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  cat(
#    "Catalytic `cat_cox_bayes_joint` - Predition Score:",
#    get_cox_partial_likelihood(
#      X = test_x,
#      time = test_time,
#      status = test_status,
#      coefs = coef(cat_cox_bayes_joint_model),
#      entry_point = test_entry_points
#    ) - MPLE_0
#  )

## ----eval = FALSE-------------------------------------------------------------
#  names(cat_cox_bayes_joint_model)

## ----eval = FALSE-------------------------------------------------------------
#  # The estimated tau parameter from the MCMC sampling `rstan::sampling`,
#  cat_cox_bayes_joint_model$tau
#  
#  # The mean estimated coefficients from the Bayesian model
#  cat_cox_bayes_joint_model$coefficients


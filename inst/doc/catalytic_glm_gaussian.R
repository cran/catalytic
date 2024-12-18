## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 7,
  fig.height = 6,
  cache = TRUE
)

## -----------------------------------------------------------------------------
library(catalytic)

set.seed(1)

n <- 20 # Number of observations
p <- 5 # Number of predictors

obs_x <- matrix(rnorm(n * (p - 1)), ncol = (p - 1)) # Observation covariates
true_coefs <- rnorm(p) # True coefficient
noise <- rnorm(n) # Noise for more response variability
obs_y <- true_coefs[1] + obs_x %*% true_coefs[-1] + noise # Observation response

obs_data <- as.data.frame(cbind(obs_x, obs_y))
names(obs_data) <- c(paste0("X", 1:(p - 1)), "Y")

# Seperate observation data into train and test data
train_idx <- sample(n, 10)
train_data <- obs_data[train_idx, ]
test_data <- obs_data[-train_idx, ]

print(dim(train_data))

## -----------------------------------------------------------------------------
# Fit a Linear regression model (GLM)
glm_model <- stats::glm(
  formula = Y ~ .,
  family = gaussian,
  data = train_data
)

predicted_y <- predict(
  glm_model,
  newdata = test_data
)

cat(
  "MLE GLM gaussian Model - Mean Square Error (Data):",
  mean((predicted_y - test_data$Y)^2)
)

cat(
  "\nMLE GLM gaussian Model - Sum Square Error (Coefficients):",
  sum((coef(glm_model) - true_coefs)^2)
)

## -----------------------------------------------------------------------------
plot(test_data$Y,
  predicted_y,
  main = "Scatter Plot of true Y vs Predicted Y",
  xlab = "true Y",
  ylab = "Predicted Y",
  pch = 19,
  col = "blue"
)

# Add a 45-degree line for reference
abline(a = 0, b = 1, col = "red", lwd = 2)

## -----------------------------------------------------------------------------
cat_init <- cat_glm_initialization(
  formula = Y ~ 1,
  family = gaussian,
  data = train_data,
  syn_size = 50,
  custom_variance = NULL,
  gaussian_known_variance = FALSE,
  x_degree = NULL,
  resample_only = FALSE,
  na_replace = stats::na.omit
)

cat_init

## -----------------------------------------------------------------------------
names(cat_init)

## -----------------------------------------------------------------------------
# The number of observations (rows) in the original dataset (`obs_data`)
cat_init$obs_size

# The information detailing the process of resampling synthetic data
cat_init$syn_x_resample_inform

## -----------------------------------------------------------------------------
cat_glm_model <- cat_glm(
  formula = Y ~ ., # Same as `~ .`
  cat_init = cat_init, # Output object from `cat_glm_initialization`
  tau = 10 # Defaults to the number of predictors / 4
)

cat_glm_model

## -----------------------------------------------------------------------------
cat_glm_predicted_y <- predict(
  cat_glm_model,
  newdata = test_data
)

cat(
  "Catalytic `cat_glm` - Mean Square Error (Data):",
  mean((cat_glm_predicted_y - test_data$Y)^2)
)

cat(
  "\nCatalytic `cat_glm` - Sum Square Error (Coefficients):",
  sum((coef(cat_glm_model) - true_coefs)^2)
)

## -----------------------------------------------------------------------------
plot(test_data$Y,
  cat_glm_predicted_y,
  main = "Scatter Plot of true Y vs Predicted Y (cat_glm)",
  xlab = "true Y",
  ylab = "Predicted Y (cat_glm)",
  pch = 19,
  col = "blue"
)

# Add a 45-degree line for reference
abline(a = 0, b = 1, col = "red", lwd = 2)

## -----------------------------------------------------------------------------
names(cat_glm_model)

## -----------------------------------------------------------------------------
# The formula used for modeling
cat_glm_model$formula

# The fitted GLM model object obtained from `stats::glm` with `tau`
cat_glm_model$coefficients

## -----------------------------------------------------------------------------
cat_glm_tune_cv <- cat_glm_tune(
  formula = Y ~ ., # Same as `~ .`
  cat_init = cat_init,
  risk_estimate_method = "cross_validation", # Default auto-select based on the data size
  discrepancy_method = "mean_square_error", # Default auto-select based on family
  tau_seq = seq(0, 5, 0.5), # Default is a numeric sequence around the number of predictors / 4
  cross_validation_fold_num = 2 # Default: 5
)

cat_glm_tune_cv

## -----------------------------------------------------------------------------
plot(cat_glm_tune_cv, legend_pos = "topright", text_pos = 3)

## -----------------------------------------------------------------------------
cat_glm_tune_boots <- cat_glm_tune(
  formula = ~., # Same as `Y ~ .`
  cat_init = cat_init,
  risk_estimate_method = "parametric_bootstrap", # Default auto-select based on the data size
  discrepancy_method = "mean_square_error", # Default auto-select based on family
  tau_0 = 2, # Default: 1
  parametric_bootstrap_iteration_times = 5, # Default: 100
)

cat_glm_tune_boots

## -----------------------------------------------------------------------------
cat_glm_tune_mallowian <- cat_glm_tune(
  formula = ~., # Same as `Y ~ .`
  cat_init = cat_init,
  risk_estimate_method = "mallowian_estimate", # Default auto-select based on the data size
  discrepancy_method = "mean_square_error", # Default auto-select based on family
)

cat_glm_tune_mallowian

## -----------------------------------------------------------------------------
cat_glm_tune_auto <- cat_glm_tune(
  formula = ~., # Same as `Y ~ .`
  cat_init = cat_init
)

cat_glm_tune_auto

## ----eval=FALSE---------------------------------------------------------------
#  cat_glm_tune_predicted_y <- predict(
#    cat_glm_tune_auto,
#    newdata = test_data
#  )
#  
#  cat(
#    "Catalytic `cat_glm_tune` - Mean Square Error (Data):",
#    mean((cat_glm_tune_predicted_y - test_data$Y)^2)
#  )
#  
#  cat(
#    "\nCatalytic `cat_glm_tune` - Sum Square Error (Coefficients):",
#    sum((coef(cat_glm_tune_auto) - true_coefs)^2)
#  )

## ----eval=FALSE---------------------------------------------------------------
#  plot(test_data$Y, cat_glm_tune_predicted_y,
#    main = "Scatter Plot of true Y vs Predicted Y (cat_glm_tune)",
#    xlab = "true Y",
#    ylab = "Predicted Y (cat_glm_tune)",
#    pch = 19,
#    col = "blue"
#  )
#  
#  # Add a 45-degree line for reference
#  abline(a = 0, b = 1, col = "red", lwd = 2)

## ----eval=FALSE---------------------------------------------------------------
#  names(cat_glm_tune_auto)

## ----eval=FALSE---------------------------------------------------------------
#  # The method used for risk estimation
#  cat_glm_tune_auto$risk_estimate_method
#  
#  # Selected optimal tau value from `tau_seq` that minimizes discrepancy error
#  cat_glm_tune_auto$tau

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_model <- cat_glm_bayes(
#    formula = Y ~ ., # Same as `~ .`
#    cat_init = cat_init,
#    tau = 50, # Default: number of predictors / 4
#    chains = 1, # Default: 4
#    iter = 100, # Default: 2000
#    warmup = 50, # Default: 1000
#    algorithm = "NUTS", # Default: NUTS
#    gaussian_variance_alpha = 1, # Default: number of predictors
#    gaussian_variance_beta = 1 # Default: number of predictors times variance of observation response
#  )
#  
#  cat_glm_bayes_model

## ----eval = FALSE-------------------------------------------------------------
#  traceplot(cat_glm_bayes_model)

## ----eval = FALSE-------------------------------------------------------------
#  traceplot(cat_glm_bayes_model, inc_warmup = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_predicted_y <- predict(
#    cat_glm_bayes_model,
#    newdata = test_data
#  )
#  
#  cat(
#    "Catalytic cat_glm_bayes - Mean Square Error (Data):",
#    mean((cat_glm_bayes_predicted_y - test_data$Y)^2)
#  )
#  
#  cat(
#    "\nCatalytic cat_glm_bayes - Sum Square Error (Coefficients):",
#    sum((coef(cat_glm_bayes_model) - true_coefs)^2)
#  )

## ----eval = FALSE-------------------------------------------------------------
#  plot(test_data$Y, cat_glm_bayes_predicted_y,
#    main = "Scatter Plot of true Y vs Predicted Y (cat_glm_bayes)",
#    xlab = "true Y",
#    ylab = "Predicted Y (cat_glm_bayes)",
#    pch = 19,
#    col = "blue"
#  )
#  
#  # Add a 45-degree line for reference
#  abline(a = 0, b = 1, col = "red", lwd = 2)

## ----eval = FALSE-------------------------------------------------------------
#  names(cat_glm_bayes_model)

## ----eval = FALSE-------------------------------------------------------------
#  # The number of Markov chains used during MCMC sampling in `rstan`.
#  cat_glm_bayes_model$chain
#  
#  # The mean estimated coefficients from the Bayesian GLM model
#  cat_glm_bayes_model$coefficients

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_model <- cat_glm_bayes_joint(
#    formula = Y ~ ., # Same as `~ .`
#    cat_init = cat_init,
#    chains = 1, # Default: 4
#    iter = 100, # Default: 2000
#    warmup = 50, # Default: 1000
#    algorithm = "NUTS", # Default: NUTS
#    tau_alpha = 2, # Default: 2
#    tau_gamma = 1 # Default: 1
#  )
#  
#  cat_glm_bayes_joint_model

## ----eval = FALSE-------------------------------------------------------------
#  traceplot(cat_glm_bayes_joint_model)

## ----eval = FALSE-------------------------------------------------------------
#  traceplot(cat_glm_bayes_joint_model, inc_warmup = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_predicted_y <- predict(
#    cat_glm_bayes_joint_model,
#    newdata = test_data
#  )
#  
#  cat(
#    "Catalytic `cat_glm_bayes_joint` - Mean Square Error (Data):",
#    mean((cat_glm_bayes_joint_predicted_y - test_data$Y)^2)
#  )
#  
#  cat(
#    "\nCatalytic `cat_glm_bayes_joint` - Sum Square Error (Coefficients):",
#    sum((coef(cat_glm_bayes_joint_model) - true_coefs)^2)
#  )

## ----eval = FALSE-------------------------------------------------------------
#  plot(test_data$Y,
#    cat_glm_bayes_joint_predicted_y,
#    main = "Scatter Plot of true Y vs Predicted Y (cat_glm_bayes_joint)",
#    xlab = "true Y",
#    ylab = "Predicted Y (cat_glm_bayes_joint)",
#    pch = 19,
#    col = "blue"
#  )
#  
#  # Add a 45-degree line for reference
#  abline(a = 0, b = 1, col = "red", lwd = 2)

## ----eval = FALSE-------------------------------------------------------------
#  names(cat_glm_bayes_joint_model)

## ----eval = FALSE-------------------------------------------------------------
#  # The estimated tau parameter from the MCMC sampling `rstan::sampling`,
#  cat_glm_bayes_joint_model$tau
#  
#  # The mean estimated coefficients from the Bayesian GLM model
#  cat_glm_bayes_joint_model$coefficients


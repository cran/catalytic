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
#  library(catalytic)
#  
#  set.seed(1)
#  
#  # Function for calculating the mean of logarithmic error between true response and estimated response
#  get_mean_logarithmic_error <- function(Y, est_Y) {
#    Y <- pmax(0.0001, pmin(0.9999, Y))
#    est_Y <- pmax(0.0001, pmin(0.9999, est_Y))
#    return(mean(Y * log(Y / est_Y) + (1 - Y) * log((1 - Y) / (1 - est_Y))))
#  }
#  
#  # Load swim dataset for binomial analysis
#  data("swim")
#  swim_data <- cbind(swim$x, swim$y)
#  
#  # Seperate observation data into train and test data
#  n <- 5 * ncol(swim$x + 1) # Size for training data
#  train_idx <- sample(1:nrow(swim_data), n)
#  train_data <- swim_data[train_idx, ]
#  test_data <- swim_data[-train_idx, ]
#  
#  dim(train_data)

## ----eval=FALSE---------------------------------------------------------------
#  # Fit a logistic regression model (GLM)
#  glm_model <- stats::glm(
#    formula = empyr1 ~ .,
#    family = binomial,
#    data = train_data
#  )
#  
#  predicted_y <- predict(
#    glm_model,
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "MLE GLM Binomial Model - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = predicted_y
#    )
#  )

## ----eval=FALSE---------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, predicted_y)
#  plot(roc_curve, main = "ROC Curve (MLE)", col = "blue", lwd = 2)

## ----eval=FALSE---------------------------------------------------------------
#  cat_init <- cat_glm_initialization(
#    formula = empyr1 ~ 1,
#    family = binomial, # Default: gaussian
#    data = train_data,
#    syn_size = 50, # Default: 4 times the number of predictors
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
#  cat_glm_model <- cat_glm(
#    formula = empyr1 ~ female + agege35, # Same as `~ female + agege35`
#    cat_init = cat_init,
#    tau = 10 # Default: number of predictors / 4
#  )
#  
#  cat_glm_model

## ----eval=FALSE---------------------------------------------------------------
#  cat_glm_predicted_y <- predict(
#    cat_glm_model,
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "Catalytic `cat_glm` - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = cat_glm_predicted_y
#    )
#  )

## ----eval=FALSE---------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, cat_glm_predicted_y)
#  plot(roc_curve, main = "ROC Curve (cat_glm)", col = "blue", lwd = 2)

## ----eval=FALSE---------------------------------------------------------------
#  names(cat_glm_model)

## ----eval=FALSE---------------------------------------------------------------
#  # The formula used for modeling
#  cat_glm_model$formula
#  
#  # The fitted GLM model object obtained from `stats::glm` with `tau`
#  cat_glm_model$coefficients

## ----eval=FALSE---------------------------------------------------------------
#  cat_glm_tune_cv <- cat_glm_tune(
#    formula = empyr1 ~ ., # Same as `~ .`
#    cat_init = cat_init,
#    risk_estimate_method = "cross_validation", # Default auto-select based on the data size
#    discrepancy_method = "logistic_deviance", # Default auto-select based on family
#    tau_seq = seq(0.1, 5.1, 0.5), # Default is a numeric sequence around the number of predictors / 4. Do not recommand to including 0 for cross validation.
#    cross_validation_fold_num = 3 # Default: 5
#  )
#  
#  cat_glm_tune_cv

## ----eval=FALSE---------------------------------------------------------------
#  plot(cat_glm_tune_cv)

## ----eval=FALSE---------------------------------------------------------------
#  cat_glm_tune_boots <- cat_glm_tune(
#    formula = ~., # Same as `empyr1 ~ .`
#    cat_init = cat_init,
#    risk_estimate_method = "parametric_bootstrap", # Default auto-select based on the data size
#    discrepancy_method = "logistic_deviance", # Default auto-select based on family
#    tau_0 = 2, # Default: number of predictors / 4
#    parametric_bootstrap_iteration_times = 10, # Default: 100
#  )
#  
#  cat_glm_tune_boots

## ----eval=FALSE---------------------------------------------------------------
#  cat_glm_tune_stein <- cat_glm_tune(
#    formula = ~., # Same as `empyr1 ~ .`
#    cat_init = cat_init,
#    risk_estimate_method = "steinian_estimate", # Default auto-select based on the data size
#    discrepancy_method = "logistic_deviance", # default  auto select based on family
#  )
#  
#  cat_glm_tune_stein

## ----eval=FALSE---------------------------------------------------------------
#  cat_glm_tune_auto <- cat_glm_tune(
#    formula = ~.,
#    cat_init = cat_init
#  )
#  
#  cat_glm_tune_auto

## ----eval=FALSE---------------------------------------------------------------
#  cat_glm_tune_predicted_y <- predict(
#    cat_glm_tune_auto,
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "Catalytic `cat_glm_tune` - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = cat_glm_tune_predicted_y
#    )
#  )

## ----eval=FALSE---------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, cat_glm_tune_predicted_y)
#  plot(roc_curve, main = "ROC Curve (cat_glm_tune)", col = "blue", lwd = 2)

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
#    formula = empyr1 ~ ., # Same as `~ .`
#    cat_init = cat_init,
#    tau = 50, # Default: number of predictors / 4
#    chains = 1, # Default: 4
#    iter = 100, # Default: 2000
#    warmup = 50, # Default: 1000
#    algorithm = "NUTS" # Default: NUTS
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
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "MLE GLM Binomial Model - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = cat_glm_bayes_predicted_y
#    )
#  )

## ----eval = FALSE-------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, cat_glm_bayes_predicted_y)
#  plot(roc_curve, main = "ROC Curve (cat_glm_bayes)", col = "blue", lwd = 2)

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
#    formula = empyr1 ~ ., # Same as `~ .`
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
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "Catlytic `cat_glm_bayes_joint` - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = cat_glm_bayes_joint_predicted_y
#    )
#  )

## ----eval = FALSE-------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, cat_glm_bayes_joint_predicted_y)
#  plot(roc_curve, main = "ROC Curve (cat_glm_bayes_joint)", col = "blue", lwd = 2)

## ----eval = FALSE-------------------------------------------------------------
#  names(cat_glm_bayes_joint_model)

## ----eval = FALSE-------------------------------------------------------------
#  # The estimated tau parameter from the MCMC sampling `rstan::sampling`,
#  cat_glm_bayes_joint_model$tau
#  
#  # The mean estimated coefficients from the Bayesian GLM model
#  cat_glm_bayes_joint_model$coefficients

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_gibbs_model <- cat_glm_bayes_joint_gibbs(
#    formula = empyr1 ~ ., # Same as `~ .`
#    cat_init = cat_init,
#    iter = 100, # Default: 1000
#    warmup = 50, # Default: 500
#    coefs_iter = 2, # Default: 5
#    tau_0 = 1, # Default: number of predictors / 4
#    tau_alpha = 2, # Default: 2
#    tau_gamma = 1, # Default: 1
#    refresh = TRUE # Default: TRUE
#  )
#  
#  cat_glm_bayes_joint_gibbs_model

## ----eval = FALSE-------------------------------------------------------------
#  traceplot(cat_glm_bayes_joint_gibbs_model)

## ----eval = FALSE-------------------------------------------------------------
#  traceplot(cat_glm_bayes_joint_gibbs_model, pars = c("female", "hsdip", "numchild"))

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_gibbs_predicted_y <- predict(
#    cat_glm_bayes_joint_gibbs_model,
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "Catalytic `cat_glm_bayes_joint_gibbs` - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = cat_glm_bayes_joint_gibbs_predicted_y
#    )
#  )

## ----eval = FALSE-------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, cat_glm_bayes_joint_gibbs_predicted_y)
#  plot(roc_curve, main = "ROC Curve (cat_glm_bayes_joint_gibbs)", col = "blue", lwd = 2)

## ----eval = FALSE-------------------------------------------------------------
#  names(cat_glm_bayes_joint_gibbs_model)

## ----eval = FALSE-------------------------------------------------------------
#  # The Date and time when sampling ended.
#  cat_glm_bayes_joint_gibbs_model$sys_time
#  
#  # The summary statistics
#  cat_glm_bayes_joint_gibbs_model$inform_df

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_tau_lower <- cat_glm_bayes_joint(
#    formula = ~., # Same as `empyr1 ~ .`
#    cat_init = cat_init,
#    binomial_tau_lower = 0.5, # Default: 0.05
#    chains = 1, # Default: 4
#    iter = 100, # Default: 2000
#    warmup = 50 # Default: 1000
#  )
#  
#  cat_glm_bayes_joint_tau_lower

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_tau_lower_predicted_y <- predict(
#    cat_glm_bayes_joint_tau_lower,
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "Catalytic `cat_glm_bayes_joint_tau_lower` - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = cat_glm_bayes_joint_tau_lower_predicted_y
#    )
#  )

## ----eval = FALSE-------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, cat_glm_bayes_joint_tau_lower_predicted_y)
#  plot(roc_curve, main = "ROC Curve (MLE)", col = "blue", lwd = 2)

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_theta <- cat_glm_bayes_joint(
#    formula = ~., # Same as `empyr1 ~ .`
#    cat_init = cat_init,
#    binomial_joint_theta = TRUE, # Default: FALSE
#    chains = 1, # Default: 4
#    iter = 100, # Default: 2000
#    warmup = 50 # Default: 1000
#  )
#  
#  cat_glm_bayes_joint_theta

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_theta_predicted_y <- predict(
#    cat_glm_bayes_joint_theta,
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "Catalytic `cat_glm_bayes_joint_theta` - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = cat_glm_bayes_joint_theta_predicted_y
#    )
#  )

## ----eval = FALSE-------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, cat_glm_bayes_joint_theta_predicted_y)
#  plot(roc_curve, main = "ROC Curve (cat_glm_bayes_joint_theta)", col = "blue", lwd = 2)

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_alpha <- cat_glm_bayes_joint(
#    formula = ~., # Same as `empyr1 ~ .`
#    cat_init = cat_init,
#    binomial_joint_theta = TRUE, # Default: FALSE
#    binomial_joint_alpha = TRUE, # Default: FALSE
#    chains = 1, # Default: 4
#    iter = 100, # Default: 2000
#    warmup = 50 # Default: 1000
#  )
#  
#  cat_glm_bayes_joint_alpha

## ----eval = FALSE-------------------------------------------------------------
#  cat_glm_bayes_joint_alpha_predicted_y <- predict(
#    cat_glm_bayes_joint_alpha,
#    newdata = test_data,
#    type = "response"
#  )
#  
#  cat(
#    "Catalytic `cat_glm_bayes_joint_alpha` - Logarithmic Error:",
#    get_mean_logarithmic_error(
#      Y = test_data$empyr1,
#      est_Y = cat_glm_bayes_joint_alpha_predicted_y
#    )
#  )

## ----eval = FALSE-------------------------------------------------------------
#  roc_curve <- pROC::roc(test_data$empyr1, cat_glm_bayes_joint_alpha_predicted_y)
#  plot(roc_curve, main = "ROC Curve (cat_glm_bayes_joint_alpha)", col = "blue", lwd = 2)


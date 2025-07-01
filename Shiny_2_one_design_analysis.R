## simulate data for a certain design.

generateGrowthModel <- function(
    design = c(1, 2, 5),
    int_mean = 0,
    slope_mean = 0,
    int_var = 1,
    slope_var = 1,
    int_slope_cov = 0,
    res_var = 1
) {
  y_names <- paste0("y", design)
  time_scores <- design - min(design)  # Align to start at 0 for slope
  
  # Factor loadings
  intercept_line <- paste0("i =~ ", paste0("1*", y_names, collapse = " + "))
  slope_line <- paste0("s =~ ", paste0(time_scores, "*", y_names, collapse = " + "))
  
  model <- c(
    intercept_line,
    slope_line,
    paste0("i ~~ ", int_var, "*i"),
    paste0("s ~~ ", slope_var, "*s"),
    paste0("i ~~ ", int_slope_cov, "*s"),
    paste0("i ~ ", int_mean, "*1"),
    paste0("s ~ ", slope_mean, "*1"),
    paste0(y_names, " ~~ ", res_var, "*", y_names),
    paste0(y_names, " ~ 0*1")
  )
  
  return(paste(model, collapse = "\n"))
}



###########
model.ana <- c(
  intercept_line,
  slope_line,
  paste0("i ~~ i"),
  paste0("s ~~ s"),
  paste0("i ~~ s"),
  paste0("i ~ 1"),
  paste0("s ~ 1")
)
############

library(lavaan)

simulateAndEstimateSlopeSE <- function(
    design = c(1, 2, 5),
    sample_size = 200,
    int_mean = 0,
    slope_mean = 0,
    int_var = 1,
    slope_var = 1,
    int_slope_cov = 0,
    res_var = 1,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  
  

  # Generate model only for observed timepoints
  model <- generateGrowthModel(
    design = design,
    int_mean = int_mean,
    slope_mean = slope_mean,
    int_var = int_var,
    slope_var = slope_var,
    int_slope_cov = int_slope_cov,
    res_var = res_var
  )
  
  # Simulate data
  data <- simulateData(model, sample.nobs = sample_size)
  
  # Fit model
  fit <- growth(model = model.ana, data = data, missing = "fiml")
  
  # Extract SE of slope mean
  pt <- parameterEstimates(fit)
  slope_se <- pt$se[pt$lhs == "s" & pt$op == "~1"]
  slope_est <- pt$est[pt$lhs == "s" & pt$op == "~1"]
  
  return(list(
    slope_est = slope_est,
    slope_se = slope_se,
    model = model,
    fit = fit
  ))
}


result <- simulateAndEstimateSlopeSE(
  design = c(1, 2 , 5),
  sample_size = 1666,
  int_mean = 39.457,
  slope_mean = 8.063,
  int_var = 28.776,
  slope_var = 8.201,
  int_slope_cov = 1.56,
  res_var = 30,
  seed = 123
)

cat("Slope estimate:", round(result$slope_est, 3), "\n")
cat("Slope SE:", round(result$slope_se, 3), "\n")

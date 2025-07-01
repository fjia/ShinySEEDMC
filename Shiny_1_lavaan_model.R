#### generate lavaan model using given parameter values

library(lavaan)

generateGrowthModel <- function(
    timepoints = 5,
    int_mean = 0,     # Intercept mean
    slope_mean = 0,   # Slope mean
    int_var = 1,      # Intercept variance
    slope_var = 1,    # Slope variance
    int_slope_cov = 0,# Covariance between intercept and slope
    res_var = 1       # Residual variance (equal across time)
) {
  # Variable names: y1, y2, ..., yT
  y_names <- paste0("y", 1:timepoints)
  
  # Loadings
  intercept_loadings <- paste0("1*", y_names, collapse = " + ")
  slope_loadings <- paste0((0:(timepoints - 1)), "*", y_names, collapse = " + ")
  
  # Latent variable definitions
  model <- c(
    paste0("i =~ ", intercept_loadings),
    paste0("s =~ ", slope_loadings),
    
    # Latent variances and covariances
    paste0("i ~~ ", int_var, "*i"),
    paste0("s ~~ ", slope_var, "*s"),
    paste0("i ~~ ", int_slope_cov, "*s"),
    
    # Latent means
    paste0("i ~ ", int_mean, "*1"),
    paste0("s ~ ", slope_mean, "*1"),
    
    # Residual variances
    paste0(y_names, " ~~ ", res_var, "*", y_names),
    
    # Indicator means fixed to 0
    paste0(y_names, " ~ 0*1")
  )
  
  lavaan_model <- paste(model, collapse = "\n")
  return(lavaan_model)
}


# Example Usage for Simulating Data

# Generate model for 6 timepoints with custom parameters
mod <- generateGrowthModel(
  timepoints = 5,
  int_mean = 39.457,
  slope_mean = 8.063,
  int_var = 28.776,
  slope_var = 8.201,
  int_slope_cov = 1.56,
  res_var = 30
)

cat(mod)  # View the lavaan model

# Simulate data
set.seed(123)
sim_data <- simulateData(mod, sample.nobs = 200)

head(sim_data)
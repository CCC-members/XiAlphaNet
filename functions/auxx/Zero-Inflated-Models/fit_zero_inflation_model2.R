# Load necessary library
library(glmmTMB)

# Read in the data (assume a CSV named "temp_data.csv")
data <- read.csv("temp_data.csv")

# Check if 'group' column is provided; if so, fit a model with random effects
if("group" %in% colnames(data)) {
  # Fit the zero-inflated Poisson model with random effects
  model <- glmmTMB(y ~ x + I(x^2)+(1|group),
                   data = data,
                   ziformula = ~ x +(1|group),
                   family = gaussian(link = "identity"))  # Use Poisson family with log link
} else {
  # Fit the zero-inflated Poisson model without random effects
  model <- glmmTMB(y ~ x + I(x^2),
                   data = data,
                   ziformula = ~ x+ I(x^2),
                   family = gaussian(link = "log"))  # Use Poisson family with log link
}

# Check if the model has converged and check summary
summary(model)

# Extract coefficients from both the conditional and zero-inflation models
cond_coefs <- summary(model)$coefficients$cond
zi_coefs <- summary(model)$coefficients$zi

# Save the coefficients to a CSV file for MATLAB to read
coef_data <- data.frame(
  param = c(rownames(cond_coefs), rownames(zi_coefs)),
  estimate = c(cond_coefs[, 1], zi_coefs[, 1]),
  std_error = c(cond_coefs[, 2], zi_coefs[, 2]),
  z_value = c(cond_coefs[, 3], zi_coefs[, 3]),
  p_value = c(cond_coefs[, 4], zi_coefs[, 4])
)

write.csv(coef_data, "model_coefficients.csv", row.names = FALSE)

# Read the x_eval values from the x_eval.csv file
x_eval_data <- read.csv("x_eval.csv")
x_eval <- x_eval_data$x  # Extract x_eval from the CSV

# If x_eval is provided, make predictions and calculate standard errors
if (!is.null(x_eval)) {
  # Create a data frame for the new evaluation points (including group as a placeholder)
  eval_data <- data.frame(
    x = x_eval,        # New x values
    x2 = x_eval^2,     # New x^2 values
    group = factor(rep(1, length(x_eval)))  # Add 'group' as a placeholder
  )
  
  # Conditional model predictions (log link)
  cond_pred <- predict(model, newdata = eval_data, type = "link", se.fit = TRUE)

  # Zero-inflation predictions (probability of zero-inflation)
  zi_pred <- predict(model, newdata = eval_data, type = "zprob", se.fit = TRUE)

  # Combine predictions and standard errors into a result data frame
  result <- data.frame(
    x_eval = x_eval,
    cond_pred = (cond_pred$fit),  # Exponentiate to get the mean rate
    cond_se = cond_pred$se.fit,      # Standard error for conditional model
    zi_pred = zi_pred$fit,           # Zero-inflation probability
    zi_se = zi_pred$se.fit          # Standard error for zero-inflation model
  )

  # Save predictions and standard errors to a CSV file for MATLAB to read
  write.csv(result, "predictions_with_se.csv", row.names = FALSE)
} else {
  cat("No x_eval data provided for predictions.\n")
}

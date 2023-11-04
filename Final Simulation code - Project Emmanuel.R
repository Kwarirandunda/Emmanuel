# Load required libraries
library(MASS)
library(mixtools)
library(ggplot2)

# Set random seed for reproducibility
set.seed(100)

# Define simulation parameters
samples <- 500  # Number of data points
common_mean <- 0  # Common mean for both components
sigma <- c(1, 2)  # Standard deviations for the two components
lambda <- c(0.8, 0.2)  # Mixing proportions

# Simulate a two-component mixture with equal means
simulated_data <- rnormmix(samples, lambda, c(common_mean, common_mean), sigma)


# Fit a single Gaussian distribution to the data
single_gaussian_fit <- fitdistr(simulated_data, "normal")
single_gaussian_fit

# Get the estimated parameters of the single Gaussian distribution
estimated_mean <- single_gaussian_fit$estimate[1]
estimated_mean

estimated_sd <- single_gaussian_fit$estimate[2]
estimated_sd

# Compute the likelihood function for the data
likelihood_function <- dnorm(simulated_data, mean = estimated_mean, sd = estimated_sd)

# Calculate the BIC and AIC
n <- length(simulated_data)
log_likelihood <- sum(log(likelihood_function))
k2 <- 2  # Number of parameters (mean and standard deviation)
bic <- -2 * log_likelihood + k2 * log(n)
aic <- -2 * log_likelihood + 2 * k2

# Print the estimated parameters, likelihood, and BIC
cat("Estimated Mean: ", estimated_mean, "\n")
cat("Estimated Standard Deviation: ", estimated_sd, "\n")
cat("Log Likelihood: ", log_likelihood, "\n")
cat("BIC: ", bic, "\n")
cat("AIC: ", aic, "\n")

# Create a histogram of the simulated data and overlay the fitted Gaussian
hist(simulated_data, breaks = 10, probability = TRUE, 
     main = "Simulated Data with fitted Gaussian")
curve(dnorm(x, mean = estimated_mean, sd = estimated_sd), 
      col = "blue", lwd = 2, add = TRUE)
################################################################################
# EM algorithm
k <- 2  # Number of components
max_iterations <- 100
tolerance <- 1e-5

# Initialize parameters
init_mean <- c(mean(simulated_data), mean(simulated_data))
init_variances <- c(var(simulated_data), var(simulated_data) + 5)
init_weights <- c(0.8, 0.2)

n <- length(simulated_data)

# Initialize log likelihood with a large negative value
log_likelihood_old <- -Inf

for (iteration in 1:max_iterations) {
  # E-step
  pdf_values <- matrix(0, n, k)
  for (i in 1:k) {
    pdf_values[, i] <- dnorm(simulated_data, mean = init_mean[i], 
                             sd = sqrt(init_variances[i]))
  }
  responsibility_test_1 <- cbind(pdf_values[, 1] * init_weights[1], 
                                 pdf_values[, 2] * init_weights[2])
  responsibility_2 <- responsibility_test_1 / rowSums(responsibility_test_1)
  
  # M-step
  Nk <- colSums(responsibility_2)
  init_weights <- Nk / n
  init_mean <- colSums(simulated_data * responsibility_2) / Nk
  init_variances <- colSums((simulated_data^2) * responsibility_2) / Nk
  
  # Calculate the log likelihood
  log_likelihood_new <- sum(log(rowSums(pdf_values * init_weights)))
  
  # Check for convergence
  if (abs(log_likelihood_new - log_likelihood_old) < tolerance) {break}
  
  log_likelihood_old <- log_likelihood_new
}
# Calculate the BIC
n <- length(simulated_data)

#log_likelihood <- sum(log(likelihood_function))
k1 <- 4  # Number of parameters (mean and standard deviation)
bic_EM <- -2 * log_likelihood_new  + k1 * log(n)

# Calculate the AIC
aic_EM <- -2 * log_likelihood_new + 2 * k1


# Estimated parameters
cat("Estimated variances: ", init_variances, "\n")
cat("Estimated weights: ", init_weights, "\n")
cat("Estimated mean: ", init_mean, "\n")
cat("Log-likelihood: ", log_likelihood_new, "\n")
cat("BIC: ", bic_EM, "\n")
cat("AIC: ", aic_EM, "\n")

###############################################################################

x <- seq(min(simulated_data), max(simulated_data), length = 1000)

# Density functions for the components
den1 <- dnorm(x, mean = init_mean[1], sd = sqrt(init_variances[1]))
den2 <- dnorm(x, mean = init_mean[2], sd = sqrt(init_variances[2]))

hist(simulated_data, prob = TRUE, breaks = 10,
     main = "Histogram and EM Fit of Two-Component Mixture", xlab = "Simulated Data")

# Plot component 1
lines(x, init_weights[1] * den1, col = "red", lwd = 2)

# Plot component 2
lines(x, init_weights[2] * den2, col = "blue", lwd = 2)

# Plot the total mixture
Mixture_density <- init_weights[1] * den1 + init_weights[2] * den2
lines(x, Mixture_density, col = "black", lwd = 2)

# Add legend
legend("topright", legend = c("Mixture", "Component 1", "Component 2"), 
       col = c("black", "red", "blue"), lty = c(1, 1, 1), lwd = c(2, 2, 2))

################################################################################
                                 #Assigning
##############################################################################
# Assign c to the greater weight between nit_weights[1] and nit_weights[2]
if (init_weights[1] > init_weights[2]) {
  c <- init_weights[1]
  a <- init_mean[1]
  b <- init_variances[1]
} else {
  c <- init_weights[2]
  a <- init_mean[1]
  b <- init_variances[2]
}

# Print the values of c, a, and b
cat("c =", c, "\n")
cat("a =", a, "\n")
cat("b =", b, "\n")


# Calculate eta based on estimated variances
estimated_eta <- init_variances[1] / init_variances[2]

#contamination_pdf_eta <- function(x, a, b, c, estimated_η)
contamination_pdf_eta <- function(x1, a1, b1, c1, estimated_η1) {
first_term <- c1 * dnorm(x1, mean = a1, sd = sqrt(b1))
second_term <- (1 - c1) * dnorm(x1, mean = a1, sd = sqrt(estimated_η1 * b1))
return(first_term + second_term)
  }

#contamination_pdf_eta(x, a, b, c, estimated_eta)

# Define the "goodness" probability function with eta
goodness_probability_eta <- function(x1, a1, b1, c1, estimated_η1) {
  num <- c1 * dnorm(x1, mean = a1, sd = sqrt(b1))
  denom <- contamination_pdf_eta(x1, a1, b1, c1, estimated_η1)
  return(num / denom)
}

#goodness_probability_eta(x, a, b, c, estimated_eta)

# Classification of points
outliers_eta <- sapply(simulated_data, function(x) {
  if (goodness_probability_eta(x, a, b, c, estimated_eta) > 0.5) {
    return("Not outlier")
  } else {
    return("Outlier")
  }
})

 #outliers_eta
 
######################################################################
 # Define custom colors for the pie chart
 custom_colors <- c("lightblue", "lightcoral")
 
 # Count the number of "Not outlier" and "Outlier" classifications
 classification_counts <- table(outliers_eta)
 
 # Create a pie chart
 pie(classification_counts, 
     labels = paste(names(classification_counts), ": ", classification_counts),
     main = "Classification of data points ",
     col = custom_colors) 
 
 # Add a legend
 legend("topright", legend = names(classification_counts), fill = custom_colors)
##################################################################################

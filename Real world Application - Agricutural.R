library(MASS)
library(mixtools)
library(ggplot2)
library(bp)
library(caret)
library(pgmm)
###############################################################################
# Load the 'bp_ghana' dataset
#data(bp_ghana)
#View(bp_ghana)
# Extract the 'SBP' variable
#bp_data <- bp_ghana$SBP


###############################################################################
#load coffee bean dataset
data("coffee")
View(coffee)
# Extract the 'Bean weight' variable
data_coffee <- coffee$`Bean Weight`

################################################################################


# Identify missing and infinite values
missing_count <- sum(is.na(data_coffee))
infinite_count <- sum(!is.finite(data_coffee))

# Print the counts
cat("Number of Missing Values: ", missing_count, "\n")
cat("Number of Infinite Values: ", infinite_count, "\n")


# Check the summary of the cleaned data
summary(data_coffee)
var(data_coffee)


# Fit a normal distribution to the cleaned data
single_gaussian_fit <- fitdistr(data_coffee, "normal")

# Calculate the BIC and AIC for the single Gaussian fit
n <- length(data_coffee)
log_likelihood <- sum(dnorm(data_coffee, mean = single_gaussian_fit$estimate[1],
                            sd = single_gaussian_fit$estimate[2], log = TRUE))
k1 <- 2  # Number of parameters (mean and standard deviation)
bic_single <- -2 * log_likelihood + k1 * log(n)
aic_single <- -2 * log_likelihood + 2 * k1

# Print the results for the single Gaussian fit
cat("Single Gaussian Fit:\n")
cat("Estimated Mean: ", single_gaussian_fit$estimate[1], "\n")
cat("Estimated Standard Deviation: ", single_gaussian_fit$estimate[2], "\n")
cat("Log Likelihood: ", log_likelihood, "\n")
cat("BIC: ", bic_single, "\n")
cat("AIC: ", aic_single, "\n")

# Create a histogram of the cleaned data and overlay the fitted Gaussian
hist(data_coffee, probability = TRUE,
     main = "Coffee bean weight histogram fitted with a single Gaussian")
curve(dnorm(x, mean = single_gaussian_fit$estimate[1], sd = single_gaussian_fit$estimate[2]), 
      col = "blue", lwd = 2, add = TRUE)
##############################################################################

# EM algorithm for Gaussian Mixture Models

#Initialization
k <- 2  # Number of components
max_iterations <- 1000
tolerance <- 1e-10

# Initialize parameters
init_weights <- c(0.8, 0.2)  # Equal mixture weights
#common_mean <- 135  # Set a common initial mean value for both components
init_mean <- c(mean(data_coffee), mean(data_coffee))
#init_mean <- rep(common_mean, k)  # Initialize means equally
init_variances <- c(var(data_coffee), var(data_coffee))

n <- length(data_coffee)
n
# Initialize log likelihood with a large negative value
log_likelihood_old <- -Inf

for (iteration in 1:max_iterations) {
  # E-step
  pdf_values <- matrix(0, n, k)
  for (i in 1:k) {
    pdf_values[, i] <- dnorm(data_coffee, mean = init_mean[i], 
                             sd = sqrt(init_variances[i]))
  }
  responsibility_test_1 <- cbind(pdf_values[, 1] * init_weights[1], 
                                 pdf_values[, 2] * init_weights[2])
  responsibility_2 <- responsibility_test_1 / rowSums(responsibility_test_1)
  
  # M-step
  Nk <- colSums(responsibility_2)
  init_weights <- Nk / n
  init_mean <- colSums(data_coffee * responsibility_2) / Nk
  init_variances <- colSums((data_coffee - init_mean[i])^2 * responsibility_2) / Nk[i]
  
  
  # Calculate the log likelihood
  log_likelihood_new <- sum(log(rowSums(pdf_values * init_weights)))
  
  # Check for convergence
  if (abs(log_likelihood_new - log_likelihood_old) < tolerance) {break}
  
  log_likelihood_old <- log_likelihood_new
  #cat("Iteration ", iteration, " - Mean 1: ", init_mean[1], " Mean 2: ", init_mean[2], "\n")
  
}
# Calculate the BIC
n <- length(data_coffee)

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

x <- seq(min(data_coffee), max(data_coffee), length = 1000)

# Density functions for the components
den1 <- dnorm(x, mean = init_mean[2], sd = sqrt(init_variances[1]))

den2 <- dnorm(x, mean = init_mean[2], sd = sqrt(init_variances[2]))

hist(data_coffee, prob = TRUE,
     main = "SBP histogram and EM Fit of Two-Component Mixture", xlab = "bp_data_DBP")

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
############################################################################

# Calculate eta based on estimated variances
estimated_eta <- init_variances[1] / init_variances[2]
estimated_eta

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
outliers_eta <- sapply(data_coffee, function(x) {
  if (goodness_probability_eta(x, a, b, c, estimated_eta) > 0.5) {
    return("Not outlier")
  } else {
    return("Outlier")
  }
})

#outliers_eta

###########################################################################
# Define custom colors for the pie chart
custom_colors <- c("lightblue", "lightcoral")

# Count the number of "Not outlier" and "Outlier" classifications
classification_counts <- table(outliers_eta)

# Create a pie chart
pie(classification_counts, 
    labels = paste(names(classification_counts), ": ", classification_counts),
    main = "Classification of coffee bean weight data points ",
    col = custom_colors) 

# Add a legend
legend("topright", legend = names(classification_counts), fill = custom_colors)

#############################################################################




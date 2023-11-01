library(mclust)
###############################################################################
# Codes to build illustrative graphs
###############################################################################
# Generate x values
x <- seq(-4, 4, length = 1000)

# Generate y values for each Gaussian distribution
y1 <- dnorm(x, mean = 0, sd = 1)
y2 <- dnorm(x, mean = 0, sd = sqrt(5))
y3 <- dnorm(x, mean = 0, sd = sqrt(10))

# Plot the standard normal distribution
plot(x, y1, type = "l", xlab = "x", ylab = "Density", 
     main = "Gaussian distributions with different variances",
     col = "blue", lwd = 2)

# Add lines for the other two distributions
lines(x, y2, col = "red", lwd = 2)
lines(x, y3, col = "green", lwd = 2)

# Add a vertical line at x=0
abline(v = 0, col = "red", lty = 2)

# Add a legend
legend("topleft", legend = c("Variance = 1", "Variance = 5", 
                             "Variance = 10", "Mean (x=0)"), 
       col = c("blue", "red", "green", "red"), lty = c(1, 1, 1, 2),
       lwd = c(2, 2, 2, 2))

##############################################################################
waiting <- faithful $ waiting
n <- length (waiting)
waiting.Mclust <- Mclust (waiting ,model ="V",G=2)

# Plot Densities:
x <- seq(from = min(waiting), to = max(waiting), length = 1000)
den1 <- dnorm(x, mean = waiting.Mclust$parameters$mean[1],
              sd = sqrt(waiting.Mclust$parameters$variance$sigmasq[1]))
den2 <- dnorm(x, mean = waiting.Mclust$parameters$mean[2],
              sd = sqrt(waiting.Mclust$parameters$variance$sigmasq[2]))
tau1 <- waiting.Mclust$parameters$pro[1]
tau2 <- waiting.Mclust$parameters$pro[2]
dens <- tau1 * den1 + tau2 * den2

# Create a histogram of the data
hist_data <- hist(waiting, plot = FALSE)

# Plot the histogram with specified breaks
hist(waiting, breaks = hist_data$breaks, probability = TRUE, col = "gray",
     xlab = "y", ylab = "Probability Density",
     main = "Probability density function for a two-component Gaussian mixture model")

# Overlay the density plot
lines(x, dens, type = "l", lwd = 2)
lines(x, tau1 * den1, col = "red")
lines(x, tau2 * den2, col = "blue")

# Add legend
legend(x = min(x), y = max(dens), legend = c("Mixture", "Component 1", "Component 2"),
       col = c("black", "red", "blue"), lty = c(1, 1, 1), lwd = c(2, 1, 1))
#####################################################################################
# Section 1.5 of Chapter 1: illustrating how the tails miss the outliers 
# Hence the need for mixture modelling 
#############################################################################
# Set the seed for reproducibility
set.seed(123)

# Number of data points
n <- 1000

# Simulate a mixture of Gaussian distributions (contaminated model)
# Parameters for the mixture
mu1 <- 0
sigma1 <- 1
mu2 <- 5
sigma2 <- 1
proportion_outliers <- 0.1

# Simulate the data
data <- c(rnorm(n * (1 - proportion_outliers), mean = mu1, sd = sigma1),
          rnorm(n * proportion_outliers, mean = mu2, sd = sigma2))

# Fit a single Gaussian to the data
fitted_mean <- mean(data)
fitted_sd <- sd(data)

# Plot a histogram of the data and the fitted Gaussian
hist(data, breaks = 30, main = "Simulation of Contaminated Data", xlab = "Value")
curve(dnorm(x, mean = fitted_mean, sd = fitted_sd) * n * diff(hist(data)$breaks)[1],
      col = "blue", lwd = 2, add = TRUE)


# Add legend
legend("topright", legend = c("Data", "Fitted Gaussian"), col = c("black", "blue"), lwd = c(1, 2))

# Display the estimated parameters
cat("Estimated Mean:", fitted_mean, "\n")
cat("Estimated Standard Deviation:", fitted_sd, "\n")

################################################################################




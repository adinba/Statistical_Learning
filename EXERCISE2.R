########################################################### STATISTICAL LEARNING MIDTERM
########################### EXERCISE 2
###########################################################

#### Set random seed 
set.seed(122)

#### Library
library(ggplot2)

######### Question 1)

# Define BIC calculation function
bic <- function(residual_sum_squares, num_params) {
  n * log(residual_sum_squares / n) + num_params * tau**2
}

# Define EBIC calculation function
ebic <- function(residual_sum_squares, num_params) {
  n * log(residual_sum_squares / n) + num_params * ln(n) + num_params*ln(p)
}

forward_backward_bic <- function(X, y, tau, criteria = bic, max_iter = 100,cv_criteria = "fix_point" ,tol = 1e-100) {
  # X: Design matrix
  # y: Response vector
  # tau: Penalty parameter for BIC
  # criteria: Objectif function
  # max_iter: Maximum iterations for convergence
  # cv_criteria : Criteria of convergence
  # tol: Convergence tolerance
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize model with no variables
  J <- integer(0)  # Start with an empty set of selected variables
  J_minus_1 <- integer(0)
  bic_best <- Inf
  converged <- FALSE
  iter <- 0
  
  while (!converged) {
    J_minus_2 = J_minus_1
    J_minus_1 = J
    
    iter <- iter + 1
    improved <- FALSE
    
    if(length(J) != 0){
      model <- lm(y ~ X[, J] - 1)
      rss <- sum(resid(model)^2)  # Residual sum of squares
      old_bic <- bic(rss, length(J) + 1)
    }
    else{
      old_bic = 1000000000000000000 # Setting the old bic very high 
    }
      
    # Forward step
    forward_bic <- rep(Inf, p)  # Default to Inf for indices not considered
    for (j in setdiff(1:p, J)) {  # Check variables not in J
      model <- lm(y ~ X[, c(J, j)] - 1)  # Fit model with J and j
      rss <- sum(resid(model)^2)  # Residual sum of squares
      forward_bic[j] <- bic(rss, length(J) + 1)
    }
    # Find the variable with the best BIC reduction
    j_add <- which.min(forward_bic)
    bic_add <- forward_bic[j_add]
    
    
    # If adding a variable improves BIC, update model
    if (bic_add < bic_best) {
      J <- c(J, j_add)
      improved <- TRUE
    }
    
    # Backward step (only if J is at least of size 2)
    if (length(J) > 1) {
      backward_bic <- rep(Inf, length(J))  # Default to Inf for indices not considered
      for (k in seq_along(J)) {  # Check variables in J
        model <- lm(y ~ X[, J[-k]] - 1)  # Fit model without J[k]
        rss <- sum(resid(model)^2)
        backward_bic[k] <- bic(rss, length(J) - 1)
      }
      # Find the variable with the best BIC reduction to remove
      k_remove <- which.min(backward_bic)
      bic_remove <- backward_bic[k_remove]
      
      # If removing a variable improves BIC, update model
      if (bic_remove < bic_best) {
        J <- J[-k_remove]
        bic_best <- bic_remove
        improved <- TRUE
      }
    }
    
    # Check for convergence
    if(cv_criteria == "fix_point"){
      if(length(J_minus_2) >= 1 && length(J_minus_2) == length(J)){
        converged <- all(J == J_minus_2)
      }
    }
    else if(cv_criteria == "close_enough"){
      #print(converged)
      #print(improved)
      #print(bic_best)
      #print(old_bic)
      #print(abs(bic_best - old_bic))
        
      converged <- !improved || abs(bic_best - old_bic) < tol
      if(bic_best == -Inf && old_bic ==-Inf){
        converged <- TRUE
      }
    }
    else if(cv_criteria == "maximum_iteration"){
      if(iter >= max_iter){
        converged <- TRUE   
      }
    }
  }

  # Final model fit (only if J is not empty)
  if (length(J) > 0) {
    final_model <- lm(y ~ X[, J] - 1)
  } else {
    final_model <- NULL
  }
  
  return(list(selected_variables = J, model = final_model, bic = bic_best))
}

# Example usage with a random dataset
n <- 100
p <- 100
X <- matrix(rnorm(n * p), n, p)
y <- X[, c(1, 3, 5)] %*% c(2, -1, 0.5) + rnorm(n)  # True model uses variables 1, 3, 5
tau <- log(n)  # Common choice for BIC penalty

# Examples of the use of the various convergence criteria.
forward_backward_bic(X,y,tau,cv_criteria = "fix_point")
forward_backward_bic(X,y,tau,cv_criteria = "close_enough",tol = 1e-4)
forward_backward_bic(X,y,tau,cv_criteria = "maximum_iteration", max_iter = 10)

# Examples of the use with ebic objectif criteria.
forward_backward_bic(X,y,tau,criteria = ebic,cv_criteria = "maximum_iteration", max_iter = 10)

######### Question 2)

simulate_observations <- function(n, p, s, sigma, amplitude_min = 0.5) {
  # n: number of observations
  # p: number of basis functions
  # s: sparsity level (number of non-zero coefficients)
  # sigma: standard deviation of noise
  # amplitude_min: minimum amplitude for non-zero coefficients
  
  # Generate random points ti in [0, 1]
  t <- sort(runif(n, 0, 1))
  
  # Construct the design matrix X with trigonometric basis functions
  X <- matrix(0, n, p)
  X[, 1] <- 1  # First basis function is constant
  
  # Generate trigonometric basis functions
  r <- (p - 1) / 2  # Number of cosine/sine pairs
  for (j in 1:r) {
    X[, 2 * j] <- cos(2 * pi * j * t)
    X[, 2 * j + 1] <- sin(2 * pi * j * t)
  }
  
  # Create sparse coefficient vector theta*
  theta_star <- rep(0, p)
  non_zero_indices <- sample(1:p, s)  # Randomly select s non-zero entries
  theta_star[non_zero_indices] <- runif(s, amplitude_min, amplitude_min + 1)  # Random amplitudes
  
  # Generate response vector y with noise
  noise <- rnorm(n, mean = 0, sd = sigma)
  y <- X %*% theta_star + noise
  
  return(list(X = X, y = y, theta_star = theta_star, t = t))
}

# Example usage
n <- 100  # Number of observations
p <- 11   # Number of basis functions (e.g., 5 cosine/sine pairs + 1 constant)
s <- 3    # Sparsity level
sigma <- 0.5  # Noise standard deviation

sim_data <- simulate_observations(n, p, s, sigma)
X <- sim_data$X  # Design matrix
y <- sim_data$y  # Observed response vector
theta_star <- sim_data$theta_star  # True sparse coefficients
t <- sim_data$t  # Observation points

print(forward_backward_bic(X,y,t))


######### Question 3)

# Simulation function for different n, p, s, tau
run_simulations <- function(n_values, p_values, s_values, sigma, criteria = bic, tau_values = c(1), tau_function_of_n = FALSE,f_tau, repetitions = 1,cv_criteria = "fix_point") {
  results <- data.frame()
  
  for (n in n_values) {
    print(n)
    for (p in p_values) {
      for (s in s_values) {
        for (tau in tau_values) {
          if(tau_function_of_n == TRUE){
            tau = f_tau(n)
          }
          for (rep in 1:repetitions) {
            # Simulate observations
            sim_data <- simulate_observations(n, p, s, sigma)
            X <- sim_data$X
            y <- sim_data$y
            theta_star <- sim_data$theta_star
            
            # Apply forward-backward BIC estimator
            result <- forward_backward_bic(X, y, tau,criteria = criteria, cv_criteria = cv_criteria)
            
            # Store results
            selected_vars <- length(result$selected_variables)
            correct_vars <- sum(result$selected_variables %in% which(theta_star != 0))
            bic_value <- result$bic
            
            results <- rbind(results, data.frame(
              n = n,
              p = p,
              s = s,
              tau = tau,
              repetition = rep,
              selected_vars = selected_vars,
              correct_vars = correct_vars,
              bic_value = bic_value
            ))
          }
        }
      }
    }
  }
  return(results)
}


########### TABLE NUMBER 1 in report
n_values <- c(10,(3:7)*8)       # Different sample sizes
p_values <- (3:6)*3   # Different numbers of basis functions
s_values <- c(3,7)       # Different levels of sparsity
tau_values <- c(log(50))  # Fixed penalty values
sigma <- 0.5  # Noise level
results <- run_simulations(n_values, p_values, s_values, sigma, tau_values)
res <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
write.csv(res, "~/3A/Statistical_Learning/Midterm/result.csv", row.names = FALSE)
###########

########### GRAPH 1 in report
n_values <- c(20,30,50,80)       # Different sample sizes
p_values <- c(20,30,50,80)   # Different numbers of basis functions
s_values <- c(10:20)       # Different levels of sparsity
sigma <- 0.5  # Noise level

f_tau <- function(n) log(n)
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau1 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) log(n)**2
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau2 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) sqrt(log(n))
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau3 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) log(n)/2
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau4 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) log(n)*2
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau5 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) n
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau6 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

# Add a column indicating the tau value in each data frame
results_tau1$tau <- "penalty = ln(n)"
results_tau2$tau <- "penalty = ln(n)**2"
results_tau3$tau <- "penalty = sqrt(ln(n))"
results_tau4$tau <- "penalty = ln(n)/2"
results_tau5$tau <- "penalty = ln(n)*2"
results_tau6$tau <- "penalty = n"

# Add the ratio column to each data frame
results_tau1$ratio_correct_selected <- with(results_tau1, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau2$ratio_correct_selected <- with(results_tau2, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau3$ratio_correct_selected <- with(results_tau3, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau4$ratio_correct_selected <- with(results_tau4, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau5$ratio_correct_selected <- with(results_tau5, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau6$ratio_correct_selected <- with(results_tau6, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))

# Combine the data frames into one
results <- rbind(results_tau1, results_tau2, results_tau3, results_tau4, results_tau5, results_tau6)

# Plot correct_vars vs selected_vars with facets for tau and grid for n and p
ggplot(results, aes(x = selected_vars, y = correct_vars, color = factor(n))) +
  geom_point(size = 2, alpha = 0.7) +  # Increase point size and adjust transparency for clarity
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Add a dashed trend line for lm
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +  # Add the f(x) = x line
  facet_grid(p ~ tau, scales = "free") +
  labs(
    title = "Correctly Selected vs. Total Selected Variables",
    x = "Total Selected Variables",
    y = "Correctly Selected Variables",
    color = "Sample Size (n)"
  ) +
  theme_minimal(base_size = 14) +  # Use a larger base font size for readability
  theme(
    panel.spacing = unit(1, "lines"),      # Increase space between panels
    strip.background = element_rect(fill = "grey90", color = "grey50"),  # Different background for facet labels
    strip.text = element_text(face = "bold", size = 12),  # Bold and increase size of facet labels
    panel.grid.major = element_line(color = "grey80"),    # Lighter grid lines
    panel.grid.minor = element_blank(),                   # Remove minor grid lines
    plot.title = element_text(face = "bold", hjust = 0.5) # Center and bold title
  )

########### GRAPH 2 in report
# Configuration 1
n_values <- sample(x = 1:1000, size = 3)       # Different random sample sizes
p_values <- sample(x = 700:1000, size = 3)   # Different random numbers of basis functions
s_values <- sample(x = 1:700, size = 5)       # Different random levels of sparsity
sigma <- 0.7  # Random noise level

# Configuration 2
n_values <- sample(x = 1:30, size = 3)       # Different random sample sizes
p_values <- sample(x = 10:20, size = 3)   # Different random numbers of basis functions
s_values <- sample(x = 1:10, size = 5)       # Different random levels of sparsity
sigma <- 0.7  # Random noise level


f_tau <- function(n) log(n)
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau1 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) log(n)**2
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau2 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) sqrt(log(n))
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau3 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) log(n)/2
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau4 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) log(n)*2
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau5 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

f_tau <- function(n) n
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE,f_tau = f_tau)
results_tau6 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

# Add a column indicating the tau value in each data frame
results_tau1$tau <- "penalty = ln(n)"
results_tau2$tau <- "penalty = ln(n)**2"
results_tau3$tau <- "penalty = sqrt(ln(n))"
results_tau4$tau <- "penalty = ln(n)/2"
results_tau5$tau <- "penalty = ln(n)*2"
results_tau6$tau <- "penalty = n"

# Add the ratio column to each data frame
results_tau1$ratio_correct_selected <- with(results_tau1, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau2$ratio_correct_selected <- with(results_tau2, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau3$ratio_correct_selected <- with(results_tau3, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau4$ratio_correct_selected <- with(results_tau4, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau5$ratio_correct_selected <- with(results_tau5, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_tau6$ratio_correct_selected <- with(results_tau6, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))

# Combine the data frames into one
results <- rbind(results_tau1, results_tau2, results_tau3, results_tau4, results_tau5, results_tau6)


# Distribution plot for the ratio of correct to selected variables by tau
ggplot(results, aes(x = ratio_correct_selected, fill = factor(tau))) +
  geom_density(alpha = 0.6) +  # Plot the distribution with transparency
  facet_wrap(~ tau, scales = "free_y") +
  labs(
    title = "Distribution of Selection Accuracy Ratio by BIC Penalty (Tau)",
    x = "Ratio of Correct to Selected Variables",
    y = "Density",
    fill = "Tau Value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


###########

########### GRAPH 3 in report
n_values <- c(10)       # Different sample sizes
p_values <- 20 + c(1:7)*60   # Different numbers of basis functions
s_values <- (1:5)*10      # Different levels of sparsity
sigma <- 0.5  # Noise level

f_tau <- function(n) log(n)
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau)
results_bic <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

results <- run_simulations(n_values, p_values, s_values, sigma, criteria = ebic, tau_function_of_n = TRUE, f_tau = f_tau)
results_ebic <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)

# Add the ratio column to each data frame
results_bic$ratio_correct_selected <- with(results_bic, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_ebic$ratio_correct_selected <- with(results_ebic, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))

# Plot for BIC ratio of correct to selected variables vs true sparsity level s
ggplot(results_bic, aes(x = s, y = ratio_correct_selected, color = factor(n))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Add a dashed trend line for lm
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +  # Add the f(x) = x line
  facet_grid(. ~ p, scales = "free_y") +
  labs(
    title = "BIC: Ratio of Correct to Selected Variables vs True Sparsity (s)",
    x = "True Sparsity (s)",
    y = "Ratio of Correct to Selected Variables",
    color = "Sample Size (n)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/ex_2_graph_4.png", width = 600)

# Plot for EBIC ratio of correct to selected variables vs true sparsity level s
ggplot(results_ebic, aes(x = s, y = ratio_correct_selected, color = factor(n))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Add a dashed trend line for lm
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +  # Add the f(x) = x line
  facet_grid(. ~ p, scales = "free_y") +
  labs(
    title = "EBIC: Ratio of Correct to Selected Variables vs True Sparsity (s)",
    x = "True Sparsity (s)",
    y = "Ratio of Correct to Selected Variables",
    color = "Sample Size (n)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/ex_2_graph_5.png", width = 600)


# Label each dataset with its respective criterion
results_bic$criterion <- "BIC"
results_ebic$criterion <- "EBIC"

# Combine the datasets into one
results_combined_criteria <- rbind(results_bic, results_ebic)

# Plot the distribution of selection accuracy ratio by BIC and EBIC criteria
ggplot(results_combined_criteria, aes(x = ratio_correct_selected, fill = criterion)) +
  geom_density(alpha = 0.6) +  # Plot the distribution with transparency
  facet_wrap(~ criterion, scales = "free_y") +
  labs(
    title = "Distribution of Selection Accuracy Ratio (BIC vs EBIC)",
    x = "Ratio of Correct to Selected Variables",
    y = "Density",
    fill = "Criterion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/ex_2_graph_6.png", width = 600)
###########

########### GRAPH 3

#################### Setting 1) n>p
n_values <- 30       # Different sample sizes
p_values <- 20   # Different numbers of basis functions
s_values <- c(10:16,19)      # Different levels of sparsity
sigma <- 1.5  # Noise level

f_tau <- function(n) log(n)

# Label each dataset with its respective convergence criterion
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "fix_point")
results_crit1 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit1$convergence_criterion <- "Fix point"

results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "close_enough")
results_crit2 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit2$convergence_criterion <- "Close enough"

results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "maximum_iteration")
results_crit3 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit3$convergence_criterion <- "Maximum iteration"

# Add the ratio column to each data frame
results_crit1$ratio_correct_selected <- with(results_crit1, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_crit2$ratio_correct_selected <- with(results_crit2, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_crit3$ratio_correct_selected <- with(results_crit3, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))

# Combine the datasets into one
results_combined <- rbind(results_crit1, results_crit2, results_crit3)

# Density plot for the ratio of correct to selected variables by convergence criterion
ggplot(results_combined, aes(x = ratio_correct_selected, fill = convergence_criterion)) +
  geom_density(alpha = 0.6) +  # Add density plot with transparency for overlapping areas
  facet_wrap(~ convergence_criterion, scales = "free_y") +
  labs(
    title = "Density of the Ratio of Number of Correct to Selected Variables",
    x = "Ratio of Correct to Selected Variables",
    y = "Density",
    fill = "Convergence Criterion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/ex_2_graph_7.png", width = 600)


#################### Setting 2) n=p
n_values <- 25       # Different sample sizes
p_values <- 25   # Different numbers of basis functions
s_values <- c(2:12)*2      # Different levels of sparsity
sigma <- 0.5  # Noise level

f_tau <- function(n) log(n)

# Label each dataset with its respective convergence criterion
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "fix_point")
results_crit1 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit1$convergence_criterion <- "Fix point"

results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "close_enough")
results_crit2 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit2$convergence_criterion <- "Close enough"

results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "maximum_iteration")
results_crit3 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit3$convergence_criterion <- "Maximum iteration"

# Add the ratio column to each data frame
results_crit1$ratio_correct_selected <- with(results_crit1, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_crit2$ratio_correct_selected <- with(results_crit2, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_crit3$ratio_correct_selected <- with(results_crit3, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))

# Combine the datasets into one
results_combined <- rbind(results_crit1, results_crit2, results_crit3)

# Density plot for the ratio of correct to selected variables by convergence criterion
ggplot(results_combined, aes(x = ratio_correct_selected, fill = convergence_criterion)) +
  geom_density(alpha = 0.6) +  # Add density plot with transparency for overlapping areas
  facet_wrap(~ convergence_criterion, scales = "free_y") +
  labs(
    title = "Density of the Ratio of Number of Correct to Selected Variables",
    x = "Ratio of Correct to Selected Variables",
    y = "Density",
    fill = "Convergence Criterion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/ex_2_graph_8.png", width = 600)

#################### Setting 1) n>p
n_values <- 15       # Different sample sizes
p_values <- 100   # Different numbers of basis functions
s_values <- (6:13)*5      # Different levels of sparsity
sigma <- 0.5  # Noise level

f_tau <- function(n) log(n)

# Label each dataset with its respective convergence criterion
results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "fix_point")
results_crit1 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit1$convergence_criterion <- "Fix point"

results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "close_enough")
results_crit2 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit2$convergence_criterion <- "Close enough"

results <- run_simulations(n_values, p_values, s_values, sigma, tau_function_of_n = TRUE, f_tau = f_tau,cv_criteria = "maximum_iteration")
results_crit3 <- aggregate(cbind(selected_vars, correct_vars, bic_value) ~ n + p + s + tau, data = results, mean)
results_crit3$convergence_criterion <- "Maximum iteration"

# Add the ratio column to each data frame
results_crit1$ratio_correct_selected <- with(results_crit1, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_crit2$ratio_correct_selected <- with(results_crit2, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))
results_crit3$ratio_correct_selected <- with(results_crit3, ifelse(selected_vars == 0, 0, correct_vars / selected_vars))

# Combine the datasets into one
results_combined <- rbind(results_crit1, results_crit2, results_crit3)

# Density plot for the ratio of correct to selected variables by convergence criterion
ggplot(results_combined, aes(x = ratio_correct_selected, fill = convergence_criterion)) +
  geom_density(alpha = 0.6) +  # Add density plot with transparency for overlapping areas
  facet_wrap(~ convergence_criterion, scales = "free_y") +
  labs(
    title = "Density of the Ratio of Number of Correct to Selected Variables",
    x = "Ratio of Correct to Selected Variables",
    y = "Density",
    fill = "Convergence Criterion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/ex_2_graph_9.png", width = 600)

###########

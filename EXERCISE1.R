########################################################### STATISTICAL LEARNING MIDTERM
########################### EXERCISE 1
###########################################################

#### Set random seed 
set.seed(122)

#### Library
library(ggplot2)

######### Question 1)
simulate_observations <- function(n, sigma, min_jump, max_jump, num_jumps = NULL, random_jumps = TRUE, seed = 123) {
  set.seed(seed)
  
  # Grid points t_i
  t <- seq(0, 1, length.out = n)
  
  # Determine the number of jumps
  if (is.null(num_jumps)) {
    num_jumps <- floor(0.1 * n)  # Default to 10% of n if not specified
  } else if (random_jumps) {
    num_jumps <- sample(1:num_jumps, 1)  # Choose a random number of jumps within the specified range
  }
  
  # Generate theta* with specified number of jumps, constrained to [0,1]
  theta_star <- rep(0, n)
  jump_indices <- sample(1:n, size = num_jumps, replace = FALSE)
  theta_star[jump_indices] <- runif(length(jump_indices), min_jump, max_jump)
  
  # Design matrix X
  X <- outer(1:n, 1:n, function(i, j) as.numeric(i >= j))
  
  # Generate observations
  noise <- rnorm(n, mean = 0, sd = sigma)
  y <- X %*% theta_star + noise
  
  # Ensure we are not jumping after we reach 1.
  #print(theta_star)
  cumulative_theta_star <- cumsum(theta_star)
  #print(cumulative_theta_star)
  theta_star <- ifelse(cumulative_theta_star > 1,0,theta_star)
  #print(theta_star)
  cumulative_theta_star <- ifelse(cumulative_theta_star>1,1,cumulative_theta_star)
  #print(cumulative_theta_star)
  
  # Ensure y is within [0,1]
  y <- pmax(0, pmin(1, y))
  
  
  list(y = y, theta_star = theta_star, X = X, jump_indices = jump_indices,t = t,cumulative_theta_star= cumulative_theta_star)
}

# Example usage
sim_data <- simulate_observations(n = 39, sigma = 0.01, min_jump = 0.45, max_jump = 0.46, num_jumps = 10, random_jumps = FALSE)
df <- data.frame(t = sim_data$t,y = sim_data$y,theta_star = sim_data$theta_star,cumulative_theta_star = sim_data$cumulative_theta_star)


# Create plots of y
plot_example <- function(sim_data, title="") {
  ggplot(sim_data, aes(x=t,y = y)) +
    geom_point(color = "blue") +
    labs(title = title, x = "t",y = "y") +
    theme_minimal()
}

# Create plots of y with true values of theta star.
plot_example_2 <- function(sim_data, title = "") {
  # Identify jump positions: indices where theta_star changes from zero to nonzero or vice versa
  jump_positions <- which(sim_data$theta_star != 0)
  
  # Create the plot
  ggplot(sim_data, aes(x = t, y = y)) +
    # Add noisy observations
    geom_point(color = "blue") +
    # Add vertical lines at the jump positions
    geom_vline(xintercept = sim_data$t[jump_positions], color = "red", linetype = "dashed", alpha = 0.7) +
    # Add piecewise constant line for cumulative_theta
    geom_step(aes(x = t, y = cumulative_theta_star), color = "green", size = 1, direction = "hv") +
    # Add labels and theme
    labs(title = title, x = "t", y = "y") +
    theme_minimal()
}



# Simulate and plot examples with different parameters
sim1 <- simulate_observations(n = 100, sigma = 0.01, min_jump = 0.1, max_jump = 0.2, num_jumps = 5, random_jumps = FALSE)
sim2 <- simulate_observations(n = 100, sigma = 0.1, min_jump = 0.1, max_jump = 0.3, num_jumps = 8, random_jumps = TRUE)
sim3 <- simulate_observations(n = 500, sigma = 0.05, min_jump = 0.1, max_jump = 0.11, num_jumps = 10, random_jumps = FALSE)
sim4 <- simulate_observations(n = 10, sigma = 0.05, min_jump = 0.2, max_jump = 0.21, num_jumps = 1, random_jumps = FALSE)
sim5 <- simulate_observations(n = 500, sigma = 0.05, min_jump = 0.3, max_jump = 0.31, num_jumps = 2, random_jumps = FALSE)
sim6 <- simulate_observations(n = 500, sigma = 0.5, min_jump = 0.01, max_jump = 0.02, num_jumps = 100, random_jumps = FALSE)
sim7 <- simulate_observations(n = 500, sigma = 0.03, min_jump = 0.7, max_jump = 0.8, num_jumps = 1, random_jumps = FALSE)
sim8 <- simulate_observations(n = 500, sigma = 0.17, min_jump = 0.1, max_jump = 0.11, num_jumps = 4, random_jumps = FALSE)
sim9 <- simulate_observations(n = 500, sigma = 0.001, min_jump = 0.05, max_jump = 0.06, num_jumps = 20, random_jumps = FALSE)
sim10 <- simulate_observations(n = 500, sigma = 0.09, min_jump = 0.01, max_jump = 0.04, num_jumps = 10, random_jumps = TRUE)
sim11 <- simulate_observations(n = 100, sigma = 0.00001, min_jump = 0.0001, max_jump = 0.0002, num_jumps = 100, random_jumps = FALSE)
sim12 <- simulate_observations(n = 500, sigma = 0.03, min_jump = 0.0009, max_jump = 0.001, num_jumps = 499, random_jumps = FALSE)
sim13 <- simulate_observations(n = 500, sigma = 1, min_jump = 0.1, max_jump = 0.11, num_jumps = 4, random_jumps = TRUE)
sim14 <- simulate_observations(n = 500, sigma = 0.04, min_jump = 0.05, max_jump = 0.06, num_jumps = 19, random_jumps = FALSE)
sim15 <- simulate_observations(n = 50, sigma = 0.09, min_jump = 0.02, max_jump = 0.03, num_jumps = 50, random_jumps = TRUE)
sim16 <- simulate_observations(n = 500, sigma = 0.4, min_jump = 1, max_jump = 1, num_jumps = 1, random_jumps = FALSE)

# Create data-frames to plot
df1 <- data.frame(t = sim1$t,y = sim1$y,theta_star = sim1$theta_star,cumulative_theta_star = sim1$cumulative_theta_star)
df2 <- data.frame(t = sim2$t,y = sim2$y,theta_star = sim2$theta_star,cumulative_theta_star = sim2$cumulative_theta_star)
df3 <- data.frame(t = sim3$t,y = sim3$y,theta_star = sim3$theta_star,cumulative_theta_star = sim3$cumulative_theta_star)
df4 <- data.frame(t = sim4$t,y = sim4$y,theta_star = sim4$theta_star,cumulative_theta_star = sim4$cumulative_theta_star)
df5 <- data.frame(t = sim5$t,y = sim5$y,theta_star = sim5$theta_star,cumulative_theta_star = sim5$cumulative_theta_star)
df6 <- data.frame(t = sim6$t,y = sim6$y,theta_star = sim6$theta_star,cumulative_theta_star = sim6$cumulative_theta_star)
df7 <- data.frame(t = sim7$t,y = sim7$y,theta_star = sim7$theta_star,cumulative_theta_star = sim7$cumulative_theta_star)
df8 <- data.frame(t = sim8$t,y = sim8$y,theta_star = sim8$theta_star,cumulative_theta_star = sim8$cumulative_theta_star)
df9 <- data.frame(t = sim9$t,y = sim9$y,theta_star = sim9$theta_star,cumulative_theta_star = sim9$cumulative_theta_star)
df10 <- data.frame(t = sim10$t,y = sim10$y,theta_star = sim10$theta_star,cumulative_theta_star = sim10$cumulative_theta_star)
df11 <- data.frame(t = sim11$t,y = sim11$y,theta_star = sim11$theta_star,cumulative_theta_star = sim11$cumulative_theta_star)
df12 <- data.frame(t = sim12$t,y = sim12$y,theta_star = sim12$theta_star,cumulative_theta_star = sim12$cumulative_theta_star)
df13 <- data.frame(t = sim13$t,y = sim13$y,theta_star = sim13$theta_star,cumulative_theta_star = sim13$cumulative_theta_star)
df14 <- data.frame(t = sim14$t,y = sim14$y,theta_star = sim14$theta_star,cumulative_theta_star = sim14$cumulative_theta_star)
df15 <- data.frame(t = sim15$t,y = sim15$y,theta_star = sim15$theta_star,cumulative_theta_star = sim15$cumulative_theta_star)
df16 <- data.frame(t = sim16$t,y = sim16$y,theta_star = sim16$theta_star,cumulative_theta_star = sim16$cumulative_theta_star)

# Plot each simulation
plot_example(df1)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_1.png", width = 600)
plot_example(df2)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2.png", width = 600)
plot_example(df3)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3.png", width = 600)
plot_example(df4)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4.png", width = 600)
plot_example(df5)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5.png", width = 600)
plot_example(df6)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_6.png", width = 600)
plot_example(df7)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_7.png", width = 600)
plot_example(df8)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_8.png", width = 600)
plot_example(df9)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_9.png", width = 600)
plot_example(df10)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_10.png", width = 600)
plot_example(df11)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_11.png", width = 600)
plot_example(df12)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_12.png", width = 600)
plot_example(df13)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_13.png", width = 600)
plot_example(df14)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_14.png", width = 600)
plot_example(df15)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_15.png", width = 600)
plot_example(df16)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_16.png", width = 600)

# Plot each simulation
plot_example_2(df1)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_1.png", width = 600)
plot_example_2(df2)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_2.png", width = 600)
plot_example_2(df3)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_3.png", width = 600)
plot_example_2(df4)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_4.png", width = 600)
plot_example_2(df5)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_5.png", width = 600)
plot_example_2(df6)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_6.png", width = 600)
plot_example_2(df7)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_7.png", width = 600)
plot_example_2(df8)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_8.png", width = 600)
plot_example_2(df9)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_9.png", width = 600)
plot_example_2(df10)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_10.png", width = 600)
plot_example_2(df11)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_11.png", width = 600)
plot_example_2(df12)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_12.png", width = 600)
plot_example_2(df13)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_13.png", width = 600)
plot_example_2(df14)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_14.png", width = 600)
plot_example_2(df15)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_15.png", width = 600)
plot_example_2(df16)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_2_16.png", width = 600)

######### Question 2)

# Function to compute N2(s) for a specific s
compute_N2_single <- function(y, s) {
  if (s == 1) {
    return(0)  # N2(1) is defined as 0
  }
  sum(y[1:(s - 1)]^2)
}

# Function to compute R2(k, j) for specific k and j
compute_R2_single <- function(y, k, j) {
  if (k > j) {
    stop("k must be less than j")
  }
  segment <- y[k:(j - 1)]
  mean_segment <- mean(segment)
  sum((segment - mean_segment)^2)
}

# Example usage
# Simulate observations
sim_data <- simulate_observations(n = 100, sigma = 0.01, min_jump = 0.2, max_jump = 0.21, num_jumps = 4, random_jumps = FALSE)
y <- sim_data$y

# Compute N2(s) and R2(k, j) for specific values
N2_s <- compute_N2_single(y, s = 5)     # Example with s = 5
R2_kj <- compute_R2_single(y, k = 3, j = 10)  # Example with k = 3, j = 10


##### Implementation of the algorithm :
compute_best_jumps <- function(y, tau, s_max) {
  n <- length(y)  # Number of observations
  
  # Initialize dynamic programming tables
  delta <- matrix(Inf, nrow = n, ncol = s_max)  # Minimum error for s jumps
  J <- matrix(0, nrow = n, ncol = s_max)       # Stores jump positions
  
  # Precompute N2 and R2 for all possible intervals
  N2 <- sapply(1:n, function(k) compute_N2_single(y, k))
  R2 <- matrix(NA, nrow = n, ncol = n + 1)
  for (k in 1:n) {
    for (j in (k + 1):(n + 1)) {
      R2[k, j] <- compute_R2_single(y, k, j)
    }
  }
  
  # Step 1: Initialize for s = 1
  for (k in 1:n) {
    delta[k, 1] <- N2[k] + R2[k, n + 1]
  }
  
  # Step 2: Dynamic programming for s >= 2
  for (s in 2:s_max) {
    print(s)
    for (k in s:n) {
      min_value <- Inf
      min_position <- 0
      for (i in s:k) {
        current_value <- delta[i - 1, s - 1] + R2[i, k + 1]
        if (current_value < min_value) {
          min_value <- current_value
          min_position <- i
        }
      }
      delta[k, s] <- min_value
      J[k, s] <- min_position
    }
  }
  
  # Step 3: Compute BIC for each s
  BIC <- numeric(s_max)
  for (s in 1:s_max) {
    BIC[s] <- delta[n, s] + tau^2 * s
  }
  
  # Step 4: Select the best s and recover jumps
  best_s <- which.min(BIC)  # Optimal number of jumps
  
  recover_jumps <- function(k, s) {
    if (s == 0) return(NULL)
    jump <- J[k, s]
    c(recover_jumps(jump - 1, s - 1), jump)
  }
  
  best_jumps <- recover_jumps(n, best_s)  # Optimal jump positions
  
  # Ensure the jump at 0 is excluded
  best_jumps <- best_jumps[best_jumps > 1]
  
  # Return the optimal number of jumps and their positions
  list(best_s = length(best_jumps), best_jumps = best_jumps)
}

res <- compute_best_jumps(y,tau = log(length(y))/length(y),s_max=length(y))

#### We can now build a function to recover the "matter" function f. 
matter_function_f_BIC <- function(jump_points, local_y_predicted) {
  n <- length(local_y_predicted)  # Total number of observations
  print(n)
  if (length(jump_points) == 0) {
    # If there are no jumps, return the global mean
    return(rep(mean(local_y_predicted), n))
  }
  
  # Ensure jump points are sorted and include boundaries
  jump_points <- sort(unique(jump_points))
  jump_points <- c(1, jump_points, n + 1)  # Include start and end boundaries
  
  # Compute the predicted function f
  f_predicted <- numeric(n)
  for (i in seq_len(length(jump_points) - 1)) {
    start <- jump_points[i]
    end <- jump_points[i + 1] - 1
    f_predicted[start:end] <- mean(local_y_predicted[start:end])  # Assign segment mean
  }
  
  return(f_predicted)
}

f_predicted <- matter_function_f_BIC(res$best_jumps,y)

######### Question 3)
###### Experiments

plot_simulated_data <- function(df, title = "") {
  # Base plot with noisy observations
  p <- ggplot(df, aes(x = t, y = y)) +
    geom_point(color = "blue", alpha = 0.7) +  # Noisy observations
    geom_step(aes(y = f_predicted), color = "red", size = 1) +  # Predicted function
    labs(title = title, x = "t", y = "y") +
    theme_minimal()
  
  # Add true function if provided
  if (!is.null(df$cumulative_theta_star)) {
    p <- p + geom_step(aes(y = cumulative_theta_star), color = "green", linetype = "dashed", size = 1)  # True function
  }
  
  return(p)
}


############ Graph 1)
res <- compute_best_jumps(df1$y,tau = log(length(df1$y)),s_max=length(df1$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df$y)
df1["predicted_f"] <- f_predicted
plot_simulated_data(df1)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_1.png", width = 600)

res <- compute_best_jumps(df2$y,tau = log(length(df2$y)),s_max=length(df2$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df2$y)
df2["predicted_f"] <- f_predicted
plot_simulated_data(df2)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_2.png", width = 600)

res <- compute_best_jumps(df3$y,tau = log(length(df3$y)),s_max=length(df3$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df3$y)
df3["predicted_f"] <- f_predicted
plot_simulated_data(df3)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_3.png", width = 600)

res <- compute_best_jumps(df4$y,tau = log(length(df4$y)),s_max=length(df4$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df4$y)
df4["predicted_f"] <- f_predicted
plot_simulated_data(df4)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_4.png", width = 600)

res <- compute_best_jumps(df5$y,tau = log(length(df5$y)),s_max=length(df5$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df5$y)
df5["predicted_f"] <- f_predicted
plot_simulated_data(df5)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_5.png", width = 600)

res <- compute_best_jumps(df6$y,tau = log(length(df6$y)),s_max=length(df6$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df6$y)
df6["predicted_f"] <- f_predicted
plot_simulated_data(df6)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_6.png", width = 600)

res <- compute_best_jumps(df7$y,tau = log(length(df7$y)),s_max=length(df7$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df7$y)
df7["predicted_f"] <- f_predicted
plot_simulated_data(df7)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_7.png", width = 600)

res <- compute_best_jumps(df8$y,tau = log(length(df8$y)),s_max=length(df8$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df8$y)
df8["predicted_f"] <- f_predicted
plot_simulated_data(df8)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_8.png", width = 600)

res <- compute_best_jumps(df9$y,tau = log(length(df9$y)),s_max=length(df9$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df9$y)
df9["predicted_f"] <- f_predicted
plot_simulated_data(df9)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_9.png", width = 600)

res <- compute_best_jumps(df10$y,tau = log(length(df10$y)),s_max=length(df10$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df10$y)
df10["predicted_f"] <- f_predicted
plot_simulated_data(df10)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_10.png", width = 600)

res <- compute_best_jumps(df11$y,tau = log(length(df11$y)),s_max=length(df11$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df11$y)
df11["predicted_f"] <- f_predicted
plot_simulated_data(df11)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_11.png", width = 600)

res <- compute_best_jumps(df12$y,tau = log(length(df12$y)),s_max=length(df12$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df12$y)
df12["predicted_f"] <- f_predicted
plot_simulated_data(df12)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_12.png", width = 600)

res <- compute_best_jumps(df13$y,tau = log(length(df13$y)),s_max=length(df13$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df13$y)
df13["predicted_f"] <- f_predicted
plot_simulated_data(df13)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_13.png", width = 600)

res <- compute_best_jumps(df14$y,tau = log(length(df14$y)),s_max=length(df14$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df14$y)
df14["predicted_f"] <- f_predicted
plot_simulated_data(df14)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_14.png", width = 600)

res <- compute_best_jumps(df15$y,tau = log(length(df15$y)),s_max=length(df15$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df15$y)
df15["predicted_f"] <- f_predicted
plot_simulated_data(df15)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_15.png", width = 600)

res <- compute_best_jumps(df16$y,tau = log(length(df16$y)),s_max=length(df16$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df16$y)
df16["predicted_f"] <- f_predicted
plot_simulated_data(df16)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_3_16.png", width = 600)

############ Graph 2)
res <- compute_best_jumps(df1$y,tau = log(length(df1$y))/length(df1$y),s_max=length(df1$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df$y)
df1["predicted_f"] <- f_predicted
plot_simulated_data(df1)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_1.png", width = 600)

res <- compute_best_jumps(df2$y,tau = log(length(df2$y))/length(df2$y),s_max=length(df2$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df2$y)
df2["predicted_f"] <- f_predicted
plot_simulated_data(df2)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_2.png", width = 600)

res <- compute_best_jumps(df3$y,tau = log(length(df3$y))/length(df3$y),s_max=length(df3$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df3$y)
df3["predicted_f"] <- f_predicted
plot_simulated_data(df3)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_3.png", width = 600)

res <- compute_best_jumps(df4$y,tau = log(length(df4$y))/length(df4$y),s_max=length(df4$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df4$y)
df4["predicted_f"] <- f_predicted
plot_simulated_data(df4)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_4.png", width = 600)

res <- compute_best_jumps(df5$y,tau = log(length(df5$y))/length(df5$y),s_max=length(df5$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df5$y)
df5["predicted_f"] <- f_predicted
plot_simulated_data(df5)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_5.png", width = 600)

res <- compute_best_jumps(df6$y,tau = log(length(df6$y))/length(df6$y),s_max=length(df6$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df6$y)
df6["predicted_f"] <- f_predicted
plot_simulated_data(df6)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_6.png", width = 600)

res <- compute_best_jumps(df7$y,tau = log(length(df7$y))/length(df7$y),s_max=length(df7$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df7$y)
df7["predicted_f"] <- f_predicted
plot_simulated_data(df7)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_7.png", width = 600)

res <- compute_best_jumps(df8$y,tau = log(length(df8$y))/length(df8$y),s_max=length(df8$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df8$y)
df8["predicted_f"] <- f_predicted
plot_simulated_data(df8)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_8.png", width = 600)

res <- compute_best_jumps(df9$y,tau = log(length(df9$y))/length(df9$y),s_max=length(df9$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df9$y)
df9["predicted_f"] <- f_predicted
plot_simulated_data(df9)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_9.png", width = 600)

res <- compute_best_jumps(df10$y,tau = log(length(df10$y))/length(df10$y),s_max=length(df10$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df10$y)
df10["predicted_f"] <- f_predicted
plot_simulated_data(df10)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_10.png", width = 600)

res <- compute_best_jumps(df11$y,tau = log(length(df11$y))/length(df11$y),s_max=length(df11$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df11$y)
df11["predicted_f"] <- f_predicted
plot_simulated_data(df11)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_11.png", width = 600)

res <- compute_best_jumps(df12$y,tau = log(length(df12$y))/length(df12$y),s_max=length(df12$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df12$y)
df12["predicted_f"] <- f_predicted
plot_simulated_data(df12)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_12.png", width = 600)

res <- compute_best_jumps(df13$y,tau = log(length(df13$y))/length(df13$y),s_max=length(df13$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df13$y)
df13["predicted_f"] <- f_predicted
plot_simulated_data(df13)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_13.png", width = 600)

res <- compute_best_jumps(df14$y,tau = log(length(df14$y))/length(df14$y),s_max=length(df14$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df14$y)
df14["predicted_f"] <- f_predicted
plot_simulated_data(df14)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_14.png", width = 600)

res <- compute_best_jumps(df15$y,tau = log(length(df15$y))/length(df15$y),s_max=length(df15$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df15$y)
df15["predicted_f"] <- f_predicted
plot_simulated_data(df15)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_15.png", width = 600)

res <- compute_best_jumps(df16$y,tau = log(length(df16$y))/length(df16$y),s_max=length(df16$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df16$y)
df16["predicted_f"] <- f_predicted
plot_simulated_data(df16)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_4_16.png", width = 600)

############ Graph 3)

tau_function_of_volatility <- function(y){
  return(var(y)/mean(y))
}

res <- compute_best_jumps(df1$y,tau = tau_function_of_volatility(df1$y) ,s_max=length(df1$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df$y)
df1["predicted_f"] <- f_predicted
plot_simulated_data(df1)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_1.png", width = 600)

res <- compute_best_jumps(df2$y,tau = tau_function_of_volatility(df2$y),s_max=length(df2$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df2$y)
df2["predicted_f"] <- f_predicted
plot_simulated_data(df2)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_2.png", width = 600)

res <- compute_best_jumps(df3$y,tau = tau_function_of_volatility(df3$y),s_max=length(df3$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df3$y)
df3["predicted_f"] <- f_predicted
plot_simulated_data(df3)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_3.png", width = 600)

res <- compute_best_jumps(df4$y,tau = tau_function_of_volatility(df4$y),s_max=length(df4$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df4$y)
df4["predicted_f"] <- f_predicted
plot_simulated_data(df4)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_4.png", width = 600)

res <- compute_best_jumps(df5$y,tau = tau_function_of_volatility(df5$y),s_max=length(df5$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df5$y)
df5["predicted_f"] <- f_predicted
plot_simulated_data(df5)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_5.png", width = 600)

res <- compute_best_jumps(df6$y,tau = tau_function_of_volatility(df6$y),s_max=length(df6$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df6$y)
df6["predicted_f"] <- f_predicted
plot_simulated_data(df6)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_6.png", width = 600)

res <- compute_best_jumps(df7$y,tau = tau_function_of_volatility(df7$y),s_max=length(df7$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df7$y)
df7["predicted_f"] <- f_predicted
plot_simulated_data(df7)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_7.png", width = 600)

res <- compute_best_jumps(df8$y,tau = tau_function_of_volatility(df8$y),s_max=length(df8$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df8$y)
df8["predicted_f"] <- f_predicted
plot_simulated_data(df8)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_8.png", width = 600)

res <- compute_best_jumps(df9$y,tau = tau_function_of_volatility(df9$y),s_max=length(df9$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df9$y)
df9["predicted_f"] <- f_predicted
plot_simulated_data(df9)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_9.png", width = 600)

res <- compute_best_jumps(df10$y,tau = tau_function_of_volatility(df10$y),s_max=length(df10$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df10$y)
df10["predicted_f"] <- f_predicted
plot_simulated_data(df10)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_10.png", width = 600)

res <- compute_best_jumps(df11$y,tau = tau_function_of_volatility(df11$y),s_max=length(df11$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df11$y)
df11["predicted_f"] <- f_predicted
plot_simulated_data(df11)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_11.png", width = 600)

res <- compute_best_jumps(df12$y,tau = tau_function_of_volatility(df12$y),s_max=length(df12$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df12$y)
df12["predicted_f"] <- f_predicted
plot_simulated_data(df12)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_12.png", width = 600)

res <- compute_best_jumps(df13$y,tau = tau_function_of_volatility(df13$y),s_max=length(df13$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df13$y)
df13["predicted_f"] <- f_predicted
plot_simulated_data(df13)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_13.png", width = 600)

res <- compute_best_jumps(df14$y,tau = tau_function_of_volatility(df14$y),s_max=length(df14$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df14$y)
df14["predicted_f"] <- f_predicted
plot_simulated_data(df14)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_14.png", width = 600)

res <- compute_best_jumps(df15$y,tau = tau_function_of_volatility(df15$y),s_max=length(df15$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df15$y)
df15["predicted_f"] <- f_predicted
plot_simulated_data(df15)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_15.png", width = 600)

res <- compute_best_jumps(df16$y,tau = tau_function_of_volatility(df16$y),s_max=length(df16$y))
f_predicted <- matter_function_f_BIC(res$best_jumps,df16$y)
df16["predicted_f"] <- f_predicted
plot_simulated_data(df16)
dev.print(device = png, file = "~/3A/Statistical_Learning/Midterm/graph_5_16.png", width = 600)

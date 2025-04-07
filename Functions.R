# R File to hold functions 
library(ggplot2)
suppressWarnings(library(dplyr)) 
library(dplyr) 
library(glue)

# Function to price an American put using a Binomial Tree
binomial_tree_american_put <- function(S0, K, T, r, sigma, M) {
  dt <- T / M  #Time Step Size 
  u <- exp(sigma * sqrt(dt)) # Upturn factor
  d <- exp(-sigma * sqrt(dt)) # Downward factor
  p <- (exp(r * dt) - d) / (u - d) # Risk-neutral probability

  # Stock price tree
  stock_tree <- matrix(0, nrow = M + 1, ncol = M + 1) # Create matrix to store stock prices
  for (j in 0:M) {
    for (i in 0:j) {
      stock_tree[i + 1, j + 1] <- S0 * u^i * d^(j - i) # Number of steps multiplied by upturn and downturn factors
    }
  }
  # Option value tree
  option_tree <- matrix(0, nrow = M + 1, ncol = M + 1) # Matrix to store option values 
  option_tree[, M + 1] <- pmax(K - stock_tree[, M + 1], 0)  # Payoff is max(K - S, 0) at maturity, applied to entire column M+1

  # Work backward through the tree
  for (j in (M - 1):0) { # Go through each column in reverse order starting at M-1  
    for (i in 0:j) { # Go though each row in the column
      continuation <- exp(-r * dt) * (p * option_tree[i + 2, j + 2] + (1 - p) * option_tree[i + 1, j + 2])
      exercise <- pmax(K - stock_tree[i + 1, j + 1], 0)
      option_tree[i + 1, j + 1] <- pmax(continuation, exercise)
    }
  }
  
  return(option_tree[1, 1])  # Option price at root
}

# Function to price an American put using Monte Carlo simulations
price_american_put_longstaff_schwartz_MC <- function(K, M, N, r, S0,sigma, polynomial) {
  
  dt <- 1/M
  discount <- exp(-r * dt)  
  set.seed(123)
  Z <- matrix(rnorm(N * M), nrow = N, ncol = M) # Vectorize Brownian Motion Simulation  
  S <- S0 * exp(sigma * sqrt(dt) * t(apply(Z, 1, cumsum)))
  
  Cash_flow <- matrix(0, nrow = N, ncol = M)
  Cash_flow[, M] <- pmax(K - S[, M], 0) 
  
  # Cash Flows at each time step 
  for (m in M:2) {
    X <- S[, m-1]
    Y <- Cash_flow[, m] * discount
    XY <- cbind(X, Y)
    XY[X > K, ] <- NA  

    if (all(is.na(XY))) {
        Cash_flow[, m-1] <- 0  # Skip regression if no in the money paths
        next
    }
    
    regression <- lm(polynomial, data = as.data.frame(XY))
    
    immediate_exercise <- pmax(K - S[, m-1], 0)
    continuation <- predict(regression, newdata = as.data.frame(X))
    
    full_step <- cbind(continuation, immediate_exercise)
    full_step[immediate_exercise == 0, ] <- NA
    
    result_vector <- ifelse(
      is.na(full_step[, 2]), 0,                       
      ifelse(full_step[, 1] > full_step[, 2], 0, full_step[, 2])
    )
    Cash_flow[, m-1] <- result_vector
  }
  
  # Discounting 
  for (i in 1:nrow(Cash_flow)) {
    for (j in 1:ncol(Cash_flow)) {
      if (Cash_flow[i, j] != 0) { 
        Cash_flow[i, j] <- Cash_flow[i, j] * round(exp(-r * j), 5)
        if (j < ncol(Cash_flow)) {
          Cash_flow[i, (j+1):ncol(Cash_flow)] <- 0
        }
        break
      }
    }
  }
  
  option_price <- mean(rowSums(Cash_flow))
  return(option_price)
}















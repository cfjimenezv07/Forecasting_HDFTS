###############################
# Numbers of factors selection
###############################

### The method of den Reijer et al. 2021 ###

k_criterion <- function(eigenvalue, plot = TRUE)
{
  n = length(eigenvalue)
  Hn = sum(1/(1:n))
  
  eigen_diff = eigenvalue[-n] - eigenvalue[-1]
  lambda_bar = 1/((2:n)*Hn)
  
  n_cross = which(eigen_diff <= lambda_bar)[1]
  
  if(n_cross > 1)
  {
    n_return = n_cross - 1
  } else 
    {
      n_return = 1
    }
  
  if(plot)
  {
    plot(eigen_diff[1:(n_cross+6)], type = "b", xlab = expression("Eigenvalue number" ~ italic(k)), ylab = expression("Eigenvalue" ~~~ lambda))
    lines(lambda_bar[1:(n_cross+6)], type ="b", col = 2)
    abline(v = n_return, lty = 2, col = 4)
    legend("topright", lty = c(1,1), lwd = c(2,2), col = c(1,2), c(expression(lambda ~~ "(eigenvalue)  "), expression(bar(lambda) ~~ "(threshold)  ")))
  }
  
  return(n_return)
}


### Eigenratio $k$ selection method ###

select_K <- function(tau, eigenvalue)
{
  
  k_max = length(eigenvalue)
  k_all = rep(0, k_max-1)
  for(k in 1:(k_max-1))
  {
    k_all[k] = (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
  }
  
  K_hat = which.min(k_all)
  
  return(K_hat)
}


######################################
# Simulations for performance testing 
######################################

library(sde)
library(ftsa)

# standard Brownian Motion [0,1]

BrownMat <- function(N, refinement)
{
  mat <- matrix(nrow = refinement, ncol = N)
  c <- 1
  while(c <= N)
  {
    vec <- BM(N = refinement-1)
    mat[,c] <- vec
    c <- c+1
  }
  return(mat)
}

phi_1 <- function(u)
{
  sin(2*pi*u)
}

phi_2 <- function(u)
{
  cos(2*pi*u)
}

phi_3 <- function(u)
{
  sin(4*pi*u)
}

sim_factor_model <- function(NT = 100, no_grid = 101, seed_number)
{
  set.seed(123 + seed_number)
  
  # generate epsilon(u)
  epsilon_all = BrownMat(N = NT+50, refinement = no_grid)
  
  # factor loading functions
  u_grids = seq(0, 1, length.out = no_grid)
  phi_all = cbind(phi_1(u_grids), phi_2(u_grids), phi_3(u_grids))
  
  # factor time series
  w_all = matrix(0, nrow = 3, ncol = NT+50)
  w_all[1,] = rnorm(NT+50, mean = 0, sd = 0.5)
  w_all[2,] = rnorm(NT+50, mean = 0, sd = 0.25)
  w_all[3,] = rnorm(NT+50, mean = 0, sd = 0.05)
  
  f_all = matrix(0, nrow = 3, ncol = NT+50)
  f_all[,1] = rbind(20, 10, 5)
  for(it in 2:(NT+50))
  {
    f_all[,it] = 0.6*f_all[,it-1] + w_all[,it]
  }
  
  # combine results
  X_all = phi_all %*% f_all + epsilon_all
  
  X_out = X_all[,51:(NT+50)]
  
  return(X_out)
}


### Test performance ###


data_sim = sim_factor_model(NT = 200, no_grid = 300, seed_number = 10)
eigen_sim = eigen(crossprod(t(data_sim))/(nrow(data_sim)*ncol(data_sim)), symmetric = TRUE)$values
  
# eigen ratio method
k_old_method = select_K(tau = 0.001, eigen_sim)
k_old_method

k_new_method = k_criterion(eigen_sim, plot = TRUE)
k_new_method


# Point forecast computation

#General function for point forecasts based on the FMP-ANOVA decomposition

#1. General function for point forecasts obtained  to the functional residuals
# after removing deterministic components from any decomposition approach

#2. Computation of the point forecasts based on the rolling window approach, 
# using the function state_forecast_curves.

#3. Computation of the point forecasts based on the expanding window approach, 
# using the function state_forecast_curves.

###############################################################################
# 1. Computation of the point forecast based on the FMP decomposition
################################################################################
library("demography")
library("ftsa")
library("vars")
source("forecast_Arima.R")

#Auxiliary function to select the functional component scores with the EVR criterion.
select_k <- function(tau, eigenvalue)
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
#############
# fixed_com: matrix of dimension (n_states X n_year) by 2*n_age with the deterministic components for both genders and all states.
# Residuals_f:  matrix of dimension (n_states X n_year) by 2*n_age with the residual components for both genders and all states.
# est_method: denotes the estimation method. lrc denotes the long-range covariance and cov denotes the covariance operator.
# fh: denotes the forecast horizon.
# prediction_method: denotes the prediction method. ARIMA or VAR.
# select_K: denotes the method to select the number of PC to be retained.
# K=6. If a pre-specified number of PC can be set by the username.
# B number of bootstrap samples, 1000 by default


state_forecast_curves <- function(fixed_com, Residuals_f,
                               est_method = c("lrc", "cov"),
                               fh = 30,  
                               prediction_method=c("ARIMA","VAR"),select_K=c("Fixed","EVR"), 
                               K=6, B = 1000){

  med_polish_resi=t(Residuals_f)
  if(est_method == "lrc"){
    # estimate long-run covariance by kernel sandwich estimator
    med_polish_resi_lrc = long_run_covariance_estimation(med_polish_resi)
  }else if(est_method == "cov"){
    # estimate empirical covariance function
    med_polish_resi_lrc = cov(t(med_polish_resi))
  }
  # perform eigen-decomposition
  med_polish_resi_eigen = eigen(med_polish_resi_lrc)
  
  if(select_K=="Fixed"){
    
    retain_component = K
    
  }else if(select_K=="EVR"){
    # determine retained number of components via eigenvalue ratio criterion
    
    lambda_val = med_polish_resi_eigen$values[which(med_polish_resi_eigen$values > 0)]
    retain_component = select_k(tau = 1/log(max(lambda_val[1], length(lambda_val))), eigenvalue = lambda_val)
    
    
  }
  
  # The total variance explained by the number of retained functional principal components.
  var_total_variations = (sum(med_polish_resi_eigen$values[1:retain_component])/sum(med_polish_resi_eigen$values))*100
  
  
  # determine a set of basis function and its scores
  med_polish_resi_basis = as.matrix(med_polish_resi_eigen$vectors[,1:retain_component])
  med_polish_resi_score = crossprod(med_polish_resi, med_polish_resi_basis)
  
  
  med_polish_resi_score_forecast = matrix(NA, retain_component, fh)
  if(prediction_method=="ARIMA"){
    # obtain forecasts of PC scores via auto.arima
    for(ik in 1:retain_component){
      dum = forecast_Arima(object=auto.arima(med_polish_resi_score[,ik]), h = fh, bootstrap = TRUE, npaths = B) 
      med_polish_resi_score_forecast[ik,] = dum$mean
    }
  }else if(prediction_method=="VAR"){
    # obtain forecasts of PC scores via VAR.
    object  <- med_polish_resi_score
    colnames(object) <- 1:dim(object)[2]
    lag=VARselect(y=object,type = "const")$selection[1]
    model_VAR <- VAR(y=object,type = "const",ic="AIC",p=lag)
    pred=predict(model_VAR,n.ahead=fh)$fcst
    for (ik in 1:retain_component) {
      pred1=pred[[ik]]
      med_polish_resi_score_forecast[ik,]=pred1[,1]
    }
    
  }
  
  # Compute the forecast for the functional residuals
  med_polish_resi_forecast = med_polish_resi_basis %*% med_polish_resi_score_forecast
  
  # add back the deterministic components
  Fixed <- t(fixed_com)[,1:fh]
  med_polish_curve_forecast <- med_polish_resi_forecast + Fixed
  
  return(list(med_polish_curve_forecast=med_polish_curve_forecast, 
              med_polish_resi_forecast=  med_polish_resi_forecast,TV = var_total_variations))
  
}


################################################################################
# 2. With rolling window approach
################################################################################
# This function computes the forecast based on a rolling window approach. The testing set is fixed 
library(foreach)
library(doParallel)
max_h = 10
ForecastC <- function(i,fixed_com,Residuals_f,est_method = c("lrc", "cov"),prediction_method=c("ARIMA","VAR"),select_K=c("Fixed","EVR"), K=6)
{
  state=(n_year * i - (n_year - 1)):(n_year * i)
  forecasted_state <- matrix(NA, nrow = n_age * 2, ncol = max_h)
  TV <- matrix(NA, nrow = 1, ncol = max_h)
  for (k in 1:max_h) {
    state_k <- state[k:(n_year - max_h - 1 + k)]
    frc <- state_forecast_curves(fixed_com = fixed_com[state_k, ],
                                Residuals_f = Residuals_f[state_k, ],
                                est_method = est_method,
                                fh = 1, prediction_method=prediction_method,select_K=select_K, K=K, B = 1000)
    forecasted_state[ , k] = frc$med_polish_curve_forecast
    TV[,k]=frc$TV

  }

  return(list(forecasted_state,TV))
}


################################################################################
# 3.  With an expanding window approach
################################################################################
# This function computes the forecast based on a expanding window approach. The testing set increases every iteration. 
library(foreach)
library(doParallel)
max_h = 10
ForecastC_expanding <- function(i,fixed_com,Residuals_f,est_method = c("lrc", "cov"),prediction_method=c("ARIMA","VAR"),select_K=c("Fixed","EVR"), K=6){
  state=(n_year * i - (n_year - 1)):(n_year * i)
  forecasted_state <- matrix(0, nrow = n_age * 2, ncol = max_h)
  TV <- matrix(NA, nrow = 1, ncol = max_h)
  for (k in 1:max_h) {
    state_k<-state[1:(n_year-max_h-1+k)]
    frc=state_forecast_curves(fixed_com=fixed_com[state_k,],
                             Residuals_f=Residuals_f[state_k,],
                             est_method = est_method,
                             fh = k ,prediction_method=prediction_method,select_K=select_K, K=K, B = 1000)
    forecasted_state[,k]=frc$med_polish_curve_forecast[,1]
    TV[,k]=frc$TV
  }
  return(list(forecasted_state,TV))
}


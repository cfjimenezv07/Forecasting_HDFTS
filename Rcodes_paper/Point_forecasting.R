#General function for point forecasts based on the FMP-ANOVA decomposition

#1. General function for point forecasts obtained  to the functional residuals
# after removing deterministic components from the FMP-ANOVA approach

#2. Computation of the point forecasts based on the rolling window approach, 
# using the function Pref_forecast_curves.

#3. Computation of the point forecasts based on the expanding window approach, 
# using the function Pref_forecast_curves.

###############################################################################
# Computation of the point forecast based on the FMP decomposition
################################################################################
library(ftsa)
library(demography)
source("forecast_Arima.R")
Pref_forecast_curves<-function(fixed_com,Residuals_f,
                                 est_method = c("lrc", "cov"),
                                 fh = 30, PI = NULL, B = 1000){
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
  # determine retained number of components via eigenvalue ratio
  ret_component = vector("numeric", length(med_polish_resi_eigen$values) - 1)
  for(ik in 1:(length(med_polish_resi_eigen$values) - 1)){
    ret_component[ik] = med_polish_resi_eigen$values[ik+1]/med_polish_resi_eigen$values[ik]
  }
  retain_component = which.min(ret_component)
  # determine 1st set of basis function and its scores
  med_polish_resi_basis = as.matrix(med_polish_resi_eigen$vectors[,1:retain_component])
  med_polish_resi_score = crossprod(med_polish_resi, med_polish_resi_basis)
  
  # obtain forecasts of PC scores via auto.arima
  med_polish_resi_score_forecast = matrix(NA, retain_component, fh)
  med_polish_resi_score_forecast_boot = array(NA, dim = c(retain_component, fh, B))
  for(ik in 1:retain_component){
    dum = forecast_Arima(object=auto.arima(med_polish_resi_score[,ik]), h = fh, bootstrap = TRUE, npaths = B)
    med_polish_resi_score_forecast[ik,] = dum$mean
    med_polish_resi_score_forecast_boot[ik,,] = t(dum$sim)
    rm(ik); rm(dum)
  }
  med_polish_resi_forecast = med_polish_resi_basis %*% med_polish_resi_score_forecast
  
  # add the fixed parts
  
  Fixed=t(fixed_com)[,1:fh]
  med_polish_curve_forecast = med_polish_resi_forecast + Fixed
  
  return(list(med_polish_curve_forecast=med_polish_curve_forecast, 
              med_polish_resi_forecast=  med_polish_resi_forecast))
  
}


################################################################################
# With rolling window approach
################################################################################
library(foreach)
library(doParallel)
max_h = 10
ForecastC <- function(i,fixed_com,Residuals_f){
  pref=(n_year * i - (n_year - 1)):(n_year * i)
  forecasted_prefecture <- matrix(0, nrow = n_age * 2, ncol = max_h)
  for (k in 1:max_h) {
    pref_k <- pref[k:(n_year - max_h - 1 + k)]
    frc <- Pref_forecast_curves(fixed_com = fixed_com[pref_k, ],
                                  Residuals_f = Residuals_f[pref_k, ],
                                  est_method = "lrc",
                                  fh = 1, PI = NULL, B = 1000)
    forecasted_prefecture[ , k] = frc$med_polish_curve_forecast
  }
  return(forecasted_prefecture)
}

################################################################################
# With an expanding window approach
################################################################################
library(foreach)
library(doParallel)
max_h = 10
ForecastC_expanding <- function(i,fixed_com,Residuals_f){
  pref=(n_year * i - (n_year - 1)):(n_year * i)
  forecasted_prefecture <- matrix(0, nrow = n_age * 2, ncol = max_h)
  for (k in 1:max_h) {
    pref_k<-pref[1:(n_year-max_h-1+k)]
    frc=Pref_forecast_curves(fixed_com=fixed_com[pref_k,],
                               Residuals_f=Residuals_f[pref_k,],
                               est_method = "lrc",
                               fh = k, PI = NULL, B = 1000)
    forecasted_prefecture[,k]=frc$med_polish_curve_forecast[,1]
  }
  return(forecasted_prefecture)
}




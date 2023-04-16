#General function for point forecasts based on the FM-ANOVA decomposition

#1. General function for point forecasts obtained  to the functional residuals
# after removing deterministic components from the FM-ANOVA approach

#2. Computation of the point forecasts based on the rolling window approach, 
# using the function Pref_forecast_curves_mean.

#3. Computation of the point forecasts based on the expanding window approach, 
# using the function Pref_forecast_curves_mean.


###############################################################################
# Computation of the point forecast based on the FM-ANOVA decomposition
################################################################################
library(ftsa)
library(demography)
source("forecast_Arima.R")
Pref_forecast_curves_mean<-function(fixed_com,Residuals_f,
                                      est_method = c("lrc", "cov"),
                                      fh = 30, PI = NULL, B = 1000){
  mean_resi=t(Residuals_f)
  if(est_method == "lrc"){
    # estimate long-run covariance by kernel sandwich estimator
    mean_resi_lrc = long_run_covariance_estimation(mean_resi)
  }else if(est_method == "cov"){
    # estimate empirical covariance function
    mean_resi_lrc = cov(t(mean_resi))
  }
  # perform eigen-decomposition
  mean_resi_eigen = eigen(mean_resi_lrc)
  # determine retained number of components via eigenvalue ratio
  ret_component = vector("numeric", length(mean_resi_eigen$values) - 1)
  for(ik in 1:(length(mean_resi_eigen$values) - 1)){
    ret_component[ik] = mean_resi_eigen$values[ik+1]/mean_resi_eigen$values[ik]
  }
  retain_component = which.min(ret_component)
  # determine 1st set of basis function and its scores
  mean_resi_basis = as.matrix(mean_resi_eigen$vectors[,1:retain_component])
  mean_resi_score = crossprod(mean_resi, mean_resi_basis)
  # obtain forecasts of PC scores via auto.arima
  mean_resi_score_forecast = matrix(NA, retain_component, fh)
  mean_resi_score_forecast_boot = array(NA, dim = c(retain_component, fh, B))
  for(ik in 1:retain_component){
    dum = forecast_Arima(auto.arima(mean_resi_score[,ik]), h = fh, bootstrap = TRUE, npaths = B)
    mean_resi_score_forecast[ik,] = dum$mean
    mean_resi_score_forecast_boot[ik,,] = t(dum$sim)
    rm(ik); rm(dum)
  }
  mean_resi_forecast = mean_resi_basis %*% mean_resi_score_forecast
  
  # add the fixed parts
  
  Fixed=t(fixed_com)[,1:fh]
  mean_curve_forecast = mean_resi_forecast + Fixed
  
  return(list(mean_curve_forecast=mean_curve_forecast, 
              mean_resi_forecast=  mean_resi_forecast))
  
}


################################################################################
# Forecast with a  rolling window approach 
################################################################################
library(foreach)
library(doParallel)
max_h = 10
ForecastC_mean <- function(i,fixed_com,Residuals_f){
  pref=(n_year*i-(n_year-1)):(n_year*i)
  forecasted_prefecture_m<-matrix(0,nrow = n_age*2,ncol=max_h)
  for (k in 1:max_h) {
    pref_k<-pref[k:(n_year-max_h-1+k)]
    frc=Pref_forecast_curves_mean(fixed_com=fixed_com[pref_k,],
                                    Residuals_f=Residuals_f[pref_k,],
                                    est_method = "lrc",
                                    fh = 1, PI = NULL, B = 1000)
    forecasted_prefecture_m[,k]=frc$mean_curve_forecast
  }
  return(forecasted_prefecture_m)
}

################################################################################
# Forecast with an  expanding window approach 
################################################################################
library(foreach)
library(doParallel)
max_h = 10
ForecastC_mean_expanding <- function(i,fixed_com,Residuals_f){
  pref=(n_year*i-(n_year-1)):(n_year*i)
  forecasted_prefecture_m<-matrix(0,nrow = n_age*2,ncol=max_h )
  for (k in 1:max_h) {
    pref_k<-pref[1:(n_year-max_h-1+k)]
    frc=Pref_forecast_curves_mean(fixed_com=fixed_com[pref_k,],
                                    Residuals_f=Residuals_f[pref_k,],
                                    est_method = "lrc",
                                    fh = k, PI = NULL, B = 1000)
    forecasted_prefecture_m[,k]=frc$mean_curve_forecast[,1]
  }
  return(forecasted_prefecture_m)
}

# Script for computation of the coverage probabilities in prediction intervals. 


####################################################################
# Load the source file with the sieve bootstrap code.
library(doMC)
source("sieve_code.R")
####################################################################


#############################################
# Function for computing the interval score
#############################################

# l: lower bound
# u: upper bound
# x: actual holdout data
# alpha: level of significance alpha = 0.2

interval_score <- function(holdout, lb, ub, alpha)
{
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score  = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  return(score)
}

############################################################################################
# Function for the computation of the pointwise coverage probability and CPD
############################################################################################
# sim_data:   denotes the functional residuals for each state and gender. A matrix of dimension N_year X n_age
# fixed_comp: denotes the deterministic components for each state and gender. A matrix of dimension N_year X n_age
# sample_number: denotes the number of curves in each group. In this case it referes to n_year
# no_boot: the number of bootstrap samples in the sieve bootstrap procedure.
# no_core: number of determined cores to carried out the computation.
# method_pred: It refers to the prediction method to be used in the sieve bootstrap procedure. Bosq represents FAR(1) prediction method to forecast curves, and Bosq is the default
# K_val: NULL by default. If K refers to the fixed number of PC.
# selection_ncomp_porder: CPV_AICC is the default choice
# percent_CPV: 0.85 (by default), for robustness check, we also consider 0.8 and 0.9

coverage_comp <- function(sim_data,fixed_comp, sample_number=NULL, no_boot = 100,
                          no_core, method_pred = "Bosq", K_val = NULL,
                          selection_ncomp_porder="CPV_ACC",
                          prediction_method = c("sieve_bootstrap", "FAR_naive_bootstrap", "FAR_block_bootstrap")
                          ,percent_CPV= 0.90)
{
  
  grid_point = seq(0, 1, length.out = nrow(sim_data))
  colnames(sim_data) = 1:sample_number
  rownames(sim_data) = grid_point
  
  # define testing and training samples
  n_val = ncol(sim_data)
  n_testing = round(n_val/5)
  n_training_ini = round(n_val*4/5)
  if((n_training_ini + n_testing) != n_val)
  {
    warning("length of training sample + testing sample != total sample")
  }
  
  parametric_pointwise_PI_80_lb = parametric_pointwise_PI_80_ub = parametric_pointwise_PI_95_lb = parametric_pointwise_PI_95_ub = matrix(NA, nrow(sim_data), n_testing)
  if(prediction_method == "sieve_bootstrap")
  {
    if(is.null(K_val))
    {
      registerDoMC(no_core)
      sieve_boot = foreach(iwk = 1:n_testing) %dopar% sieve_bootstrap(fun_dat = t(sim_data[,1:(n_training_ini - 1 + iwk)]),
                                                                      k_observed_fun = nrow(t(sim_data[,1:(n_training_ini - 1 + iwk)])) - 1,
                                                                      grid = grid_point, ncomp_porder_selection = selection_ncomp_porder,
                                                                      CPV_percent = percent_CPV, pred_method = method_pred, B = no_boot
                                                                      ,PI_alpha = c(0.8, 0.95), semi_metric = "deriv", 
                                                                      q_order = 2,VAR_type = "none")
    }          
    else
    {
      registerDoMC(no_core)
      sieve_boot = foreach(iwk = 1:n_testing) %dopar% sieve_bootstrap(fun_dat = t(sim_data[,1:(n_training_ini - 1 + iwk)]),
                                                                      k_observed_fun = nrow(t(sim_data[,1:(n_training_ini - 1 + iwk)])) - 1,
                                                                      grid = grid_point, ncomp_porder_selection = selection_ncomp_porder,
                                                                      CPV_percent = percent_CPV, pred_method = method_pred, B = no_boot
                                                                      ,PI_alpha = c(0.8, 0.95), semi_metric = "deriv", 
                                                                      q_order = 2,VAR_type = "none")
    }
    for(ik in 1:n_testing)
    {
      # pointwise
      
      parametric_pointwise_PI_80_lb[,ik] = sieve_boot[[ik]]$parametric_pointwise_out_fore_80[,1]+fixed_comp[,1]
      parametric_pointwise_PI_80_ub[,ik] = sieve_boot[[ik]]$parametric_pointwise_out_fore_80[,2]+fixed_comp[,2]
      
      parametric_pointwise_PI_95_lb[,ik] = sieve_boot[[ik]]$parametric_pointwise_out_fore_95[,1]+fixed_comp[,1]
      parametric_pointwise_PI_95_ub[,ik] = sieve_boot[[ik]]$parametric_pointwise_out_fore_95[,2]+fixed_comp[,2]
      
      nonparametric_pointwise_PI_80_lb[,ik] = sieve_boot[[ik]]$nonparametric_pointwise_out_fore_80[,1]+fixed_comp[,1]
      nonparametric_pointwise_PI_80_ub[,ik] = sieve_boot[[ik]]$nonparametric_pointwise_out_fore_80[,2]+fixed_comp[,2]
      
      nonparametric_pointwise_PI_95_lb[,ik] = sieve_boot[[ik]]$nonparametric_pointwise_out_fore_95[,1]+fixed_comp[,1]
      nonparametric_pointwise_PI_95_ub[,ik] = sieve_boot[[ik]]$nonparametric_pointwise_out_fore_95[,2]+fixed_comp[,2]
      
      print(paste("sieve_", ik, sep="")); rm(ik)
    }
    rm(sieve_boot)
  }
  if(prediction_method == "FAR_naive_bootstrap"|prediction_method == "FAR_block_bootstrap")
  {
    registerDoMC(no_core)
    # far_boot = foreach(iwk = 1:n_testing) %dopar% far_boot(sim_dat = sim_data[,1:(n_training_ini - 1 + iwk)], bootrep = no_boot)
    if(prediction_method == "FAR_naive_bootstrap")
    {
      far_boot = foreach(iwk = 1:n_testing) %dopar% naive_bootstrap(simu_data = sim_data[,1:(n_training_ini - 1 + iwk)], bootrep = no_boot)
    }
    if(prediction_method == "FAR_block_bootstrap")
    {
      far_boot = foreach(iwk = 1:n_testing) %dopar% block_bootstrap(simu_data = sim_data[,1:(n_training_ini - 1 + iwk)], bootrep = no_boot)
    }
    for(ik in 1:n_testing)
    {
      # pointwise
      
      parametric_pointwise_PI_80_lb[,ik] = apply(far_boot[[ik]], 1, quantile, 0.1)
      parametric_pointwise_PI_80_ub[,ik] = apply(far_boot[[ik]], 1, quantile, 0.9)
      
      parametric_pointwise_PI_95_lb[,ik] = apply(far_boot[[ik]], 1, quantile, 0.025)
      parametric_pointwise_PI_95_ub[,ik] = apply(far_boot[[ik]], 1, quantile, 0.975)
      
      # uniform
      
      far_boot_sd = apply(far_boot[[ik]], 1, sd)
      far_boot_mean = apply(far_boot[[ik]], 1, mean)
      far_boot_scale = t(scale(t(far_boot[[ik]]), center = TRUE, scale = TRUE))
      
      far_boot_max = vector("numeric", no_boot)
      for(iw in 1:no_boot)
      {
        far_boot_max[iw] = max(abs((far_boot_scale)[,iw]))
      }
      far_boot_max_critical_0.8 = quantile(far_boot_max, 0.8)
      far_boot_max_critical_0.95 = quantile(far_boot_max, 0.95)
      
      rm(far_boot_sd); rm(far_boot_mean); rm(far_boot_scale); rm(far_boot_max); rm(far_boot_max_critical_0.8); rm(far_boot_max_critical_0.95)
      print(paste("far_", ik, sep="")); rm(ik)
    }
    rm(far_boot)
  }
  
  # compute pointwise coverage probability
  
  test_data = as.matrix(sim_data[,(n_training_ini + 1):n_val])+as.matrix(fixed_comp[,(n_training_ini + 1):n_val])
  
  parametric_lb_ind_80 = which(parametric_pointwise_PI_80_lb >= test_data)
  parametric_ub_ind_80 = which(parametric_pointwise_PI_80_ub <= test_data)
  parametric_pointwise_coverage_80 = 1 - length(union(parametric_lb_ind_80, parametric_ub_ind_80))/length(test_data)
  
  parametric_lb_ind_95 = which(parametric_pointwise_PI_95_lb >= test_data)
  parametric_ub_ind_95 = which(parametric_pointwise_PI_95_ub <= test_data)
  parametric_pointwise_coverage_95 = 1 - length(union(parametric_lb_ind_95, parametric_ub_ind_95))/length(test_data)
  
  # compute pointwise interval scores
  
  parametric_pointwise_score_80 = mean(interval_score(holdout = test_data, lb = parametric_pointwise_PI_80_lb,
                                                      ub = parametric_pointwise_PI_80_ub, alpha = 0.2))
  
  parametric_pointwise_score_95 = mean(interval_score(holdout = test_data, lb = parametric_pointwise_PI_95_lb,
                                                      ub = parametric_pointwise_PI_95_ub, alpha = 0.05))
  
 
  rm(parametric_lb_ind_80); rm(parametric_ub_ind_80); rm(parametric_lb_ind_95); rm(parametric_ub_ind_95)
  
  
  CPD_pointwise_0.80=abs(parametric_pointwise_coverage_80-0.80)
  CPD_pointwise_0.95=abs(parametric_pointwise_coverage_95-0.95)
  result = matrix(c(parametric_pointwise_coverage_80, parametric_pointwise_coverage_95,
                    parametric_pointwise_score_80, parametric_pointwise_score_95,CPD_pointwise_0.80
                    , CPD_pointwise_0.95), nrow = 1, byrow = TRUE)
  colnames(result) = c("Pointwise 80% coverage", "Pointwise 95% coverage",
                       "80% interval score", "95% interval score",
                       "CPD Pointwise 80% coverage","CPD Pointwise 95% coverage")
  
  rm(test_data); rm(ik)
  return(result)
}

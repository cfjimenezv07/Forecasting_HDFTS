##################
# Sieve bootstrap
##################

source("nonparametric_fof_regression.R")
source("method.FPE.R")
# install.packages(c("sde", "sandwich", "vars", "ftsa", "doMC", "fda.usc", "far"))
require(sde)
require(ftsa)
require(vars)
require(sandwich)
require(doMC)
require(fda.usc)
require(far)

# fun_dat: data matrix (n by p)
# k_observed_fun: any number between 1 and n, where n is sample size; customarily k = 1 (if computational time is not that long)
# grid: grid points in a curve
# p_m_order: give a pre-specified number of components and VAR order
# ncomp_porder_selection: CPV_AICC is the default choice
# CPV_percent: 0.85 (by default), for robustness check, we also consider 0.8 and 0.9
# VAR_type: default is none
# pred_method: Bosq represents FAR(1) prediction method to forecast curves, and Bosq is the default
# PI_alpha: 80% and 95% prediction intervals by default
# semi_metric: only useful when pred_method = "FV"
# q_order: only useful when pred_method = "FV"

sieve_bootstrap <- function(fun_dat, k_observed_fun, grid, p_m_order, ncomp_porder_selection, CPV_percent,
                            VAR_type = "none", pred_method = c("Bosq", "FV", "farforecast"),
                            B, PI_alpha = c(0.8, 0.95), semi_metric = "deriv", q_order = 2)
{
  #########
  # Step 1
  #########
  
  n_row = nrow(fun_dat)
  n_col = ncol(fun_dat)
  
  # de-centering functional data
  
  center_fun_dat = scale(fun_dat, center = TRUE, scale = FALSE)
  rownames(center_fun_dat)<-1:n_row
  
  # computing the mean of functional data
  
  X_bar = colMeans(fun_dat)
  
  # determine p_order
  
  if(missing(grid))
  {
    if(is.null(colnames(fun_dat)))
    {
      grid = seq(1, n_col, length.out = n_col)
    }
    else
    {
      grid = seq(min(as.numeric(colnames(fun_dat))),
                 max(as.numeric(colnames(fun_dat))),
                 length.out = n_col)
    }
  }
  
  # determine m_component based on FPE of Aue et al. (2015)
  
  fun_dat_object = fts(grid, t(center_fun_dat))
  
  if(missing(p_m_order))
  {
    if(ncomp_porder_selection == "FPE")
    {
      ncomp_porder = method.FPE(object = fun_dat_object, D = 6, var_type = VAR_type, Pmax = 5)
    }
    if(ncomp_porder_selection == "CPV_AICC")
    {
      ncomp = head(which(cumsum(ftsm(y = fun_dat_object, order = min(nrow(center_fun_dat), ncol(center_fun_dat)))$varprop)>= CPV_percent), 1)
      ftsm_object_scores = ftsm(y = fun_dat_object, order = ncomp, mean = FALSE)$coeff
      
      VAR_AICC = vector("numeric", 9)
      if(ncomp == 1)
      {
        for(VAR_order in 1:9)
        {
          ftsm_object_resi = na.omit(residuals(ar(x = ftsm_object_scores, aic = FALSE, order.max = VAR_order)))
          sigma = det(crossprod(ftsm_object_resi)/n_row)
          VAR_AICC[VAR_order] = n_row * log(sigma) + n_row * (n_row * ncomp + VAR_order * (ncomp^2))/(n_row - ncomp * (VAR_order + 1) - 1)
          rm(sigma); rm(ftsm_object_resi)
        }
      }
      if(ncomp > 1)
      {
        for(VAR_order in 1:9)
        {
          ftsm_object_resi = residuals(VAR(y = ftsm_object_scores, p = VAR_order, type = "none"))
          sigma = det(crossprod(ftsm_object_resi)/n_row)
          VAR_AICC[VAR_order] = n_row * log(sigma) + n_row * (n_row * ncomp + VAR_order * (ncomp^2))/(n_row - ncomp * (VAR_order + 1) - 1)
          rm(sigma); rm(ftsm_object_resi)
        }
      }
      ncomp_porder = c(which.min(VAR_AICC[is.finite(VAR_AICC)][VAR_AICC[is.finite(VAR_AICC)]>0]), ncomp)
      rm(VAR_AICC); rm(ncomp); rm(ftsm_object_scores)
    }
  }
  else
  {
    ncomp_porder = p_m_order
  }
  
  # think about removing p_order = 0
  
  p_order = ncomp_porder[1]
  m_component = ncomp_porder[1]
  rm(ncomp_porder)
  
  # compute the matrix variance (X %*% t(X))
  
  variance_est = var(center_fun_dat)
  
  # compute eigen-decomposition
  
  eigen_decomp = eigen(variance_est)
  eigen_decomp_vector = as.matrix(eigen_decomp$vectors[,1:m_component])
  eigen_decomp_score  = center_fun_dat %*% eigen_decomp_vector
  colnames(eigen_decomp_score) = 1:ncol(eigen_decomp_score)
  
  # reconstruction of original curves
  
  reconstruction = eigen_decomp_score %*% t(eigen_decomp_vector[,1:m_component])
  reconstruction_mean = (sweep(reconstruction, 2, X_bar, FUN = "+"))
  reconstruction_err = fun_dat - reconstruction_mean
  
  # determine p_order
  
  rm(variance_est); rm(eigen_decomp); rm(reconstruction); rm(reconstruction_mean)
  
  #########
  # Step 2
  #########
  
  if(m_component == 1)
  {
    A_est = matrix(ar(eigen_decomp_score, aic = FALSE, order = p_order)$ar, p_order, 1)
    residual_e = matrix(ar(eigen_decomp_score, aic = FALSE, order = p_order)$resid[-(1:p_order)],,1)
  }
  if(m_component > 1)
  {
    VAR_forward = VAR(eigen_decomp_score, p = p_order, type = VAR_type)
    A_est = matrix(NA, p_order * m_component, m_component)
    for(iw in 1:m_component)
    {
      A_est[,iw] = VAR_forward$varresult[[iw]]$coef
    }
    residual_e = resid(VAR_forward)
    rownames(residual_e) = (p_order + 1):n_row
    rm(VAR_forward)
  }
  colnames(A_est) = colnames(eigen_decomp_score) = colnames(residual_e) = paste("score", 1:m_component, sep="_")
  A_est_list = list()
  for(j in 1:p_order)
  {
    A_est_list[[j]] = A_est[((j-1)*m_component+1):(j*m_component),]
  }
  
  # eigen_decomp_score[141,] %*% A_est_list[[2]] + eigen_decomp_score[142,] %*% A_est_list[[1]] + eigen_decomp_score[35,] %*% A_est_list[[1]](checked)
  
  # obtain centered residuals later for bootstrapping
  
  residual_e_centered = scale(residual_e, center = TRUE, scale = FALSE)
  residual_e_centered_boot = array(NA, dim = c(nrow(residual_e_centered) * 2, ncol(residual_e_centered), B))
  for(iw in 1:B)
  {
    residual_e_centered_boot[,,iw] = residual_e_centered[sample(1:nrow(residual_e_centered), size = nrow(residual_e_centered) * 2, replace = TRUE),]
  }
  
  if(p_order > k_observed_fun)
  {
    series_len = length(1:(n_row - p_order))
    ell = p_order - k_observed_fun
    
    series_predictor_forward_predict = matrix(NA, p_order, m_component)
    ik = 1
    while(ik <= p_order)
    {
      series_predictor_forward_predict[ik,] = eigen_decomp_score[(series_len+ik),]
      ik = ik + 1
    }
    colnames(series_predictor_forward_predict) = paste("score", 1:m_component, sep="_")
    forward_predict_seq = as.matrix(series_predictor_forward_predict[rev(1:p_order),])
    
    epsilon_star = matrix(0, ell, m_component)
    for(ij in 1:p_order)
    {
      epsilon_star[1,] = epsilon_star[1,] + forward_predict_seq[ij,] %*% A_est_list[[ij]]
    }
    
    ij = 1
    while(ij < ell)
    {
      if(m_component == 1)
      {
        epsilon_new = matrix(c(epsilon_star[rev(1:ij)], forward_predict_seq)[1:p_order], ncol = 1)
      }
      else
      {
        epsilon_new = rbind(epsilon_star[rev(1:ij),], forward_predict_seq)[1:p_order,]
      }
      epsilon_star[ij+1,] = epsilon_new[ij,] %*% A_est_list[[ij]]
      ij = ij + 1
    }
    colnames(epsilon_star) = paste("predicted_score", 1:m_component, sep="_")
    rownames(epsilon_star) = paste("predicted_observation", 1:ell, sep="_")
    
    residual_e_centered_boot_ell = array(NA, dim = c(nrow(epsilon_star), ncol(epsilon_star), B))
    for(iw in 1:B)
    {
      residual_e_centered_boot_ell[,,iw] = residual_e_centered[sample(1:nrow(residual_e_centered), size = nrow(epsilon_star), replace = TRUE),]
    }
    
    epsilon_star_boot = array(NA, dim = c(nrow(epsilon_star), ncol(epsilon_star), B))
    for(ij in 1:B)
    {
      epsilon_star_boot[,,ij] = residual_e_centered_boot_ell[,,ij] + epsilon_star
    }
    rm(series_len)
  }
  
  #########
  # Step 3
  #########
  
  if(m_component == 1)
  {
    B_est = matrix(ar(eigen_decomp_score[rev(1:n_row),], aic = FALSE, order = p_order)$ar, p_order, 1)
  }
  if(m_component > 1)
  {
    VAR_backward = VAR(y = eigen_decomp_score[rev(1:n_row),], p = p_order, type = VAR_type)
    B_est = matrix(NA, p_order * m_component, m_component)
    for(iw in 1:m_component)
    {
      B_est[,iw] = VAR_backward$varresult[[iw]]$coef
    }
  }
  colnames(B_est) = paste("score", 1:m_component, sep = "_")
  B_est_list = list()
  for(j in 1:p_order)
  {
    B_est_list[[j]] = B_est[((j-1)*m_component+1):(j*m_component),]
  }
  
  # eigen_decomp_score[143,] %*% B_est_list[[2]] + eigen_decomp_score[142,] %*% B_est_list[[1]] (checked)
  
  #########
  # Step 4
  #########
  
  if(p_order <= k_observed_fun)
  {
    boot_score = array(NA, dim = c(n_row, m_component, B))
    for(ik in 1:B)
    {
      boot_score[n_row:(n_row - k_observed_fun + 1),,ik] = as.matrix(eigen_decomp_score[n_row:(n_row - k_observed_fun+1),])
    }
    n_row_new = nrow(boot_score)
  }
  if(p_order > k_observed_fun)
  {
    boot_score = array(NA, dim = c(n_row + (p_order - k_observed_fun), m_component, B))
    boot_score[(n_row+1):(n_row+(p_order - k_observed_fun)),,] = epsilon_star_boot
    for(ik in 1:B)
    {
      boot_score[n_row:(n_row - k_observed_fun + 1),,ik] = as.matrix(eigen_decomp_score[n_row:(n_row - k_observed_fun + 1),])
    }
    n_row_new = nrow(boot_score)
  }
  
  for(ik in 1:B)
  {
    for(ij in 1:(n_row - k_observed_fun))
    {
      mean_term_dum = vector("numeric", m_component)
      for(j in 1:p_order)
      {
        # mean_term_dum = mean_term_dum + boot_score[(n_row - k_observed_fun - ij + j),,ik] %*% B_est_list[[j]]
        mean_term_dum = mean_term_dum + boot_score[(n_row - k_observed_fun + j),,ik] %*% B_est_list[[j]]
      }
      mean_term = mean_term_dum
      
      truncation_order_est <- function(truncation_order)
      {
        psi = list()
        for(j in 1:p_order)
        {
          if(m_component == 1)
          {
            dum_mat_sum = 0
            for(u in 1:j)
            {
              dum_mat = try(psi[[j-u]] %*% A_est_list[[u]], silent = TRUE)
              if(any(class(dum_mat) == "try-error"))
              {
                dum_mat = diag(m_component) %*% A_est_list[[u]]
              }
              dum_mat_sum = dum_mat_sum + dum_mat
            }
            psi[[j]] = dum_mat_sum
            rm(dum_mat)
          }
          else
          {
            dum_mat_sum = matrix(0, nrow(A_est_list[[1]]), ncol(A_est_list[[1]]))
            for(u in 1:j)
            {
              dum_mat = try(psi[[j-u]] %*% A_est_list[[u]], silent = TRUE)
              if(any(class(dum_mat) == "try-error"))
              {
                dum_mat = diag(m_component) %*% A_est_list[[u]]
              }
              dum_mat_sum = dum_mat_sum + dum_mat
            }
            psi[[j]] = dum_mat_sum
            rm(dum_mat)
          }
        }
        
        if(truncation_order < (p_order + 1))
        {
          warning("Please supply a larger truncation order.")
        }
        if(m_component == 1)
        {
          for(j in (p_order+1):truncation_order)
          {
            dum_mat_sum = 0
            for(u in 1:p_order)
            {
              dum_mat = psi[[j-u]] %*% A_est_list[[u]]
              dum_mat_sum = dum_mat_sum + dum_mat
            }
            psi[[j]] = dum_mat_sum
            rm(dum_mat)
          }
        }
        else
        {
          for(j in (p_order+1):truncation_order)
          {
            dum_mat_sum = matrix(0, nrow(A_est_list[[1]]), ncol(A_est_list[[1]]))
            for(u in 1:p_order)
            {
              dum_mat = psi[[j-u]] %*% A_est_list[[u]]
              dum_mat_sum = dum_mat_sum + dum_mat
            }
            psi[[j]] = dum_mat_sum
            rm(dum_mat)
          }
        }
        backward_mat_recon = diag(m_component) %*% var(residual_e_centered) %*% diag(m_component)
        for(iw in 1:truncation_order)
        {
          backward_mat_recon = backward_mat_recon + t(psi[[iw]]) %*% var(residual_e_centered) %*% psi[[iw]]
        }
        return(list(psi = psi, err = norm(backward_mat_recon - var(eigen_decomp_score), type = "F")))
      }
      truncation_err = vector("numeric", (19 - p_order))
      for(iwj in (p_order + 2):20)
      {
        truncation_err[iwj - (p_order + 1)] = truncation_order_est(truncation_order = iwj)$err
      }
      s_value_1 = 20
      if(any(abs(diff(truncation_err)) <= 1*(10^-5)))
      {
        s_value_1 = which(abs(diff(truncation_err)) <= 1*(10^-5))[1]
      }
      else if(any(abs(diff(truncation_err)/truncation_err[2:length(truncation_err)]) <= 1*(10^-3)))
      {
        s_value_1 = which.min(abs(diff(truncation_err)/truncation_err[2:length(truncation_err)])) + 1
      }
      s_value_2 = which.min(truncation_err) + (p_order + 1)
      s_value = min(s_value_1 + (p_order + 1), s_value_2)
      psi = truncation_order_est(truncation_order = s_value)$psi
      
      gamma_u = list()
      for(u in (-p_order + 10):(-1 + 10))
      {
        if(m_component == 1)
        {
          dum_mat_sum = 0
          for(j in 0:(p_order + (u - 10)))
          {
            dum_mat = try(B_est_list[[p_order - j]] %*% psi[[p_order + u - j - 10]], silent = TRUE)
            if(any(class(dum_mat) == "try-error"))
            {
              dum_mat = B_est_list[[p_order - j]] %*% diag(m_component)
            }
            dum_mat_sum = dum_mat_sum + dum_mat
          }
          gamma_u[[u]] = -1 * dum_mat_sum
          rm(dum_mat)
        }
        else
        {
          dum_mat_sum = matrix(0, nrow(A_est_list[[1]]), ncol(A_est_list[[1]]))
          for(j in 0:(p_order + (u - 10)))
          {
            dum_mat = try(B_est_list[[p_order - j]] %*% psi[[p_order + u - j - 10]], silent = TRUE)
            if(any(class(dum_mat) == "try-error"))
            {
              dum_mat = B_est_list[[p_order - j]] %*% diag(m_component)
            }
            dum_mat_sum = dum_mat_sum + dum_mat
          }
          gamma_u[[u]] = -1 * dum_mat_sum
          rm(dum_mat)
        }
      }
      
      if(m_component == 1)
      {
        dum_mat_sum = 0
        for(j in 1:p_order)
        {
          dum_mat = B_est_list[[j]] %*% psi[[j]]
          dum_mat_sum = dum_mat_sum + dum_mat
        }
        gamma_u[[10]] = diag(m_component) - dum_mat_sum
      }
      else
      {
        dum_mat_sum = matrix(0, nrow(A_est_list[[1]]), ncol(A_est_list[[1]]))
        for(j in 1:p_order)
        {
          dum_mat = B_est_list[[j]] %*% psi[[j]]
          dum_mat_sum = dum_mat_sum + dum_mat
        }
        gamma_u[[10]] = diag(m_component) - dum_mat_sum
      }
      
      if(m_component == 1)
      {
        for(u in 1:(s_value - p_order))
        {
          dum_mat_sum = 0
          for(j in 1:p_order)
          {
            dum_mat = B_est_list[[j]] %*% psi[[u+j]]
            dum_mat_sum = dum_mat_sum + dum_mat
          }
          gamma_u[[10+u]] = psi[[u]] - dum_mat_sum
        }
      }
      else
      {
        for(u in 1:(s_value - p_order))
        {
          dum_mat_sum = matrix(0, nrow(A_est_list[[1]]), ncol(A_est_list[[1]]))
          for(j in 1:p_order)
          {
            dum_mat = B_est_list[[j]] %*% psi[[u+j]]
            dum_mat_sum = dum_mat_sum + dum_mat
          }
          gamma_u[[10+u]] = psi[[u]] - dum_mat_sum
        }
      }
      gamma_u_final = list()
      if(all(!sapply(gamma_u, is.null)))
      {
        gamma_u_final = gamma_u
      }
      else
      {
        gamma_u_final = gamma_u[-which(sapply(gamma_u, is.null))]
      }
      
      if(p_order <= k_observed_fun)
      {
        # n_last = n_row_new - p_order + ij
        n_last = n_row_new + p_order - ij
      }
      else
      {
        if(p_order <= m_component)
        {
          n_last = max(n_row_new - p_order - ij, 0)
        }
        else
        {
          n_last = n_row_new - ij - p_order - (p_order - m_component)
        }
      }
      dum_mat_sum = vector("numeric", m_component)
      if(ij == 1)
      {
        if((n_last - 1) != nrow(residual_e))
        {
          warnings("n_last does not equal to the nrow(residual_e)")
          for(iwk in 1:(s_value + 1))
          {
            dum_mat_sum = dum_mat_sum + residual_e_centered_boot[(nrow(residual_e) * 2 + 2 - ij) - iwk,,ik] %*% gamma_u_final[[iwk]]
          }
        }
        else
        {
          for(iwk in 1:(s_value + 1))
          {
            dum_mat_sum = dum_mat_sum + residual_e_centered_boot[nrow(residual_e) + n_last - iwk,,ik] %*% gamma_u_final[[iwk]]
          }
        }
      }
      else
      {
        if(n_last <= nrow(residual_e))
        {
          for(iwk in 1:(s_value + 1))
          {
            dum_mat_sum = dum_mat_sum + residual_e_centered_boot[nrow(residual_e) + n_last - iwk,,ik] %*% gamma_u_final[[iwk]]
          }
        }
        else
        {
          for(iwk in 1:(s_value + 1))
          {
            dum_mat_sum = dum_mat_sum + residual_e_centered_boot[nrow(residual_e) * 2 - iwk,,ik] %*% gamma_u_final[[iwk]]
          }
        }
      }
      boot_v = dum_mat_sum
      boot_score[(n_row - k_observed_fun - ij + 1),,ik] = mean_term + boot_v
      rm(mean_term); rm(boot_v); rm(n_last); rm(ij)
    }
    dup_index = duplicated(boot_score[,,ik])
    if(any(which(which(dup_index) > (s_value + 2))))
    {
      index_remove = which(dup_index)[(which(dup_index) > (s_value + 2))]
      n_index_remove = length(index_remove)
      if(m_component == 1)
      {
        boot_score[,,ik] = c(rep(boot_score[1,,ik], n_index_remove),
                             boot_score[-index_remove,,ik])
      }
      else
      {
        boot_score[,,ik] = rbind(matrix(rep(boot_score[1,,ik], n_index_remove),
                                        n_index_remove, m_component, byrow = TRUE),
                                 boot_score[-index_remove,,ik])
      }
      rm(index_remove); rm(n_index_remove)
    }
    print(ik); rm(ik); rm(dup_index)
  }
  
  #########
  # Step 5
  #########
  
  # bootstrap samples (used to produce \widehat{X}_star and \widehat{X})
  
  U_centered = scale(reconstruction_err, center = TRUE, scale = FALSE)[1:(n_row - k_observed_fun),]
  X_boot = array(NA, dim = c(n_col, n_row, B))
  
  # check boot_score
  
  if((n_row - k_observed_fun) == 1)
  {
    for(ik in 1:B)
    {
      X_boot[,1:(n_row - k_observed_fun),ik] = eigen_decomp_vector %*% as.matrix(boot_score[1:(n_row - k_observed_fun),,ik]) + U_centered + X_bar
      X_boot[,((n_row - k_observed_fun + 1): n_row),ik] = t(fun_dat[(n_row - k_observed_fun + 1): n_row,])
    }
  }
  else
  {
    for(ik in 1:B)
    {
      X_boot[,1:(n_row - k_observed_fun),ik] = eigen_decomp_vector %*% t(boot_score[1:(n_row - k_observed_fun),,ik]) + t(U_centered[sample(1:(n_row - k_observed_fun), replace = TRUE),]) + X_bar
      X_boot[,((n_row - k_observed_fun + 1): n_row),ik] = t(fun_dat[(n_row - k_observed_fun + 1): n_row,])
    }
  }
  
  # bootstrap one-step-ahead prediction intervals
  
  e_boot_val = as.matrix(residual_e_centered[sample(1:nrow(residual_e_centered), size = B, replace = TRUE),])
  U_boot_err = as.matrix(reconstruction_err[sample(1:n_row, size = B, replace = TRUE),])
  predict_boot = matrix(NA, n_col, B)
  for(iw in 1:B)
  {
    mean_term = 0
    for(ij in 1:p_order)
    {
      mean_term = mean_term + eigen_decomp_score[(n_row + 1 - ij),] %*% A_est_list[[ij]]
    }
    predict_boot[,iw] = X_bar + eigen_decomp_vector %*% t(mean_term + e_boot_val[iw,]) + U_boot_err[iw,]
    rm(mean_term); rm(iw)
  }
  colnames(predict_boot) = 1:B
  
  # produce one-step-ahead prediction based on bootstrap samples
  
  X_boot_fore = matrix(NA, n_col, B)
  for(iw in 1:B)
  {
    X_boot_mat = X_boot[,,iw]
    colnames(X_boot_mat) = rownames(fun_dat)
    rownames(X_boot_mat) = colnames(fun_dat)
    # X_boot_mat_centered = sweep(X_boot_mat, 1, X_bar, FUN = "-")
    X_boot_mat_centered = t(scale(t(X_boot_mat), center = TRUE, scale = FALSE))
    if(pred_method == "Bosq")
    {
      # implementing FAR(1)
      
      X_boot_far = as.fdata(object = X_boot_mat_centered, name = "X")
      Last_observed = as.fdata(object = X_boot_mat_centered[,ncol(X_boot_mat)], name = "X")
      
      model1.cv = far.cv(data = X_boot_far, y = "X", x = NULL, ncv = 10,
                         center = FALSE, na.rm = FALSE, joined = FALSE)
      k1 = model1.cv$minL2[1]
      
      model1 = far(data = X_boot_far, y = "X", kn = k1, center = FALSE, na.rm = FALSE)
      X_boot_fore[,iw] = X_bar + predict.far(model1, newdata = Last_observed, na.rm = FALSE)$X
      rm(model1); rm(model1.cv); rm(k1)
    }
    if(pred_method == "farforecast")
    {
      X_boot_fore[,iw] = farforecast(object = fts(1:n_col, X_boot_mat), h = 1, Dmax_value = 10, Pmax_value = 5)$point_fore$y
    }
    if(pred_method == "FV")
    {
      # implementing NFR
      
      X_boot_fore[,iw] = X_bar + as.numeric(ffunopare.knn.gcv(RESPONSES = t(X_boot_mat_centered)[2:ncol(X_boot_mat_centered),],
                                                              CURVES = t(X_boot_mat_centered)[1:(ncol(X_boot_mat_centered) - 1),],
                                                              PRED = t(X_boot_mat_centered)[ncol(X_boot_mat_centered),],
                                                              semimetric = semi_metric, q = q_order)$Predicted.values)
    }
    rm(X_boot_mat_centered); rm(X_boot_mat)
    print(iw); rm(iw)
  }
  colnames(X_boot_fore) = 1:B
  
  # compute prediction error
  
  prediction_err_boot_pre = prediction_err_boot = predict_boot - X_boot_fore
  
  # one-step-ahead prediction intervals based on historical observations (without standardization)
  
  if(pred_method == "Bosq")
  {
    out_fore_mat_mean = colMeans(fun_dat)
    out_fore_mat_centered = t(scale(fun_dat, center = TRUE, scale = FALSE))
    
    out_fore_far = as.fdata(object = out_fore_mat_centered, name = "X")
    out_fore_last_observed  = as.fdata(object = out_fore_mat_centered[,ncol(out_fore_mat_centered)], name = "X")
    
    model1.cv = far.cv(data = out_fore_far, y = "X", x = NULL, ncv = 10,
                       center = FALSE, na.rm = FALSE, joined = FALSE)
    k1 = model1.cv$minL2[1]
    
    model1 = far(data = out_fore_far, y = "X", kn = k1, center = FALSE, na.rm = FALSE)
    out_fore_actual = out_fore_mat_mean + predict.far(model1, newdata = out_fore_last_observed, na.rm = FALSE)$X
    rm(model1); rm(model1.cv); rm(k1)
  }
  if(pred_method == "farforecast")
  {
    out_fore_actual = farforecast(object = fts(1:n_col, t(fun_dat)), h = 1, Dmax_value = 10, Pmax_value = 5)$point_fore$y
  }
  if(pred_method == "FV")
  {
    out_fore_actual = X_bar + as.numeric(ffunopare.knn.gcv(RESPONSES = center_fun_dat[2:nrow(fun_dat),],
                                                           CURVES = center_fun_dat[1:(nrow(fun_dat) - 1),],
                                                           PRED = center_fun_dat[nrow(fun_dat),],
                                                           semimetric = semi_metric, q = q_order)$Predicted.values)
  }
  rm(center_fun_dat); rm(X_bar)
  
  # parametric approach
  
  prediction_err_sd = apply(prediction_err_boot, 1, sd)
  
  select_pointwise_tuning_parameter <- function(tuning_parameter, nominal_coverage)
  {
    nB = ncol(prediction_err_boot)
    lb = -1 * tuning_parameter * prediction_err_sd
    ub = tuning_parameter * prediction_err_sd
    out_number = vector("numeric", nB)
    for(ik in 1:nB)
    {
      holdout = prediction_err_boot[,ik]
      out_number[ik] = length(union(which(lb > holdout), which(ub < holdout)))
    }
    empirical_coverage = 1 - (sum(out_number)/length(prediction_err_boot))
    return(abs(empirical_coverage - nominal_coverage))
  }
  
  select_uniform_tuning_parameter <- function(tuning_parameter, nominal_coverage)
  {
    nB = ncol(prediction_err_boot)
    lb = -1 * tuning_parameter * prediction_err_sd
    ub = tuning_parameter * prediction_err_sd
    out_number = vector("numeric", nB)
    for(ik in 1:nB)
    {
      holdout = prediction_err_boot[,ik]
      out_number[ik] = ifelse(any(union(which(lb > holdout), which(ub < holdout))), 1, 0)
    }
    empirical_coverage = 1 - sum(out_number)/nB
    return(abs(empirical_coverage - nominal_coverage))
  }
  
  pointwise_tuning_parameter_PI_alpha_1 = optimize(f = select_pointwise_tuning_parameter, interval = c(0.1, 10), nominal_coverage = PI_alpha[1])$minimum
  pointwise_tuning_parameter_PI_alpha_2 = optimize(f = select_pointwise_tuning_parameter, interval = c(0.1, 10), nominal_coverage = PI_alpha[2])$minimum
  
  uniform_tuning_parameter_PI_alpha_1 = optimize(f = select_uniform_tuning_parameter, interval = c(0.1, 10), nominal_coverage = PI_alpha[1])$minimum
  uniform_tuning_parameter_PI_alpha_2 = optimize(f = select_uniform_tuning_parameter, interval = c(0.1, 10), nominal_coverage = PI_alpha[2])$minimum
  
  tuning_parameters = c(pointwise_tuning_parameter_PI_alpha_1, pointwise_tuning_parameter_PI_alpha_2,
                        uniform_tuning_parameter_PI_alpha_1, uniform_tuning_parameter_PI_alpha_2)
  
  parametric_pointwise_out_fore_80 = cbind(out_fore_actual - prediction_err_sd * pointwise_tuning_parameter_PI_alpha_1,
                                           out_fore_actual + prediction_err_sd * pointwise_tuning_parameter_PI_alpha_1)
  
  parametric_pointwise_out_fore_95 = cbind(out_fore_actual - prediction_err_sd * pointwise_tuning_parameter_PI_alpha_2,
                                           out_fore_actual + prediction_err_sd * pointwise_tuning_parameter_PI_alpha_2)
  
  parametric_uniform_out_fore_80 = cbind(out_fore_actual - prediction_err_sd * uniform_tuning_parameter_PI_alpha_1,
                                         out_fore_actual + prediction_err_sd * uniform_tuning_parameter_PI_alpha_1)
  
  parametric_uniform_out_fore_95 = cbind(out_fore_actual - prediction_err_sd * uniform_tuning_parameter_PI_alpha_2,
                                         out_fore_actual + prediction_err_sd * uniform_tuning_parameter_PI_alpha_2)
  
  # nonparametric bootstrap
  
  rm(prediction_err_boot)
  prediction_err_boot = t(scale(t(prediction_err_boot_pre), center = TRUE, scale = FALSE))
  
  # boot_fore_samples = apply(prediction_err_boot, 2, "+", out_fore_actual)
  nonparametric_pointwise_out_fore_80 = apply(t(apply(prediction_err_boot, 1, quantile, c((1 - PI_alpha[1])/2, (1 + PI_alpha[1])/2))), 2, "+", out_fore_actual)
  nonparametric_pointwise_out_fore_95 = apply(t(apply(prediction_err_boot, 1, quantile, c((1 - PI_alpha[2])/2, (1 + PI_alpha[2])/2))), 2, "+", out_fore_actual)
  
  # uniform prediction interval: one-step-ahead prediction intervals baased on historical observations (with standardization)
  
  pointwise_sd = apply(prediction_err_boot, 1, sd)
  standardized_prediction_err_boot = matrix(NA, n_col, B)
  for(iw in 1:B)
  {
    standardized_prediction_err_boot[,iw] = prediction_err_boot[,iw]/pointwise_sd
  }
  uniform_prediction_err = vector("numeric", B)
  for(iw in 1:B)
  {
    uniform_prediction_err[iw] = max(abs(standardized_prediction_err_boot[,iw]))
  }
  nonparametric_uniform_out_fore_80 = cbind(out_fore_actual - quantile(uniform_prediction_err, PI_alpha[1], na.rm = TRUE) * pointwise_sd,
                                            out_fore_actual + quantile(uniform_prediction_err, PI_alpha[1], na.rm = TRUE) * pointwise_sd)
  
  nonparametric_uniform_out_fore_95 = cbind(out_fore_actual - quantile(uniform_prediction_err, PI_alpha[2], na.rm = TRUE) * pointwise_sd,
                                            out_fore_actual + quantile(uniform_prediction_err, PI_alpha[2], na.rm = TRUE) * pointwise_sd)
  
  return(list(X_star = predict_boot, X_boot = X_boot, X_boot_fore = X_boot_fore,
              prediction_err_boot = prediction_err_boot,
              p_order = p_order, m_component = m_component,
              out_fore_actual = out_fore_actual,
              out_fore_boot = prediction_err_boot + matrix(rep(out_fore_actual, B), nrow(prediction_err_boot), B),
              prediction_err_sd = prediction_err_sd,
              tuning_parameters = tuning_parameters,
              parametric_pointwise_out_fore_80 = parametric_pointwise_out_fore_80,
              parametric_pointwise_out_fore_95 = parametric_pointwise_out_fore_95,
              parametric_uniform_out_fore_80 = parametric_uniform_out_fore_80,
              parametric_uniform_out_fore_95 = parametric_uniform_out_fore_95,
              nonparametric_pointwise_out_fore_80 = nonparametric_pointwise_out_fore_80,
              nonparametric_pointwise_out_fore_95 = nonparametric_pointwise_out_fore_95,
              nonparametric_uniform_out_fore_80 = nonparametric_uniform_out_fore_80,
              nonparametric_uniform_out_fore_95 = nonparametric_uniform_out_fore_95))
}

#############################################
# Function for computing the pointwise score
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

#######################################################
# DGP for simulating stationary functional time series
#######################################################

# FAR_psi_1

funKernel = function(ref, k)
{
  Mat = matrix(nrow=ref, ncol=ref)
  for(i in 1:ref)
  {
    for(j in 1:ref)
    {
      Mat[i,j] = k*exp(.5*((i/ref)^2+(j/ref)^2))
    }
  }
  return(Mat)
}

funIntegral = function(ref, Mat, X)
{
  Mat = Mat%*%X
  return(Mat/ref)
}

funARMat_b = function(N, refinement, b_val, para_val)
{
  Mat = funKernel(refinement, k = para_val)
  firstmat <- matrix(nrow = refinement, ncol = N+50)
  firstmat[,1] = BM(N = refinement - 1)
  firstmat[,2] = BM(N = refinement - 1)
  for(i in 3:(N+50))
  {
    firstmat[,i] = funIntegral(refinement, Mat, firstmat[,(i-1)]) + b_val * firstmat[,(i-2)] + BM(N = refinement-1)
  }
  finalMat = firstmat[1:refinement,51:(N+50)]
  return(finalMat)
}

funARMAMat_b <- function(N, refinement, b_val, para_val, c_val)
{
  Mat = funKernel(refinement, k = para_val)
  error = matrix(nrow = refinement, ncol = N + N)
  for(i in 1:(N + N)) error[,i] = BM(N = refinement - 1)
  
  firstmat <- matrix(nrow = refinement, ncol = N + N)
  firstmat[,1] = error[,1]
  firstmat[,2] = error[,2]
  for(i in 3:(N + N))
  {
    firstmat[,i] = funIntegral(refinement, Mat, firstmat[,(i-1)]) + b_val * firstmat[,(i-2)] + error[,i] + c_val * error[,(i-1)]
  }
  finalMat = firstmat[1:refinement,(N + 1):(N + N)]
  return(finalMat)
}

###########################
# one-step-ahead forecasts
###########################

coverage_sim <- function(DGP = c("FAR", "FARMA"), seed_number, sample_number, no_boot = 100,
                         no_core, method_pred, b_val_sim, norm_val, c_val_sim, K_val = NULL,
                         selection_ncomp_porder, percent_CPV,
                         prediction_method = c("sieve_bootstrap", "FAR_naive_bootstrap", "FAR_block_bootstrap"))
{
  DGP = match.arg(DGP)
  prediction_method = match.arg(prediction_method)
  
  # set random seed
  
  set.seed(123 + seed_number)
  
  # simulate data
  
  if(DGP == "FAR")
  {
    sim_data = funARMat_b(N = sample_number, refinement = 101, b_val = b_val_sim,
                          para_val = norm_val)
  }
  if(DGP == "FARMA")
  {
    sim_data = funARMAMat_b(N = sample_number, refinement = 101, b_val = b_val_sim,
                            para_val = norm_val, c_val = c_val_sim)
  }
  grid_point = seq(0, 1, length.out = nrow(sim_data))
  colnames(sim_data) = 1:sample_number
  rownames(sim_data) = grid_point
  
  # define testing and training samples
  n_val = ncol(sim_data)
  n_testing = trunc(n_val/5)
  n_training_ini = trunc(n_val*4/5)
  if((n_training_ini + n_testing) != n_val)
  {
    warning("length of training sample + testing sample != total sample")
  }
  
  parametric_pointwise_PI_80_lb = parametric_pointwise_PI_80_ub = parametric_pointwise_PI_95_lb = parametric_pointwise_PI_95_ub =
    parametric_uniform_PI_80_lb = parametric_uniform_PI_80_ub = parametric_uniform_PI_95_lb = parametric_uniform_PI_95_ub =
    nonparametric_pointwise_PI_80_lb = nonparametric_pointwise_PI_80_ub = nonparametric_pointwise_PI_95_lb = nonparametric_pointwise_PI_95_ub =
    nonparametric_uniform_PI_80_lb = nonparametric_uniform_PI_80_ub = nonparametric_uniform_PI_95_lb = nonparametric_uniform_PI_95_ub = matrix(NA, nrow(sim_data), n_testing)
  if(prediction_method == "sieve_bootstrap")
  {
    if(is.null(K_val))
    {
      registerDoMC(no_core)
      sieve_boot = foreach(iwk = 1:n_testing) %dopar% sieve_bootstrap(fun_dat = t(sim_data[,1:(n_training_ini - 1 + iwk)]),
                                                                      k_observed_fun = nrow(t(sim_data[,1:(n_training_ini - 1 + iwk)])) - 1,
                                                                      grid = grid_point, ncomp_porder_selection = selection_ncomp_porder,
                                                                      CPV_percent = percent_CPV, pred_method = method_pred, B = no_boot)
    }
    else
    {
      registerDoMC(no_core)
      sieve_boot = foreach(iwk = 1:n_testing) %dopar% sieve_bootstrap(fun_dat = t(sim_data[,1:(n_training_ini - 1 + iwk)]),
                                                                      k_observed_fun = K_val,
                                                                      grid = grid_point, ncomp_porder_selection = selection_ncomp_porder,
                                                                      CPV_percent = percent_CPV, pred_method = method_pred, B = no_boot)
    }
    for(ik in 1:n_testing)
    {
      # pointwise
      
      parametric_pointwise_PI_80_lb[,ik] = sieve_boot[[ik]]$parametric_pointwise_out_fore_80[,1]
      parametric_pointwise_PI_80_ub[,ik] = sieve_boot[[ik]]$parametric_pointwise_out_fore_80[,2]
      
      parametric_pointwise_PI_95_lb[,ik] = sieve_boot[[ik]]$parametric_pointwise_out_fore_95[,1]
      parametric_pointwise_PI_95_ub[,ik] = sieve_boot[[ik]]$parametric_pointwise_out_fore_95[,2]
      
      nonparametric_pointwise_PI_80_lb[,ik] = sieve_boot[[ik]]$nonparametric_pointwise_out_fore_80[,1]
      nonparametric_pointwise_PI_80_ub[,ik] = sieve_boot[[ik]]$nonparametric_pointwise_out_fore_80[,2]
      
      nonparametric_pointwise_PI_95_lb[,ik] = sieve_boot[[ik]]$nonparametric_pointwise_out_fore_95[,1]
      nonparametric_pointwise_PI_95_ub[,ik] = sieve_boot[[ik]]$nonparametric_pointwise_out_fore_95[,2]
      
      # uniform
      
      parametric_uniform_PI_80_lb[,ik] = sieve_boot[[ik]]$parametric_uniform_out_fore_80[,1]
      parametric_uniform_PI_80_ub[,ik] = sieve_boot[[ik]]$parametric_uniform_out_fore_80[,2]
      
      parametric_uniform_PI_95_lb[,ik] = sieve_boot[[ik]]$parametric_uniform_out_fore_95[,1]
      parametric_uniform_PI_95_ub[,ik] = sieve_boot[[ik]]$parametric_uniform_out_fore_95[,2]
      
      nonparametric_uniform_PI_80_lb[,ik] = sieve_boot[[ik]]$nonparametric_uniform_out_fore_80[,1]
      nonparametric_uniform_PI_80_ub[,ik] = sieve_boot[[ik]]$nonparametric_uniform_out_fore_80[,2]
      
      nonparametric_uniform_PI_95_lb[,ik] = sieve_boot[[ik]]$nonparametric_uniform_out_fore_95[,1]
      nonparametric_uniform_PI_95_ub[,ik] = sieve_boot[[ik]]$nonparametric_uniform_out_fore_95[,2]
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
      
      parametric_uniform_PI_80_lb[,ik] = far_boot_mean - far_boot_max_critical_0.8 * far_boot_sd
      parametric_uniform_PI_80_ub[,ik] = far_boot_mean + far_boot_max_critical_0.8 * far_boot_sd
      
      parametric_uniform_PI_95_lb[,ik] = far_boot_mean - far_boot_max_critical_0.95 * far_boot_sd
      parametric_uniform_PI_95_ub[,ik] = far_boot_mean + far_boot_max_critical_0.95 * far_boot_sd
      rm(far_boot_sd); rm(far_boot_mean); rm(far_boot_scale); rm(far_boot_max); rm(far_boot_max_critical_0.8); rm(far_boot_max_critical_0.95)
      print(paste("far_", ik, sep="")); rm(ik)
    }
    rm(far_boot)
  }
  
  # compute pointwise coverage probability
  
  test_data = as.matrix(sim_data[,(n_training_ini + 1):n_val])
  
  parametric_lb_ind_80 = which(parametric_pointwise_PI_80_lb >= test_data)
  parametric_ub_ind_80 = which(parametric_pointwise_PI_80_ub <= test_data)
  parametric_pointwise_coverage_80 = 1 - length(union(parametric_lb_ind_80, parametric_ub_ind_80))/length(test_data)
  
  parametric_lb_ind_95 = which(parametric_pointwise_PI_95_lb >= test_data)
  parametric_ub_ind_95 = which(parametric_pointwise_PI_95_ub <= test_data)
  parametric_pointwise_coverage_95 = 1 - length(union(parametric_lb_ind_95, parametric_ub_ind_95))/length(test_data)
  
  nonparametric_lb_ind_80 = which(nonparametric_pointwise_PI_80_lb >= test_data)
  nonparametric_ub_ind_80 = which(nonparametric_pointwise_PI_80_ub <= test_data)
  nonparametric_pointwise_coverage_80 = 1 - length(union(nonparametric_lb_ind_80, nonparametric_ub_ind_80))/length(test_data)
  
  nonparametric_lb_ind_95 = which(nonparametric_pointwise_PI_95_lb >= test_data)
  nonparametric_ub_ind_95 = which(nonparametric_pointwise_PI_95_ub <= test_data)
  nonparametric_pointwise_coverage_95 = 1 - length(union(nonparametric_lb_ind_95, nonparametric_ub_ind_95))/length(test_data)
  
  # compute pointwise interval scores
  
  parametric_pointwise_score_80 = mean(interval_score(holdout = test_data, lb = parametric_pointwise_PI_80_lb,
                                                      ub = parametric_pointwise_PI_80_ub, alpha = 0.2))
  
  parametric_pointwise_score_95 = mean(interval_score(holdout = test_data, lb = parametric_pointwise_PI_95_lb,
                                                      ub = parametric_pointwise_PI_95_ub, alpha = 0.05))
  
  nonparametric_pointwise_score_80 = mean(interval_score(holdout = test_data, lb = nonparametric_pointwise_PI_80_lb,
                                                         ub = nonparametric_pointwise_PI_80_ub, alpha = 0.2))
  
  nonparametric_pointwise_score_95 = mean(interval_score(holdout = test_data, lb = nonparametric_pointwise_PI_95_lb,
                                                         ub = nonparametric_pointwise_PI_95_ub, alpha = 0.05))
  rm(parametric_lb_ind_80); rm(parametric_ub_ind_80); rm(parametric_lb_ind_95); rm(parametric_ub_ind_95)
  rm(nonparametric_lb_ind_80); rm(nonparametric_ub_ind_80); rm(nonparametric_lb_ind_95); rm(nonparametric_ub_ind_95)
  
  # compute uniform coverage probability
  
  parametric_count_80 = parametric_count_95 = nonparametric_count_80 = nonparametric_count_95 = vector("numeric", n_testing)
  for(ik in 1:n_testing)
  {
    parametric_count_80[ik] = ifelse(any(c(any(which(parametric_uniform_PI_80_lb[,ik] >= test_data[,ik])),
                                           any(which(parametric_uniform_PI_80_ub[,ik] <= test_data[,ik])))), 1, 0)
    
    parametric_count_95[ik] = ifelse(any(c(any(which(parametric_uniform_PI_95_lb[,ik] >= test_data[,ik])),
                                           any(which(parametric_uniform_PI_95_ub[,ik] <= test_data[,ik])))), 1, 0)
    
    nonparametric_count_80[ik] = ifelse(any(c(any(which(nonparametric_uniform_PI_80_lb[,ik] >= test_data[,ik])),
                                              any(which(nonparametric_uniform_PI_80_ub[,ik] <= test_data[,ik])))), 1, 0)
    
    nonparametric_count_95[ik] = ifelse(any(c(any(which(nonparametric_uniform_PI_95_lb[,ik] >= test_data[,ik])),
                                              any(which(nonparametric_uniform_PI_95_ub[,ik] <= test_data[,ik])))), 1, 0)
  }
  parametric_uniform_coverage_80 = 1 - sum(parametric_count_80)/ncol(test_data)
  parametric_uniform_coverage_95 = 1 - sum(parametric_count_95)/ncol(test_data)
  
  nonparametric_uniform_coverage_80 = 1 - sum(nonparametric_count_80)/ncol(test_data)
  nonparametric_uniform_coverage_95 = 1 - sum(nonparametric_count_95)/ncol(test_data)
  rm(parametric_count_80); rm(parametric_count_95); rm(nonparametric_count_80); rm(nonparametric_count_95)
  
  result = matrix(c(parametric_pointwise_coverage_80, parametric_pointwise_coverage_95,
                    parametric_pointwise_score_80, parametric_pointwise_score_95,
                    parametric_uniform_coverage_80, parametric_uniform_coverage_95,
                    
                    nonparametric_pointwise_coverage_80, nonparametric_pointwise_coverage_95,
                    nonparametric_pointwise_score_80, nonparametric_pointwise_score_95,
                    nonparametric_uniform_coverage_80, nonparametric_uniform_coverage_95), nrow = 2, byrow = TRUE)
  colnames(result) = c("Pointwise 80% coverage", "Pointwise 95% coverage",
                       "80% interval score", "95% interval score",
                       "Uniform 80% coverage", "Uniform 95% coverage")
  
  rm(test_data); rm(ik)
  return(result)
}

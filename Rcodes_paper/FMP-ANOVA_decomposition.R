
packages <- c("generics", "demography", "forecast","fda","fdaoutlier", "rlist", "mrfDepth")
## Now load or install&load all
package_check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
###########################################################################
# Two way median polish 
###########################################################################
Two_way_median_polish <- function(Y, row_partition_index,column_partition_index) {
  ### 1. when Y is a matrix, it represents univariate functional observations 
  ## with the number of rows = n (sample size) and the number of cols = the grid length
  ### 2. partition_index for the rows. The indexes must contain the amount of replicates per row
  ### 3. partition_index for the columns. The indexes must contain the amount of replicates per columns
  ### It shows the possible categories of the effect and the length of partition_index >=2
  if (class(Y) != "matrix" && class(Y) != "array") {
    stop("Y should be functional data")
  } else if (length(class(Y)) == 2 && 
             class(Y)[1] == "matrix" && 
             class(Y)[2] == "array") {
    
    sample_size <- 2*nrow(Y)
    grid_size <- ncol(Y)/2
  } 
  
  if (length(row_partition_index) < 2|length(column_partition_index) < 2) {
    stop("partition_index is not a qualified partition")
  }
  
  Functional_grand_effect <- rep(0, length.out = grid_size)
  Functional_row_effect <- matrix(0, nrow = length(row_partition_index), ncol = grid_size)
  Functional_column_effect <- matrix(0, nrow = length(column_partition_index), ncol = grid_size)
  new_Functional_row_effect <- Functional_row_effect + 1
  Y_iter<-Y
  Iteration <- 1
  while (sum(new_Functional_row_effect) != 0) {
    cat("Iteration=", Iteration,"\n")
    cat("Sum =", sum(new_Functional_row_effect),"\n")
    
    partition_row_Y_f <- function(Y_iter){
      partition_Y_1 <- lapply(1:length(column_partition_index), function(k) {
        Y_iter[,column_partition_index[[k]] ]})
      partition_by_row=list()
      for (h in 1:length(partition_Y_1)) {
        YY=partition_Y_1[[h]]
        partition_by_row[[h]] <- lapply(1:length(row_partition_index), function(k) {
          YY[row_partition_index[[k]], ]})
      }
      Final_partition_by_row<-list()
      for (k in 1:length(partition_by_row[[1]])) {
        PBR_1=partition_by_row[[1]]
        PBR_2=partition_by_row[[2]]
        Final_partition_by_row[[k]]=rbind(PBR_1[[k]],PBR_2[[k]])
      }
      return(Final_partition_by_row)
    }
    partition_row_Y<-partition_row_Y_f(Y_iter)
    
    #Compute the median in each row with the curves in the row. Depth based selection
    get_median<-function(data){
      if(nrow(data)==2){
        med_curve=colMeans(data)
      }else{
        fb=functional_boxplot(data,depth_method = "mbd")
        med=fb$median_curve
        med_curve=data[med,]
      }
      return(med_curve)
    }
    Median_at_each_row<-lapply(partition_row_Y,get_median)
    
    
    #Get the functional median of the row functional medians
    Medians_medians_row<-get_median(list.rbind( Median_at_each_row))
    Functional_grand_effect<-Functional_grand_effect+Medians_medians_row
    
    #Get the functional row effect
    new_Functional_row_effect<-t(sapply(Median_at_each_row, function(k) {k -Medians_medians_row}))
    #sum(new_Functional_row_effect)
    Functional_row_effect<-Functional_row_effect+new_Functional_row_effect
    
    # Substract the row-median in each row from the rest of the curves in the row
    Substraction<-function(data){
      med<-get_median(data)
      new_data=matrix(0,nrow = dim(data)[1],ncol = dim(data)[2])
      for (i in 1:dim(data)[1]) {
        new_data[i,]=data[i,]-med
      }
      return(new_data)
    }
    partition_row_Y<-lapply(partition_row_Y,Substraction)
    
    #Match again Y_iter
    f<-function(partition_row_Y){
      #combine all dataset
      New_data_1<-list()
      New_data_2<-list()
      for (i in 1:length(partition_row_Y)) {
        ND<-partition_row_Y[[i]]
        New_data_1[[i]]=ND[1:n_year,]
        New_data_2[[i]]=ND[(n_year+1):(2*n_year),]
      }
      P_1<-list.rbind(New_data_1)
      P_2<-list.rbind(New_data_2)
      Y_matched<-cbind(P_1,P_2)
      return(Y_matched)
    }
    Y_iter<-f(partition_row_Y)
    
    # Now by columns
    partition_column_Y_f<-function(Y_iter){
      P_1=Y_iter[,1:n_age]
      P_2=Y_iter[,(n_age+1):(2*n_age)]
      Partition_col=list( P_1,P_2)
      return(Partition_col)
    } 
    partition_column_Y<-partition_column_Y_f(Y_iter)
    #get the median at each column
    Median_at_each_column<-lapply(partition_column_Y,get_median)
    
    #Get the functional median of the row functional medians
    Medians_medians_col<-get_median(list.rbind( Median_at_each_column))
    
    #Add the FGE
    Functional_grand_effect=Functional_grand_effect+Medians_medians_col
    
    #Get the functional column effect
    new_Functional_column_effect<-t(sapply(Median_at_each_column, function(k) {k -Medians_medians_col}))
    Functional_column_effect<-Functional_column_effect+new_Functional_column_effect
    
    #data
    partition_column_Y<-lapply(partition_column_Y,Substraction)
    
    Y_iter<-list.cbind(partition_column_Y)
    sum(new_Functional_row_effect)
    Iteration <- Iteration + 1
    
  }
  

  return(list(grand_effect = Functional_grand_effect,
              row_effect = Functional_row_effect,
              col_effect = Functional_column_effect))
}

###########################################################################
# Two way median polish for the time-varying components
###########################################################################
# 1. The following code uses the results from the two way functional median polish
# to obtain the functional residuals. 
#2. The code also proves that the deterministic and time-varying components sum up to the data.

Two_way_Residuals<-function(Y,Both,n_prefectures,n_year,n_age){
  # Y should be the matrix data that contains both male and female population stacked by columns
  # Both is the result after running the function Two_way_median_polish
  Two_FGE=Both$grand_effect
  Two_FRE=Both$row_effect
  Two_FCE=Both$col_effect
  all_male=Y[,1:n_age]
  all_female=Y[,(n_age+1):(2*n_age)]
  
  residuals_b1<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  residuals_b2<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  for (i in 1:nrow(residuals_b1)) {
    residuals_b1[i,]=all_male[i,]-Two_FGE
    residuals_b2[i,]=all_female[i,]-Two_FGE
  }
  
  #remove the column
  residuals_b1c<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  residuals_b2c<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  for (i in 1:nrow(residuals_b1c)) {
    residuals_b1c[i,]=residuals_b1[i,]-Two_FCE[1,]
    residuals_b2c[i,]=residuals_b2[i,]-Two_FCE[2,]
  }
  
  #Remove the row
  residuals_b1r<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  residuals_b2r<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  for (j in 1:n_prefectures) {
    residuals_b1r[(n_year*j-(n_year-1)):(n_year*j),]=residuals_b1c[(n_year*j-(n_year-1)):(n_year*j),]-t(replicate(n_year,Two_FRE[j,]))
    residuals_b2r[(n_year*j-(n_year-1)):(n_year*j),]=residuals_b2c[(n_year*j-(n_year-1)):(n_year*j),]-t(replicate(n_year,Two_FRE[j,]))
  }
  
  #Proof Residuals
  data_c1<-matrix(0,nrow=n_prefectures,ncol=n_age)
  data_c2<-matrix(0,nrow=n_prefectures,ncol=n_age)
  for (j in 1:n_prefectures) {
    data_c1[j,]=Two_FGE+Two_FRE[j,]+Two_FCE[1,]
    data_c2[j,]=Two_FGE+Two_FRE[j,]+Two_FCE[2,]
  }
  
  fixed_data_1<-matrix(0,nrow=(n_year*n_prefectures),ncol=n_age)
  fixed_data_2<-matrix(0,nrow=(n_year*n_prefectures),ncol=n_age)
  for (i in 1:n_prefectures) {
    fixed_data_1[(n_year*i-(n_year-1)):(n_year*i),]=do.call(rbind, replicate(n_year, data_c1[i,], simplify=FALSE))
    fixed_data_2[(n_year*i-(n_year-1)):(n_year*i),]=do.call(rbind, replicate(n_year, data_c2[i,], simplify=FALSE))
  }
  FD<-cbind(fixed_data_1,fixed_data_2)
  Recovered_data_1<-fixed_data_1+residuals_b1r
  Recovered_data_2<-fixed_data_2+residuals_b2r
  RD<-cbind(Recovered_data_1,Recovered_data_2)
  
  
  
  #reconstruction proof
  R1=all(round(all_male,4)==round(Recovered_data_1,4))
  R2=all(round(all_female,4)==round(Recovered_data_2,4))
  R<-c(R1,R2)
  
  return(list(residuals1= residuals_b1r,residuals2=residuals_b2r,rd=RD,R=R,Fixed_comp=FD))
  
}

###############################################################################
# Point  forecasts 
################################################################################
library(ftsa)
library(demography)
source("forecast_Arima.R")
Pref_forecasted_curves<-function(fixed_com,Residuals_f,
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
    dum = forecast_Arima(auto.arima(med_polish_resi_score[,ik]), h = fh, bootstrap = TRUE, npaths = B)
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


# R code for all the computations for the US mortality database with the proposed 
#Functional median polish ANOVA (FMP-ANOVA) approach.

# This script can be easily incorporated to the other two data sets (France and Japan)
# by simply changing the data set entries and the data. The data can easily loaded by
# changing the folder name (France or Japan)e.g., ./dataset_entries/France/

###########################################################################################
# 1. First install all the R packages and libraries required for the project
###########################################################################################


packages <- c("generics", "demography", "forecast","fda","fdaoutlier", 
              "rlist", "mrfDepth","ftsa","rainbow")

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

###########################################################################################
# 2.  Set a working directory
###########################################################################################
# This directory must constain the folder with all the codes for this project
setwd("~/My Drive/Fall 2023/STAT 397/PhD project 4/Revisions from JCGS/Final_Rcodes_editor/")
# Define an extra folder within the working directory where the results will be saved. 
# dirl.p is the path to that folder.
dirl.p <- "./Test_Results/"
##########################################################################################################
# 3. Load dataset entries and define the row (by states) and column partitions (by gender)
# Note: If you want to change to another dataset, please change the word USA onwards. (e.g. France/Japan) 
##########################################################################################################
dataset_entries<-readRDS("./dataset_entries/USA/dataset_entries.rds")
year = dataset_entries[[1]]
n_year = dataset_entries[[2]]
age = dataset_entries[[3]]
n_age = dataset_entries[[4]]
n_states=dataset_entries[[5]]

# Row partition
part_list = list()
for(ik in 1:n_states) {
  part_list[[ik]] = (n_year*ik-(n_year-1)):(n_year*ik)
}

#Column partition
n_populations=2
part_list_c = list()
for(ik in 1:n_populations) {
  part_list_c[[ik]] = (n_age*ik-(n_age-1)):(n_age*ik)
}
names_states <- readRDS("./names_states/USA/names_states.rds")
################################################################################
# This function removes the zeros in the dataset by replacing the values with the
# one in the instant inemdiatly before. This function is proposed by (Hallin et. al. 2023b)
################################################################################
remove_zeroes <- function(data_raw,n_states,n_year,n_age) {
  N       <- n_states
  T       <- n_year
  age_max <- n_age
  
  data_nozero <- list()
  for(i in 1:N) {
    data_nozero[[i]] <- matrix(data_raw[[i]],age_max,T)
  }
  for(i in 1:N) {
    for(j in 2:age_max) {#For j=1, the two zeroes have been removed before
      for(k in 1:T){ 
        if(data_nozero[[i]][j,k]==-Inf){
          data_nozero[[i]][j,k] <- data_nozero[[i]][j-1,k]
        }
      }
    }
  }
  return(data_nozero)
}


################################################################################
# 4.  Load the datasets
################################################################################

USA_male   <- readRDS("./datasets/USA/USA_male.rds")
USA_female <- readRDS("./datasets/USA/USA_female.rds")
# Remove the zeros in the datasets
USA_male   <- remove_zeroes(data_raw=USA_male,n_states,n_year,n_age) 
USA_female <- remove_zeroes(data_raw=USA_female,n_states,n_year,n_age)

# column bind the datasets for FMP-ANOVA decomposition
all_male    <- t(list.cbind(USA_male))
all_female  <- t(list.cbind(USA_female))

###########################################################################################
# 5. Apply the functional median polish decomposition and obtained the functional residuals
###########################################################################################
# 5.1 column bind all the dataset for both genders.
Y <- cbind(all_male,all_female)
# Note:  Y is a rectangular matrix with dimension (n_states x n_year) by 2* n_age
# 5.2 Apply the FMP-ANOVA decomposition and obtain the functional residuals
Residuals<- ftsa::Two_way_Residuals(Y,n_states,year,age,n_populations)
#. Functional residuals for the male population
Res1 <- Residuals$residuals1
# Functional residuals for the female population
Res2 <- Residuals$residuals2
# 5.3 Residuals for both populations.
Residuals_ <- cbind(Res1,Res2)
# Reconstructed data
RR <- Residuals$rd #Matrix with the original data reconstructed from the FMP decomposition
#  It's the proof of the reconstruction of the residuals. 
Residuals$R #The result should be a vector with two entries TRUE, TRUE.
#Indicating that after adding both deterministic and residual components the FTS are recovered.
# 5.4 deterministic components to be added up after forecasting
Fixed_part <- Residuals$Fixed_comp 

#############################################################################################
# 6. Computation of the point forecasts
#############################################################################################
source("New_Point_forecast.R") # This auxiliary file contain the functions for computing the 
# as well as their respective forecast errors. 

# There are several alternatives to compute the point forecasts. The implementation is in parallel
# In particular here, EVR denotes using the eigenvalue ratio criterion, K6 indicates that the
# retained PC are equal to 6. The following code uses the decomposition of the covariance operator and ARIMA as
# forecasting method for the functional scores. R denotes rolling and E expanding.


# 6.1 With a rolling window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
R_forecasted_FMP_cov_ARIMA_USA_EVR <- foreach(i = 1:n_states, .packages = c("ftsa")) %dopar% ForecastC(i,fixed_com=Fixed_part,Residuals_f=Residuals_,est_method = "cov",prediction_method="ARIMA",select_K="EVR", K=6)
R_forecasted_FMP_cov_ARIMA_USA_K6 <- foreach(i = 1:n_states, .packages = c("ftsa")) %dopar% ForecastC(i,fixed_com=Fixed_part,Residuals_f=Residuals_,est_method = "cov",prediction_method="ARIMA",select_K="Fixed", K=6)
stopCluster(cl)

# 6.2  With an expanding window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
E_forecasted_FMP_cov_ARIMA_USA_EVR <- foreach(i = 1:n_states, .packages = c("ftsa")) %dopar% ForecastC_expanding(i,fixed_com=Fixed_part,Residuals_f=Residuals_,est_method = "cov",prediction_method="ARIMA",select_K="EVR", K=6)
E_forecasted_FMP_cov_ARIMA_USA_K6 <- foreach(i = 1:n_states, .packages = c("ftsa")) %dopar% ForecastC_expanding(i,fixed_com=Fixed_part,Residuals_f=Residuals_,est_method = "cov",prediction_method="ARIMA",select_K="Fixed", K=6)
stopCluster(cl)

# The result from each computation is a list with n_states lists within. 
# For each state, the list contains two sub-lists: the first one is the 
# forecasted curves, the second is the percentage of retained variance. 

# It is recommended to save this results that contain the forecast curves either 
# for the rolling window or expanding window approach. Results will be saves in dirl.p defined before.
saveRDS(R_forecasted_FMP_cov_ARIMA_USA_EVR,paste0(dirl.p,"R_forecasted_FMP_cov_ARIMA_USA_3.rds"))
saveRDS(R_forecasted_FMP_cov_ARIMA_USA_K6,paste0(dirl.p,"R_forecasted_FMP_cov_ARIMA_USA_4.rds"))

###########################################################################################################
#7. Compute the error for the point forecasts obtained in the rolling window approach. 
#  A complete similar approach can be used to compute forecast errors from the expanding window approach.
###########################################################################################################
# This auxiliary file contains all the functions to compute the forcast errors
source(".forecast_errors.R")
# 7.1 Load the forecasted curves obtained either in the step 6.1 or 6.2. 
# e.g. Here it is loaded the results from the rolling window approach and with the EVR criterion.
forecasted_FMP_cov_ARIMA_USA <- R_forecasted_FMP_cov_ARIMA_USA_EVR
# Define the forecast horizon.
max_h=10
# Define the empty matrices where the errors will be saved.
error_MAPE_male<-matrix(0,nrow=max_h,ncol = n_states)
error_RMSPE_male<-matrix(0,nrow=max_h,ncol = n_states)
error_MAPE_female<-matrix(0,nrow=max_h,ncol = n_states)
error_RMSPE_female<-matrix(0,nrow=max_h,ncol = n_states)
for (i in 1:n_states) {
  forecasted_state_male=forecasted_FMP_cov_ARIMA_USA[[i]][[1]][1:n_age,]
  forecasted_state_female=forecasted_FMP_cov_ARIMA_USA[[i]][[1]][(n_age+1):(2*n_age),]
  True_pop_male=t(male[(n_year*i-(max_h-1)):(n_year*i),])
  True_pop_female=t(female[(n_year*i-(max_h-1)):(n_year*i),])
  for (j in 1:max_h) {
    error_MAPE_male[j,i]<-mape(forecast=forecasted_state_male[,j], true=True_pop_male[,j])
    error_RMSPE_male[j,i]<-rmspe(forecast=forecasted_state_male[,j], true=True_pop_male[,j])
    error_MAPE_female[j,i]<-mape(forecast=forecasted_state_female[,j], true=True_pop_female[,j])
    error_RMSPE_female[j,i]<-rmspe(forecast=forecasted_state_female[,j], true=True_pop_female[,j])
    
  }
  
}

# save all the errors in a matrix. This errirs are per each state at each forecast horizon 1:10
All_errors_rolling <- list(error_MAPE_male,error_MAPE_female,
                         error_RMSPE_male,error_RMSPE_female)
# Save the results
saveRDS(All_errors_rolling,paste0(dirl.p,"All_errors_R_FMP_cov_ARIMA_USA_EVR.rds"))

# 7.2 Average the results across the forecast horizon.
Errors_mean_rolling<-matrix(0,nrow <- n_states,ncol=length(All_errors_rolling))
Errors_sd_rolling<-matrix(0,nrow   <- n_states,ncol=length(All_errors_rolling))
for (i in 1:length(All_errors_rolling)) {
  error=All_errors_rolling[[i]]
  Errors_mean_rolling[,i] <- apply( error,2,mean)
  Errors_sd_rolling[,i]   <- apply( error,2,sd)
}
colnames(Errors_mean_rolling) <- c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
colnames(Errors_sd_rolling)   <- c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
rownames(Errors_mean_rolling) <- names_states
rownames(Errors_sd_rolling)   <- names_states

# Save the results of the forecast errors accross forecast horizon.
saveRDS(Errors_mean_rolling,paste0(dirl.p,"Errors_R_FMP_cov_ARIMA_USA_EVR.rds"))

#####################################################################################
#8.  Compute interval forecasts and coverage probability
#####################################################################################
# Split residuals per states. Each sub matrix will be of dimension n_year x n_age

male_states_res <- lapply(1:length(part_list), 
                            function(k){Res1[part_list[[k]], ]})
female_states_res <- lapply(1:length(part_list), 
                              function(k){Res2[part_list[[k]], ]})

# Split the fixed components by states. Each sub matrix will be of dimension n_year x n_age
male_states_fixed <- lapply(1:length(part_list), 
                              function(k){Fixed_part[part_list[[k]],1:n_age ]})

female_states_fixed <- lapply(1:length(part_list), 
                                function(k){Fixed_part[part_list[[k]],(n_age+1):(2*n_age) ]})

# This auxiliary function contains the functions and auxiliary files to compute the coverage probability.
source("Compute_coverages.R")

Emp_cov_male_0.85<-matrix(0,nrow=n_states,ncol = max_h)
Emp_cov_female_0.85<-matrix(0,nrow=n_states,ncol = max_h)
for (i in 1:n_states) {
  Emp_cov_male_0.85[i,]<-coverage_comp(sim_data=t(male_states_res[[i]]),fixed_comp=t(male_states_fixed[[i]]), sample_number=n_year, no_boot = 100,
                                       no_core=(detectCores()-2),K_val = NULL,prediction_method = "sieve_bootstrap"
                                       ,selection_ncomp_porder="CPV_AICC",method_pred = "Bosq",percent_CPV=0.85)
  Emp_cov_female_0.85[i,]<-coverage_comp(sim_data=t(female_states_res[[i]]),fixed_comp=t(female_states_fixed[[i]]), sample_number=n_year, no_boot = 100,
                                         no_core=(detectCores()-2),K_val = NULL,prediction_method = "sieve_bootstrap"
                                         ,selection_ncomp_porder="CPV_AICC",method_pred = "Bosq",percent_CPV=0.85)
}

Emp_cov_0.85<-cbind(Emp_cov_male_0.85,Emp_cov_female_0.85)
colnames(Emp_cov_0.85) = c("Pointwise 80% coverage", "Pointwise 95% coverage",
                           "80% interval score", "95% interval score",
                           "CPD Pointwise 80% coverage","CPD Pointwise 95% coverage",
                           "Pointwise 80% coverage", "Pointwise 95% coverage",
                           "80% interval score", "95% interval score",
                           "CPD Pointwise 80% coverage","CPD Pointwise 95% coverage")

rownames(Emp_cov_0.85)<-names_states
# R code for all the computations for the USA mortality database with the  
#Functional mean ANOVA (FM-ANOVA) approach.

#1. This scripts can be easily incorporated to the other two data sets (France and Japan)
# by simply changing the data set entries and the data.
# changing the folder name (France or Japan). e.g., ./dataset_entries/France/
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
################################################################################
# Data set entries
################################################################################
dataset_entries<-readRDS("./dataset_entries/USA/dataset_entries.rds")
year = dataset_entries[[1]]
n_year = dataset_entries[[2]]
age = dataset_entries[[3]]
n_age = dataset_entries[[4]]
n_prefectures=dataset_entries[[5]]

# Row partition
part_list = list()
for(ik in 1:n_prefectures) {
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
# Data set 
################################################################################
USA_male <- readRDS("./datasets/USA/USA_male.rds")
USA_female <- readRDS("./datasets/USA/USA_female.rds")
all_male<-t(list.cbind(USA_male))
all_female<-t(list.cbind(USA_female))
########################################################################
# Comparison with the functional mean ANOVA approach (FM-ANOVA)
########################################################################
source("FM-ANOVA_decomposition.R")
# This function computes the functional mean ANOVA decomposition based on means
FANOVA_means<-FANOVA(data_pop1=t(all_male),data_pop2=t(all_female),n_year
                     ,n_prefectures,n_age,n_populations=2,row_par=part_list)
#This function computes the functional residuals after removing the deterministic components
# obtained from the FANOVA function.
Residuals_means<-Two_way_Residuals_means(FANOVA_means,data_pop1=t(all_male),data_pop2=t(all_female)
                                         ,n_prefectures,n_year,n_age)
Res1_means=Residuals_means$residuals1_mean
Res2_means=Residuals_means$residuals2_mean
Residuals_mean<-cbind(Res1_means,Res2_means)
# Reconstructed data
RR<-Residuals_means$rd #Matrix with the original data reconstructed from the FMP decomposition
#  It's the proof of the reconstruction of the residuals. 
Residuals_means$R #The result should be a vector with two entries TRUE, TRUE.
#Indicating that after adding both deterministic and time-varying components the FTS are recovered.
Fixed_part_means<-Residuals_means$Fixed_comp_mean # deterministic components to be added up after forecasting

################################################################################
# Computation of the point forecasts based on functional mean ANOVA (FM-ANOVA)
################################################################################
source("Point_forecasting_FANOVA.R")
# With a rolling window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
forecasted_curves_triangular_means_rolling <- foreach(i = 1:n_prefectures, .packages = c("ftsa")) %dopar% ForecastC_mean(i,fixed_com=Fixed_part_means,Residuals_f=Residuals_mean)
stopCluster(cl)

# With an expanding window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
forecasted_curves_triangular_means_expanding <- foreach(i = 1:n_prefectures, .packages = c("ftsa")) %dopar% ForecastC_mean_expanding(i,fixed_com=Fixed_part_means,Residuals_f=Residuals_mean)
stopCluster(cl)

##############################################################################################################
#Compute the error based on FM-ANOVA approach, with the forecasts obtained in the rolling window approach
##############################################################################################################
source("forecast_errors.R")
max_h=10
error_MAPE_means_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_RMSPE_means_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_MAPE_means_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_RMSPE_means_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
for (i in 1:n_prefectures) {
  forecasted_pref_male=forecasted_curves_triangular_means_rolling[[i]][1:n_age,]
  forecasted_pref_female=forecasted_curves_triangular_means_rolling[[i]][(n_age+1):(2*n_age),]
  True_pop_male=t(male[(n_prefectures*i-(max_h-1)):(n_prefectures*i),])
  True_pop_female=t(female[(n_prefectures*i-(max_h-1)):(n_prefectures*i),])
  for (j in 1:max_h) {
    error_MAPE_means_male_1[j,i]<-mape(forecast=forecasted_pref_male[,j], true=True_pop_male[,j])
    error_RMSPE_means_male_1[j,i]<-rmspe(forecast=forecasted_pref_male[,j], true=True_pop_male[,j])
    error_MAPE_means_female_1[j,i]<-mape(forecast=forecasted_pref_female[,j], true=True_pop_female[,j])
    error_RMSPE_means_female_1[j,i]<-rmspe(forecast=forecasted_pref_female[,j], true=True_pop_female[,j])
    
  }
  
}

All_errors_mean_rolling<-list(error_MAPE_means_male,error_MAPE_means_female,
                                  error_RMSPE_means_male,error_RMSPE_means_female)
Errors_mean_basedmeans_rolling<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_mean_rolling))
Errors_sd_basedmeans_rolling<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_mean_rolling))
for (i in 1:length(All_errors_mean_rolling)) {
  error=All_errors_mean_rolling[[i]]
  Errors_mean_rolling[,i]=apply( error,2,mean)
  Errors_sd_rolling[,i]=apply( error,2,sd)
}
colnames(Errors_mean_rolling)<-c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
colnames(Errors_sd_rolling)<-c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
rownames(Errors_mean_rolling)<-names_states
rownames(Errors_sd_rolling)<-names_states

##############################################################################################################
#Compute the error based on FM-ANOVA approach, with the forecasts obtained in the expanding window approach
##############################################################################################################
source("forecast_errors.R")
max_h=10
error_MAPE_means_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_RMSPE_means_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_MAPE_means_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_RMSPE_means_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
for (i in 1:n_prefectures) {
  forecasted_pref_male=forecasted_curves_triangular_means_expanding[[i]][1:n_age,]
  forecasted_pref_female=forecasted_curves_triangular_means_expanding[[i]][(n_age+1):(2*n_age),]
  True_pop_male=t(male[(n_prefectures*i-(max_h-1)):(n_prefectures*i),])
  True_pop_female=t(female[(n_prefectures*i-(max_h-1)):(n_prefectures*i),])
  for (j in 1:max_h) {
    error_MAPE_means_male_1[j,i]<-mape(forecast=forecasted_pref_male[,j], true=True_pop_male[,j])
    error_RMSPE_means_male_1[j,i]<-rmspe(forecast=forecasted_pref_male[,j], true=True_pop_male[,j])
    error_MAPE_means_female_1[j,i]<-mape(forecast=forecasted_pref_female[,j], true=True_pop_female[,j])
    error_RMSPE_means_female_1[j,i]<-rmspe(forecast=forecasted_pref_female[,j], true=True_pop_female[,j])
    
  }
  
}

All_errors_mean_expanding<-list(error_MAPE_means_male,error_MAPE_means_female,
                              error_RMSPE_means_male,error_RMSPE_means_female)
Errors_mean_basedmeans_expanding<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_mean_expanding))
Errors_sd_basedmeans_expanding<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_mean_expanding))
for (i in 1:length(All_errors_mean_expanding)) {
  error=All_errors_mean_expanding[[i]]
  Errors_mean_expanding[,i]=apply( error,2,mean)
  Errors_sd_expanding[,i]=apply( error,2,sd)
}
colnames(Errors_mean_expanding)<-c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
colnames(Errors_sd_expanding)<-c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
rownames(Errors_mean_expanding)<-names_states
rownames(Errors_sd_expanding)<-names_states

########################################################################
# Compute interval forecasts and coverage probability
########################################################################

# Split residuals per prefecture
male_prefecture_res_means<-lapply(1:length(part_list), 
                                  function(k){Res1_means[part_list[[k]], ]})
female_prefecture_res_means<-lapply(1:length(part_list), 
                                    function(k){Res2_means[part_list[[k]], ]})

#split the fixed components by prefecture
male_prefecture_fixed_means<-lapply(1:length(part_list), 
                                    function(k){Fixed_part_means[part_list[[k]],1:n_age ]})

female_prefecture_fixed_means<-lapply(1:length(part_list), 
                                      function(k){Fixed_part_means[part_list[[k]],(n_age+1):(2*n_age) ]})

source("Compute_coverages.R")
Emp_cov_male_0.85_means<-matrix(0,nrow=n_prefectures,ncol = max_h)
Emp_cov_female_0.85_means<-matrix(0,nrow=n_prefectures,ncol = max_h)
for (i in 1:n_prefectures) {
  Emp_cov_male_0.85_means[i,]<-coverage_comp(sim_data=t(male_prefecture_res_means[[i]]),fixed_comp=t(male_prefecture_fixed_means[[i]]), sample_number=n_year, no_boot = 100,
                                             no_core=(detectCores()-2),K_val = NULL,prediction_method = "sieve_bootstrap"
                                             ,selection_ncomp_porder="CPV_AICC",method_pred = "Bosq",percent_CPV=0.85)
  Emp_cov_female_0.85_means[i,]<-coverage_comp(sim_data=t(female_prefecture_res_means[[i]]),fixed_comp=t(female_prefecture_fixed_means[[i]]), sample_number=n_year, no_boot = 100,
                                               no_core=(detectCores()-2),K_val = NULL,prediction_method = "sieve_bootstrap"
                                               ,selection_ncomp_porder="CPV_AICC",method_pred = "Bosq",percent_CPV=0.85)
}

Emp_cov_0.85_means<-cbind(Emp_cov_male_0.85_means,Emp_cov_female_0.85_means)
colnames(Emp_cov_0.85_means) = c("Pointwise 80% coverage", "Pointwise 95% coverage",
                                   "80% interval score", "95% interval score",
                                   "CPD Pointwise 80% coverage","CPD Pointwise 95% coverage",
                                   "Pointwise 80% coverage", "Pointwise 95% coverage",
                                   "80% interval score", "95% interval score",
                                   "CPD Pointwise 80% coverage","CPD Pointwise 95% coverage")

rownames(Emp_cov_0.85_means)<-names_departments

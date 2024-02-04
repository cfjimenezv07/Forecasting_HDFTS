# Replicate Tavakoli et. al. 2023 (HNT) and Gao et al. 2019 (GSY) methods for our datasets. 
# These Rcodes are public posted in https://zenodo.org/record/7408999 from HNT23, it includes GSY19.
# We adapted to our paper cases.


# This script can be easily incorporated to the other two data sets (France and Japan)
# by simply changing the data set entries and the data. The data can easily loaded by
# changing the folder name (France or Japan) e.g., ./dataset_entries/France/

###########################################################################################
# 1. First install all the R packages and libraries required for the project
###########################################################################################
packages <- c("generics", "demography", "forecast","fda","fdaoutlier", 
              "rlist", "mrfDepth","ftsa","rainbow", 
              "foreach", "doParallel", "vars", "fda.usc", "far",
              "gdata", "MTS", "RSpectra", "tictoc", "VAR.etp")

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
# one in the instant inemdiatly before. This function is proposed by (Tavakoli et. al. 2023)
################################################################################
remove_zeroes <- function(data_raw,n_states,n_year,n_age) {
  N       <- n_states
  t.max       <- n_year
  age_max <- n_age
  
  data_nozero <- list()
  for(i in 1:N) {
    data_nozero[[i]] <- matrix(data_raw[[i]],age_max,t.max)
  }
  for(i in 1:N) {
    for(j in 2:age_max) {#For j=1, the two zeroes have been removed before
      for(k in 1:t.max){ 
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

################################################################################
# auxiliary functions for the computation of the HNT and GSY methods
################################################################################
source('./MITS_class.R') 
source("aux_HNT.R")


################################################################################
# For the data set: USA female
################################################################################

order <- 3
r <- 3
p <- 9
h_max <- 10 # forecast horizon
years <-  year
tmp <- init_data_new(data=USA_female,p,Names=names_states)
data_female <- tmp[[1]]  # datasets for Gao2019
mixed_fts_female <- tmp[[2]]
# basis splines
basis <- create.bspline.basis(c(0, 1), nbasis = p, norder = 4)
args <- seq(0, 1, length=n_age)


model_female <- list()
training_set_female <- list()
# Training model for Gao 2019 approach
model_female <- list()
for(t in 1:h_max){ 
  training_set_female[[years[n_year+t-h_max]]] <- trim_data_new(data_female,t)[[1]]
  model_female[[years[n_year+t-h_max]]] <- hdfpca(training_set_female[[years[n_year+t-h_max]]], order, q = sqrt(dim(training_set_female[[years[n_year+t-h_max]]][[1]])[2]),r) 
}


predictions_GAO_female <- list()
predictions_HNT_female <- list()
for(h in 1:h_max){
  predictions_GAO_female[[h]] <- array(0,dim=c(n_states,1,n_age))
  predictions_HNT_female[[h]] <- array(0,dim=c(n_states,1,n_age))
}

r_hat_female <- list()
cp_female <- list()
for(t in 1:h_max){
  cp_female[[years[n_year+t-h_max]]] <- mixed_fts_female$copy()
  cp_female[[years[n_year+t-h_max]]]$trimSampleSize(n_year+t-h_max)
  r_hat_female[[years[n_year+t-h_max]]] <-   est_r_abc(cp_female[[years[n_year+t-h_max]]],h_max)
  print(paste('estimated r for year ',years[n_year+t-h_max],': ',r_hat_female[[years[n_year+t-h_max]]]))
}


for(h in 1:h_max){
  predictions_HNT_female[[h]][,1,] <- predict_NTH_rknown_v2(cp_female[[years[n_year+h-h_max]]],h,p,basis,args,r_hat_female[[years[n_year+h-h_max]]])
  
}


for(t in 1:h_max){ 
  tmp <- forecast.hdfpca(model_female[[years[n_year+t-h_max]]],h=h_max,B=1)
  for(i in 1:n_states){
    predictions_GAO_female[[t]][i,,] <- tmp$forecast[[i]][,t]
  }
}


# It is recommended to save the results in the folder /Test_Results

saveRDS(predictions_HNT_female,paste0(dirl.p,"predictions_HNT_female_USA.rds"))
saveRDS(predictions_GAO_female,paste0(dirl.p,"predictions_GAO_female_USA.rds"))

################################################################################
# For the dataset: USA_male
################################################################################


tmp <- init_data_new(data=USA_male,p,Names=names_states)
data_male <- tmp[[1]]  # datasets for Gao2019
mixed_fts_male <- tmp[[2]]
# basis splines
basis <- create.bspline.basis(c(0, 1), nbasis = p, norder = 4)
args <- seq(0, 1, length=n_age)
model_male <- list()
training_set_male <- list()


# Training model for Gao 2019
model_male <- list()
for(t in 1:h_max){ 
  training_set_male[[years[n_year+t-h_max]]] <- trim_data_new(data_male,t)[[1]]
  model_male[[years[n_year+t-h_max]]] <- hdfpca(training_set_male[[years[n_year+t-h_max]]], order, q = sqrt(dim(training_set_female[[years[n_year+t-h_max]]][[1]])[2]),r) 
}


predictions_GAO_male <- list()
predictions_HNT_male <- list()
for(h in 1:h_max){
  predictions_GAO_male[[h]] <- array(0,dim=c(n_states,1,n_age))
  predictions_HNT_male[[h]] <- array(0,dim=c(n_states,1,n_age))
}

r_hat_male <- list()
cp_male <- list()
for(t in 1:h_max){
  cp_male[[years[n_year+t-h_max]]] <- mixed_fts_male$copy()
  cp_male[[years[n_year+t-h_max]]]$trimSampleSize(n_year+t-h_max)
  r_hat_male[[years[n_year+t-h_max]]] <-   est_r_abc(cp_male[[years[n_year+t-h_max]]],h_max)
  print(paste('estimated r for year ',years[n_year+t-h_max],': ',r_hat_male[[years[n_year+t-h_max]]]))
}


for(h in 1:h_max){
  predictions_HNT_male[[h]][,1,] <- predict_NTH_rknown_v2(cp_male[[years[n_year+h-h_max]]],h,p,basis,args,r_hat_male[[years[n_year+h-h_max]]])
  
}


for(t in 1:h_max){ 
  tmp <- forecast.hdfpca(model_male[[years[n_year+t-h_max]]],h=h_max,B=1)
  for(i in 1:n_states){
    predictions_GAO_male[[t]][i,,] <- tmp$forecast[[i]][,t]
  }
}

# It is recommended to save the results in the folder /Test_Results

saveRDS(predictions_HNT_male,paste0(dirl.p,"predictions_HNT_male_USA.rds"))
saveRDS(predictions_GAO_male,paste0(dirl.p,"predictions_GAO_male_USA.rds"))

##############################################################################################################
#Compute the point forecast errors obtained in the rolling window approach for both GSY and HNT methods.
##############################################################################################################
# Load the predictions previously computed
predictions_GAO_male_USA   <- readRDS("./Test_Results/predictions_GAO_male_USA.rds")
predictions_GAO_female_USA <- readRDS("./Test_Results/predictions_GAO_female_USA.rds")
predictions_HNT_male_USA   <- readRDS("./Test_Results/predictions_HNT_male_USA.rds")
predictions_HNT_female_USA <- readRDS("./Test_Results/predictions_HNT_female_USA.rds")


source("./forecast_errors.R")
GAO_error_MAPE_male<-matrix(0,nrow=h_max,ncol = n_states)
GAO_error_RMSPE_male<-matrix(0,nrow=h_max,ncol = n_states)
GAO_error_MAPE_female<-matrix(0,nrow=h_max,ncol = n_states)
GAO_error_RMSPE_female<-matrix(0,nrow=h_max,ncol = n_states)
HNT_error_MAPE_male<-matrix(0,nrow=h_max,ncol = n_states)
HNT_error_RMSPE_male<-matrix(0,nrow=h_max,ncol = n_states)
HNT_error_MAPE_female<-matrix(0,nrow=h_max,ncol = n_states)
HNT_error_RMSPE_female<-matrix(0,nrow=h_max,ncol = n_states)
for (i in 1:n_states) {
  for (j in 1:h_max) {
    GAO_forecasted_pref_male=predictions_GAO_male_USA[[j]][i,1,]
    GAO_forecasted_pref_female=predictions_GAO_female_USA[[j]][i,1,]
    HNT_forecasted_pref_male=predictions_HNT_male_USA[[j]][i,1,]
    HNT_forecasted_pref_female=predictions_HNT_female_USA[[j]][i,1,]
    True_pop_male= data_male[[i]][,n_year+j-h_max]
    True_pop_female=data_female[[i]][,n_year+j-h_max]
    GAO_error_MAPE_male[j,i]<-mape(forecast=GAO_forecasted_pref_male, true=True_pop_male)
    GAO_error_RMSPE_male[j,i]<-rmspe(forecast=GAO_forecasted_pref_male, true=True_pop_male)
    GAO_error_MAPE_female[j,i]<-mape(forecast=GAO_forecasted_pref_female, true=True_pop_female)
    GAO_error_RMSPE_female[j,i]<-rmspe(forecast=GAO_forecasted_pref_female, true=True_pop_female)
    HNT_error_MAPE_male[j,i]<-mape(forecast=HNT_forecasted_pref_male, true=True_pop_male)
    HNT_error_RMSPE_male[j,i]<-rmspe(forecast=HNT_forecasted_pref_male, true=True_pop_male)
    HNT_error_MAPE_female[j,i]<-mape(forecast=HNT_forecasted_pref_female, true=True_pop_female)
    HNT_error_RMSPE_female[j,i]<-rmspe(forecast=HNT_forecasted_pref_female, true=True_pop_female)
  }
  
}

All_errors_rolling<-list(GAO_error_MAPE_male,GAO_error_MAPE_female,
                         GAO_error_RMSPE_male,GAO_error_RMSPE_female,
                         HNT_error_MAPE_male,HNT_error_MAPE_female,
                         HNT_error_RMSPE_male,HNT_error_RMSPE_female)

Errors_mean_rolling<-matrix(0,nrow=n_states,ncol=length(All_errors_rolling))
Errors_sd_rolling<-matrix(0,nrow=n_states,ncol=length(All_errors_rolling))
for (i in 1:length(All_errors_rolling)) {
  error=All_errors_rolling[[i]]
  Errors_mean_rolling[,i]=apply( error,2,mean)
  Errors_sd_rolling[,i]=apply( error,2,sd)
}
colnames(Errors_mean_rolling)<-c("Gao_MAPE_male","Gao_MAPE_female","Gao_RMSPE_male","Gao_RMSPE_female",
                                 "HNT_MAPE_male","HNT_MAPE_female","HNT_RMSPE_male","HNT_RMSPE_female")
colnames(Errors_sd_rolling)<-c("Gao_MAPE_male","Gao_MAPE_female","Gao_RMSPE_male","Gao_RMSPE_female",
                               "HNT_MAPE_male","HNT_MAPE_female","HNT_RMSPE_male","HNT_RMSPE_female")
rownames(Errors_mean_rolling)<-names_states
rownames(Errors_sd_rolling)<-names_states




# It is recommended to save the results in the folder /Test_Results
# 1. For GSY19 
saveRDS(Errors_mean_rolling[,1:4],paste0(dirl.p,"Errors_GAO_USA.rds"))
# 2. For HNT23 
saveRDS(Errors_mean_rolling[,5:8],paste0(dirl.p,"Errors_HNT_USA.rds"))









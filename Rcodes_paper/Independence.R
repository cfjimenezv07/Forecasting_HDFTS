# R code for the computations of the naive approach (independence)
# this R code computes point forecasts assuming independece across genders and states

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
# Data set 
################################################################################
USA_male   <- readRDS("./datasets/USA/USA_male.rds")
USA_female <- readRDS("./datasets/USA/USA_female.rds")
# Remove the zeros in the datasets
USA_male   <- remove_zeroes(data_raw=USA_male,n_states,n_year,n_age) 
USA_female <- remove_zeroes(data_raw=USA_female,n_states,n_year,n_age)

# column bind the datasets for later decomposing per state and gender
all_male    <- t(list.cbind(USA_male))
all_female  <- t(list.cbind(USA_female))


#Data per state per population
states_male<-list()
states_female<-list()
for (i in 1:n_states) {
  states_male[[i]]=all_male[(n_year*i-(n_year-1)):(n_year*i),]
  states_female[[i]]=all_female[(n_year*i-(n_year-1)):(n_year*i),]
}


##########################################################################################
# Naive Forecasting for 10 years ahead with the function forecast.ftsm from the Ftsa Rcpkg.
##########################################################################################
naive_forecast_male_10<-list()
naive_forecast_female_10<-list()
h_max=10
for (j in 1:n_states) {
  data_male   <- t(states_male[[j]])[,1:(n_year-h_max)]
  data_female <- t(states_female[[j]])[,1:(n_year-h_max)]
  colnames(data_male)<-1:(n_year-h_max)
  colnames(data_female)<-1:(n_year-h_max)
  data_male   <-rainbow::fts(x=1:n_age, y=data_male, xname='age',yname='USA')
  data_female <-rainbow::fts(x=1:n_age, y=data_female, xname='age',yname='USA')
  ftsm_male<-ftsm(data_male)
  ftsm_female<-ftsm(data_female)
  naive_forecast_male_10[[j]]<-forecast.ftsm(ftsm_male)$mean$y
  naive_forecast_female_10[[j]]<-forecast.ftsm(ftsm_female)$mean$y
}


################################################################################
# Compute the errors for the naive approach assuming independence
################################################################################
source("./forecast_errors.R")
error_MAPE_male_naive<-matrix(0,nrow=h_max,ncol = n_states)
error_RMSPE_male_naive<-matrix(0,nrow=h_max,ncol = n_states)
error_MAPE_female_naive<-matrix(0,nrow=h_max,ncol = n_states)
error_RMSPE_female_naive<-matrix(0,nrow=h_max,ncol = n_states)

for (i in 1:n_states) {
  forecasted_pref_male_naive=naive_forecast_male_10[[i]]
  forecasted_pref_female_naive=naive_forecast_female_10[[i]]
  True_pop_male_naive=t(states_male[[i]])
  True_pop_female_naive=t(states_female[[i]])
  for (j in 1:h_max) {
    error_MAPE_male_naive[j,i]<-mape(forecast=forecasted_pref_male_naive[,j], true=True_pop_male_naive[,j])
    error_RMSPE_male_naive[j,i]<-rmspe(forecast=forecasted_pref_male_naive[,j], true=True_pop_male_naive[,j])
    error_MAPE_female_naive[j,i]<-mape(forecast=forecasted_pref_female_naive[,j], true=True_pop_female_naive[,j])
    error_RMSPE_female_naive[j,i]<-rmspe(forecast=forecasted_pref_female_naive[,j], true=True_pop_female_naive[,j])
    
  }
  
}


All_errors_naive<-list(error_MAPE_male_naive,error_MAPE_female_naive
                       ,error_RMSPE_male_naive,error_RMSPE_female_naive)

Errors_mean_naive<-matrix(0,nrow=n_states,ncol=length(All_errors_naive))
Errors_sd_naive<-matrix(0,nrow=n_states,ncol=length(All_errors_naive))
for (i in 1:length(All_errors_naive)) {
  error=All_errors_naive[[i]]
  Errors_mean_naive[,i]=apply( error,2,mean)
  Errors_sd_naive[,i]=apply( error,2,sd)
}
colnames(Errors_mean_naive)<-c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
colnames(Errors_sd_naive)<-c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
rownames(Errors_mean_naive)<-names_states
rownames(Errors_sd_naive)<-names_states

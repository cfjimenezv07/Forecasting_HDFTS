# R code for the computations of the naive approach (independence)
# this R code computes point forecasts assuming independece across genders and states
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

#Data per prefecture per population
prefectures_male<-list()
prefectures_female<-list()
for (i in 1:n_prefectures) {
  prefectures_male[[i]]=all_male[(n_year*i-(n_year-1)):(n_year*i),]
  prefectures_female[[i]]=all_female[(n_year*i-(n_year-1)):(n_year*i),]
}


################################################################################
# Naive Forecasting for 10 years ahead
################################################################################
naive_forecast_male_10<-list()
naive_forecast_female_10<-list()
h_max=10
for (j in 1:n_prefectures) {
  data_male   <- t(prefectures_male[[j]])[,1:(n_year-h_max)]
  data_female <- t(prefectures_female[[j]])[,1:(n_year-h_max)]
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
source("forecast_errors.R")
error_MAPE_male_naive<-matrix(0,nrow=h_max,ncol = n_prefectures)
error_RMSPE_male_naive<-matrix(0,nrow=h_max,ncol = n_prefectures)
error_MAPE_female_naive<-matrix(0,nrow=h_max,ncol = n_prefectures)
error_RMSPE_female_naive<-matrix(0,nrow=h_max,ncol = n_prefectures)

for (i in 1:n_prefectures) {
  forecasted_pref_male_naive=naive_forecast_male_10[[i]]
  forecasted_pref_female_naive=naive_forecast_female_10[[i]]
  True_pop_male_naive=t(prefectures_male[[i]])
  True_pop_female_naive=t(prefectures_female[[i]])
  for (j in 1:h_max) {
    error_MAPE_male_naive[j,i]<-mape(forecast=forecasted_pref_male_naive[,j], true=True_pop_male_naive[,j])
    error_RMSPE_male_naive[j,i]<-rmspe(forecast=forecasted_pref_male_naive[,j], true=True_pop_male_naive[,j])
    error_MAPE_female_naive[j,i]<-mape(forecast=forecasted_pref_female_naive[,j], true=True_pop_female_naive[,j])
    error_RMSPE_female_naive[j,i]<-rmspe(forecast=forecasted_pref_female_naive[,j], true=True_pop_female_naive[,j])
    
  }
  
}


All_errors_naive<-list(error_MAPE_male_naive,error_MAPE_female_naive
                       ,error_RMSPE_male_naive,error_RMSPE_female_naive)

Errors_mean_naive<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_naive))
Errors_sd_naive<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_naive))
for (i in 1:length(All_errors_naive)) {
  error=All_errors_naive[[i]]
  Errors_mean_naive[,i]=apply( error,2,mean)
  Errors_sd_naive[,i]=apply( error,2,sd)
}
colnames(Errors_mean_naive)<-c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
colnames(Errors_sd_naive)<-c("MAPE_male","MAPE_female","RMSPE_male","RMSPE_female")
rownames(Errors_mean_naive)<-names_states
rownames(Errors_sd_naive)<-names_states

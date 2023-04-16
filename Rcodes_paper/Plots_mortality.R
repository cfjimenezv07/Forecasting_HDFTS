#Plots

packages <- c("generics", "demography", "forecast","fda","fdaoutlier", "rlist", "mrfDepth","ftsa","rainbow")
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




#Plots for US 
################################################################################
year = 1959:2020
n_year = length(year)
age = 0:100
n_age = length(age)
n_prefectures=51
USA_male <- readRDS("~/My Drive/Spring 2023/STAT 397/PhD project 4/Rcodes/USA_male.rds")
USA_female <- readRDS("~/My Drive/Spring 2023/STAT 397/PhD project 4/Rcodes/USA_female.rds")
all_male<-t(list.cbind(USA_male))
all_female<-t(list.cbind(USA_female))


# Plot mortality for the US

Male<-matrix(0,nrow = n_age,ncol=n_year)
Female<-matrix(0,nrow = n_age,ncol=n_year)
for (j in 1:n_year) {
  Male[,j]<-colMeans(all_male[seq(from=j,to=n_year*n_prefectures,by=n_year),]) 
  Female[,j]<-colMeans(all_female[seq(from=j,to=n_year*n_prefectures,by=n_year),])
}


 colnames(Male)<-1:dim(Male)[2]
  curves_male<-rainbow::fts(x=1:n_age, y=Male , xname='Age',yname='Log death rate')
  plot.fds(curves_male,main="USA: male death rates (1959-2020)")
  colnames(Female)<-1:dim(Female)[2]
  curves_female<-rainbow::fts(x=1:n_age, y=Female , xname='Age',yname='Log death rate')
  plot.fds(curves_female,main="USA: female death rates (1959-2020)")

  
  #Plots for US 
  year = 1959:2020
  n_year = length(year)
  age = 0:100
  n_age = length(age)
  n_prefectures=51
  USA_male <- readRDS("~/My Drive/Spring 2023/STAT 397/PhD project 4/Rcodes/USA_male.rds")
  USA_female <- readRDS("~/My Drive/Spring 2023/STAT 397/PhD project 4/Rcodes/USA_female.rds")
  all_male<-t(list.cbind(USA_male))
  all_female<-t(list.cbind(USA_female))
  
  
  # Plot mortality for the US
  
  Male<-matrix(0,nrow = n_age,ncol=n_year)
  Female<-matrix(0,nrow = n_age,ncol=n_year)
  for (j in 1:n_year) {
    Male[,j]<-colMeans(all_male[seq(from=j,to=n_year*n_prefectures,by=n_year),]) 
    Female[,j]<-colMeans(all_female[seq(from=j,to=n_year*n_prefectures,by=n_year),])
  }
  
  library(refund.shiny)
  
  colnames(Male)<-1:dim(Male)[2]
  curves_male<-rainbow::fts(x=1:n_age, y=Male , xname='Age',yname='Log death rate')
  savepdf("~/My Drive/Spring 2023/STAT 397/PhD project 4/Plots_paper/USA_mortality_male", width = 12, height = 10, toplines = 0.8)
  plot.fds(curves_male,main="USA: male death rates (1959-2020)")
  
  
  colnames(Female)<-1:dim(Female)[2]
  curves_female<-rainbow::fts(x=1:n_age, y=Female , xname='Age',yname='Log death rate')
  plot.fds(curves_female,main="USA: female death rates (1959-2020)")
  
  
  
  
################################################################################
# For France
  year = 1968:2021
  n_year = length(year)
  age = 0:100
  n_age = length(age)
  n_prefectures=95
  
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
  
  # load the data
  France_male <- readRDS("~/My Drive/Spring 2023/STAT 397/PhD project 4/Rcodes/France_male.rds")
  France_female <- readRDS("~/My Drive/Spring 2023/STAT 397/PhD project 4/Rcodes/France_female.rds")
  all_male<-t(list.cbind(France_male))
  all_female<-t(list.cbind(France_female))
  
  # Plot mortality for the US
  
  Male<-matrix(0,nrow = n_age,ncol=n_year)
  Female<-matrix(0,nrow = n_age,ncol=n_year)
  for (j in 1:n_year) {
    Male[,j]<-colMeans(all_male[seq(from=j,to=n_year*n_prefectures,by=n_year),]) 
    Female[,j]<-colMeans(all_female[seq(from=j,to=n_year*n_prefectures,by=n_year),])
  }
  
  
  colnames(Male)<-1:dim(Male)[2]
  curves_male<-rainbow::fts(x=1:n_age, y=Male , xname='Age',yname='Log death rate')
  plot.fds(curves_male,main="France: male death rates (1968-2021)")
  colnames(Female)<-1:dim(Female)[2]
  curves_female<-rainbow::fts(x=1:n_age, y=Female , xname='Age',yname='Log death rate')
  plot.fds(curves_female,main="France: female death rates (1968-2021)")
  
  
  ################################################################################ 
  #Plots for Japan
  ################################################################################
  year = 1975:2020
  n_year = length(year)
  age = 0:98
  n_age = length(age)
  n_prefectures=47
  library(readr)
  all_male <- as.matrix(read_csv("all_male"))
  all_female <- as.matrix(read_csv("all_female"))

  
  
  # Plot mortality for the US
  
  Male<-matrix(0,nrow = n_age,ncol=n_year)
  Female<-matrix(0,nrow = n_age,ncol=n_year)
  for (j in 1:n_year) {
    Male[,j]<-colMeans(all_male[seq(from=j,to=n_year*n_prefectures,by=n_year),]) 
    Female[,j]<-colMeans(all_female[seq(from=j,to=n_year*n_prefectures,by=n_year),])
  }
  
  
  colnames(Male)<-1:dim(Male)[2]
  curves_male<-rainbow::fts(x=1:n_age, y=Male , xname='Age',yname='Log death rate')
  plot.fds(curves_male,main="Japan: male death rates (1975-2020)")
  colnames(Female)<-1:dim(Female)[2]
  curves_female<-rainbow::fts(x=1:n_age, y=Female , xname='Age',yname='Log death rate')
  plot.fds(curves_female,main="Japan: female death rates (1975-2020)")
  
  ################################################################################ 
  #Plots for Australia
  ################################################################################
  year = 1971:2021
  n_year = length(year)
  age = 0:99
  n_age = length(age)
  n_prefectures=8
  AUS_male <- readRDS("~/My Drive/Spring 2023/STAT 397/PhD project 4/Rcodes/AUS_male.rds")
  AUS_female <- readRDS("~/My Drive/Spring 2023/STAT 397/PhD project 4/Rcodes/AUS_female.rds")
  all_male<-t(list.cbind(AUS_male))
  all_female<-t(list.cbind(AUS_female))
  
  
  # Plot mortality for the US
  
  Male<-matrix(0,nrow = n_age,ncol=n_year)
  Female<-matrix(0,nrow = n_age,ncol=n_year)
  for (j in 1:n_year) {
    Male[,j]<-colMeans(all_male[seq(from=j,to=n_year*n_prefectures,by=n_year),]) 
    Female[,j]<-colMeans(all_female[seq(from=j,to=n_year*n_prefectures,by=n_year),])
  }
  
  
  colnames(Male)<-1:dim(Male)[2]
  curves_male<-rainbow::fts(x=1:n_age, y=Male , xname='Age',yname='Log death rate')
  plot.fds(curves_male,main="Australia: male death rates (1971-2021)")
  colnames(Female)<-1:dim(Female)[2]
  curves_female<-rainbow::fts(x=1:n_age, y=Female , xname='Age',yname='Log death rate')
  plot.fds(curves_female,main="Australia: female death rates (1971-2021)")
  
  
  
  
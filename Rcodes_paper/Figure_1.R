# Code to reproduce the Figure 1 in the paper. 

################################################################
# Set a working directory
################################################################
setwd("~/My Drive/Fall 2023/STAT 397/PhD project 4/Revisions from JCGS/Final_Rcodes_editor")
# All plots from this code will be saved in the folder plots2paper in the working directory
rm(list = ls())

################################################################

# Function to save the plots.
savepdf <- function(file, width=16, height=10)
{
  fname <- paste(file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}

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
USA_male <- readRDS("./datasets/USA/USA_male.rds")
USA_female <- readRDS("./datasets/USA/USA_female.rds")
all_male<-t(list.cbind(USA_male))
all_female<-t(list.cbind(USA_female))


# Plot mortality for the US

Male<-matrix(0,nrow = n_age,ncol=n_year)
Female<-matrix(0,nrow = n_age,ncol=n_year)
for (j in 1:n_year) {
  Male[,j]<-colMeans(all_male[seq(from=j,to=n_year*n_states,by=n_year),]) 
  Female[,j]<-colMeans(all_female[seq(from=j,to=n_year*n_states,by=n_year),])
}

savepdf(paste0("./plots2paper/", 
               "Fig_1a"))
colnames(Male)<-1:dim(Male)[2]
curves_male<-rainbow::fts(x=1:n_age, y=Male , xname='Age',yname='Log death rate')
plot.fds(curves_male,main="USA: male death rates (1959-2020)",ylim=c(-4,0))
dev.off()

savepdf(paste0("./plots2paper/", 
               "Fig_1b"))
colnames(Female)<-1:dim(Female)[2]
curves_female<-rainbow::fts(x=1:n_age, y=Female , xname='Age',yname='Log death rate')
plot.fds(curves_female,main="USA: female death rates (1959-2020)",ylim=c(-4,0))
dev.off()

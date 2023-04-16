# Two-way functional mean ANOVA (FM-ANOVA)
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

###############################################################################
#Functional grand effect based on means
###############################################################################
mu_hat<-function(data_pop1,data_pop2,n_year
                 ,n_prefectures,n_age,n_populations=2,row_par=part_list){
  all_data<-array(0,dim=c(n_prefectures,n_populations,n_year,n_age))
  Mu_hat<-rep(0,n_age)
  for (i in 1:n_prefectures) {
    all_data[i,1,,]=t(data_pop1[,(n_year*i-(n_year-1)):(n_year*i)])
    all_data[i,2,,]=t(data_pop2[,(n_year*i-(n_year-1)):(n_year*i)])
    for (j in 1:n_populations) {
      for (k in 1:n_year) {
        Mu_hat=Mu_hat+all_data[i,j,k,]
      }
    }
  }
  final_GE<-1/(n_year*n_populations*n_prefectures)*Mu_hat
  return(final_GE)
}

###############################################################################
#Functional row effect based on means
###############################################################################
FRE_means<-function(data_pop1=t(all_male),data_pop2=t(all_female),n_year
                    ,n_prefectures,n_age,n_populations=2,row_par=part_list){
  Mu_hat<-mu_hat(data_pop1=t(all_male),data_pop2=t(all_female),n_year
                 ,n_prefectures,n_age,n_populations=2,row_par=part_list)
  all_data<-array(0,dim=c(n_prefectures,n_populations,n_year,n_age))
  FRE_means<-matrix(0,nrow=n_age,ncol=n_prefectures)
  for (i in 1:n_prefectures) {
    alpha_hat<-rep(0,n_age)
    all_data[i,1,,]=t(data_pop1[,(n_year*i-(n_year-1)):(n_year*i)])
    all_data[i,2,,]=t(data_pop2[,(n_year*i-(n_year-1)):(n_year*i)])
    for (j in 1:n_populations) {
      for (k in 1:n_year) {
        alpha_hat=alpha_hat+all_data[i,j,k,]
      }
    }
    FRE_means[,i]<-(1/(n_year*2))*alpha_hat-Mu_hat
  }
  
  return(FRE_means=FRE_means)
  
}

###############################################################################
#Functional column effect based on means
###############################################################################
FCE_means<-function(data_pop1,data_pop2,n_year
                    ,n_prefectures,n_age,n_populations=2,row_par=part_list){
  Mu_hat<-mu_hat(data_pop1=t(all_male),data_pop2=t(all_female),n_year
                 ,n_prefectures,n_age,n_populations=2,row_par=part_list)
  all_data<-array(0,dim=c(n_prefectures,n_populations,n_year,n_age))
  FCE_means<-matrix(0,nrow=n_age,ncol=n_populations)
  for (j in 1:n_populations) {
    beta_hat<-rep(0,n_age)
    for (i in 1:n_prefectures) {
      all_data[i,1,,]=t(data_pop1[,(n_year*i-(n_year-1)):(n_year*i)])
      all_data[i,2,,]=t(data_pop2[,(n_year*i-(n_year-1)):(n_year*i)])
      for (k in 1:n_year) {
        beta_hat=beta_hat+all_data[i,j,k,]
      }
    }
    FCE_means[,j]<-(1/(n_year*n_prefectures))*beta_hat-Mu_hat
  }
  return(FCE_means=FCE_means)
}

###############################################################################
#Two-way Functional Anova
###############################################################################
#This is a general function for the Functional ANOVA based on means. It provides
#the functional grand effect, the functional row and column effects.
FANOVA<-function(data_pop1,data_pop2,n_year
                 ,n_prefectures,n_age,n_populations=2,row_par=part_list)
{
  
  #functional grand effect
  FGE_means<-mu_hat(data_pop1=t(all_male),data_pop2=t(all_female),n_year
                    ,n_prefectures,n_age,n_populations=2,row_par=part_list)
  
  #functional row effect
  FRE_mean<-FRE_means(data_pop1=t(all_male),data_pop2=t(all_female),n_year
                      ,n_prefectures,n_age,n_populations=2,row_par=part_list)
  
  #Functional column effect
  FCE_mean<-FCE_means(data_pop1=t(all_male),data_pop2=t(all_female),n_year
                      ,n_prefectures,n_age,n_populations=2,row_par=part_list)
  
  return(list(FGE_mean=FGE_means,FRE_mean=t(FRE_mean),FCE_mean=t(FCE_mean)))
}

###############################################################################
#Residuals obtained based on means
###############################################################################
Two_way_Residuals_means<-function(FANOVA_means,data_pop1,data_pop2
                                  ,n_prefectures,n_year,n_age){
  
  # FANOVA_means is the result after running the function FANOVA
  Two_FGE=FANOVA_means$FGE_mean
  Two_FRE=FANOVA_means$FRE_mean
  Two_FCE=FANOVA_means$FCE_mean
  
  all_male=t(data_pop1)
  all_female=t(data_pop2)
  
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
  residuals_b1r<-matrix(0,nrow=(n_prefectures*n_year),ncol = n_age)
  residuals_b2r<-matrix(0,nrow=(n_prefectures*n_year),ncol = n_age)
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
  
  return(list(residuals1_mean= residuals_b1r,residuals2_mean=residuals_b2r,rd_mean=RD,R_mean=R,Fixed_comp_mean=FD))
  
}



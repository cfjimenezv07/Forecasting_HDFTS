# Code to reproduce the Figure 5  in the paper. 

# All the data required to reproduce this figure can be recomputed from the codes 
# FMP_ANOVA/FM_ANOVA. Comparison with VAR method when using K=6 PC, can be obtained from the same R codes
# by simply changing ARIMA to VAR in the point forecast section (Section 6.)


################################################################
# Set a working directory
################################################################
setwd("~/My Drive/Fall 2023/STAT 397/PhD project 4/Revisions from JCGS/Final_Rcodes_editor")
# All plots from this code will be saved in the folder plots2paper in the working directory
rm(list = ls())
library(tidyverse)
library(ggplot2)
################################################################

# Function to save the plots.
savefig <- function(filename, height = 10, width = (1 + sqrt(5))/2*height, dpi = 300) {
  ggsave(filename, height = height/2.54, width = width/2.54, dpi = dpi)
}


################################################################
# Load the datasets
################################################################
#USA
Errors_FMP_cov_ARIMA_USA <- readRDS("./Results_Figure_A1/Errors_FMP_cov_ARIMA_USA_3.rds")
Errors_FMP_cov_VAR_USA <- readRDS("./Results_Figure_A1/Errors_FMP_cov_VAR_USA.rds")

Errors_FM_cov_ARIMA_USA <- readRDS("./Results_Figure_A1/Errors_FM_cov_ARIMA_USA_3.rds")
Errors_FM_cov_VAR_USA <- readRDS("./Results_Figure_A1/Errors_FM_cov_VAR_USA.rds")


#####################################################################################################################################################

#France
Errors_FMP_cov_ARIMA_France <- readRDS("./Results_Figure_A1/Errors_FMP_cov_ARIMA_France_3.rds")
Errors_FMP_cov_VAR_France <- readRDS("./Results_Figure_A1/Errors_FMP_cov_VAR_France.rds")

Errors_FM_cov_ARIMA_France <- readRDS("./Results_Figure_A1/Errors_FM_cov_ARIMA_France_3.rds")
Errors_FM_cov_VAR_France <- readRDS("./Results_Figure_A1/Errors_FM_cov_VAR_France.rds")


# Japan
Errors_FMP_cov_ARIMA_Japan <- readRDS("./Results_Figure_A1/Errors_FMP_cov_ARIMA_Japan_3.rds")
Errors_FMP_cov_VAR_Japan <- readRDS("./Results_Figure_A1/Errors_FMP_cov_VAR_Japan.rds")

Errors_FM_cov_ARIMA_Japan <- readRDS("./Results_Figure_A1/Errors_FM_cov_ARIMA_Japan_3.rds")
Errors_FM_cov_VAR_Japan <- readRDS("./Results_Figure_A1/Errors_FMP_cov_VAR_Japan.rds")

to_long <- function(x, male = T, method) {
  temp <- rbind.data.frame(as.matrix(x[, 1]), as.matrix(x[, -1]))
  temp$gender <- ifelse(male, "Male", "Female")
  temp$metric <- rep(c("MAPE", "RMSPE"), each = nrow(temp)/2)
  temp$method <- method
  return(temp)
}

df_toplot <- function(x, country, method){
  female <- 2*(1:2)
  long_male <- to_long(x[, -female], method = method)
  long_female <- to_long(x[, female], F, method)
  out <- rbind.data.frame(long_male, long_female)
  out$country <- country
  out
}



france.df <- list(Errors_FMP_cov_ARIMA_France,Errors_FMP_cov_VAR_France, Errors_FM_cov_ARIMA_France,Errors_FM_cov_VAR_France)
france.country <- as.list(rep("France",4))
france.methods <- list("FMP ARIMA","FMP VAR", "FM ARIMA", "FM VAR")
france <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, france.df, france.country, france.methods, SIMPLIFY = F))




us.df <- list(Errors_FMP_cov_ARIMA_USA,Errors_FMP_cov_VAR_USA, Errors_FM_cov_ARIMA_USA,Errors_FM_cov_VAR_USA)
us.country <- as.list(rep("USA", 4))
us.methods <- list("FMP ARIMA","FMP VAR", "FM ARIMA", "FM VAR")
us <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, us.df, us.country, us.methods, SIMPLIFY = F))




japan.df <- list(Errors_FMP_cov_ARIMA_Japan,Errors_FMP_cov_VAR_Japan, Errors_FM_cov_ARIMA_Japan,Errors_FM_cov_VAR_Japan)
japan.country <- as.list(rep("Japan", 4))
japan.methods <- list("FMP ARIMA","FMP VAR", "FM ARIMA", "FM VAR")
japan <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, japan.df, japan.country, japan.methods, SIMPLIFY = F))



gen_plot <- function(df, legend = T) {
  p <-df %>% 
    mutate(method = factor(method, levels = c("FMP ARIMA","FMP VAR", "FM ARIMA", "FM VAR"), ordered = T),
           metric = factor(metric, levels = c("RMSPE", "MAPE"), ordered = T)
    ) %>% 
    ggplot(aes(x = method, y = V1, fill = gender)) +
    geom_boxplot(show.legend = legend, outlier.shape = NA) +
    labs(x = "", y = "") +
    theme(legend.position = c(0.1,0.85))
  co <- 1
  for (metric in c("RMSPE", "MAPE")) {
    (p +
       ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = co))
    # savefig(paste0("~/My Drive/Fall 2023/STAT 397/PhD project 4/Revisions from JCGS/plots2paper/",filename, metric, ".pdf"))
    co <- co + 1
  }
  return(p)
}


us.metrics <- gen_plot(us)
france.metrics <- gen_plot(france, F)
japan.metrics <- gen_plot(japan, F)


# plots USA
us.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 5.5) 
savefig(paste0("./plots2paper/", 
               "Fig_5a",  ".pdf"))

us.metrics <- gen_plot(us,F)
us.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 5.5) 

savefig(paste0("./plots2paper/", 
               "Fig_5b",  ".pdf"))


#Deleting the top caption from France and Japan

#France
france.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 60) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_5c",  ".pdf"))

france.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 60) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_5d",  ".pdf"))


#Japan
japan.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 75) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_5e",  ".pdf"))


japan.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 75) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_5f",  ".pdf"))


rm(list = ls())

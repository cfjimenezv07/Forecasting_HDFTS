# Code to reproduce the Figure 6  in the paper. 

# All the data required to reproduce this figure can be recomputed from the codes 
# FMP_ANOVA/FM_ANOVA. Computation of the forecast for independent case can be 
# obtained in the Rcode Independence.R 


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

us.means <- readRDS("./Results_Figure_A2/Errors_FMP_cov_ARIMA_USA_4.rds")
us.basedmeans <- readRDS("./Results_Figure_A2/Errors_FM_cov_ARIMA_USA_4.rds")
us.naive <- readRDS("./Results_Figure_A2/Errors_mean_naive_USA.rds")


france.means <- readRDS("./Results_Figure_A2/Errors_FMP_cov_ARIMA_France_4.rds")
france.basedmeans <- readRDS("./Results_Figure_A2/Errors_FM_cov_ARIMA_France_4.rds")
france.naive <- readRDS("./Results_Figure_A2Errors_mean_naive_France.rds")


japan.means <- readRDS("./Results_Figure_A2/Errors_FMP_cov_ARIMA_Japan_4.rds")
japan.basedmeans <- readRDS("./Results_Figure_A2/Errors_FM_cov_ARIMA_Japan_4.rds")
japan.naive <- readRDS("./Results_Figure_A2/Errors_mean_naive_japan.rds")

column_name_fixer <- function(x) {
  str_remove(x, "(GSY|HNT)_")
}


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



france.df <- list(france.means, france.basedmeans,  france.naive)
france.country <- as.list(rep("France", 3))
france.methods <- list("FMP-ANOVA", "FM-ANOVA", "Independent")
france <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, france.df, france.country, france.methods, SIMPLIFY = F))




us.df <- list(us.means, us.basedmeans, us.naive)
us.country <- as.list(rep("USA", 3))
us.methods <- list("FMP-ANOVA", "FM-ANOVA", "Independent")
us <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, us.df, us.country, us.methods, SIMPLIFY = F))




japan.df <- list(japan.means, japan.basedmeans,  japan.naive)
japan.country <- as.list(rep("Japan", 3))
japan.methods <- list("FMP-ANOVA", "FM-ANOVA", "Independent")
japan <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, japan.df, japan.country, japan.methods, SIMPLIFY = F))



gen_plot <- function(df,  legend = T) {
  p <-df %>% 
    mutate(method = factor(method, levels = c("FMP-ANOVA", "FM-ANOVA",
                                              "Independent"), ordered = T),
           metric = factor(metric, levels = c("RMSPE", "MAPE"), ordered = T)
    ) %>% 
    ggplot(aes(x = method, y = V1, fill = gender)) +
    geom_boxplot(show.legend = legend, outlier.shape = NA) +
    labs(x = "", y = "")+
    theme(legend.position = c(0.10,0.85))
  co <- 1
  for (metric in c("RMSPE", "MAPE")) {
    (p +
       ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = co))
    # savefig(paste0("./Plots/",filename, metric, ".pdf"))
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
  ylim(0, 25) 
savefig(paste0("./plots2paper/", 
               "Fig_6a",  ".pdf"))

us.metrics <- gen_plot(us,F)
us.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 25) 

savefig(paste0("./plots2paper/", 
               "Fig_6b",  ".pdf"))


#Deleting the top caption from France and Japan

#France
france.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 35) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_6c",  ".pdf"))

france.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 35) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_6d",  ".pdf"))


#Japan
japan.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 35) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_6e",  ".pdf"))


japan.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 24) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_6f",  ".pdf"))


rm(list = ls())
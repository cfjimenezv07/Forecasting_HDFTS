# Code to reproduce the Figure 2 in the paper. Based on EVR and static FPCA

# All the data required to reproduce this figure can be recomputed from the codes FMP_ANOVA
# and FM_ANOVA, The results for HNT and GSY methods can be obtained from HNT23_GSY19.R

################################################################
# Set a working directory
################################################################
setwd("~/My Drive/Fall 2023/STAT 397/PhD project 4/Revisions from JCGS/Final_Rcodes_editor")
# All plots from this code will be saved in the folder plots2paper in the working directory
rm(list = ls())
library(tidyverse)

################################################################

# Function to save the plots.
savefig <- function(filename, height = 10, width = (1 + sqrt(5))/2*height, dpi = 300) {
  ggsave(filename, height = height/2.54, width = width/2.54, dpi = dpi)
}
################################################################
# Load the datasets
################################################################

us.means <- readRDS("./Results_Figure_2/Errors_FMP_cov_ARIMA_USA_4.rds")
us.basedmeans <- readRDS("./Results_Figure_2/Errors_FM_cov_ARIMA_USA_4.rds")
us.gao <- readRDS("./Results_Figure_2/Errors_GAO_USA.rds")
us.tnh <- readRDS("./Results_Figure_2/Errors_TNH_USA.rds")

france.means <- readRDS("./Results_Figure_2/Errors_FMP_cov_ARIMA_France_4.rds")
france.basedmeans <- readRDS("./Results_Figure_2/Errors_FM_cov_ARIMA_France_4.rds")
france.gao <- readRDS("./Results_Figure_2/Errors_GAO_France.rds")
france.tnh <- readRDS("./Results_Figure_2/Errors_TNH_France.rds")

japan.means <- readRDS("./Results_Figure_2/Errors_FMP_cov_ARIMA_Japan_4.rds")
japan.basedmeans <- readRDS("./Results_Figure_2/Errors_FM_cov_ARIMA_Japan_4.rds")
japan.gao <- readRDS("./Results_Figure_2/Errors_GAO_Japan.rds")
japan.tnh <- readRDS("./Results_Figure_2/Errors_TNH_Japan.rds")


column_name_fixer <- function(x) {
  str_remove(x, "(GSY|HNT)_")
}


savefig <- function(filename, height = 10, width = (1 + sqrt(5))/2*height, dpi = 300) {
  ggsave(filename, height = height/2.54, width = width/2.54, dpi = dpi)
}

colnames(france.gao) <- sapply(colnames(france.gao), column_name_fixer)

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

#standardizing column names
colnames(france.gao) <- sapply(colnames(france.gao), column_name_fixer)
colnames(france.tnh) <- sapply(colnames(france.tnh), column_name_fixer)

france.df <- list(france.means, france.basedmeans,  france.gao, france.tnh)
france.country <- as.list(rep("France", 4))
france.methods <- list("FMP-ANOVA", "FM-ANOVA", "GSY", "HNT")
france <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, france.df, france.country, france.methods, SIMPLIFY = F))

#standardizing column names
colnames(us.gao) <- sapply(colnames(us.gao), column_name_fixer)
colnames(us.tnh) <- sapply(colnames(us.tnh), column_name_fixer)


us.df <- list(us.means, us.basedmeans,  us.gao, us.tnh)
us.country <- as.list(rep("USA", 4))
us.methods <- list("FMP-ANOVA", "FM-ANOVA",  "GSY", "HNT")
us <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, us.df, us.country, us.methods, SIMPLIFY = F))


#standardizing column names
colnames(japan.gao) <- sapply(colnames(japan.gao), column_name_fixer)
colnames(japan.tnh) <- sapply(colnames(japan.tnh), column_name_fixer)

japan.df <- list(japan.means, japan.basedmeans,  japan.gao, japan.tnh)
japan.country <- as.list(rep("Japan", 4))
japan.methods <- list("FMP-ANOVA", "FM-ANOVA", "GSY", "HNT")
japan <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, japan.df, japan.country, japan.methods, SIMPLIFY = F))



gen_plot <- function(df,  legend = T) {
  p <-df %>% 
    mutate(method = factor(method, levels = c("FMP-ANOVA", "FM-ANOVA", "HNT",
                                              "GSY"), ordered = T),
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
  ylim(0, 10) 
savefig(paste0("./plots2paper/", 
               "Fig_2a",  ".pdf"))

us.metrics <- gen_plot(us,F)
us.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 10) 

savefig(paste0("./plots2paper/", 
               "Fig_2b",  ".pdf"))


#Deleting the top caption from France and Japan

#France
france.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 10) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_2c",  ".pdf"))

france.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 10) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_2d",  ".pdf"))


#Japan
japan.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 10) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_2e",  ".pdf"))


japan.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 10) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_2f",  ".pdf"))


rm(list = ls())

# Code to reproduce the Figure 3 in the paper. 

# All the data required to reproduce this figure can be recomputed from the codes FMP_ANOVA


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
us.base <- readRDS("./Results_Figure_3/USA/Emp_cov_USA.rds")
france.base <- readRDS("./Results_Figure_3/France/Emp_cov_France.rds")
japan.base <- readRDS("./Results_Figure_3/Japan/Emp_cov_Japan.rds")

df_toplot <- function(x, lb, ub, type, country) {
  male <- rbind.data.frame(as.matrix(x[, lb]), as.matrix(x[, ub]))
  male$gender <- "Male"
  male$type <- rep(type, each = nrow(x))
  female <- rbind.data.frame(as.matrix(x[, lb+6]), as.matrix(x[, ub+6]))
  female$gender <- "Female"
  female$type <- rep(type, each = nrow(x))
  out <- rbind.data.frame(male, female)
  out$country <- country
  return(out)
} 

list_maker <- function(x) {as.list(rep(x, 3))}

dfs <- list(us.base, france.base, japan.base)
country <- list("USA", "France", "Japan")
types <- list(c("80%", "95%"), c("80%", "95%"), c("80%", "95%"))

pw <- mapply(df_toplot, x = dfs, lb = list_maker(1), ub = list_maker(2), 
             type = types, country = country, SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

is <- mapply(df_toplot, x = dfs, lb = list_maker(3), ub = list_maker(4), 
             type = types, country = country, SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

cpd <- mapply(df_toplot, x = dfs, lb = list_maker(5), ub = list_maker(6), 
              type = types, country = country, SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

gen_plot <- function(df, legend = T) {
  p <- df %>% 
    mutate(country = factor(country, levels = c("USA", "France", "Japan"))) %>% 
    ggplot(aes(country, V1, fill = type)) +
    geom_boxplot(show.legend = legend) +
    scale_fill_manual(values = c("darkblue", "darkgreen"))+
    labs(x = "", y = "")+
    theme(legend.position = c(0.1,0.8))
  co <- 1
  for (gender in c("Female", "Male")) {
    (p +
       ggforce::facet_wrap_paginate(~gender, nrow = 1, ncol = 1, page = co))
    # if (legend) {
    #   temp <- temp +
    #     ggforce::facet_wrap_paginate(~gender, nrow = 1, ncol = 1, page = co) +
    #     theme(strip.background = element_blank(), strip.text = element_blank())
    # }
    # ggsave(paste0("./Plot/",filename, gender, ".eps"))
    co <- co + 1
  }
  return(p)
}

pw.plot <- gen_plot(pw)
is.plot <- gen_plot(is, F)
cpd.plot <- gen_plot(cpd, F)

# plots PW
pw.plot +
  ggforce::facet_wrap_paginate(~gender, ncol = 1, nrow = 1, page = 1) +
  ylim(0.6, 1.2) 

savefig(paste0("./plots2paper/", 
               "Fig_3a",  ".pdf"))

pw.plot <- gen_plot(pw,F)
pw.plot +
  ggforce::facet_wrap_paginate(~gender, ncol = 1, nrow = 1, page = 2) +
  ylim(0.6, 1.2)  

savefig(paste0("./plots2paper/", 
               "Fig_3b",  ".pdf"))


#Deleting the top caption from is.plot and CPD

#cpd.plot
cpd.plot +
  ggforce::facet_wrap_paginate(~gender, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 0.25) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_3c",  ".pdf"))


cpd.plot +
  ggforce::facet_wrap_paginate(~gender, ncol = 1, nrow = 1, page = 2) +
  ylim(0,0.25) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_3d",  ".pdf"))
#is.plot
is.plot +
  ggforce::facet_wrap_paginate(~gender, ncol = 1, nrow = 1, page = 1,) +
  ylim(0.2, 1.6) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_3e",  ".pdf"))

is.plot +
  ggforce::facet_wrap_paginate(~gender, ncol = 1, nrow = 1, page = 2) +
  ylim(0.2, 1.6) +
  theme(strip.background = element_blank(), strip.text = element_blank())

savefig(paste0("./plots2paper/", 
               "Fig_3f",  ".pdf"))


rm(list = ls())


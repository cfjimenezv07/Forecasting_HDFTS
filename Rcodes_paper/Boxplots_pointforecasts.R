# This code produces the plots of Figure 2.

library(tidyverse)

france.means <- readRDS("./results_for_the plots/France/Errors_mean_1_France.rds")
france.basedmeans <- readRDS("./results_for_the plots/France/Errors_mean_basedmeans_1_France.rds") 
france.naive <- readRDS("./results_for_the plots/France/Errors_mean_naive_France.rds")

us.means <- readRDS("./results_for_the plots/USA/Errors_mean_1_USA.rds")
us.basedmeans <- readRDS("./results_for_the plots/USA/Errors_mean_basedmeans_1_USA.rds")
us.naive <- readRDS("./results_for_the plots/USA/Errors_mean_naive_USA.rds")

japan.means <- readRDS("./results_for_the plots/Japan/Errors_mean_1_Japan.rds")
japan.basedmeans <- readRDS("./results_for_the plots/Japan/Errors_mean_basedmeans_1_Japan.rds")
japan.naive <- readRDS("./results_for_the plots/Japan/Errors_mean_naive_japan.rds")




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

france.df <- list(france.means, france.basedmeans, france.naive)
france.country <- as.list(rep("France", 3))
france.methods <- list("FMP-ANOVA", "FM-ANOVA", "Independence")
france <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, france.df, france.country, france.methods, SIMPLIFY = F))

us.df <- list(us.means, us.basedmeans, us.naive)
us.country <- as.list(rep("USA", 3))
us.methods <- list("FMP-ANOVA", "FM-ANOVA", "Independence")
us <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, us.df, us.country, us.methods, SIMPLIFY = F))

japan.df <- list(japan.means, japan.basedmeans, japan.naive)
japan.country <- as.list(rep("Japan", 3))
japan.methods <- list("FMP-ANOVA", "FM-ANOVA", "Independence")
japan <- 
  do.call(rbind.data.frame,
          mapply(df_toplot, japan.df, japan.country, japan.methods, SIMPLIFY = F))



gen_plot <- function(df, filename) {
  p <-df %>% 
    mutate(method = factor(method, levels = c("FMP-ANOVA", "FM-ANOVA", "Independence"), ordered = T),
           metric = factor(metric, levels = c("RMSPE", "MAPE"), ordered = T)
    ) %>% 
    ggplot(aes(x = method, y = V1, fill = gender)) +
    geom_boxplot(show.legend = F, outlier.shape = NA) +
    labs(x = "", y = "")+
    theme(legend.position = c(0.10,0.85))
  co <- 1
  for (metric in c("RMSPE", "MAPE")) {
    (p +
      ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = co))
    ggsave(paste0("./Plots/",filename, metric, ".eps"))
    co <- co + 1
  }
  return(p)
}


us.metrics <- gen_plot(us, "PF_comparisons_USA")
france.metrics <- gen_plot(france, "PF_comparisons_France")
japan.metrics <- gen_plot(japan, "PF_comparisons_Japan")


#Deleting the top caption from France and Japan

#France
france.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  ylim(0, 45) +
  theme(strip.background = element_blank(), strip.text = element_blank())

ggsave(paste0("./Plots/", 
              "PF_comparisons_France", "RMSPE", ".eps"))

france.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  ylim(0, 45) +
  theme(strip.background = element_blank(), strip.text = element_blank())

ggsave(paste0("./Plots/", 
              "PF_comparisons_France", "MAPE", ".eps"))

#Japan
japan.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 1,) +
  theme(strip.background = element_blank(), strip.text = element_blank())

ggsave(paste0("./Plots/", 
              "PF_comparisons_Japan", "RMSPE", ".eps"))

japan.metrics +
  ggforce::facet_wrap_paginate(~metric, ncol = 1, nrow = 1, page = 2) +
  theme(strip.background = element_blank(), strip.text = element_blank())

ggsave(paste0("./Plots/", 
              "PF_comparisons_Japan", "MAPE", ".eps"))


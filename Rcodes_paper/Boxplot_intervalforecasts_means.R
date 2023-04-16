#This R code produces the plots of the Figure 4. For the computations of interval 
#coverages based on  FM-ANOVA approach

library(ggplot2)
library(tidyverse)

us.fanova <- readRDS("./results_for the plots/USA/Emp_cov_USA_means.rds") 
france.fanova <- readRDS("./results_for the plots/France/Emp_cov_France_means.rds")
japan.fanova <- readRDS("./results_for the plots/Japan/Emp_cov_Japan_means.rds")

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

dfs <- list(us.base.means, france.base.means, japan.base.means)
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

gen_plot <- function(df, filename, remove = F) {
  p <- df %>% 
    mutate(country = factor(country, levels = c("USA", "France", "Japan"))) %>% 
    ggplot(aes(country, V1, fill = type)) +
    geom_boxplot(show.legend = F) +
    labs(x = "", y = "") +
    scale_fill_manual(values = c("darkblue", "darkgreen"))
  co <- 1
  for (gender in c("FEMALE", "MALE")) {
    temp <- p +
      ggforce::facet_wrap_paginate(~gender, nrow = 1, ncol = 1, page = co)
    if (remove) {
      temp <- temp +
        ggforce::facet_wrap_paginate(~gender, nrow = 1, ncol = 1, page = co) +
        theme(strip.background = element_blank(), strip.text = element_blank())
    }
    ggsave(paste0("./Plots/",filename, gender, ".eps"))
    co <- co + 1
  }
  return(p)
}

pw.plot <- gen_plot(pw, "PointWiseMeans")

is.plot <- gen_plot(is, "IntervalScoreMeans", T)

cpd.plot <- gen_plot(cpd, "CPDMeans", T)

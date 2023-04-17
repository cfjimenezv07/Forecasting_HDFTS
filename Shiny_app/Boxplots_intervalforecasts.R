library(ggplot2)
library(tidyverse)

us.base <- readRDS("USA/Emp_cov_USA.rds")
france.base <- readRDS("France/Emp_cov_France.rds")
japan.base <- readRDS("Japan/Emp_cov_Japan.rds")

us.fanova <- readRDS("USA/Emp_cov_USA_means.rds") 
france.fanova <- readRDS("France/Emp_cov_France_means.rds")
japan.fanova <- readRDS("Japan/Emp_cov_Japan_means.rds")

df_toplot <- function(x, lb, ub, type, country, method) {
  male <- rbind.data.frame(as.matrix(x[, lb]), as.matrix(x[, ub]))
  male$gender <- "Male"
  male$type <- rep(type, each = nrow(x))
  female <- rbind.data.frame(as.matrix(x[, lb+6]), as.matrix(x[, ub+6]))
  female$gender <- "Female"
  female$type <- rep(type, each = nrow(x))
  out <- rbind.data.frame(male, female)
  out$country <- country
  out$state <- gsub("\\d+", "", rownames(out))
  out$method <- method
  return(out)
} 

list_maker <- function(x) {as.list(rep(x, 3))}

dfs.fmp <- list(us.base, france.base, japan.base)
dfs.fanova <- list(us.fanova, france.fanova, japan.fanova)
country <- list("USA", "France", "Japan")
types <- list(c("80%", "95%"), c("80%", "95%"), c("80%", "95%"))

pw.fmp <- mapply(df_toplot, x = dfs.fmp, lb = list_maker(1), ub = list_maker(2), 
             type = types, country = country, method = list_maker("FMP"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

pw.fanova <- mapply(df_toplot, x = dfs.fanova, lb = list_maker(1), ub = list_maker(2), 
                 type = types, country = country, method = list_maker("FANOVA"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

pw <- rbind.data.frame(pw.fmp, pw.fanova)


is.fmp <- mapply(df_toplot, x = dfs.fmp, lb = list_maker(3), ub = list_maker(4), 
                 type = types, country = country, method = list_maker("FMP"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

is.fanova <- mapply(df_toplot, x = dfs.fanova, lb = list_maker(3), ub = list_maker(4), 
                 type = types, country = country, method = list_maker("FANOVA"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

is <- rbind.data.frame(is.fmp, is.fanova)

cpd.fmp <- mapply(df_toplot, x = dfs.fmp, lb = list_maker(5), ub = list_maker(6), 
                 type = types, country = country, method = list_maker("FMP"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

cpd.fanova <- mapply(df_toplot, x = dfs.fanova, lb = list_maker(5), ub = list_maker(6), 
                    type = types, country = country, method = list_maker("FANOVA"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

cpd <- rbind.data.frame(cpd.fmp, cpd.fanova)


saveRDS(pw, "intervaldfs/pw.Rds")
saveRDS(is, "intervaldfs/is.Rds")
saveRDS(cpd, "intervaldfs/cpd.Rds")

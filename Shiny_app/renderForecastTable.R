#point forecast errors per country per method per gender
us.errors <- readRDS("tabledfs/us.Rds")
france.errors <- readRDS("tabledfs/france.Rds")
japan.errors <- readRDS("tabledfs/japan.Rds")
australia.errors <- readRDS("tabledfs/australia.Rds")

dflist <- list(USA = us.errors, France = france.errors, 
               Japan = japan.errors, Australia = australia.errors)

point.forecast.table <- function(df, .gender, .metric, .state, ltable = T) {
  df <- dflist[[df]]
  temp <- df %>% 
    filter(gender == .gender, metric == .metric)
  l.df <- temp %>% 
    filter(method == "FMP") %>% 
    select(.dots = .state) %>% 
    rename(FMP = .dots) 
    if(ltable) {
      l.df <- l.df %>% 
        mutate(h = 1:nrow(.)) %>% 
        select(h, everything())
    }
  r.df <- temp %>% 
    filter(method != "FMP") %>% 
    select(.dots = .state) %>% 
    rename(FANOVA = .dots)
  out <- cbind(
    l.df,
    r.df
  )
  colnames(out) <- paste(.metric, colnames(out))
  out
}

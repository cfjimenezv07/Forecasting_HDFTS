#point forecast errors per country per method per gender


point.forecast.table <- function(df, .gender, .metric, .state, type, ltable = T) {
  
  type <- ifelse(type == "EVR", 4, 3)
  us.errors <- readRDS(paste0("tabledfs/us", type, ".Rds"))
  france.errors <- readRDS(paste0("tabledfs/france", type, ".Rds"))
  japan.errors <- readRDS(paste0("tabledfs/japan", type, ".Rds"))

  dflist <- list(USA = us.errors, France = france.errors, Japan = japan.errors)

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
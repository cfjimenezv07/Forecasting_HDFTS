pw <- readRDS("intervaldfs/pw.Rds")
is <- readRDS("intervaldfs/is.Rds")
cpd <- readRDS("intervaldfs/cpd.Rds")

df.list <- list(pw, cpd, is)
list_maker <- function(x) {as.list(rep(x, 3))}
kind.list <- list("Pointwise coverage probability",
                  "CPD", "Interval score")

individualMeasurement <- function(df, .state, .gender, kind) {
  temp <- df %>% 
    filter(state == .state, gender == .gender) %>% 
    mutate(type = paste(method, type)) %>% 
    select(V1, type) %>% 
    pivot_wider(names_from = type, values_from = V1) %>% 
    as.data.frame()
  rownames(temp) <- kind
  temp
}

measurementsTableGenerator <- function(state, gender) {
  temp <- mapply(individualMeasurement, df = df.list, .state = list_maker(state),
                 .gender = list_maker(gender), kind = kind.list, 
                  SIMPLIFY = F) %>% 
    do.call(rbind.data.frame, .)
  round(temp, 4)
}


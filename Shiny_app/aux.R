#turning messy lists into treatable dataframes

gen.lists <- function(idx) {
  us.fmp <<- readRDS(paste0("USA/All_errors_FMP_cov_ARIMA_USA_", idx, ".rds"))
  us.fm <<- readRDS(paste0("USA/All_errors_FM_cov_ARIMA_USA_", idx, ".rds")) 
  us.names <<- readRDS("names/names_states.rds")

  france.fmp <<- readRDS(paste0("France/All_errors_FMP_cov_ARIMA_France_", idx, ".rds"))
  france.fm <<- readRDS(paste0("France/All_errors_FM_cov_ARIMA_France_", idx, ".rds"))
  france.names <<- readRDS("names/names_departments.rds")

  japan.fmp <<- readRDS(paste0("Japan/All_errors_FMP_cov_ARIMA_Japan_", idx, ".rds"))
  japan.fm <<- readRDS(paste0("Japan/All_errors_FM_cov_ARIMA_Japan_", idx, ".rds"))
  japan.names <<- readRDS("names/names_prefectures.rds")
}

join.dataframe <- function(lst, cnames, method) {
  cnames <- sapply(cnames, function(x) x)
  temp <- do.call(rbind.data.frame, lst)
  colnames(temp) <- cnames
  err.type <- matrix(rep(c("MAPE", "RMSPE"), each = 10), ncol = 1)
  temp$metric <- as.character(rbind(err.type, err.type))
  temp$gender <- rep(c("Male", "Female"), each = 20)
  temp$method <- method
  temp
}

join.methods <- function(fmp, fm, cnames) {
  fmp <- join.dataframe(fmp, cnames, "FMP")
  fm <- join.dataframe(fm, cnames, "FM")
  rbind.data.frame(fmp, fm)
}

gen.data <- function(idx) {
  gen.lists(idx)
  us.full <- join.methods(us.fmp, us.fm, us.names)
  saveRDS(us.full, paste0("tabledfs/us", idx, ".Rds"))
  france.full <- join.methods(france.fmp, france.fm, france.names)
  saveRDS(france.full, paste0("tabledfs/france", idx, ".Rds"))
  japan.full <- join.methods(japan.fmp, japan.fm, japan.names)
  saveRDS(japan.full, paste0("tabledfs/japan", idx, ".Rds"))
}

gen.data(3)
gen.data(4)

#rainbow plots
# us.male <- readRDS("USA/USA_male.rds")
# alabama <- us.male[["Alabama"]]
# us.male.forecasted <- readRDS("USA/forecasted_FM_cov_ARIMA_USA_3.rds")
# us.male.forecasted.fmp <- readRDS("USA/forecasted_FMP_cov_ARIMA_USA_3.rds")
# alabama.forecasts <- us.male.forecasted[[1]][[1]][1:101, ]
# alabama.forecasts.fmp <- us.male.forecasted.fmp[[1]][[1]][1:101, ]


# m <- matrix(c(1, 1, 2, 2), nrow = 2, byrow = T)
# layout(m)
# matplot(alabama[, 1:52],type='l',col='gray86'
#                   ,main=paste("State","Alabama"),xlab = 'Age'
#                   ,ylab='Log mortality rate')

# colnames(alabama.forecasts)<-1:10
# Res_forcasted_curves_male<-rainbow::fts(x=0:100, y=alabama.forecasts, 
#                                         xname='Age',yname='USA Mr')
# rainbow::lines.fds(Res_forcasted_curves_male)

# matplot(alabama[, 1:52],type='l',col='gray86'
#                   ,main=paste("State","Alabama"),xlab = 'Age'
#                   ,ylab='Log mortality rate')

# colnames(alabama.forecasts.fmp)<-1:10
# Res_forcasted_curves_male.fmp<-rainbow::fts(x=0:100, y=alabama.forecasts.fmp, 
#                                         xname='Age',yname='USA Mr')
# rainbow::lines.fds(Res_forcasted_curves_male.fmp)

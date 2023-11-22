us.male <- readRDS("USA/USA_male.rds")
us.female <- readRDS("USA/USA_female.rds")
us.names.vec <- sapply(readRDS("names/names_states.rds"), 
                   function(x) x)
# us.forecasted <- readRDS("USA/forecasted_curves_triangular_USA.rds")

france.male <- readRDS("France/France_male.rds")
france.female <- readRDS("France/France_female.rds")
france.names.vec <- sapply(readRDS("names/names_departments.rds"), 
                       function(x) x)
#france.forecasted <- readRDS("France/forecasted_curves_triangular_France.rds")

japan.male <- readRDS("Japan/Japan_male.rds")
japan.female <- readRDS("Japan/Japan_female.rds")
japan.names.vec <- sapply(readRDS("names/names_prefectures.rds"), 
                      function(x) x)
#japan.forecasted <- readRDS("Japan/forecasted_curves_triangular_Japan.rds")

name.mapper <- function(x, name) {
  which.max(x == name)
}

male.base <- list(USA = us.male, France = france.male, Japan = japan.male)
female.base <- list(USA = us.female, France = france.female, Japan = japan.female)
gender.base <- list(Male = male.base, Female = female.base)
# forecast.list <- list(USA = us.forecasted, France = france.forecasted, 
#                       Japan = japan.forecasted)
names.list <- list(USA = us.names.vec, France = france.names.vec, 
                   Japan = japan.names.vec)
limits.gender <- list(USA = 1:101, France = 1:101, Japan = 1:99)
limits.forecast <- list(USA = 1:52, France = 1:44, Japan = 1:36)
titles <- list(Male = list(USA = "USA: male death rates (1959-2020)",
                   France = "France: male death rates (1968-2021)",
                   Japan = "Japan: male death rates (1975-2020)"),
               Female = list(USA = "USA: female death rates (1959-2020)",
                             France = "France: female death rates (1968-2021)",
                             Japan = "Japan: female death rates (1975-2020)"))

generate_forecast_list <- function(type, forecast) {
  type <- ifelse(type == "EVR", 4, 3)
  us.forecasted <- readRDS(paste0("USA/forecasted_", forecast, "_cov_ARIMA_USA_", type, ".rds"))
  france.forecasted <- readRDS(paste0("France/forecasted_", forecast, "_cov_ARIMA_France_", type, ".rds"))
  japan.forecasted <- readRDS(paste0("Japan/forecasted_", forecast, "_cov_ARIMA_Japan_", type, ".rds"))
  list(USA = us.forecasted, France = france.forecasted, 
                      Japan = japan.forecasted)
}

rainbow.generator <- function(country, state, type, forecast, gender = 'Male', plot.title = "", xlabel = F) {
  gender.limits <- limits.gender[[country]]
  forecast.limits <- limits.forecast[[country]]
  country.data <- gender.base[[gender]][[country]]
  state.data <- country.data[[state]]
  gender.factor <- ifelse(gender == "Male", 1, -1)
  forecast.idx <- name.mapper(names.list[[country]], state) 
  forecast.list <- generate_forecast_list(type, forecast)
  forecast.data <- forecast.list[[country]][[forecast.idx]][[1]][gender.factor*gender.limits, ]
  matplot(state.data[, forecast.limits],type='l',col='gray86', 
          xlab = ifelse(xlabel, 'Age', ""), ylab='Log death rate',
          main = paste0(titles[[gender]][[country]], "\n", plot.title))
  
  colnames(forecast.data)<-1:10
  Res_forcasted_curves<-rainbow::fts(x=gender.limits - 1, y=forecast.data, 
                                          xname='Age',yname='')
  rainbow::lines.fds(Res_forcasted_curves)
  
}


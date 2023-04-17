#turning messy lists into treatable dataframes

us.fmp <- readRDS("USA/PFE_USA.rds")
us.fanova <- readRDS("USA/PFE_USA_basedmeans.rds")
us.names <- readRDS("names/names_states.rds")

france.fmp <- readRDS("France/PFE_France.rds")
france.fanova <- readRDS("France/PFE_France_means.rds")
france.names <- readRDS("names/names_departments.rds")

japan.fmp <- readRDS("Japan/PFE_Japan.rds")
japan.fanova <- readRDS("Japan/PFE_Japan_means.rds")
japan.names <- readRDS("names/names_prefectures.rds")

australia.fmp <- readRDS("Australia/PFE_AUS.rds")
australia.fanova <- readRDS("Australia/PFE_AUS_means.rds")
australia.names <- readRDS("names/names_regions.rds")

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

join.methods <- function(fmp, fanova, cnames) {
  fmp <- join.dataframe(fmp, cnames, "FMP")
  fanova <- join.dataframe(fanova, cnames, "FANOVA")
  rbind.data.frame(fmp, fanova)
}

us.full <- join.methods(us.fmp, us.fanova, us.names)
saveRDS(us.full, "tabledfs/us.Rds")
france.full <- join.methods(france.fmp, france.fanova, france.names)
saveRDS(france.full, "tabledfs/france.Rds")
japan.full <- join.methods(japan.fmp, japan.fanova, japan.names)
saveRDS(japan.full, "tabledfs/japan.Rds")
australia.full <- join.methods(australia.fmp, australia.fanova, australia.names)
saveRDS(australia.full, "tabledfs/australia.Rds")

#rainbow plots
us.male <- readRDS("USA/USA_male.rds")
alabama <- us.male[["Alabama"]]
us.male.forecasted <- readRDS("USA/forecasted_curves_triangular_USA.rds")
alabama.forecasts <- us.male.forecasted[[1]][1:101, ]

matplot(alabama[, 1:52],type='l',col='gray86'
                  ,main=paste("State","Alabama"),xlab = 'Age'
                  ,ylab='Log mortality rate')

          colnames(alabama.forecasts)<-1:10
          Res_forcasted_curves_male<-rainbow::fts(x=0:100, y=alabama.forecasts, 
                                                  xname='Age',yname='USA Mr')
          rainbow::lines.fds(Res_forcasted_curves_male)

#generating datasets for interval forecast measures

us.base <- readRDS("USA/Emp_cov_USA.rds")
france.base <- readRDS("France/Emp_cov_France.rds")
japan.base <- readRDS("Japan/Emp_cov_Japan.rds")
australia.base <- readRDS("Australia/Emp_cov_AUS.rds")

us.fanova <- readRDS("USA/Emp_cov_USA_means.rds") 
france.fanova <- readRDS("France/Emp_cov_France_means.rds")
japan.fanova <- readRDS("Japan/Emp_cov_Japan_means.rds")
australia.fanova <- readRDS("Australia/Emp_cov_AUS_means.rds")

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

list_maker <- function(x) {as.list(rep(x, 4))}

dfs.fmp <- list(us.base, france.base, japan.base, australia.base)
dfs.fanova <- list(us.fanova, france.fanova, japan.fanova, australia.fanova)
country <- list("USA", "France", "Japan", "Australia")
types <- list(c("80%", "95%"), c("80%", "95%"), c("80%", "95%"), c("80%", "95%"))

pw.fmp <- mapply(df_toplot, x = dfs.fmp, lb = list_maker(1), ub = list_maker(2), 
                 type = types, country = country, method = list_maker("FMP"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

pw.fanova <- mapply(df_toplot, x = dfs.fanova, lb = list_maker(1), ub = list_maker(2), 
                    type = types, country = country, method = list_maker("FANOVA"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

pw <- rbind.data.frame(pw.fmp, pw.fanova) %>% 
  mutate(state = if_else(state == "Souther Australia", "Southern Australia", state))


is.fmp <- mapply(df_toplot, x = dfs.fmp, lb = list_maker(3), ub = list_maker(4), 
                 type = types, country = country, method = list_maker("FMP"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

is.fanova <- mapply(df_toplot, x = dfs.fanova, lb = list_maker(3), ub = list_maker(4), 
                    type = types, country = country, method = list_maker("FANOVA"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

is <- rbind.data.frame(is.fmp, is.fanova) %>% 
  mutate(state = if_else(state == "Souther Australia", "Southern Australia", state))

cpd.fmp <- mapply(df_toplot, x = dfs.fmp, lb = list_maker(5), ub = list_maker(6), 
                  type = types, country = country, method = list_maker("FMP"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

cpd.fanova <- mapply(df_toplot, x = dfs.fanova, lb = list_maker(5), ub = list_maker(6), 
                     type = types, country = country, method = list_maker("FANOVA"), SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .)

cpd <- rbind.data.frame(cpd.fmp, cpd.fanova) %>% 
  mutate(state = if_else(state == "Souther Australia", "Southern Australia", state))


saveRDS(pw, "intervaldfs/pw.Rds")
saveRDS(is, "intervaldfs/is.Rds")
saveRDS(cpd, "intervaldfs/cpd.Rds")


#adequating maps for proper visualization

us.map <- st_read("shp/USA/cb_2018_us_state_5m.shp") %>% 
  st_transform("+proj=longlat +datum=WGS84") %>% 
  filter(!(NAME %in% c("Puerto Rico", "American Samoa", "United States Virgin Islands", 
                       "Guam", "Commonwealth of the Northern Mariana Islands")))

france.map <- st_read("shp/France/FRA_adm2.shp")
japan.map <- st_read("shp/Japan/JPN_adm1.shp")

australia.map <- st_read("shp/Australia/STE_2021_AUST_GDA2020.shp") %>% 
  st_transform("+proj=longlat +datum=WGS84")

#usa
#checking compatibility between states names
sum(us.names.vec %in% us.map$NAME) == length(us.names.vec) #all contained

#fixing alaska
temp.geometry <- us.map[us.map$NAME == "Alaska", ]$geometry #selecting geometry of Alaska

polys <- temp.geometry %>%
  st_cast("POLYGON") #turning the multipolygon into different single polygons

real.alaska <- polys %>%
  lapply(st_coordinates) %>%
  sapply(function(df) all(df[, "X"] < 0)) %>% #filtering polygons whose lat is negative
  polys[.] %>% #selecting those polygons
  lapply(function(x) list(st_coordinates(x)[, 1:2])) %>% #extracting lat and lng
  st_multipolygon %>% #turning them into a single multipolygon
  st_sfc #making the polygon sfc instead of sfg

st_crs(real.alaska) <- st_crs(temp.geometry) #asigning crs


us.map[us.map$NAME == "Alaska", ]$geometry <- real.alaska

us.map <- us.map %>%
  select(NAME)

saveRDS(us.map, "shp/USA/usamap.Rds")

#france
france.map <- france.map %>%
  rename(NAME = NAME_2) %>%
  select(NAME)

france.map$NAME[c(17, 16, 27, 41)] <- c("Val-D'Oise", "Seine-St-Denis", "Côtes d'Armor", "Corse")

contained <- france.names.vec %in% france.map$NAME
no.accent.france <- iconv(france.map$NAME,to="ASCII//TRANSLIT")
for (i in france.names.vec[!contained]) {
  france.map$NAME[no.accent.france == iconv(i, to="ASCII//TRANSLIT")] <- i
}

saveRDS(france.map, "shp/France/francemap.Rds")

#japan
japan.map <- japan.map %>%
  rename(NAME = NAME_1) %>%
  select(NAME)
japan.map$NAME[c(13, 27)] <- c("Hyogo", "Nagasaki")
contained <- japan.names.vec %in% japan.map$NAME
japan.names.vec[!contained]

japan.map$geometry <- st_simplify(japan.map$geometry,
                                  dTolerance = 1000)

saveRDS(japan.map, "shp/Japan/japanmap.Rds")

#australia
australia.names.vec <- readRDS("names/names_regions.rds")
get.state.name <- Vectorize(function(x) {
  if (x == "South Australia") {
    return("Southern Australia")
  }
  idx <- which.max(tolower(x) == tolower(australia.names.vec))
  return(australia.names.vec[idx])
})

australia.map <- australia.map %>%
  filter(STE_NAME21 %in% STE_NAME21[-(9:10)]) %>% 
  mutate(STE_NAME21 = get.state.name(STE_NAME21)) %>% 
  rename(NAME = STE_NAME21) %>% 
  select(NAME)
  
australia.names.vec %in% australia.map$NAME
australia.map$geometry <- st_simplify(australia.map$geometry, preserveTopology = FALSE, 
                                      dTolerance = 1000)
saveRDS(australia.map, "shp/Australia/australiamap.Rds")

# f <- readRDS("names/names_regions.rds")
# f[5] <- "Southern Australia"
# saveRDS(f, "names/names_regions.rds")
# plot(australia.map$geometry)

add.centroids <- function(df, crs, path) {
  df <- cbind(df, st_centroid(st_transform(df$geometry, crs)) %>% 
                    st_transform("+proj=longlat +datum=WGS84") %>% 
                    st_coordinates)
  saveRDS(df, path)
  df
}

plotter <- function(df) {
  plot(df$geometry)
  with(df, points(X, Y, col = "red", pch = 20))
}

add.centroids(us.map, 3857, "shp/USA/usamap.Rds") %>% plotter

add.centroids(france.map, 2154, "shp/France/francemap.Rds") %>% plotter

add.centroids(japan.map, 6677, "shp/Japan/japanmap.Rds") %>% plotter

add.centroids(australia.map, 3112, "shp/Australia/australiamap.Rds") %>% plotter

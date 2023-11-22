library(leaflet) 
library(sf)

us.map <- readRDS("shp/USA/usamap.Rds")
france.map <- readRDS("shp/France/francemap.Rds")
japan.map <- readRDS("shp/Japan/japanmap.Rds")

map.list <- list(USA = us.map, France = france.map, 
                 Japan = japan.map)

fill.color <- c("gray", "gold")

render.map <- function(country, state) {

  temp.map.df <- map.list[[country]] %>%
  mutate(
    NAME == if_else(rep(country == "JAPAN", nrow(.)), gsub("\\t", "", NAME), NAME),
    id = (NAME == state) + 1)

  fill <- fill.color[temp.map.df$id]
  selected.coords <- as.numeric(temp.map.df[temp.map.df$NAME == state, c("X", "Y")])[1:2]
  
  leaflet(temp.map.df) %>%
    addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
    addPolygons(
      fillColor = fill, color = "black", fillOpacity = 0.5, weight = 1) %>% 
    addPopups(lng = selected.coords[1], lat = selected.coords[2], 
              popup = paste0("<b>", state, "</b>"))
}

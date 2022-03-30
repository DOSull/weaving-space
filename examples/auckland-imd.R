source("weaving-space-utils.R")
source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weave-map.R")

library(sf)        # vector spatial data
library(tmap)      # thematic maps
library(dplyr)     # data wrangling

ak <- st_read("data/imd-auckland-2018.gpkg")
ak %>% plot()

weave_unit <- get_weave_unit(spacing = 50, aspect = 0.8, type = "twill", n = 3,
                             strands = "abc|defg", crs = 2193)
weave_unit %>% plot_unit()

weave <- weave_layer(weave_unit, ak, angle = 30)
# weave[, 1] %>% plot()

lyrs <- split(weave, as.factor(weave$strand))

names(weave)
tmap_mode("view")
tm_shape(lyrs$a) + 
  tm_fill(col = "Rank_Emplo", palette = "-BrBG", title = "Employment", style = "cont") +
  tm_shape(lyrs$b) + 
  tm_fill(col = "Rank_Incom", palette = "-PiYG", title = "Income", style = "cont") +
  tm_shape(lyrs$c) + 
  tm_fill(col = "Rank_Crime", palette = "-RdBu", title = "Crime", style = "cont") +
  tm_shape(lyrs$d) + 
  tm_fill(col = "Rank_Housi", palette = "-RdGy", title = "Housing", style = "cont") +
  tm_shape(lyrs$e) + 
  tm_fill(col = "Rank_Healt", palette = "-PRGn", title = "Health", style = "cont") +
  tm_shape(lyrs$f) + 
  tm_fill(col = "Rank_Educa", palette = "-RdYlBu", title = "Education", style = "cont") +
  tm_shape(lyrs$g) + 
  tm_fill(col = "Rank_Acces", palette = "-PuOr", title = "Access", style = "cont")
  


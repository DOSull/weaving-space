source("weaving-space-utils.R")
source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weave-map.R")

library(sf)        # vector spatial data
library(tmap)      # thematic maps
library(dplyr)     # data wrangling

ak <- st_read("data/imd-auckland-2018.gpkg")
ak %>% plot()

weave_unit <- get_biaxial_weave_unit(spacing = 50, aspect = 0.8, 
                                     strands = "abc|defg", crs = 2193)
weave_unit$primitive %>% plot()

weave <- weave_layer(weave_unit, ak, angle = 30)
weave[, 1] %>% plot()

lyrs <- split(weave, as.factor(weave$strand))

names(weave)

tm_shape(lyrs$a) + 
  tm_fill(col = "Decile_Emp", palette = "-BrBG", title = "Employment", n = 3, style = "quantile") +
  tm_shape(lyrs$b) + 
  tm_fill(col = "Decile_Inc", palette = "-PiYG", title = "Income", n = 3, style = "quantile") +
  tm_shape(lyrs$c) + 
  tm_fill(col = "Decile_Cri", palette = "-RdBu", title = "Crime", n = 3, style = "quantile") +
  tm_shape(lyrs$d) + 
  tm_fill(col = "Decile_Hou", palette = "-RdGy", title = "Housing", n = 3, style = "quantile") +
  tm_shape(lyrs$e) + 
  tm_fill(col = "Decile_Hea", palette = "-PRGn", title = "Health", n = 3, style = "quantile") +
  tm_shape(lyrs$f) + 
  tm_fill(col = "Decile_Edu", palette = "-RdYlBu", title = "Education", n = 3, style = "quantile") +
  tm_shape(lyrs$g) + 
  tm_fill(col = "Decile_Acc", palette = "-PuOr", title = "Access", n = 3, style = "quantile")
  

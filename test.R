source("weaving-space-utils.R")
source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weave-map.R")

library(sf)        # vector spatial data
library(tmap)      # thematic maps
library(dplyr)     # data wrangling

region <- st_read("data/vax-auckland-20211006.gpkg")

tu <- matrix(c(0, 1, 1, 0, 1, 1,
               1, 1, 0, 1, 0, 1,
               1, 0, 1, 0, 1, 0,
               0, 1, 0, 0, 0, 1,
               1, 0, 1, 0, 1, 0,
               1, 1, 0, 1, 0, 1), 6, 6, byrow = TRUE)

unit1 <- get_biaxial_weave_unit(type = "this", tie_up = tu, 
                                spacing = 100, aspect = 0.75,
                                strands = "a|c")

weave1 <- weave_layer(unit1, region)

tmap_mode("view")
tm_shape(weave1) + 
  tm_fill(col = "strand", palette = "Paired") +
  tm_layout(title = "test weave")


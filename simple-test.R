source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weaving-space-utils.R")
source("render-weave-grids.R")
source("weave-map.R")

ak <- st_read("data/vax-auckland-20211006.gpkg")

bi_weave <- get_biaxial_weave_unit(spacing = 250, aspect = 0.8, crs = 2193)
tri_weave <- get_triaxial_weave_unit(spacing = 250, type = "cube", strands = "a|b|c",
                                     aspect = 0.5, crs = 2193)

w1 <- weave_layer(bi_weave, ak, angle = 30)
w2 <- weave_layer(tri_weave, ak)

tmap_options(check.and.fix = TRUE)
tm_shape(w1) + 
  tm_fill(col = "strand", palette = "Set1")

tm_shape(w2) +
  tm_fill(col = "strand", palette = "Set1")

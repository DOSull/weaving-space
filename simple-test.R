source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weaving-space-utils.R")
source("render-weave-grids.R")
source("weave-map.R")

ak <- st_read("data/vax-auckland-20211006.gpkg")

bi_weave <- get_biaxial_weave_unit(spacing = 250, aspect = 0.8, crs = 2193)
w1 <- weave_layer(bi_weave, ak, angle = 30)

tmap_options(check.and.fix = TRUE)
tm_shape(w1) + 
  tm_fill(col = "strand", palette = "Set1")


tri_weave <- get_triaxial_weave_unit(spacing = 100, type = "cube", strands = "ab-|cd-|ef-",
                                     aspect = 1, crs = 2193)
w2 <- weave_layer(tri_weave, ak, angle = 15)
tm_shape(w2) + 
  tm_fill(col = "strand", style = "cat")

layers <- w2 %>% split(as.factor(w2$strand))
tm_shape(layers$a, name = "Pākehā") +
  tm_fill(col = "pEuropean", palette = "Greys", title = "%Pākehā") +
  tm_shape(layers$b, name = "Māori") + 
  tm_fill(col = "pMaori", palette = "Reds", title = "%Māori") +
  tm_shape(layers$c, name = "Pasifika") +
  tm_fill(col = "pPacific", palette = "Purples", title = "%Pasifika") +
  tm_shape(layers$d, name = "Asian") +
  tm_fill(col = "pAsian", palette = "Greens", title = "%Asian") +
  tm_shape(layers$e, name = "Dose 1") +
  tm_fill(col = "dose1_uptake", palette = "magma", title = "Dose 1") +
  tm_shape(layers$f, name = "Dose 2") + 
  tm_fill(col = "dose2_uptake", palette = "magma", title = "Dose 2")

write_weave_layers(w2, ak, "other-data/w2.gpkg")


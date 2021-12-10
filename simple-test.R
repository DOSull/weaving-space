library(tmap)

source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weaving-space-utils.R")
source("render-weave-grids.R")
source("weave-map.R")

ak <- st_read("data/vax-auckland-20211006.gpkg")

get_weave_unit(spacing = 250, aspect = 0.8) %>% plot_unit()
get_weave_unit(spacing = 250, aspect = 0.8, strands = "abc|de") %>% plot_unit()
get_weave_unit(spacing = 250, aspect = 0.8, strands = "a-(bc)|d-e") %>% plot_unit()
get_weave_unit(spacing = 250, aspect = 0.8, type = "twill", strands = "a-b|c-d") %>% plot_unit()
get_weave_unit(spacing = 250, aspect = 0.8, type = "basket", strands = "ab|cd") %>% plot_unit()

# weird not weave weave
bi_weave <- get_weave_unit(strands = "(ab)|(cd)",
                           spacing = 200, aspect = 0.7, margin = -0.2, 
                           crs = 2193)
tmap_mode("plot")
bi_weave$weave_unit %>% tm_shape() + tm_fill(col = "strand", alpha = 0.75)
w1 <- weave_layer(bi_weave, ak, angle = 30)

layers <- w1 %>% split(as.factor(w1$strand))
tmap_mode("view")
tmap_options(check.and.fix = TRUE)
tm_shape(layers$a, name = "Pākehā") +
  tm_fill(col = "pEuropean", palette = "Greys", style = "cont", 
          title = "%Pākehā", alpha = 0.65) +
  tm_shape(layers$b, name = "Māori") + 
  tm_fill(col = "pMaori", palette = "YlOrRd", style = "cont",
          title = "%Māori", alpha = 0.65) +
  tm_shape(layers$c, name = "Pasifika") +
  tm_fill(col = "pPacific", palette = "RdPu", style = "cont",
          title = "%Pasifika", alpha = 0.65) +
  tm_shape(layers$d, name = "Asian") +
  tm_fill(col = "pAsian", palette = "YlGnBu", style = "cont",
          title = "%Asian", alpha = 0.65)
  

tmap_mode("plot")
tri_weave <- get_weave_unit(spacing = 200, type = "cube", strands = "a-b|c-d|e-f",
                            aspect = 0.8, crs = 2193)
tm_shape(tri_weave$weave_unit) +
  tm_fill(col = "strand") +
  tm_shape(tri_weave$tile) +
  tm_borders(col = "red")

w2 <- weave_layer(tri_weave, ak, angle = 15) %>%
  st_make_valid()
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


library(tmap)
source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weaving-space-utils.R")
source("render-weave-grids.R")
source("weave-map.R")
sf_use_s2(FALSE)

bc <- st_read("data/data-for-dos-test.gpkg")

unit <- get_weave_unit(type = "twill", 
                       spacing = 36000, 
                       strands = "(abcd)-|(efgh)-", 
                       crs = st_crs(bc))
w_bc <- weave_layer(unit, bc, angle = 30)
fabric <- w_bc %>% 
  st_union() %>%
  st_buffer(0.1) %>%
  st_buffer(-0.1)
squares <- st_difference(st_combine(bc) %>% 
                           st_buffer(-1), 
                         fabric)

sst_write(w_bc, "data/weave_bc.gpkg", delete_dsn = TRUE)
st_write(fabric, "data/fabric.gpkg", delete_dsn = TRUE)
st_write(squares, "data/squares.gpkg", delete_dsn = TRUE)

unit2 <- get_weave_unit(type = "basket", n = 8, strands = "abcd----|efgh----", spacing = 9000)
w_bc2 <- weave_layer(unit2, bc, angle = 30)
fabric2 <- w_bc2 %>% st_union()
squares2 <- st_difference(bc, fabric2)

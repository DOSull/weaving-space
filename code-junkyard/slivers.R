library(tmap)
source("biaxial-weave-units.R")
source("triaxial-weave-units.R")
source("weaving-space-utils.R")
source("render-weave-grids.R")
source("weave-map.R")
sf_use_s2(FALSE)

bc <- st_read("data/data-for-dos-test.gpkg")
bc_bb <- bc %>% 
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_sf()
unit <- get_weave_unit(type = "basket", n = 8, spacing = 9000, 
                       strands = "abcd----|efgh----", crs = st_crs(bc))
w_bc <- weave_layer(unit, bc_bb, angle = 30)
fabric <- w_bc %>% 
  st_buffer(.1) %>%
  st_union() %>%
  st_buffer(-.1) %>%
  st_sf()
squares <- st_difference(bc_bb, fabric) %>% 
  st_cast("POLYGON") %>% 
  st_sf() %>% 
  slice(-1) %>% 
  st_filter(bc)
pts <- squares %>%
  st_centroid() %>%
  st_filter(bc, .predicate = st_within)
squares <- squares %>%
  st_intersection(bc)
w_bc <- w_bc %>% 
  st_intersection(bc)
tm_shape(squares) + 
  tm_fill(col = "pink") + 
  tm_shape(pts) + 
  tm_dots() + 
  tm_shape(w_bc) + 
  tm_fill(col = "LHA_Name") +
  tm_layout(legend.outside = TRUE)


st_write(w_bc, "data/weave_bc.gpkg", delete_dsn = TRUE)
st_write(fabric, "data/fabric.gpkg", delete_dsn = TRUE)
st_write(squares, "data/squares.gpkg", delete_dsn = TRUE)

unit2 <- get_weave_unit(type = "basket", n = 8, strands = "abcd----|efgh----", spacing = 9000)
w_bc2 <- weave_layer(unit2, bc, angle = 30)
fabric2 <- w_bc2 %>% st_union()
squares2 <- st_difference(bc, fabric2)

unit <- get_weave_unit(type = "basket", n = 8, spacing = 9000, 
                       strands = "abcd----|efgh----", crs = st_crs(bc))
hole <- unit$tile %>% 
  st_difference(st_union(unit$weave_unit)) %>% 
  st_sf() %>%
  mutate(strand = "hole")
hole_centre <- hole %>%
  st_buffer(-12000) %>%
  mutate(strand = "centre")
unit$weave_unit <- unit$weave_unit %>% 
  bind_rows(hole, hole_centre)
unit %>% plot_unit()
